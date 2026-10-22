# Compile expressions to Julia functions.
#
# The expression DAG is turned into a Julia `Expr` and compiled with
# RuntimeGeneratedFunctions, so the result can be called right away (no world-age
# problem) and runs at the speed of hand-written code. The generated code follows the
# semantics of calling an expression (see `evaluate`):
#   * exact constants are Julia integers and rationals, converted to the arguments'
#     floating point type when there is one (so BigFloat keeps its precision);
#     `F = _liftfloat(...)` is a compile-time constant, so for Float64 they fold;
#   * integer powers are literal powers (`x^3`), negative ones use `inv`;
#   * a negative base with a fractional power gives the principal complex value.
# Subexpressions that occur more than once are computed once, in local variables.

"""
    lambdify(ex, vars...)
    lambdify(ex, vars::Union{Tuple,AbstractVector})
    lambdify(ex)

Compile `ex` into a Julia function of the symbols `vars`, in that order. With no `vars`,
`ex` may have at most one free symbol. `ex` may also be an array or tuple of
expressions, and the function then returns an array of the same shape (or a tuple),
with common subexpressions shared between the entries.

The function computes what calling the expression computes, at the speed of compiled
Julia code: `lambdify(x^2 + 1)(3) === 10`, and `lambdify(sin(x) * y, x, y)(1.0, 2.0)`.
Arguments may be any number type (`BigFloat`, dual numbers, intervals, ...).

`lambdify_expr` gives the generated Julia expression.
"""
lambdify(ex, vars::Sym...) = _lambdify(ex, collect(Sym, vars))
lambdify(ex, vars::Union{Tuple,AbstractVector}) = _lambdify(ex, Sym[Sym(v) for v in vars])
function lambdify(ex)
    vars = free_symbols(ex)
    length(vars) <= 1 || throw(ArgumentError(
        "the expression has free symbols $(join(vars, ", ")); give their order, e.g. lambdify(ex, x, y)"))
    _lambdify(ex, vars)
end

_lambdify(ex, vars::Vector{Sym}) = RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, lambdify_expr(ex, vars))

"""
    lambdify_expr(ex, vars)

The Julia expression (an anonymous function) that [`lambdify`](@ref) compiles.
"""
function lambdify_expr(ex, vars::Union{Tuple,AbstractVector})
    vars = Sym[Sym(v) for v in vars]
    all(v -> v.kind == SYMBOL, vars) || throw(ArgumentError("variables must be symbols"))
    args = [Symbol(:arg, i) for i in eachindex(vars)]
    cg = _Codegen(Dict{Sym,Symbol}(v => a for (v, a) in zip(vars, args)))
    outs = _outputs(ex)
    for o in outs
        _count!(cg, o)
    end
    vals = Any[_gen(cg, o) for o in outs]
    result = ex isa Tuple ? Expr(:tuple, vals...) :
             ex isa AbstractArray ? :(reshape($(Expr(:vect, vals...)), $(size(ex)))) :
             only(vals)
    lift = isempty(args) ? :(F = nothing) :
           :(F = $(GlobalRef(SymJ, :_liftfloat))(promote_type($([:(typeof($a)) for a in args]...))))
    Expr(:->, Expr(:tuple, args...), Expr(:block, lift, cg.stmts..., result))
end

_outputs(ex::Union{Tuple,AbstractArray}) = Sym[Sym(e) for e in ex]
_outputs(ex::Poly) = Sym[Sym(ex)]
_outputs(ex) = Sym[Sym(ex)]

mutable struct _Codegen
    vars::Dict{Sym,Symbol}
    counts::Dict{Sym,Int}
    names::Dict{Sym,Symbol}
    stmts::Vector{Any}
end
_Codegen(vars) = _Codegen(vars, Dict{Sym,Int}(), Dict{Sym,Symbol}(), Any[])

_children(ex::Sym) = ex.kind == MUL ? vcat(ex.args, ex.exps) : ex.args

# how often each node is used; children of a repeated node are counted once
function _count!(cg::_Codegen, ex::Sym)
    n = get(cg.counts, ex, 0)
    cg.counts[ex] = n + 1
    n == 0 && foreach(c -> _count!(cg, c), _children(ex))
end

function _gen(cg::_Codegen, ex::Sym)
    name = get(cg.names, ex, nothing)
    name === nothing || return name
    code = _code(cg, ex)
    if get(cg.counts, ex, 0) > 1 && ex.kind in (ADD, MUL, POW, FN)
        name = Symbol(:t, length(cg.names) + 1)
        push!(cg.stmts, :($name = $code))
        cg.names[ex] = name
        return name
    end
    code
end

const _CNUM = GlobalRef(SymJ, :_cnum)
@inline _cnum(::Nothing, v) = v
@inline _cnum(::Type{F}, v) where {F} = convert(F, v)

_numcode(n::SymNumber) = :($_CNUM(F, $(_jlnum(n))))

function _code(cg::_Codegen, ex::Sym)
    k = ex.kind
    k == NUMBER && return _numcode(ex.num)
    if k == SYMBOL
        a = get(cg.vars, ex, nothing)
        a === nothing && throw(ArgumentError("free symbol $ex is not among the variables"))
        return a
    elseif k == CONSTANT
        ex.name === :pi && return :($_CNUM(F, $(π)))
        ex.name === :E && return :($_CNUM(F, $(ℯ)))
        throw(ArgumentError("no numeric value for constant $ex"))
    elseif k == ADD
        terms = Any[_scaled(ex.coefs[i], _gen(cg, ex.args[i])) for i in eachindex(ex.args)]
        iszero(ex.num) || push!(terms, _numcode(ex.num))
        return length(terms) == 1 ? terms[1] : Expr(:call, :+, terms...)
    elseif k == MUL
        factors = Any[_powcode(cg, ex.args[i], ex.exps[i]) for i in eachindex(ex.args)]
        _isN1(ex.num) || pushfirst!(factors, _numcode(ex.num))
        return length(factors) == 1 ? factors[1] : Expr(:call, :*, factors...)
    elseif k == POW
        return _powcode(cg, ex.args[1], ex.args[2])
    end
    rule = get(FUNCTION_RULES, ex.name, nothing)
    (rule === nothing || rule.julia === nothing) &&
        throw(ArgumentError("$(ex.name) has no numeric method; register one with `register_function!`"))
    Expr(:call, _fnref(rule.julia), Any[_intarg(a) ? _toint(a.num) : _gen(cg, a) for a in ex.args]...)   # see _intarg
end

_scaled(c::SymNumber, t) = _isN1(c) ? t : c == NM1 && isexact(c) ? :(-$t) : :($(_numcode(c)) * $t)

# a reference to a named function where possible (readable code), else the object
function _fnref(f)
    m, n = parentmodule(f), nameof(f)
    isdefined(m, n) && getfield(m, n) === f ? GlobalRef(m, n) : f
end

const _SQRT = GlobalRef(SymJ, :_cg_sqrt)
const _RPOW = GlobalRef(SymJ, :_cg_pow)
@inline _cg_sqrt(b) = b isa Real && b < 0 ? sqrt(complex(b)) : sqrt(b)
@inline _cg_pow(b, e) = b isa Real && e isa Real && b < 0 && !isinteger(e) ? complex(b)^e : b^e

function _powcode(cg::_Codegen, b::Sym, e::Sym)
    _isone(e) && return _gen(cg, b)
    b === E && return :(exp($(_gen(cg, e))))
    bc = _gen(cg, b)
    if e.kind == NUMBER && isexact(e.num)
        if _isint(e.num)
            k = _toint(e.num)
            return k >= 0 ? :($bc^$k) : k == -1 ? :(inv($bc)) : :(inv($bc)^$(-k))
        end
        e.num == SymNumber(1 // 2) && return :($_SQRT($bc))
        e.num == SymNumber(-1 // 2) && return :(inv($_SQRT($bc)))
        ec = :($_CNUM(F === nothing ? Float64 : F, $(_jlnum(e.num))))
    else
        ec = _gen(cg, e)
    end
    # a positive constant base needs no complex branch
    positive = (b.kind == NUMBER && sign(b.num) > 0) || b.kind == CONSTANT
    positive ? :($bc^$ec) : :($_RPOW($bc, $ec))
end
