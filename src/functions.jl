# Unevaluated function applications (kind FN) and the rules that simplify them.
#
# A function is identified by its name. Behaviour is looked up in `FUNCTION_RULES`:
#   julia: a Julia function for numeric evaluation (and TermInterface.operation)
#   eval:  args -> Union{Nothing,Sym}; return a simplified value or `nothing`
#   diff:  (args, i) -> partial derivative with respect to the i-th argument
# Names with no rule (user functions) stay unevaluated.

"""
    SymFunction(name)

A symbolic function. Calling it builds an unevaluated application:

```julia
f = SymFunction(:f)     # or `@syms f()`
f(x, 2)                 # f(x, 2)
```
"""
struct SymFunction <: AbstractSymbolic
    name::Symbol
end
(f::SymFunction)(args...) = apply_function(f.name, Sym[Sym(a) for a in args])
Base.show(io::IO, f::SymFunction) = print(io, f.name)
Base.show(io::IO, ::MIME"text/plain", f::SymFunction) = show(io, f)

# f' (f'', …): the derivative of a function of one argument; f'(x) is diff(f(x), x) and
# f'(0) the derivative at 0 (for initial conditions, `y'(0) => 1`)
struct SymDerivative <: AbstractSymbolic
    f::SymFunction
    k::Int
end
Base.adjoint(f::SymFunction) = SymDerivative(f, 1)
Base.adjoint(d::SymDerivative) = SymDerivative(d.f, d.k + 1)
function (d::SymDerivative)(a)
    a = Sym(a)
    t = a.kind == SYMBOL ? a : Sym(:x)
    e = d.f(t)
    for _ in 1:d.k
        e = _mkfn_node(:Derivative, Sym[e, t])
    end
    a.kind == SYMBOL ? e : apply_function(:Derivative, Sym[e.args[1], t, a])
end
Base.show(io::IO, d::SymDerivative) = print(io, d.f.name, "'"^d.k)
Base.show(io::IO, ::MIME"text/plain", d::SymDerivative) = show(io, d)

struct FunctionRule
    julia::Union{Nothing,Function}
    eval::Union{Nothing,Function}
    diff::Union{Nothing,Function}
end

const FUNCTION_RULES = Dict{Symbol,FunctionRule}()

"""
    register_function!(name; julia=nothing, eval=nothing, diff=nothing)

Attach behaviour to the function `name`: a Julia function used for numeric evaluation,
a simplification rule `args -> Union{Nothing,Sym}`, and a partial-derivative rule
`(args, i) -> Union{Nothing,Sym}` (`nothing` for a derivative with no closed form,
which is left unevaluated).
"""
function register_function!(name::Symbol; julia = nothing, eval = nothing, diff = nothing)
    FUNCTION_RULES[name] = FunctionRule(julia, eval, diff)
    name
end

"""
    apply_function(name::Symbol, args::Vector{Sym})

Apply the function `name`, simplifying with its registered rule if any.
"""
function apply_function(name::Symbol, args::Vector{Sym})
    rule = get(FUNCTION_RULES, name, nothing)
    if rule !== nothing
        # Float arguments are values: evaluate first, before symbolic rules such as
        # asin(-u) = -asin(u), which do not hold for Julia's principal values on a
        # branch cut (asin(-17.0 + 0im) is not -asin(17.0 + 0im)).
        if rule.julia !== nothing && _float_args(args)
            v = _try_numeric(rule.julia, args)
            v === nothing && (v = _try_numeric_complex(rule.julia, args))
            v === nothing || return v
        end
        if rule.eval !== nothing
            r = rule.eval(args)
            r === nothing || return r::Sym
        end
    end
    _mkfn_node(name, args)
end

_float_args(args) = all(a -> a.kind == NUMBER, args) && any(a -> isfloat(a.num), args)

# outside the real domain (asin(17.0)): the principal complex value
function _try_numeric_complex(f, args)
    try
        return Sym(f(map(a -> complex(Float64(a.num)), args)...))
    catch err
        err isa Union{DomainError,MethodError,InexactError} || rethrow()
        return nothing
    end
end

function _try_numeric(f, args)
    try
        return Sym(f(map(a -> Float64(a.num), args)...))
    catch err
        err isa Union{DomainError,MethodError,InexactError} || rethrow()   # e.g. factorial(2.5)
        return nothing
    end
end

_fn_operation(name::Symbol) =
    (r = get(FUNCTION_RULES, name, nothing); r === nothing || r.julia === nothing ? SymFunction(name) : r.julia)

## helpers for rules

# Would `-x` print more simply than `x`? Exactly one of `x`, `-x` answers true (or both
# false), which makes it usable to canonicalize odd and even functions.
function _could_extract_minus(x::Sym)
    k = x.kind
    (k == NUMBER || k == MUL) && return sign(x.num) < 0
    k == ADD && return iszero(x.num) ? sign(x.coefs[1]) < 0 : sign(x.num) < 0
    false
end

# x == r*π for an exact rational r
function _pi_multiple(x::Sym)
    iszero(x) && isexact(x.num) && return Rational{BigInt}(0)
    x === PI && return Rational{BigInt}(1)
    if x.kind == MUL && length(x.args) == 1 && x.args[1] === PI && _isone(x.exps[1]) && isexact(x.num)
        return _rat(x.num)
    end
    nothing
end

# sin(rπ) when it has a simple closed form
function _sinpi(r::Rational{BigInt})
    r = mod(r, 2)
    neg = r >= 1
    neg && (r -= 1)
    r > 1 // 2 && (r = 1 - r)
    v = r == 0 ? ZERO :
        r == 1 // 6 ? HALF :
        r == 1 // 4 ? sqrt(Sym(2)) / 2 :
        r == 1 // 3 ? sqrt(Sym(3)) / 2 :
        r == 1 // 2 ? ONE : nothing
    v === nothing && return nothing
    neg ? -v : v
end

_isfn(x::Sym, name::Symbol) = x.kind == FN && x.name === name && length(x.args) == 1

_odd(name, a) = _could_extract_minus(a) ? -apply_function(name, Sym[-a]) : nothing
_even(name, a) = _could_extract_minus(a) ? apply_function(name, Sym[-a]) : nothing

function _unary(f)
    args -> length(args) == 1 ? f(args[1]) : throw(ArgumentError("expected one argument, got $(length(args))"))
end

## elementary functions

register_function!(:sin; julia = sin,
    eval = _unary(a -> begin
        r = _pi_multiple(a)
        r !== nothing && (v = _sinpi(r); v === nothing || return v)
        _isfn(a, :asin) && return a.args[1]
        _odd(:sin, a)
    end),
    diff = (args, i) -> cos(args[1]))

register_function!(:cos; julia = cos,
    eval = _unary(a -> begin
        r = _pi_multiple(a)
        r !== nothing && (v = _sinpi(r + 1 // 2); v === nothing || return v)
        _isfn(a, :acos) && return a.args[1]
        _even(:cos, a)
    end),
    diff = (args, i) -> -sin(args[1]))

register_function!(:tan; julia = tan,
    eval = _unary(a -> begin
        r = _pi_multiple(a)
        if r !== nothing
            s, c = _sinpi(r), _sinpi(r + 1 // 2)
            s !== nothing && c !== nothing && !iszero(c) && return s / c
        end
        _isfn(a, :atan) && return a.args[1]
        _odd(:tan, a)
    end),
    diff = (args, i) -> 1 + tan(args[1])^2)

# exact values of asin and atan at the standard points (odd: the negatives follow);
# acos(v) = π/2 - asin(v)
const _ASIN_VALUES = Pair{Sym,Sym}[ZERO => ZERO, HALF => PI / 6, sqrt(Sym(2)) / 2 => PI / 4,
                                   sqrt(Sym(3)) / 2 => PI / 3, ONE => PI / 2]
const _ATAN_VALUES = Pair{Sym,Sym}[ZERO => ZERO, sqrt(Sym(3)) / 3 => PI / 6, ONE => PI / 4,
                                   sqrt(Sym(3)) => PI / 3]
function _table_value(table, a::Sym)
    for (v, θ) in table
        _eq(a, v) && return θ
        _eq(a, -v) && return -θ
    end
    nothing
end

register_function!(:asin; julia = asin,
    eval = _unary(a -> something(_table_value(_ASIN_VALUES, a), _odd(:asin, a), Some(nothing))),
    diff = (args, i) -> inv(sqrt(1 - args[1]^2)))

register_function!(:acos; julia = acos,
    eval = _unary(a -> (θ = _table_value(_ASIN_VALUES, a); θ === nothing ? nothing : PI / 2 - θ)),
    diff = (args, i) -> -inv(sqrt(1 - args[1]^2)))

register_function!(:atan; julia = atan,
    eval = _unary(a -> something(_table_value(_ATAN_VALUES, a), _odd(:atan, a), Some(nothing))),
    diff = (args, i) -> inv(1 + args[1]^2))

register_function!(:sinh; julia = sinh,
    eval = _unary(a -> iszero(a) ? ZERO : _odd(:sinh, a)),
    diff = (args, i) -> cosh(args[1]))

register_function!(:cosh; julia = cosh,
    eval = _unary(a -> iszero(a) ? ONE : _even(:cosh, a)),
    diff = (args, i) -> sinh(args[1]))

register_function!(:tanh; julia = tanh,
    eval = _unary(a -> iszero(a) ? ZERO : _odd(:tanh, a)),
    diff = (args, i) -> 1 - tanh(args[1])^2)

register_function!(:asinh; julia = asinh,
    eval = _unary(a -> iszero(a) ? ZERO : _odd(:asinh, a)),
    diff = (args, i) -> inv(sqrt(args[1]^2 + 1)))

register_function!(:acosh; julia = acosh,
    eval = _unary(a -> isone(a) ? ZERO : nothing),
    diff = (args, i) -> inv(sqrt(args[1] - 1) * sqrt(args[1] + 1)))

register_function!(:atanh; julia = atanh,
    eval = _unary(a -> iszero(a) ? ZERO : _odd(:atanh, a)),
    diff = (args, i) -> inv(1 - args[1]^2))

register_function!(:abs; julia = abs,
    eval = _unary(a -> begin
        a.kind == NUMBER && return Sym(abs(a.num))
        _isfn(a, :abs) && return a
        _even(:abs, a)
    end),
    diff = (args, i) -> args[1] / abs(args[1]))

# log(x) and log(b, x) (base first, as in Base.log); log(b, x) is kept unevaluated
register_function!(:log; julia = log,
    eval = function (args)
        if length(args) == 1
            a = args[1]
            isone(a) && return ZERO
            a === E && return ONE
            if a.kind == POW && a.args[1] === E && a.args[2].kind == NUMBER
                return a.args[2]                 # log(e^r) == r for real r
            end
            return nothing
        end
        length(args) == 2 || throw(ArgumentError("log takes one or two arguments"))
        b, x = args
        isone(x) && return ZERO
        _eq(b, x) && return ONE
        b === E && return apply_function(:log, Sym[x])
        nothing
    end,
    diff = function (args, i)
        length(args) == 1 && return inv(args[1])
        b, x = args
        i == 2 ? inv(x * log(b)) : -log(x) / (b * log(b)^2)
    end)

## reciprocal and inverse-reciprocal functions
#
# These stay as their own functions (sec(x), not 1/cos(x)). Their exact values come
# from the rules of the function they are defined through, accepted only when that
# function no longer appears: sec(π/3) = 1/cos(π/3) = 2 and sec(acos(x)) = 1/x, while
# sec(x) stays sec(x). A pole (sec(π/2)) stays unevaluated.

_hasfn(ex::Sym, name::Symbol) =
    (ex.kind == FN && ex.name === name) || any(a -> _hasfn(a, name), ex.args) || any(a -> _hasfn(a, name), ex.exps)

# f(a) = num(a) / den(a), with num and den among the existing functions (num may be ONE)
function _ratio_rule(name::Symbol, num::Union{Nothing,Symbol}, den::Symbol, parity::Symbol)
    _unary(function (a)
        d = apply_function(den, Sym[a])
        n = num === nothing ? ONE : apply_function(num, Sym[a])
        if !_hasfn(d, den) && (num === nothing || !_hasfn(n, num))
            iszero(d) || return n / d
            return nothing                                # a pole
        end
        parity === :odd ? _odd(name, a) : parity === :even ? _even(name, a) : nothing
    end)
end

# f(a) = g(1/a), with g an existing inverse function
function _invrecip_rule(name::Symbol, g::Symbol, parity::Symbol, atzero)
    _unary(function (a)
        iszero(a) && return atzero
        r = apply_function(g, Sym[inv(a)])
        _hasfn(r, g) || return r
        parity === :odd ? _odd(name, a) : nothing
    end)
end

register_function!(:sec; julia = sec, eval = _ratio_rule(:sec, nothing, :cos, :even),
    diff = (args, i) -> sec(args[1]) * tan(args[1]))
register_function!(:csc; julia = csc, eval = _ratio_rule(:csc, nothing, :sin, :odd),
    diff = (args, i) -> -csc(args[1]) * cot(args[1]))
register_function!(:cot; julia = cot, eval = _ratio_rule(:cot, :cos, :sin, :odd),
    diff = (args, i) -> -1 - cot(args[1])^2)
register_function!(:sech; julia = sech, eval = _ratio_rule(:sech, nothing, :cosh, :even),
    diff = (args, i) -> -sech(args[1]) * tanh(args[1]))
register_function!(:csch; julia = csch, eval = _ratio_rule(:csch, nothing, :sinh, :odd),
    diff = (args, i) -> -csch(args[1]) * coth(args[1]))
register_function!(:coth; julia = coth, eval = _ratio_rule(:coth, :cosh, :sinh, :odd),
    diff = (args, i) -> 1 - coth(args[1])^2)

# as in Base: asec(x) = acos(1/x), acot(x) = atan(1/x), asech(x) = acosh(1/x), ...
register_function!(:asec; julia = asec, eval = _invrecip_rule(:asec, :acos, :none, nothing),
    diff = (args, i) -> inv(args[1]^2 * sqrt(1 - inv(args[1]^2))))
register_function!(:acsc; julia = acsc, eval = _invrecip_rule(:acsc, :asin, :odd, nothing),
    diff = (args, i) -> -inv(args[1]^2 * sqrt(1 - inv(args[1]^2))))
register_function!(:acot; julia = acot, eval = _invrecip_rule(:acot, :atan, :odd, PI / 2),
    diff = (args, i) -> -inv(1 + args[1]^2))
register_function!(:asech; julia = asech, eval = _invrecip_rule(:asech, :acosh, :none, nothing),
    diff = (args, i) -> -inv(args[1]^2 * sqrt(inv(args[1]) - 1) * sqrt(inv(args[1]) + 1)))
register_function!(:acsch; julia = acsch, eval = _invrecip_rule(:acsch, :asinh, :odd, nothing),
    diff = (args, i) -> -inv(args[1]^2 * sqrt(1 + inv(args[1]^2))))
register_function!(:acoth; julia = acoth, eval = _invrecip_rule(:acoth, :atanh, :odd, nothing),
    diff = (args, i) -> inv(1 - args[1]^2))

## real cube root: differs from x^(1/3) (the principal root) for negative x
register_function!(:cbrt; julia = cbrt,
    eval = _unary(function (a)
        (a.kind == NUMBER && isexact(a.num)) || return _odd(:cbrt, a)
        sign(a.num) < 0 ? -pow(-a, Sym(SymNumber(1 // 3))) : pow(a, Sym(SymNumber(1 // 3)))
    end),
    diff = (args, i) -> inv(3 * cbrt(args[1])^2))

## piecewise constant functions. Their derivatives are left unevaluated, since they
## are zero only away from the jumps.

# the value of a constant expression, when safely away from `avoid(v)` (a jump)
function _constant_float(a::Sym, avoid)
    (a.kind == NUMBER || !_hassymbol(a)) || return nothing
    v = try
        N(a)
    catch
        return nothing
    end
    v isa Real && abs(v - avoid(v)) > 1e-10 * max(1, abs(v)) ? v : nothing
end

register_function!(:sign; julia = sign,
    eval = _unary(function (a)
        a.kind == NUMBER && isexact(a.num) && return Sym(sign(a.num))
        _isfn(a, :sign) && return a
        v = _constant_float(a, _ -> 0.0)
        v === nothing ? _odd(:sign, a) : Sym(Int(sign(v)))
    end))

# the jumps to stay away from: integers, or for round the half-integers
for (f, g, jump) in ((:floor, floor, round), (:ceil, ceil, round), (:trunc, trunc, round),
                     (:round, round, v -> floor(v) + 0.5))
    register_function!(f; julia = g,
        eval = _unary(function (a)
            a.kind == NUMBER && isexact(a.num) && return Sym(SymNumber(g(BigInt, _rat(a.num))))
            _isfn(a, f) && return a
            v = _constant_float(a, jump)
            v === nothing ? nothing : Sym(SymNumber(g(BigInt, v)))
        end))
end

# Signatures for extending a function of `n` arguments to expressions: each argument
# a `Sym` or a `Number`, at least one a `Sym`. Used by the package extensions.
_symbolic_signatures(n::Int) =
    [Tuple(T) for T in Iterators.product(ntuple(_ -> (Sym, Number), n)...) if any(==(Sym), T)]

for f in (:sin, :cos, :tan, :asin, :acos, :atan, :sinh, :cosh, :tanh, :asinh, :acosh, :atanh, :abs, :log,
          :sec, :csc, :cot, :asec, :acsc, :acot, :sech, :csch, :coth, :asech, :acsch, :acoth,
          :cbrt, :sign, :floor, :ceil, :trunc, :round)
    @eval Base.$f(x::Sym) = apply_function($(QuoteNode(f)), Sym[x])
end

# functions with an existing symbolic form
Base.log2(x::Sym) = log(Sym(2), x)
Base.log10(x::Sym) = log(Sym(10), x)
Base.log1p(x::Sym) = log(1 + x)
Base.exp2(x::Sym) = Sym(2)^x
Base.exp10(x::Sym) = Sym(10)^x
Base.expm1(x::Sym) = exp(x) - 1
Base.sinpi(x::Sym) = sin(PI * x)
Base.cospi(x::Sym) = cos(PI * x)
Base.sincos(x::Sym) = (sin(x), cos(x))
Base.log(b::Sym, x::Sym) = apply_function(:log, Sym[b, x])
Base.log(b::Number, x::Sym) = log(Sym(b), x)
Base.log(b::Sym, x::Number) = log(b, Sym(x))
