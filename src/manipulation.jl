## traversal

"""
    free_symbols(ex)

The symbols appearing in `ex`, sorted by name.
"""
function free_symbols(ex::Sym)
    s = Set{Sym}()
    _free!(s, ex)
    sort!(collect(s); by = x -> x.name)
end
free_symbols(::Number) = Sym[]

function _free!(s, ex::Sym)
    ex.kind == SYMBOL && return push!(s, ex)
    foreach(a -> _free!(s, a), ex.args)
    foreach(a -> _free!(s, a), ex.exps)
    s
end

"""
    has(ex, x)

Does `x` occur as a subexpression of `ex`? (Matches whole nodes only.)
"""
has(ex::Sym, x::Sym) = _eq(ex, x) || any(a -> has(a, x), ex.args) || any(a -> has(a, x), ex.exps)

## substitution and calling

"""
    subs(ex, x => v, ...)
    subs(ex, d::AbstractDict)

Replace subexpressions of `ex` (whole nodes only) and re-simplify. Expressions are
callable, so `ex(x => 1, y => 2)` is the same as `subs(ex, x => 1, y => 2)`; see
[`Sym`](@ref) calling for the single-value forms.
"""
subs(ex::Sym, ps::Pair...) = _subs(ex, Dict{Sym,Sym}(Sym(p.first) => Sym(p.second) for p in ps))
subs(ex::Sym, d::AbstractDict) = _subs(ex, Dict{Sym,Sym}(Sym(k) => Sym(v) for (k, v) in d))
subs(ex::Number, args...) = Sym(ex)

function _subs(ex::Sym, d::Dict{Sym,Sym})
    v = get(d, ex, nothing)
    v === nothing || return v
    k = ex.kind
    if k == ADD
        return _add_many(Iterators.flatten(((Sym(ex.num),),
            (_scale(_subs(ex.args[i], d), ex.coefs[i]) for i in eachindex(ex.args)))))
    elseif k == MUL
        r = Sym(ex.num)
        for i in eachindex(ex.args)
            r = mul(r, pow(_subs(ex.args[i], d), _subs(ex.exps[i], d)))
        end
        return r
    elseif k == POW
        return pow(_subs(ex.args[1], d), _subs(ex.args[2], d))
    elseif k == FN
        return apply_function(ex.name, Sym[_subs(a, d) for a in ex.args])
    end
    ex
end

## calling

# One method handles every call form, so dispatch has no ambiguities to resolve:
#   ex()            -> the value of a constant expression as a Julia number
#                      (`Sym(2)() === 2`, `Sym(π)() === π`); an error if ex has a
#                      free symbol. Generic code (e.g. GruntzLimit) reads numeric
#                      leaves this way.
#   ex(x => a, ...) -> substitution, returns a Sym
#   ex(a::Sym)      -> substitution into the only free symbol, returns a Sym
#   ex(a)           -> numeric evaluation at the only free symbol, returns what Julia
#                      arithmetic on `a` returns (Int, Float64, BigFloat, a dual number, ...)
function (ex::Sym)(args...)
    isempty(args) && return _constant_value(ex)
    all(a -> a isa Pair, args) && return subs(ex, args...)
    length(args) == 1 || throw(ArgumentError(
        "call an expression with one value (it has at most one free symbol) or with pairs `x => value`"))
    a = only(args)
    a isa Sym ? _subs1(ex, a) : evaluate(ex, a)
end

function _constant_value(ex::Sym)
    isempty(free_symbols(ex)) ||
        throw(ArgumentError("$ex has free symbols $(join(free_symbols(ex), ", ")); give it a value, e.g. ex(2) or ex(x => 2)"))
    evaluate(ex, 0)      # the argument is unused, and an Int keeps exact values exact
end

function _subs1(ex::Sym, a::Sym)
    vars = free_symbols(ex)
    isempty(vars) && return ex
    length(vars) == 1 || _too_many_symbols(ex)
    subs(ex, only(vars) => a)
end

_too_many_symbols(ex) = throw(ArgumentError(
    "$ex has free symbols $(join(free_symbols(ex), ", ")); call it with pairs, e.g. ex(x => 1, y => 2)"))

## numeric evaluation with Julia arithmetic

"""
    evaluate(ex, a)

Evaluate `ex`, which may have at most one free symbol, at `a` using Julia arithmetic on
`a`. This is what calling an expression on a non-symbolic value does: if `u = x^2` then
`u(3) === 9`, `u(3.0) === 9.0` and `u(1//2) === 1//4`. Because the arithmetic is generic,
`a` may be any number type (big floats, dual numbers, intervals, ...).

- Exact coefficients are used as Julia integers and rationals, so exact inputs give
  exact results where Julia arithmetic would (`(x/3)(3) == 1//1`).
- When `a` is a floating point type (real or complex), constants such as `1/3`,
  `sqrt(2)` and `π` are computed in that type, so `BigFloat` inputs keep their precision.
- Negative integer powers use `inv` (`(1/x)(2) === 0.5`), and a negative base with a
  fractional power gives the principal complex value (`sqrt(x)(-4) == 2im`).
- Functions need a Julia method (registered with `julia = ...`); user functions without
  one throw an `ArgumentError`.
"""
function evaluate(ex::Sym, a)
    var = Ref{Union{Nothing,Sym}}(nothing)
    _ev(ex, a, _liftfloat(typeof(a)), var)
end

_liftfloat(::Type{T}) where {T<:AbstractFloat} = T
_liftfloat(::Type{Complex{T}}) where {T<:AbstractFloat} = T
_liftfloat(::Type) = nothing

function _jlnum(n::SymNumber)
    n.tag == T_FLOAT && return n.f
    n.tag == T_SMALL && return n.d == 1 ? n.n : n.n // n.d
    r = n.big::Rational{BigInt}
    isone(denominator(r)) ? numerator(r) : r
end

_leaf(n::SymNumber, F) = F === nothing ? _jlnum(n) : convert(F, _jlnum(n))
_constant(c, F) = F === nothing ? c : convert(F, c)

function _ev(ex::Sym, a, F, var)
    k = ex.kind
    k == NUMBER && return _leaf(ex.num, F)
    if k == SYMBOL
        if var[] === nothing
            var[] = ex
        elseif !_eq(var[]::Sym, ex)
            throw(ArgumentError("expression has more than one free symbol; call it with pairs, e.g. ex(x => 1, y => 2)"))
        end
        return a
    elseif k == CONSTANT
        ex.name === :pi && return _constant(π, F)
        ex.name === :E && return _constant(ℯ, F)
        throw(ArgumentError("no numeric value for constant $ex"))
    elseif k == ADD
        s = _leaf(ex.num, F)
        for i in eachindex(ex.args)
            t = _ev(ex.args[i], a, F, var)
            s += _isN1(ex.coefs[i]) ? t : _leaf(ex.coefs[i], F) * t
        end
        return s
    elseif k == MUL
        p = _leaf(ex.num, F)
        for i in eachindex(ex.args)
            p *= _evpow(ex.args[i], ex.exps[i], a, F, var)
        end
        return p
    elseif k == POW
        return _evpow(ex.args[1], ex.args[2], a, F, var)
    end
    rule = get(FUNCTION_RULES, ex.name, nothing)
    (rule === nothing || rule.julia === nothing) &&
        throw(ArgumentError("$(ex.name) has no numeric method; register one with `register_function!`"))
    rule.julia(map(x -> _ev(x, a, F, var), ex.args)...)
end

function _evpow(b::Sym, e::Sym, a, F, var)
    b === E && return exp(_ev(e, a, F, var))      # exp(x) is stored as ℯ^x
    bv = _ev(b, a, F, var)
    if e.kind == NUMBER && isexact(e.num)
        if _isint(e.num)
            k = _toint(e.num)
            return k >= 0 ? bv^k : inv(bv)^(-k)
        end
        if e.num == SymNumber(1 // 2)
            return bv isa Real && bv < 0 ? sqrt(complex(bv)) : sqrt(bv)
        end
        ev = convert(F === nothing ? Float64 : F, _jlnum(e.num))
    else
        ev = _ev(e, a, F, var)
    end
    bv isa Real && ev isa Real && bv < 0 && !isinteger(ev) && return complex(bv)^ev
    bv^ev
end

## expansion

"""
    expand(ex)

Distribute products and positive integer powers over sums. Expressions that are
polynomials with numeric coefficients in their symbols are expanded through [`Poly`](@ref).
"""
function expand(ex::Sym)
    if ex.kind == MUL || ex.kind == POW || ex.kind == ADD
        r = _expand_poly(ex)
        r === nothing || return r
    end
    _expand_generic(ex)
end

function _expand_poly(ex::Sym)
    gs = Set{Sym}()
    _numpoly_symbols!(gs, ex) || return nothing
    length(gs) > 16 && return nothing             # keep at least 4 bits per exponent
    g = sort!(collect(gs); lt = _lt)
    try
        return Sym(Poly{SymNumber}(ex, g))
    catch err
        err isa OverflowError || rethrow()
        return nothing                            # degree too large to pack
    end
end

# is ex a polynomial with numeric coefficients in its symbols? (collects the symbols)
function _numpoly_symbols!(gs::Set{Sym}, ex::Sym)
    k = ex.kind
    k == NUMBER && return true
    k == SYMBOL && (push!(gs, ex); return true)
    k == ADD && return all(t -> _numpoly_symbols!(gs, t), ex.args)
    k == MUL && return all(i -> _nonneg_int(ex.exps[i]) && _numpoly_symbols!(gs, ex.args[i]), eachindex(ex.args))
    k == POW && return _nonneg_int(ex.args[2]) && _numpoly_symbols!(gs, ex.args[1])
    false
end

function _expand_generic(ex::Sym)
    k = ex.kind
    if k == ADD
        return _add_many(Iterators.flatten(((Sym(ex.num),),
            (_scale(expand(ex.args[i]), ex.coefs[i]) for i in eachindex(ex.args)))))
    elseif k == MUL
        r = Sym(ex.num)
        for i in eachindex(ex.args)
            r = _expand_mul(r, expand(_mkpow(ex.args[i], ex.exps[i])))
        end
        return r
    elseif k == POW
        b, e = expand(ex.args[1]), expand(ex.args[2])
        if b.kind == ADD && e.kind == NUMBER && _isint(e.num)
            n = _toint(e.num)
            n > 0 && return _expand_pow(b, n)
            return pow(_expand_pow(b, -n), MINUS_ONE)
        end
        return pow(b, e)
    elseif k == FN
        return apply_function(ex.name, Sym[expand(a) for a in ex.args])
    end
    ex
end
expand(x::Number) = Sym(x)

function _expand_mul(a::Sym, b::Sym)
    (a.kind == ADD || b.kind == ADD) || return mul(a, b)
    ta, ca = _termlist(a)
    tb, cb = _termlist(b)
    d = Dict{Sym,SymNumber}()
    sizehint!(d, length(ta) * length(tb))
    c = N0
    for i in eachindex(ta), j in eachindex(tb)
        c = _accum_sym!(d, c, mul(ta[i], tb[j]), ca[i] * cb[j])
    end
    _add_from_dict(c, d)
end

# a sum as parallel vectors of keys and coefficients (constant term keyed by ONE)
function _termlist(x::Sym)
    if x.kind == ADD
        iszero(x.num) && return (x.args, x.coefs)
        return (_insert(x.args, 1, ONE), _insert(x.coefs, 1, x.num))
    end
    x.kind == NUMBER && return (Sym[ONE], SymNumber[x.num])
    c, t = _split_coef(x)
    (Sym[t], SymNumber[c])
end

# (t_1 + ... + t_k)^n by the multinomial theorem. Prefix products are shared down the
# recursion and powers of each term are tabulated, so each monomial costs about one `mul`.
function _expand_pow(b::Sym, n::Int)
    n == 1 && return b
    ts, cs = _termlist(b)
    k = length(ts)
    k == 1 && return _scale(pow(ts[1], Sym(n)), cs[1]^n)
    tpow = [Sym[pow(ts[i], Sym(a)) for a in 1:n] for i in 1:k]
    cpow = [SymNumber[cs[i]^a for a in 1:n] for i in 1:k]
    d = Dict{Sym,SymNumber}()
    sizehint!(d, binomial(n + k - 1, k - 1))
    c = _multinomial!(d, N0, tpow, cpow, 1, n, N1, ONE)
    _add_from_dict(c, d)
end

function _multinomial!(d, c, tpow, cpow, i, rem, coef, mono)
    if i == length(tpow) || rem == 0
        rem > 0 && (coef *= cpow[i][rem]; mono = mul(mono, tpow[i][rem]))
        return _accum_sym!(d, c, mono, coef)
    end
    for a in 0:rem
        ca = a == 0 ? coef : coef * _binom(rem, a) * cpow[i][a]
        ma = a == 0 ? mono : mul(mono, tpow[i][a])
        c = _multinomial!(d, c, tpow, cpow, i + 1, rem - a, ca, ma)
    end
    c
end

_binom(n::Int, k::Int) = n <= 60 ? SymNumber(binomial(n, k)) : SymNumber(binomial(big(n), k))

## numeric evaluation

"""
    N(ex)

Evaluate `ex` numerically as a `Float64` (or `ComplexF64` when the value is not real).
"""
function N(ex::Sym)
    z = _evalc(ex)
    iszero(imag(z)) ? real(z) : z
end
N(x::Number) = x

function _evalc(ex::Sym)::ComplexF64
    k = ex.kind
    k == NUMBER && return Float64(ex.num)
    if k == CONSTANT
        ex.name === :pi && return π
        ex.name === :E && return ℯ
        throw(ArgumentError("no numeric value for constant $ex"))
    end
    k == SYMBOL && throw(ArgumentError("cannot evaluate free symbol $ex numerically"))
    if k == ADD
        s = complex(Float64(ex.num))
        for i in eachindex(ex.args)
            s += Float64(ex.coefs[i]) * _evalc(ex.args[i])
        end
        return s
    elseif k == MUL
        p = complex(Float64(ex.num))
        for i in eachindex(ex.args)
            p *= _cpow(ex.args[i], ex.exps[i])
        end
        return p
    elseif k == POW
        return _cpow(ex.args[1], ex.args[2])
    end
    rule = get(FUNCTION_RULES, ex.name, nothing)
    (rule === nothing || rule.julia === nothing) && throw(ArgumentError("no numeric method for $(ex.name)"))
    rule.julia(map(_evalc, ex.args)...)
end

function _cpow(b::Sym, e::Sym)
    z = _evalc(b)
    e.kind == NUMBER && _isint(e.num) && return z^_toint(e.num)
    w = _evalc(e)
    iszero(imag(z)) && real(z) >= 0 && iszero(imag(w)) && return complex(real(z)^real(w))
    z^w
end

function Base.Float64(ex::Sym)
    v = N(ex)
    v isa Real || throw(InexactError(:Float64, Float64, ex))
    v
end

## convenience constructors

"""
    @syms x y f()

Create symbols (and, with `f()`, a [`SymFunction`](@ref)) and bind them to the names.
"""
macro syms(xs...)
    ex = Expr(:block)
    names = Symbol[]
    for x in xs
        if x isa Symbol
            push!(ex.args, :($(esc(x)) = Sym($(QuoteNode(x)))))
            push!(names, x)
        elseif x isa Expr && x.head === :call && length(x.args) == 1 && x.args[1] isa Symbol
            f = x.args[1]
            push!(ex.args, :($(esc(f)) = SymFunction($(QuoteNode(f)))))
            push!(names, f)
        else
            throw(ArgumentError("@syms: unsupported argument $x"))
        end
    end
    push!(ex.args, Expr(:tuple, map(esc, names)...))
    ex
end

"""
    symbols("x y z")

A tuple of symbols named by a space- or comma-separated string.
"""
symbols(s::AbstractString) = Tuple(Sym(Symbol(n)) for n in split(s, r"[\s,]+"; keepempty = false))
