# Series expansions: Taylor, Laurent and Puiseux series, at a point or at infinity, with
# logarithmic terms.
#
# `series(ex, x, x0, n)` substitutes x = x0 + t (x = 1/t at ∞, x = -1/t at -∞) and
# expands in t by truncated power series arithmetic (_PS): exponents are rationals,
# coefficients are expressions free of t, except that log(t) is carried in them as the
# symbol _LOGT, so x^x = 1 + log(x) x + … is a series in x with coefficients in log(x).
# Every operation tracks the order to which its result is correct, so cancellation or
# division only lowers that order; `series` reruns with more terms until the result is
# correct to the order asked for.

const _T = Sym(Symbol("_t"))         # the expansion variable
const _LOGT = Sym(Symbol("_logt"))   # log(_T), kept in coefficients
const _INF = 1 // 0                  # the order of an exact series

# a truncated series Σ c_i t^(e_i) + O(t^o): exponents increasing, coefficients nonzero
struct _PS
    e::Vector{Rational{Int}}
    c::Vector{Sym}
    o::Rational{Int}
end

# raised when an operation needs more terms than it was given (a division by a series
# whose leading term is not known yet): `series` retries with more
struct _SeriesPrecision <: Exception end

_isexactser(a::_PS) = a.o == _INF
_val(a::_PS) = isempty(a.e) ? a.o : a.e[1]
_coefnorm(c::Sym) = c.kind == NUMBER ? c : cancel(expand(c))

# build from terms, combining equal exponents, dropping zeros and terms past the order
function _ps(es, cs, o::Rational{Int})
    d = Dict{Rational{Int},Sym}()
    for (e, c) in zip(es, cs)
        e < o || continue
        d[e] = get(d, e, ZERO) + c
    end
    ks = sort!([k for (k, v) in d if !iszero(_coefnorm(v))])
    _PS(ks, Sym[_coefnorm(d[k]) for k in ks], o)
end
_psconst(c::Sym) = iszero(c) ? _PS(Rational{Int}[], Sym[], _INF) : _PS([0 // 1], [c], _INF)
const _PST = _PS([1 // 1], [ONE], _INF)

_psadd(a::_PS, b::_PS) = _ps(vcat(a.e, b.e), vcat(a.c, b.c), min(a.o, b.o))
_psscale(a::_PS, k::Sym) = iszero(k) ? _psconst(ZERO) : _ps(a.e, [c * k for c in a.c], a.o)
_psshift(a::_PS, s::Rational{Int}) = _PS(a.e .+ s, a.c, a.o + s)          # times t^s

function _psmul(a::_PS, b::_PS)
    o = min(a.o + _val(b), b.o + _val(a))
    es, cs = Rational{Int}[], Sym[]
    for (ea, ca) in zip(a.e, a.c), (eb, cb) in zip(b.e, b.c)
        e = ea + eb
        e < o || continue
        push!(es, e)
        push!(cs, ca * cb)
    end
    _ps(es, cs, o)
end

# the leading term c0 t^v0 and the rest u, with a = c0 t^v0 (1 + u), u of positive order
function _pslead(a::_PS)
    isempty(a.e) && throw(_SeriesPrecision())       # zero so far: need more terms
    v0, c0 = a.e[1], a.c[1]
    u = _PS(a.e[2:end] .- v0, [c / c0 for c in a.c[2:end]], a.o - v0)
    v0, c0, u
end

# Σ_k coef(k) u^k for u of positive valuation, to the order u allows (at most N: an
# exact u gives an infinite sum)
function _psgeom(u::_PS, coef, N::Rational{Int})
    o = min(u.o, N)
    vu = _val(u)
    acc = _psconst(coef(0))
    isempty(u.e) && return _PS(acc.e, acc.c, min(acc.o, o))
    p = _psconst(ONE)
    k = 0
    while true
        k += 1
        k * vu < o || break
        p = _psmul(p, u)
        acc = _psadd(acc, _psscale(p, coef(k)))
    end
    _PS(acc.e, acc.c, min(acc.o, o))
end

_psinv(a::_PS, N::Rational{Int}) = ((v0, c0, u) = _pslead(a);
    _psshift(_psscale(_psgeom(u, k -> isodd(k) ? -ONE : ONE, N + v0), inv(c0)), -v0))

# a^r, r a rational number or a constant expression (then a's valuation must be 0)
function _pspow(a::_PS, r::Sym, N::Rational{Int})
    if r.kind == NUMBER && _isint(r.num) && isexact(r.num)
        k = _toint(r.num)
        k == 0 && return _psconst(ONE)
        b = k > 0 ? a : _psinv(a, N)
        out = b
        for _ in 2:abs(k)
            out = _psmul(out, b)
        end
        return out
    end
    v0, c0, u = _pslead(a)
    s = if iszero(v0)
        0 // 1
    elseif r.kind == NUMBER && isexact(r.num)
        v0 * Rational{Int}(_numerator(r.num), _denominator(r.num))
    else
        throw(ArgumentError("series: a power with a symbolic exponent, $r, of a vanishing or infinite base"))
    end
    binom(k) = (b = ONE; for i in 0:(k - 1); b = b * (r - i) / (i + 1); end; b)
    _psshift(_psscale(_psgeom(u, binom, N - s), c0^r), s)
end

# split into the negative, constant and positive parts
function _pssplit(a::_PS)
    neg = [i for i in eachindex(a.e) if a.e[i] < 0]
    i0 = findfirst(==(0), a.e)
    pos = [i for i in eachindex(a.e) if a.e[i] > 0]
    !isempty(neg), i0 === nothing ? ZERO : a.c[i0], _PS(a.e[pos], a.c[pos], a.o)
end

function _psexp(a::_PS, N::Rational{Int})
    hasneg, a0, h = _pssplit(a)
    hasneg && throw(ArgumentError("series: exp of an expression that is infinite at the point (an essential singularity)"))
    a.o > 0 || throw(_SeriesPrecision())
    # a0 = p + q log(t): exp(a0) = exp(p) t^q, q rational
    q = cancel(diff(a0, _LOGT))
    has(q, _LOGT) && throw(ArgumentError("series: exp of a power of log at the point"))
    p = cancel(a0 - q * _LOGT)
    (q.kind == NUMBER && isexact(q.num)) || iszero(q) ||
        throw(ArgumentError("series: x^$q with a symbolic exponent at the point"))
    s = iszero(q) ? 0 // 1 : Rational{Int}(_numerator(q.num), _denominator(q.num))
    fact(k) = (f = ONE; for i in 1:k; f = f / i; end; f)
    _psshift(_psscale(_psgeom(h, fact, N - s), exp(p)), s)
end

function _pslog(a::_PS, N::Rational{Int})
    v0, c0, u = _pslead(a)
    has(c0, _LOGT) && throw(ArgumentError("series: log of an expression with a logarithmic leading term"))
    head = _psadd(_psconst(log(c0)), _psscale(_psconst(_LOGT), Sym(SymNumber(v0))))
    tail = _psgeom(u, k -> k == 0 ? ZERO : (isodd(k) ? ONE : -ONE) / k, N)
    _psadd(head, tail)
end

# f(a) for a function f with a derivative rule, analytic at a's constant term a0:
# Σ f^(k)(a0)/k! (a - a0)^k
function _psapply(name::Symbol, args::Vector{Sym}, N::Rational{Int})
    i = findall(a -> has(a, _T), args)
    length(i) == 1 || throw(ArgumentError("series: $name with several arguments depending on the variable"))
    i = only(i)
    a = _ser(args[i], N)
    hasneg, a0, h = _pssplit(a)
    hasneg && throw(ArgumentError("series: $(apply_function(name, args)) has no series here (its argument is infinite)"))
    has(a0, _LOGT) && throw(ArgumentError("series: $name of a logarithm at the point"))
    a.o > 0 || throw(_SeriesPrecision())
    y = Sym(Symbol("_y"))
    fy = apply_function(name, Sym[j == i ? y : args[j] for j in eachindex(args)])
    vh = _val(h)
    K = vh == _INF ? 0 : max(0, floor(Int, (min(a.o, N) - 1 // 10^9) / vh))
    coef = Sym[]
    d = fy
    for k in 0:K
        ck = try
            subs(d, y => a0)
        catch err
            _isdivzero(err) || rethrow()
            throw(ArgumentError("series: $(apply_function(name, args)) is singular at the point"))
        end
        _finite_coef(ck) || throw(ArgumentError("series: $(apply_function(name, args)) is singular at the point"))
        _hasfn(ck, :Derivative) && throw(ArgumentError("series: no derivative rule for $name"))
        push!(coef, ck)
        k < K && (d = diff(d, y) / (k + 1))
    end
    _psgeom(h, k -> k <= K ? coef[k+1] : ZERO, N)
end

# a coefficient that is an unevaluated Derivative, or that does not evaluate to a
# finite number when it is a constant, signals that f is not analytic at the point
function _finite_coef(c::Sym)
    _hassymbol(c) && return true
    v = try N(c) catch; return true end
    isfinite(v)
end

# the series of ex (in _T) correct to order N where possible
function _ser(ex::Sym, N::Rational{Int})
    has(ex, _T) || return _psconst(ex)
    _eq(ex, _T) && return _PST
    k = ex.kind
    if k == ADD
        acc = _psconst(Sym(ex.num))
        for i in eachindex(ex.args)
            acc = _psadd(acc, _psscale(_ser(ex.args[i], N), Sym(ex.coefs[i])))
        end
        return acc
    elseif k == MUL
        acc = _psconst(Sym(ex.num))
        for i in eachindex(ex.args)
            acc = _psmul(acc, _serpow(ex.args[i], ex.exps[i], N))
        end
        return acc
    elseif k == POW
        return _serpow(ex.args[1], ex.args[2], N)
    elseif k == FN
        ex.name === :log && length(ex.args) == 1 && return _pslog(_ser(ex.args[1], N), N)
        # through sin and cos, so that poles come out of a division (Laurent terms)
        r = get(_SER_REWRITE, ex.name, nothing)
        r !== nothing && length(ex.args) == 1 && return _ser(r(only(ex.args)), N)
        return _psapply(ex.name, ex.args, N)
    end
    throw(ArgumentError("series: cannot expand $ex"))
end

const _SER_REWRITE = Dict{Symbol,Function}(
    :tan => u -> sin(u) / cos(u), :cot => u -> cos(u) / sin(u),
    :sec => u -> 1 / cos(u), :csc => u -> 1 / sin(u),
    :tanh => u -> sinh(u) / cosh(u), :coth => u -> cosh(u) / sinh(u),
    :sech => u -> 1 / cosh(u), :csch => u -> 1 / sinh(u),
    :cbrt => u -> u^(1 // 3))                       # near the point: the principal root

function _serpow(b::Sym, e::Sym, N::Rational{Int})
    if has(e, _T)                                   # b^e = exp(e log b)
        lb = _eq(b, E) ? _psconst(ONE) : _pslog(_ser(b, N), N)
        return _psexp(_psmul(_ser(e, N), lb), N)
    end
    _eq(b, E) && return _psconst(pow(b, e))
    _pspow(_ser(b, N), e, N)
end

## the user-facing series

"""
    SymSeries

A truncated series in `x` about `x0`: `Σ cᵢ (x - x0)^eᵢ + O((x - x0)^n)`, or in `1/x` at
`x0 = ±Inf`. Returned by [`series`](@ref). `Sym(s)` is the truncated expression (the
`O` term dropped); calling `s` calls that expression. Series about the same point add,
subtract, multiply, divide and take powers, with the order of the error tracked, and
`diff(s, x)` differentiates term by term. Not exported.
"""
struct SymSeries
    x::Sym
    x0::Any                  # a Sym, or ±Inf
    ps::_PS                  # in t (x = x0 + t, x = 1/t or x = -1/t), log(t) as _LOGT
    n::Rational{Int}         # the order asked for (bounds the terms of an exact series)
end

"""
    series(ex, x, x0 = 0, n = 6)

The series of `ex` in `x` about `x0`, with the error `O((x - x0)^n)`: the terms of order
below `n` are exact. At `x0 = Inf` (or `-Inf`) it is a series in `1/x`, with error
`O((1/x)^n)`. Laurent series (negative powers), Puiseux series (fractional powers) and
logarithmic terms (`log(x)` in the coefficients, as for `x^x` at 0) are included.

The result is a [`SymSeries`](@ref); `Sym(s)` gives the truncated expression.

```julia
series(sin(x), x)              # x - x^3/6 + x^5/120 + O(x^6)
series(1/sin(x), x, 0, 3)      # 1/x + x/6 + O(x^3)
series(sqrt(x + x^2), x, 0, 3) # sqrt(x) + x^(3/2)/2 - x^(5/2)/8 + O(x^3)
series(exp(x), x, 1, 3)        # ℯ + ℯ*(x - 1) + ℯ*(x - 1)^2/2 + O((x - 1)^3)
series(sqrt(x^2 + x), x, Inf, 3)   # x + 1/2 - 1/(8x) + 1/(16x^2) + O(1/x^3)
series(x^x, x, 0, 2)           # 1 + log(x)*x + O(x^2)
```

An expression with no series at the point (an essential singularity such as `exp(1/x)`
at 0, or a function without a derivative rule) gives an `ArgumentError`.
"""
function series(ex, x::Sym, x0 = 0, n::Real = 6)
    _isnamed(x) || throw(ArgumentError("series: $x is not a symbol"))
    ex = _sinc_elementary(Sym(ex))
    sub = x0 isa Real && isinf(x0) ? (x0 > 0 ? inv(_T) : -inv(_T)) : Sym(x0) + _T
    ext = subs(ex, x => sub)
    order = Rational{Int}(rationalize(n))
    N = order
    for _ in 1:12
        ps = try
            _ser(ext, N)
        catch err
            err isa _SeriesPrecision || rethrow()
            N += 2
            continue
        end
        if ps.o >= order
            keep = findall(<(order), ps.e)
            # exact (a polynomial, say) if nothing was cut off
            o = ps.o == _INF && length(keep) == length(ps.e) ? _INF : order
            return SymSeries(x, x0 isa Real && isinf(x0) ? x0 : Sym(x0), _PS(ps.e[keep], ps.c[keep], o), order)
        end
        N += max(order - ps.o, 1 // 1) + 1
    end
    throw(ArgumentError("series: could not determine the series of $ex to order $n (too much cancellation)"))
end

# t and log(t) in terms of x
function _tvar(s::SymSeries)
    s.x0 isa Real && isinf(s.x0) && return s.x0 > 0 ? inv(s.x) : -inv(s.x)
    s.x - s.x0
end
# log(t) in terms of x: log(x - x0), or -log(x) at ∞ and -log(-x) at -∞
_logt(s::SymSeries) = s.x0 isa Real && isinf(s.x0) ? -log(s.x0 > 0 ? s.x : -s.x) : log(_tvar(s))
_coef(s::SymSeries, c::Sym) = subs(c, _LOGT => _logt(s))
_term(s::SymSeries, e::Rational{Int}, c::Sym) = _coef(s, c) * _tvar(s)^Sym(SymNumber(e))

Sym(s::SymSeries) = _add_many(_term(s, e, c) for (e, c) in zip(s.ps.e, s.ps.c))
(s::SymSeries)(args...; kws...) = Sym(s)(args...; kws...)

# a term as printed: the coefficient, then the power of (x - x0), not multiplied out
# (multiplying would distribute the coefficient over x - x0)
function _termstring(s::SymSeries, e::Rational{Int}, c::Sym)
    atinf = s.x0 isa Real && isinf(s.x0)
    (iszero(e) || atinf || s.x0 isa Sym && iszero(s.x0)) && return string(_term(s, e, c))
    p = _tvar(s)^Sym(SymNumber(e))
    ps = e == 1 ? "(" * string(_tvar(s)) * ")" : string(p)
    k = _coef(s, c)
    isone(k) && return ps
    _eq(k, -ONE) && return "-" * ps
    ks = string(k)
    (k.kind == ADD) ? "(" * ks * ")*" * ps : ks * "*" * ps
end

function Base.show(io::IO, s::SymSeries)
    parts = String[]
    for (k, (e, c)) in enumerate(zip(s.ps.e, s.ps.c))
        str = _termstring(s, e, c)
        if k == 1
            push!(parts, str)
        elseif startswith(str, "-")
            push!(parts, " - " * str[2:end])
        else
            push!(parts, " + " * str)
        end
    end
    isempty(parts) && push!(parts, "0")
    if s.ps.o != _INF                                # at ±∞, in powers of 1/x either way
        ot = s.x0 isa Real && isinf(s.x0) ? inv(s.x) : _tvar(s)
        push!(parts, " + O(" * (s.ps.o == 1 && ot.kind == ADD ? "(" * string(ot) * ")" : string(ot^Sym(SymNumber(s.ps.o)))) * ")")
    end
    print(io, join(parts))
end
Base.show(io::IO, ::MIME"text/plain", s::SymSeries) = show(io, s)

# arithmetic: series about the same point, or a series and a constant
_same(a::SymSeries, b::SymSeries) = (_eq(a.x, b.x) && (a.x0 isa Sym && b.x0 isa Sym ? _eq(a.x0, b.x0) : isequal(a.x0, b.x0))) ||
    throw(ArgumentError("series about different points or in different variables"))
_lift(s::SymSeries, c) = (c = Sym(c); has(c, s.x) && throw(ArgumentError("add a series, not an expression in $(s.x)")); _psconst(c))
_with(a::SymSeries, ps::_PS, n = a.n) = SymSeries(a.x, a.x0, ps, n)
Base.:+(a::SymSeries, b::SymSeries) = (_same(a, b); _with(a, _psadd(a.ps, b.ps), min(a.n, b.n)))
Base.:*(a::SymSeries, b::SymSeries) = (_same(a, b); _with(a, _psmul(a.ps, b.ps), min(a.n, b.n)))
Base.:-(a::SymSeries) = _with(a, _psscale(a.ps, -ONE))
Base.:-(a::SymSeries, b::SymSeries) = a + (-b)
Base.inv(a::SymSeries) = _with(a, _psinv(a.ps, a.n))
Base.:/(a::SymSeries, b::SymSeries) = a * inv(b)
Base.:^(a::SymSeries, r::Union{Integer,Rational}) = _with(a, _pspow(a.ps, Sym(r), a.n))
for op in (:+, :*, :-, :/)
    @eval Base.$op(a::SymSeries, c::Union{Number,Sym}) = $op(a, _with(a, _lift(a, c)))
    @eval Base.$op(c::Union{Number,Sym}, a::SymSeries) = $op(_with(a, _lift(a, c)), a)
end
Base.:(==)(a::SymSeries, b::SymSeries) = _eq(a.x, b.x) && isequal(a.x0, b.x0) && a.ps.e == b.ps.e &&
    a.ps.o == b.ps.o && all(_eq.(a.ps.c, b.ps.c))

"""
    diff(s::SymSeries, x)

The series of the derivative, term by term; its order is one less (at a finite point).
"""
function Base.diff(s::SymSeries, x::Sym)
    _eq(x, s.x) || return _with(s, _ps(s.ps.e, [diff(c, x) for c in s.ps.c], s.ps.o))
    # d/dt (c t^e) = (e c + dc/dlog t) t^(e-1)
    d = _ps(s.ps.e .- 1, [Sym(SymNumber(e)) * c + diff(c, _LOGT) for (e, c) in zip(s.ps.e, s.ps.c)], s.ps.o - 1)
    s.x0 isa Real && isinf(s.x0) || return _with(s, d, s.n - 1)
    # dt/dx: t = 1/x gives -t^2, t = -1/x gives t^2
    _with(s, _psscale(_psshift(d, 2 // 1), s.x0 > 0 ? -ONE : ONE), s.n + 1)
end
