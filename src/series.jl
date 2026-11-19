# Series expansions: Taylor, Laurent and Puiseux series, at a point or at infinity, with
# logarithmic terms.
#
# `series(ex, x, x0, n)` substitutes x = x0 + t (x = 1/t at ∞, x = -1/t at -∞) and
# expands in t by truncated power series arithmetic (_PS): exponents are rationals or
# rationals plus an irrational real constant (_Ex: x + x^√2), coefficients are expressions
# free of t, except that log(t) is carried in them as the symbol _LOGT, so
# x^x = 1 + log(x) x + … is a series in x with coefficients in log(x).
# Every operation tracks the order to which its result is correct, so cancellation or
# division only lowers that order; `series` reruns with more terms until the result is
# correct to the order asked for.

const _T = Sym(Symbol("_t"))         # the expansion variable
const _LOGT = Sym(Symbol("_logt"))   # log(_T), kept in coefficients
const _TPOS = ScopedValue(false)     # is _T positive (a series at ±∞)

# An exponent: a rational q plus an exact real constant that is not rational (√2, π,
# log 3/log 5) or 0, so x + x^√2 and 1/(1 + x^π) expand. The constants are interned (k
# indexes _IRR, 0 for none), so an exponent is plain bits and the rational exponents, the
# common case, cost no more than Rational{Int}. Exponents with the same constant compare by
# q; others numerically (Float64, then BigFloat when that does not decide).
struct _Ex <: Real
    q::Rational{Int}
    k::Int
end
const _IRR = Sym[]                      # the interned constants, their Float64 values
const _IRRV = Float64[]
const _IRRI = Dict{Sym,Int}()
const _IRR_LOCK = ReentrantLock()
function _intern(s::Sym)::Int
    k = lock(() -> get(_IRRI, s, 0), _IRR_LOCK)
    k == 0 || return k
    v = Float64(Base.inferencebarrier(N)(s))::Float64
    lock(_IRR_LOCK) do
        get!(_IRRI, s) do
            push!(_IRR, s)
            push!(_IRRV, v)
            length(_IRR)
        end
    end
end
_irr(k::Int) = lock(() -> _IRR[k], _IRR_LOCK)
_irrv(k::Int) = lock(() -> _IRRV[k], _IRR_LOCK)
# q + s for a constant s (its rational part moved to q)
function _Ex(q::Rational{Int}, s::Sym)
    (iszero(s) || isinf(q)) && return _Ex(q, 0)
    s.kind == NUMBER && return _Ex(q + Rational{Int}(_numerator(s.num), _denominator(s.num)), 0)
    if s.kind == ADD && !iszero(s.num)
        q += Rational{Int}(_numerator(s.num), _denominator(s.num))
        s = s - Sym(s.num)
    end
    _Ex(q, _intern(s))
end
_Ex(q::Rational{Int}) = _Ex(q, 0)
_Ex(q::Integer) = _Ex(Rational{Int}(q), 0)
_Ex(q::Rational) = _Ex(Rational{Int}(q), 0)
_Ex(e::_Ex) = e
Base.convert(::Type{_Ex}, x::Union{Integer,Rational}) = _Ex(x)
Base.promote_rule(::Type{_Ex}, ::Type{<:Union{Integer,Rational}}) = _Ex
const _INF = _Ex(1 // 0)             # the order of an exact series

_isrational(e::_Ex) = e.k == 0
_irrpart(e::_Ex) = e.k == 0 ? ZERO : _irr(e.k)
# the exponent a constant Sym is (an exact rational or an exact real constant), or nothing
function _ex_of(r::Sym)
    r.kind == NUMBER && return isexact(r.num) ? _Ex(Rational{Int}(_numerator(r.num), _denominator(r.num))) : nothing
    (_hassymbol(r) || has(r, IM) || _hasfn(r, :Integral)) && return nothing
    v = try N(r) catch; return nothing end
    (v isa Real && isfinite(v)) || return nothing
    _Ex(0 // 1, r)
end
_sym(e::_Ex) = e.k == 0 ? Sym(SymNumber(e.q)) : Sym(SymNumber(e.q)) + _irr(e.k)
Base.Float64(e::_Ex) = e.k == 0 ? Float64(e.q) : Float64(e.q) + _irrv(e.k)
Base.:+(a::_Ex, b::_Ex) = a.k == 0 && b.k == 0 ? _Ex(a.q + b.q, 0) : _Ex(a.q + b.q, _irrpart(a) + _irrpart(b))
Base.:-(a::_Ex) = a.k == 0 ? _Ex(-a.q, 0) : _Ex(-a.q, -_irr(a.k))
Base.:-(a::_Ex, b::_Ex) = a + (-b)
Base.:*(a::_Ex, b::_Ex) = a.k == 0 && b.k == 0 ? _Ex(a.q * b.q, 0) : _Ex(0 // 1, _sym(a) * _sym(b))
Base.:(==)(a::_Ex, b::_Ex) = a.q == b.q && a.k == b.k
Base.isequal(a::_Ex, b::_Ex) = a == b
Base.hash(e::_Ex, h::UInt) = e.k == 0 ? hash(e.q, h) : hash(e.k, hash(e.q, h))
@inline Base.isless(a::_Ex, b::_Ex) = a.k == b.k ? a.q < b.q : _isless_irr(a, b)   # (∞ has k = 0)
# (the slow path behind inference barriers: inferring through _evalf_big and N made isless
# infer as Any, and every sort of exponents dynamic)
@noinline function _isless_irr(a::_Ex, b::_Ex)::Bool
    (isinf(a.q) || isinf(b.q)) && return a.q < b.q
    d = Float64(a) - Float64(b)
    abs(d) > 1e-9 * (1 + abs(Float64(a)) + abs(Float64(b))) && return d < 0
    w = Base.inferencebarrier(_evalf_big)(_sym(a) - _sym(b), 256)
    w isa Real || throw(ArgumentError("series: the exponents $(_sym(a)) and $(_sym(b)) cannot be ordered"))
    abs(w) > big(2.0)^-200 && return w < 0
    throw(ArgumentError("series: whether the exponents $(_sym(a)) and $(_sym(b)) are equal is not known"))
end
@inline Base.:<(a::_Ex, b::_Ex) = isless(a, b)
Base.:<=(a::_Ex, b::_Ex) = a == b || isless(a, b)
Base.isinteger(e::_Ex) = e.k == 0 && isinteger(e.q)
Base.iszero(e::_Ex) = e.k == 0 && iszero(e.q)
Base.isinf(e::_Ex) = isinf(e.q)
Base.zero(::Type{_Ex}) = _Ex(0 // 1, 0)
Base.one(::Type{_Ex}) = _Ex(1 // 1, 0)
# a rational exponent as a Rational{Int} (code that reads a series' rational exponents)
Base.Rational{Int}(e::_Ex) = e.k == 0 ? e.q : throw(InexactError(:Rational, Rational{Int}, _sym(e)))
Base.Int(e::_Ex) = Int(Rational{Int}(e))
Base.show(io::IO, e::_Ex) = print(io, e.k == 0 ? e.q : _sym(e))
# not a whole number of t: t^e takes the side t > 0 (a fractional or irrational power)
_isfrac(e::_Ex) = e.k != 0 || denominator(e.q) > 1

# a truncated series Σ c_i t^(e_i) + O(t^o): exponents increasing, coefficients nonzero
struct _PS
    e::Vector{_Ex}
    c::Vector{Sym}
    o::_Ex
end

# raised when an operation needs more terms than it was given (a division by a series
# whose leading term is not known yet): `series` retries with more
struct _SeriesPrecision <: Exception end

_isexactser(a::_PS) = a.o == _INF
_val(a::_PS) = isempty(a.e) ? a.o : a.e[1]
_coefnorm(c::Sym) = c.kind == NUMBER ? c : cancel(expand(c))

# build from terms, combining equal exponents, dropping zeros and terms past the order
function _ps(es::AbstractVector, cs::AbstractVector, o::_Ex)
    idx = sort!([i for i in eachindex(es) if es[i] < o]; by = i -> es[i])   # (stable)
    ks, vs = _Ex[], Sym[]
    for i in idx
        if !isempty(ks) && ks[end] == es[i]
            vs[end] += cs[i]
        else
            push!(ks, es[i]); push!(vs, cs[i])
        end
    end
    # each coefficient normalized once (it was twice: to drop the zeros, then to keep)
    keep = Bool[]
    for j in eachindex(vs)
        vs[j] = _coefnorm(vs[j])
        push!(keep, !iszero(vs[j]))
    end
    _PS(ks[keep], vs[keep], o)
end
_psconst(c::Sym) = iszero(c) ? _PS(_Ex[], Sym[], _INF) : _PS([0 // 1], [c], _INF)
const _PST = _PS([1 // 1], [ONE], _INF)

_psadd(a::_PS, b::_PS) = _ps(vcat(a.e, b.e), vcat(a.c, b.c), min(a.o, b.o))
_psscale(a::_PS, k::Sym) = iszero(k) ? _psconst(ZERO) : _ps(a.e, [c * k for c in a.c], a.o)
_psshift(a::_PS, s::Real) = _PS(a.e .+ s, a.c, a.o + s)           # times t^s

function _psmul(a::_PS, b::_PS)
    o = min(a.o + _val(b), b.o + _val(a))
    es, cs = _Ex[], Sym[]
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
function _psgeom(u::_PS, coef, N::_Ex)
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

_psinv(a::_PS, N::_Ex) = ((v0, c0, u) = _pslead(a);
    _psshift(_psscale(_psgeom(u, k -> isodd(k) ? -ONE : ONE, N + v0), inv(c0)), -v0))

# a with the terms of exponent o or more dropped (then known to order o); a itself when
# there are none, so an exact series stays exact
function _pstrunc(a::_PS, o::_Ex)
    all(<(o), a.e) && return a
    keep = findall(<(o), a.e)
    _PS(a.e[keep], a.c[keep], min(a.o, o))
end

# a^r, r a rational number or a constant expression (then a's valuation must be 0)
function _pspow(a::_PS, r::Sym, N::_Ex)
    if r.kind == NUMBER && _isint(r.num) && isexact(r.num)
        k = _toint(r.num)
        k == 0 && return _psconst(ONE)
        b = k > 0 ? a : _psinv(a, N)
        # terms past what the remaining factors can bring below N dropped as they appear
        # (an exact (1 + x²)^10000 was expanded in full; SymPy #6753)
        vb = _val(b)
        out = _pstrunc(b, N - (abs(k) - 1) * vb)
        for j in 2:abs(k)
            out = _pstrunc(_psmul(out, b), N - (abs(k) - j) * vb)
        end
        return out
    end
    v0, c0, u = _pslead(a)
    s = if iszero(v0)
        _Ex(0 // 1)
    elseif (er = _ex_of(r)) !== nothing
        # (an irrational constant exponent, x^√2: SymPy #21358)
        v0 * er
    elseif r.kind == NUMBER && isfloat(r.num) && isfinite(Float64(r.num)) &&
           denominator(rationalize(Int, Float64(r.num))) <= 1000
        # a float exponent: the powers of x by the rational it is (1.2 as 6/5), the
        # coefficients in floats ((x + sin x)^1.2; SymPy #6726)
        v0 * rationalize(Int, Float64(r.num))
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

function _psexp(a::_PS, N::_Ex)
    hasneg, a0, h = _pssplit(a)
    hasneg && throw(ArgumentError("series: exp of an expression that is infinite at the point (an essential singularity)"))
    a.o > 0 || throw(_SeriesPrecision())
    # a0 = p + q log(t): exp(a0) = exp(p) t^q, q rational
    q = cancel(diff(a0, _LOGT))
    has(q, _LOGT) && throw(ArgumentError("series: exp of a power of log at the point"))
    p = cancel(a0 - q * _LOGT)
    s = iszero(q) ? _Ex(0 // 1) : _ex_of(q)
    s === nothing && throw(ArgumentError("series: x^$q with a symbolic exponent at the point"))
    fact(k) = (f = ONE; for i in 1:k; f = f / i; end; f)
    _psshift(_psscale(_psgeom(h, fact, N - s), exp(p)), s)
end

function _pslog(a::_PS, N::_Ex)
    v0, c0, u = _pslead(a)
    has(c0, _LOGT) && throw(ArgumentError("series: log of an expression with a logarithmic leading term"))
    head = _psadd(_psconst(log(c0)), _psscale(_psconst(_LOGT), _sym(v0)))
    tail = _psgeom(u, k -> k == 0 ? ZERO : (isodd(k) ? ONE : -ONE) / k, N)
    _psadd(head, tail)
end

# f(a) for a function f with a derivative rule, analytic at a's constant term a0:
# Σ f^(k)(a0)/k! (a - a0)^k
function _psapply(name::Symbol, args::Vector{Sym}, N::_Ex)
    i = findall(a -> has(a, _T), args)
    # (a derivative of an unknown function: no series is known; "several arguments" pointed
    # away from that; SymPy #11313)
    if name === :Derivative && !isempty(args) && args[1].kind == FN && !haskey(FUNCTION_RULES, args[1].name)
        throw(ArgumentError("series: $(args[1].name) is an unknown function; no series of its derivative is known"))
    end
    # a piecewise function: the piece that holds just past the point (t > 0, the series'
    # side), between it and the first breakpoint after it (SymPy #4315)
    if name === :piecewise
        fx = apply_function(name, args)
        bps = _breakpoints(fx, _T)
        (bps === nothing || any(_hassymbol, bps)) &&
            throw(ArgumentError("series: $fx: where its pieces change could not be decided"))
        pos = [real(_evalc(b)) for b in bps]          # (N here is the order)
        pos = filter(>(0), pos)
        p = isempty(pos) ? ONE : Sym(SymNumber(rationalize(BigInt, minimum(pos) / 2; tol = minimum(pos) / 8)))
        v = _pieces_at(fx, _T, p)
        (v === nothing || _haspieces(v, _T)) && throw(ArgumentError("series: $fx: no piece holds past the point"))
        return _ser(v, N)
    end
    length(i) == 1 || throw(ArgumentError("series: $name with several arguments depending on the variable"))
    i = only(i)
    a = _ser(args[i], N)
    hasneg, a0, h = _pssplit(a)
    if hasneg && name in (:atan, :acot) && length(args) == 1
        # an infinite argument: atan(a) = ±π/2 - atan(1/a), the sign of a's leading term
        # (SymPy #27109); acot(a) = atan(1/a)
        v0, c0, _ = _pslead(a)
        w = _psinv(a, N)
        at = _psgeom(w, k -> isodd(k) ? Sym((-1)^((k - 1) ÷ 2) // k) : ZERO, N)
        name === :acot && return at
        # the sign of c0 t^v0: t > 0 at ±∞ (x = ±1/t); at a point, only for an even or
        # fractional power (atan(1/x) at 0 is ±π/2 - x by side)
        (_TPOS[] || _isfrac(v0) || iseven(numerator(v0.q))) ||
            throw(ArgumentError("series: $(apply_function(name, args)) has different series on the two sides of the point"))
        sg = has(c0, _LOGT) ? nothing : _csign(c0)
        sg === nothing && throw(ArgumentError("series: atan at the point: the sign of $c0 is not known"))
        return _psadd(_psconst(sg * PI / 2), _psscale(at, MINUS_ONE))
    end
    if hasneg && name === :loggamma && length(args) == 1         # Stirling's series
        _, c0, _ = _pslead(a)
        (_TPOS[] && !has(c0, _LOGT) && _csign(c0) == 1) ||
            throw(ArgumentError("series: $(apply_function(name, args)): only for an argument tending to +∞"))
        return _psloggamma(a, N)
    end
    hasneg && throw(ArgumentError("series: $(apply_function(name, args)) has no series here (its argument is infinite)"))
    has(a0, _LOGT) && throw(ArgumentError("series: $name of a logarithm at the point"))
    # Bessel functions at 0: J_ν and I_ν of non-integer order are x^ν times a series in
    # x², not a Taylor series (SymPy #25678); Y_ν and K_ν from them; Y_n, K_n have logs
    if name in (:besselj, :besseli, :bessely, :besselk) && length(args) == 2 && i == 2 && iszero(a0)
        ν = args[1]
        if ν.kind == NUMBER && isexact(ν.num) && !_isint(ν.num)
            return _psbessel(name, ν, a, N)
        elseif name in (:bessely, :besselk)
            throw(ArgumentError("series: $(apply_function(name, args)) at 0 (logarithmic terms) is not supported"))
        end
    end
    # polylogarithms: Σ zᵏ/kˢ at 0 (the derivative rule divides by z), and at their branch
    # point z = 1 not a Taylor series (SymPy #25159)
    if (iszero(a0) || _eq(a0, ONE)) && (name === :li && length(args) == 2 && i == 2 ||
                        length(args) == 1 && name in (:li0, :li1, :li2, :li3, :li4, :li5, :li6))
        o = name === :li ? args[1] : Sym(parse(Int, string(name)[3:end]))
        iszero(a0) && return _psgeom(a, k -> k == 0 ? ZERO : inv(Sym(k)^o), N)
        (o.kind == NUMBER && isexact(o.num)) ||
            throw(ArgumentError("series: $(apply_function(name, args)) at its branch point 1 with a symbolic order"))
        return _pspolylog(o, a, N)
    end
    # Lambert W at 0: Σ (-n)^(n-1)/n! zⁿ (its derivative rule is 0/0 there; SymPy #7259)
    if name === :lambertw && length(args) == 1 && iszero(a0)
        return _psgeom(a, k -> k == 0 ? ZERO : Sym(-k)^(k - 1) / Sym(factorial(big(k))), N)
    end
    # E₁ and Eᵥ at 0 (SymPy #6057, #6725): E₁(u) = -γ - log u + Σ_{k≥1} (-1)^(k+1) u^k/(k k!);
    # Eᵥ(u) = Γ(1 - ν) u^(ν-1) - Σ_{k≥0} (-u)^k/(k! (k + 1 - ν)) for ν not an integer; for
    # an integer n ≥ 1 the k = n - 1 term is (-u)^(n-1)/(n-1)! (ψ(n) - log u)
    if name === :expint && iszero(a0) && (length(args) == 1 || (i == 2 && args[1].kind == NUMBER && isexact(args[1].num)))
        ν = length(args) == 1 ? ONE : args[1]
        u = a
        fact(k) = Sym(factorial(big(k)))
        if _isint(ν.num) && ν.num >= N1
            n = _toint(ν.num)
            logpart = _psadd(_psconst(apply_function(:digamma, Sym[Sym(n)])), _psscale(_pslog(u, N), MINUS_ONE))
            head = _psmul(_psscale(_pspow(u, Sym(n - 1), N), Sym(-1)^(n - 1) / fact(n - 1)), logpart)
            tail = _psgeom(u, k -> k == n - 1 ? ZERO : -Sym(-1)^k / (fact(k) * (k + 1 - n)), N)
            return _psadd(head, tail)
        elseif !_isint(ν.num)
            head = _psscale(_pspow(u, ν - 1, N), apply_function(:gamma, Sym[1 - ν]))
            tail = _psgeom(u, k -> -Sym(-1)^k / (fact(k) * (k + 1 - ν)), N)
            return _psadd(head, tail)
        end
    end
    # floor and ceil: constant near the point, floor(a₀) where a₀ is not an integer; at an
    # integer n, by the way the argument moves for x - x0 > 0: floor is n when it
    # increases, n - 1 when it decreases (and ceil n + 1, n; SymPy #4156)
    if name in (:floor, :ceil) && length(args) == 1
        if !(a0.kind == NUMBER && isexact(a0.num) && _isint(a0.num))
            v = apply_function(name, Sym[a0])
            _hasfn(v, name) && throw(ArgumentError("series: $(apply_function(name, args)): whether $a0 is an integer is not known"))
            return _psconst(v)
        end
        _, c1, _ = _pslead(h)
        sg = has(c1, _LOGT) ? nothing : _csign(c1)
        sg === nothing && throw(ArgumentError("series: $(apply_function(name, args)) at an integer: the sign of $c1 is not known"))
        return _psconst(name === :floor ? (sg > 0 ? a0 : a0 - 1) : (sg > 0 ? a0 + 1 : a0))
    end
    # inverse functions at their branch points ±1: through square roots, whose series
    # have half-integer powers (acosh(1 + z) = √(2z)(1 - z/12 + …); SymPy #25175)
    if length(args) == 1 && (_eq(a0, ONE) || _eq(a0, MINUS_ONE))
        w = args[1]
        # an argument passing 1 upwards along the real line (1 + x): acos = -i·acosh,
        # the value above the cut (N's, as Julia's acos(2 + 0im)), and real for x < 0; the
        # half-angle form takes the value below the cut there (SymPy #4004)
        up = false
        if _eq(a0, ONE) && name in (:acos, :asin)
            _, c1, _ = _pslead(h)
            up = !has(c1, _LOGT) && _csign(c1) == 1
        end
        r = if name === :acos
            up ? -IM * 2asinh(sqrt((w - 1) / 2)) :
            _eq(a0, ONE) ? 2asin(sqrt((1 - w) / 2)) : PI - 2asin(sqrt((1 + w) / 2))
        elseif name === :asin
            up ? PI / 2 + IM * 2asinh(sqrt((w - 1) / 2)) :
            _eq(a0, ONE) ? PI / 2 - 2asin(sqrt((1 - w) / 2)) : -PI / 2 + 2asin(sqrt((1 + w) / 2))
        elseif name === :acosh && _eq(a0, ONE)
            2asinh(sqrt((w - 1) / 2))
        elseif name === :atanh
            (log(1 + w) - log(1 - w)) / 2
        elseif name === :acoth
            (log(w + 1) - log(w - 1)) / 2
        else
            nothing
        end
        r === nothing || return _ser(r, N)
    end
    # inverse functions at a point inside their branch cut on the real line: in a form
    # analytic along the line on the side whose values N gives (Julia's, at u + 0im), so
    # the derivatives agree with the value (the derivative of acos at 2 through √(1 - 4)
    # took the other side: acos(2 + x) had +i x/√3 for -i x/√3)
    if length(args) == 1 && name in (:acos, :asin, :asec, :acsc, :asech) && !_hassymbol(a0)
        v = try _evalc(a0) catch; nothing end
        v isa Complex && iszero(imag(v)) && (v = real(v))
        w = args[1]
        r = if !(v isa Real)
            nothing
        elseif name === :acos && v > 1
            -IM * acosh(w)
        elseif name === :asin && v > 1
            PI / 2 + IM * acosh(w)
        elseif name === :asec && -1 < v < 0
            PI + IM * acosh(-1 / w)
        elseif name === :acsc && -1 < v < 0
            -PI / 2 - IM * acosh(-1 / w)
        elseif name === :asech && abs(v) > 1
            -IM * acos(1 / w)
        else
            nothing
        end
        r === nothing || return _ser(r, N)
    end
    # Γ, ψ and ψ⁽ᵏ⁾ at a pole -m (m = 0, 1, …): moved to the regular point 1 by the
    # recurrences, Γ(w) = Γ(w + m + 1)/(w (w + 1) ⋯ (w + m)) and ψ⁽ᵏ⁾(w) = ψ⁽ᵏ⁾(w + m + 1)
    # - (-1)^k k! Σ_j 1/(w + j)^(k+1), a Laurent series (SymPy #23635)
    if a0.kind == NUMBER && isexact(a0.num) && _isint(a0.num) && a0.num <= N0 && (length(args) == 1 && name in (:gamma, :digamma, :trigamma) ||
                           name === :polygamma && i == 2 && args[1].kind == NUMBER && isexact(args[1].num) && _isint(args[1].num) && args[1].num >= N0)
        m = -_toint(a0.num)
        w = args[i]
        r = if name === :gamma
            apply_function(:gamma, Sym[w + m + 1]) / prod(Sym[w + j for j in 0:m])
        else
            k = name === :digamma ? 0 : name === :trigamma ? 1 : _toint(args[1].num)
            f = k == 0 ? apply_function(:digamma, Sym[w + m + 1]) : apply_function(:polygamma, Sym[Sym(k), w + m + 1])
            f - Sym(isodd(k) ? -factorial(big(k)) : factorial(big(k))) * _add_many(inv((w + j)^(k + 1)) for j in 0:m)
        end
        return _ser(r, N)
    end
    a.o > 0 || throw(_SeriesPrecision())
    y = Sym(Symbol("_y"))
    fy = apply_function(name, Sym[j == i ? y : args[j] for j in eachindex(args)])
    vh = _val(h)
    mo = min(a.o, N)
    K = vh == _INF ? 0 : _isrational(mo) && _isrational(vh) ? max(0, floor(Int, (mo.q - 1 // 10^9) / vh.q)) :
                         max(0, floor(Int, (Float64(mo) - 1e-9) / Float64(vh)))
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
        # (a built-in function without a derivative rule; an undefined function's Taylor
        # series is in its derivatives at the point, f(0) + f'(0) x + …; SymPy #3978)
        _hasfn(ck, :Derivative) && (haskey(FUNCTION_RULES, name) || length(args) > 1) &&
            throw(ArgumentError("series: no derivative rule for $name"))
        push!(coef, ck)
        k < K && (d = diff(d, y) / (k + 1))
    end
    _psgeom(h, k -> k <= K ? coef[k+1] : ZERO, N)
end

# Li_s(a) for a → 1, with μ = -log a → 0: Γ(1 - s) μ^(s-1) + Σ ζ(s - k) (-μ)^k/k! for s
# not a positive integer; for s = n ≥ 1 the terms k = n - 1 are instead
# (-μ)^(n-1)/(n-1)! (H_(n-1) - log μ)
function _pspolylog(s::Sym, a::_PS, N::_Ex)
    μ = _psscale(_pslog(a, N), MINUS_ONE)
    _val(μ) > 0 || throw(_SeriesPrecision())
    m = _psscale(μ, MINUS_ONE)
    fact(k) = Sym(factorial(big(k)))
    if _isint(s.num) && s.num > N0
        n = _toint(s.num)
        H = _add_many(Sym(1 // j) for j in 1:(n - 1))
        reg = _psgeom(m, k -> (k == n - 1 ? H : apply_function(:zeta, Sym[s - k])) / fact(k), N)
        lg = _psscale(_psmul(_pspow(m, Sym(n - 1), N), _pslog(μ, N)), -inv(fact(n - 1)))
        return _psadd(reg, lg)
    end
    sing = _psscale(_pspow(μ, s - 1, N), apply_function(:gamma, Sym[1 - s]))
    _psadd(sing, _psgeom(m, k -> apply_function(:zeta, Sym[s - k]) / fact(k), N))
end

# J_ν(u) = Σ (-1)^k (u/2)^(2k+ν)/(k! Γ(k+ν+1)), I_ν alike without the signs, for u of
# positive valuation and ν a non-integer rational; Y_ν = (J_ν cos νπ - J_{-ν})/sin νπ,
# K_ν = π/2 (I_{-ν} - I_ν)/sin νπ
function _psbessel(name::Symbol, ν::Sym, u::_PS, N::_Ex)
    function JI(ν, alt)
        s = _val(u) * _ex_of(ν)                                                   # of u^ν
        g = apply_function(:gamma, Sym[ν + 1])
        c(k) = (alt && isodd(k) ? -ONE : ONE) /
               (Sym(factorial(big(k))) * g * prod((ν + j for j in 1:k); init = ONE))
        w = _psscale(_psmul(u, u), Sym(1 // 4))
        _psmul(_pspow(_psscale(u, Sym(1 // 2)), ν, N), _psgeom(w, c, N - s))
    end
    alt = name in (:besselj, :bessely)
    name in (:besselj, :besseli) && return JI(ν, alt)
    sn = sin(ν * PI)
    name === :bessely && return _psscale(_psadd(_psscale(JI(ν, true), cos(ν * PI)), _psscale(JI(-ν, true), MINUS_ONE)), inv(sn))
    _psscale(_psadd(JI(-ν, false), _psscale(JI(ν, false), MINUS_ONE)), PI / (2sn))
end

# a coefficient that is an unevaluated Derivative, or that does not evaluate to a
# finite number when it is a constant, signals that f is not analytic at the point
function _finite_coef(c::Sym)
    _hassymbol(c) && return true
    v = try N(c) catch; return true end
    isfinite(v)
end

# the series of ex (in _T) correct to order N where possible
function _ser(ex::Sym, N::_Ex)
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
        # asech and acsch where their argument tends to 0: log singularities, through
        # log((1 + √(1 ∓ u²))/u) (series(asech(x), x, 0, n) raised "singular"; SymPy #10885)
        if ex.name in (:asech, :acsch) && length(ex.args) == 1
            u = only(ex.args)
            us = _ser(u, N)
            if iszero(_pssplit(us)[2]) && !_pssplit(us)[1]
                sg = ex.name === :asech ? -1 : 1
                return _ser(log(1 + sqrt(1 + sg * u^2)) - log(u), N)
            end
        end
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

function _serpow(b::Sym, e::Sym, N::_Ex)
    if has(e, _T)                                   # b^e = exp(e log b)
        lb = _eq(b, E) ? _psconst(ONE) : _pslog(_ser(b, N), N)
        return _psexp(_psmul(_ser(e, N), lb), N)
    end
    _eq(b, E) && return _psconst(pow(b, e))
    _pspow(_ser(b, N), e, N)
end

## asymptotic expansions at ±∞: terms with a factor e^G, e^G cos H or e^G sin H
#
# At ±∞ an expression may have no power series in t: erf(x) = 1 - e^(-x²)(1/(√π x) - …)
# is 1 plus an exponentially small part, and sin(x) oscillates. Such an expansion is a sum
# of pieces e^G trig(H) p: G and H exact sums of negative powers of t (the parts of an
# exponent or of an argument that are infinite at the point), trig one of 1, cos and sin,
# and p a series. The piece with no factor is the ordinary series. Pieces multiply by
# adding the G's and by the product formulas for cos and sin, so exp(x²) erfc(x) is a
# power series again and cos(x)² + sin(x)² is 1. (SymPy #27109, #26257, #26207, #26210,
# #26096)

struct _Piece
    G::_PS          # exact, negative exponents; empty: no exponential factor
    k::Symbol       # :one, :cos or :sin
    H::_PS          # exact, negative exponents; empty for :one
    p::_PS
end

struct _AS
    plain::_PS
    pcs::Vector{_Piece}
end

const _PS0 = _PS(_Ex[], Sym[], _INF)
_asplain(p::_PS) = _AS(p, _Piece[])
_asconst(c) = _asplain(_psconst(Sym(c)))
_isexact0(p::_PS) = isempty(p.e) && p.o == _INF
_pseq(a::_PS, b::_PS) = a.e == b.e && a.o == b.o && all(_eq.(a.c, b.c))
_allpieces(a::_AS) = vcat([(_PS0, :one, _PS0, a.plain)], [(q.G, q.k, q.H, q.p) for q in a.pcs])

# collect pieces, normalized: cos 0 = 1, sin 0 = 0, cos(-H) = cos H, sin(-H) = -sin H
# (H's leading coefficient made positive when its sign is known), equal factors merged
function _asbuild(f)
    plain = Ref(_PS0)
    pcs = _Piece[]
    function add(G::_PS, k::Symbol, H::_PS, p::_PS)
        _isexact0(p) && return
        if k !== :one && isempty(H.e)
            k === :sin && return
            k = :one
        end
        if k !== :one && _csign(H.c[1]) == -1
            H = _psscale(H, MINUS_ONE)
            k === :sin && (p = _psscale(p, MINUS_ONE))
        end
        k === :one && (H = _PS0)
        if k === :one && isempty(G.e)
            plain[] = _psadd(plain[], p)
            return
        end
        for (i, q) in pairs(pcs)
            if q.k === k && _pseq(q.G, G) && _pseq(q.H, H)
                pcs[i] = _Piece(G, k, H, _psadd(q.p, p))
                return
            end
        end
        push!(pcs, _Piece(G, k, H, p))
    end
    f(add)
    _AS(plain[], filter(q -> !_isexact0(q.p), pcs))
end

_asadd(a::_AS, b::_AS) = _asbuild(add -> foreach(t -> add(t...), vcat(_allpieces(a), _allpieces(b))))
_asscale(a::_AS, c::Sym) = _asbuild(add -> foreach(((G, k, H, p),) -> add(G, k, H, _psscale(p, c)), _allpieces(a)))

function _asmul(a::_AS, b::_AS)
    _asbuild() do add
        for (G1, k1, H1, p1) in _allpieces(a), (G2, k2, H2, p2) in _allpieces(b)
            (_isexact0(p1) || _isexact0(p2)) && continue
            p, G = _psmul(p1, p2), _psadd(G1, G2)
            if k1 === :one
                add(G, k2, H2, p)
            elseif k2 === :one
                add(G, k1, H1, p)
            else
                h, mh = _psscale(p, Sym(1 // 2)), _psscale(p, Sym(-1 // 2))
                D, S = _psadd(H1, _psscale(H2, MINUS_ONE)), _psadd(H1, H2)
                if k1 === :cos && k2 === :cos          # cos A cos B = (cos(A-B) + cos(A+B))/2
                    add(G, :cos, D, h); add(G, :cos, S, h)
                elseif k1 === :sin && k2 === :sin      # sin A sin B = (cos(A-B) - cos(A+B))/2
                    add(G, :cos, D, h); add(G, :cos, S, mh)
                elseif k1 === :sin                     # sin A cos B = (sin(A+B) + sin(A-B))/2
                    add(G, :sin, S, h); add(G, :sin, D, h)
                else                                   # cos A sin B = (sin(A+B) - sin(A-B))/2
                    add(G, :sin, S, h); add(G, :sin, D, mh)
                end
            end
        end
    end
end

function _asinv(a::_AS, N::_Ex)
    isempty(a.pcs) && return _asplain(_psinv(a.plain, N))
    (length(a.pcs) == 1 && _isexact0(a.plain) && a.pcs[1].k === :one) ||
        throw(ArgumentError("series: the reciprocal of a sum with exponential or oscillating factors is not supported"))
    q = a.pcs[1]
    _AS(_PS0, [_Piece(_psscale(q.G, MINUS_ONE), :one, _PS0, _psinv(q.p, N))])
end

function _aspowk(a::_AS, k::Integer, N::_Ex)
    k == 0 && return _asconst(ONE)
    b = k > 0 ? a : _asinv(a, N)
    out = b
    for _ in 2:abs(k)
        out = _asmul(out, b)
    end
    out
end

# e^L: the negative powers of L stay in the factor e^G, the rest is a series
function _asexp(L::_PS, N::_Ex)
    neg = findall(<(0), L.e)
    isempty(neg) && return _asplain(_psexp(L, N))
    L.o > 0 || throw(_SeriesPrecision())
    rest = setdiff(eachindex(L.e), neg)
    G = _PS(L.e[neg], L.c[neg], _INF)
    _AS(_PS0, [_Piece(G, :one, _PS0, _psexp(_PS(L.e[rest], L.c[rest], L.o), N))])
end

# cos(u), sin(u) for u = H + r0 + h: cos(H + R) = cos H cos R - sin H sin R and so on
function _astrig(name::Symbol, u::_PS, N::_Ex)
    u.o > 0 || throw(_SeriesPrecision())
    neg = findall(<(0), u.e)
    H = _PS(u.e[neg], u.c[neg], _INF)
    _, r0, h = _pssplit(u)
    has(r0, _LOGT) && throw(ArgumentError("series: $name of a logarithm at the point"))
    fac(k) = Sym(factorial(big(k)))
    C = _psgeom(h, k -> iseven(k) ? (isodd(k ÷ 2) ? -ONE : ONE) / fac(k) : ZERO, N)
    S = _psgeom(h, k -> isodd(k) ? (isodd((k - 1) ÷ 2) ? -ONE : ONE) / fac(k) : ZERO, N)
    cR = _psadd(_psscale(C, cos(r0)), _psscale(S, -sin(r0)))
    sR = _psadd(_psscale(C, sin(r0)), _psscale(S, cos(r0)))
    _asbuild() do add
        if name === :cos
            add(_PS0, :cos, H, cR); add(_PS0, :sin, H, _psscale(sR, MINUS_ONE))
        else
            add(_PS0, :sin, H, cR); add(_PS0, :cos, H, sR)
        end
    end
end

# Stirling's series: log Γ(u) = (u - 1/2) log u - u + log(2π)/2 + Σ B_2k/(2k(2k-1)) u^(1-2k)
function _psloggamma(u::_PS, N::_Ex)
    w = _psinv(u, N)
    L = _psadd(_psmul(_psadd(u, _psconst(Sym(-1 // 2))), _pslog(u, N)), _psscale(u, MINUS_ONE))
    tail = _psgeom(w, k -> isodd(k) ? (m = (k + 1) ÷ 2; Sym(SymNumber(_bernoulli(2m))) / (2m * (2m - 1))) : ZERO, N)
    _psadd(_psadd(L, _psconst(log(2PI) / 2)), tail)
end

# the series in e^G … pieces of ex at ±∞ (t > 0), to order N where possible
function _aser(ex::Sym, N::_Ex)
    has(ex, _T) || return _asconst(ex)
    err0 = try
        return _asplain(_ser(ex, N))
    catch err
        err isa ArgumentError || rethrow()
        err
    end
    k = ex.kind
    if k == ADD
        acc = _asconst(Sym(ex.num))
        for i in eachindex(ex.args)
            acc = _asadd(acc, _asscale(_aser(ex.args[i], N), Sym(ex.coefs[i])))
        end
        return acc
    elseif k == MUL
        acc = _asconst(Sym(ex.num))
        for i in eachindex(ex.args)
            acc = _asmul(acc, _aspow(ex.args[i], ex.exps[i], N))
        end
        return acc
    elseif k == POW
        return _aspow(ex.args[1], ex.args[2], N)
    elseif k == FN
        ex.name in (:sin, :cos) && length(ex.args) == 1 && return _astrig(ex.name, _ser(only(ex.args), N), N)
        r = _asfn(ex.name, ex.args, N)
        r === nothing || return r
    end
    throw(err0)
end

function _aspow(b::Sym, e::Sym, N::_Ex)
    try
        return _asplain(_serpow(b, e, N))
    catch err
        err isa ArgumentError || rethrow()
        if _eq(b, E)
            return _asexp(_ser(e, N), N)
        elseif has(e, _T)
            return _asexp(_psmul(_ser(e, N), _pslog(_ser(b, N), N)), N)
        elseif e.kind == NUMBER && isexact(e.num) && _isint(e.num)
            return _aspowk(_aser(b, N), _toint(e.num), N)
        end
        rethrow()
    end
end

# special functions of an argument u → ±∞ (nothing: not handled here)
function _asfn(name::Symbol, args::Vector{Sym}, N::_Ex)
    name in (:erf, :erfc, :erfi, :erfcx, :expinti, :expint, :sinint, :cosint,
             :besselj, :bessely, :besseli, :besselk, :gamma, :airyai, :airybi) || return nothing
    any(a -> has(a, _T), args[1:end-1]) && return nothing
    name in (:gamma, :airyai, :airybi) && length(args) != 1 && return nothing
    u = _ser(args[end], N)
    any(<(0), u.e) || return nothing                 # finite at the point
    _, c0, _ = _pslead(u)
    sg = has(c0, _LOGT) ? nothing : _csign(c0)
    sg === nothing && throw(ArgumentError("series: $(apply_function(name, args)) at the point: the sign of $c0 is not known"))
    w = _psinv(u, N)                                 # 1/u, of positive order
    W(coef) = _psgeom(w, coef, N)
    fac(k) = Sym(factorial(big(k)))
    dfac(k) = Sym(prod(big(j) for j in k:-2:1; init = big(1)))      # k!!, (-1)!! = 1
    sgn(j) = isodd(j) ? -ONE : ONE
    # √π e^(u²) erfc(u) for u → +∞: Σ (-1)^j (2j-1)!!/2^j u^(-2j-1); without the signs for erfi
    A(alt) = W(k -> isodd(k) ? (j = (k - 1) ÷ 2; (alt ? sgn(j) : ONE) * dfac(2j - 1) / Sym(big(2)^j) / sqrt(PI)) : ZERO)
    E2(s) = _asexp(_psscale(_psmul(u, u), Sym(s)), N)                # e^(±u²)
    if name in (:erf, :erfc)
        # erfc(u) = (1 - sg) + e^(-u²) A(u) on either side (A is odd); erf = 1 - erfc
        ec = _asmul(E2(-1), _asplain(A(true)))
        name === :erfc && return _asadd(_asconst(1 - sg), ec)
        return _asadd(_asconst(sg), _asscale(ec, MINUS_ONE))
    elseif name === :erfi
        return _asmul(E2(1), _asplain(A(false)))     # u real: no constant term
    elseif name === :erfcx                           # e^(u²) erfc(u)
        sg > 0 && return _asplain(A(true))
        return _asadd(_asscale(E2(1), Sym(2)), _asplain(A(true)))
    elseif name === :expinti                         # Ei(u) = e^u Σ_{k≥1} (k-1)!/u^k, u → ±∞
        return _asmul(_asexp(u, N), _asplain(W(k -> k >= 1 ? fac(k - 1) : ZERO)))
    elseif name === :expint                          # E_ν(u) = e^(-u) Σ_{k≥1} (-1)^(k-1) (ν)_(k-1)/u^k
        sg > 0 || return nothing
        ν = length(args) == 2 ? args[1] : ONE
        poch(m) = (p = ONE; for j in 0:(m - 1); p = p * (ν + j); end; expand(p))
        return _asmul(_asexp(_psscale(u, MINUS_ONE), N), _asplain(W(k -> k >= 1 ? sgn(k - 1) * poch(k - 1) : ZERO)))
    elseif name in (:sinint, :cosint)
        name === :cosint && sg < 0 && return nothing  # Ci at -∞ is complex
        f = W(k -> isodd(k) ? sgn((k - 1) ÷ 2) * fac(k - 1) : ZERO)              # Σ (-1)^j (2j)!/u^(2j+1)
        g = W(k -> iseven(k) && k >= 2 ? sgn((k - 2) ÷ 2) * fac(k - 1) : ZERO)   # Σ (-1)^j (2j+1)!/u^(2j+2)
        cu, su = _astrig(:cos, u, N), _astrig(:sin, u, N)
        name === :sinint &&                           # Si(u) = ±π/2 - cos(u) f - sin(u) g
            return _asadd(_asconst(sg * PI / 2), _asscale(_asadd(_asmul(cu, _asplain(f)), _asmul(su, _asplain(g))), MINUS_ONE))
        return _asadd(_asmul(su, _asplain(f)), _asscale(_asmul(cu, _asplain(g)), MINUS_ONE))   # Ci = sin(u) f - cos(u) g
    elseif name in (:besselj, :bessely, :besseli, :besselk)
        sg > 0 || return nothing
        ν = args[1]
        # a_k(ν) = (4ν² - 1)(4ν² - 9)⋯(4ν² - (2k-1)²)/(k! 8^k)
        a(k) = (p = ONE; for j in 1:k; p = p * (4ν^2 - (2j - 1)^2); end; expand(p / (fac(k) * Sym(big(8)^k))))
        r = _pspow(u, Sym(-1 // 2), N)               # u^(-1/2)
        if name === :besseli                         # e^u/√(2πu) Σ (-1)^k a_k/u^k
            return _asmul(_asexp(u, N), _asplain(_psscale(_psmul(r, W(k -> sgn(k) * a(k))), ONE / (sqrt(Sym(2)) * sqrt(PI)))))
        elseif name === :besselk                     # √(π/(2u)) e^(-u) Σ a_k/u^k
            return _asmul(_asexp(_psscale(u, MINUS_ONE), N), _asplain(_psscale(_psmul(r, W(a)), sqrt(PI) / sqrt(Sym(2)))))
        end
        # √(2/(πu)) (cos ω P - sin ω Q) for J, (sin ω P + cos ω Q) for Y, ω = u - νπ/2 - π/4
        P = _asplain(W(k -> iseven(k) ? sgn(k ÷ 2) * a(k) : ZERO))
        Q = _asplain(W(k -> isodd(k) ? sgn((k - 1) ÷ 2) * a(k) : ZERO))
        ω = _psadd(u, _psconst(-ν * PI / 2 - PI / 4))
        cω, sω = _astrig(:cos, ω, N), _astrig(:sin, ω, N)
        body = name === :besselj ? _asadd(_asmul(cω, P), _asscale(_asmul(sω, Q), MINUS_ONE)) :
                                   _asadd(_asmul(sω, P), _asmul(cω, Q))
        return _asmul(_asplain(_psscale(r, sqrt(Sym(2)) / sqrt(PI))), body)
    elseif name in (:airyai, :airybi)
        # u → +∞, ζ = 2u^(3/2)/3: Ai = e^(-ζ)/(2√π u^(1/4)) Σ (-1)^k c_k/ζ^k and
        # Bi = e^ζ/(√π u^(1/4)) Σ c_k/ζ^k, c_k = (6k)!/(16^k 54^k (3k)! (2k)!) (SymPy #23563)
        sg > 0 || return nothing
        ck(k) = Sym(factorial(big(6k)) // (big(16)^k * big(54)^k * factorial(big(3k)) * factorial(big(2k))))
        zi = _psscale(_pspow(u, Sym(-3 // 2), N), Sym(3 // 2))                  # 1/ζ
        ai = name === :airyai
        S = _psgeom(zi, k -> ai && isodd(k) ? -ck(k) : ck(k), N)
        pre = _psscale(_pspow(u, Sym(-1 // 4), N), ai ? inv(2 * sqrt(PI)) : inv(sqrt(PI)))
        ζ = _psscale(_pspow(u, Sym(3 // 2), N), Sym(ai ? -2 // 3 : 2 // 3))
        return _asmul(_asexp(ζ, N), _asplain(_psmul(pre, S)))
    elseif name === :gamma                           # Stirling's series
        sg > 0 || return nothing
        g = _asexp(_psadd(_psloggamma(u, N), _psconst(-log(2PI) / 2)), N)
        return _asscale(g, sqrt(2PI))
    end
    nothing
end

## the user-facing series

"""
    SymSeries

A truncated series in `x` about `x0`: `Σ cᵢ (x - x0)^eᵢ + O((x - x0)^n)`, or in `1/x` at
`x0 = ±Inf`. At `±Inf` it may also have terms with a factor `exp(g)`, `exp(g) cos(h)` or
`exp(g) sin(h)`, each times its own series in `1/x`, as in `erf(x) = 1 - exp(-x^2)*(…)`. Returned by [`series`](@ref). `Sym(s)` is the truncated expression (the
`O` term dropped); calling `s` calls that expression. Series about the same point add,
subtract, multiply, divide and take powers, with the order of the error tracked, and
`diff(s, x)` differentiates term by term. Not exported.
"""
struct SymSeries
    x::Sym
    x0::Any                  # a Sym, or ±Inf
    ps::_PS                  # in t (x = x0 + t, x = 1/t or x = -1/t), log(t) as _LOGT
    n::_Ex         # the order asked for (bounds the terms of an exact series)
    pcs::Vector{_Piece}      # at ±∞: terms with a factor e^g, e^g cos h or e^g sin h
    pre::Sym                 # a factor t^α of the whole (α not a rational number: x^π, x^a), or 1
end
SymSeries(x, x0, ps::_PS, n) = SymSeries(x, x0, ps, n, _Piece[], ONE)
SymSeries(x, x0, ps::_PS, n, pcs::Vector{_Piece}) = SymSeries(x, x0, ps, n, pcs, ONE)

"""
    series(ex, x, x0 = 0, n = 6)

The series of `ex` in `x` about `x0`, with the error `O((x - x0)^n)`: the terms of order
below `n` are exact. At `x0 = Inf` (or `-Inf`) it is a series in `1/x`, with error
`O((1/x)^n)`. Laurent series (negative powers), Puiseux series (fractional powers) and
logarithmic terms (`log(x)` in the coefficients, as for `x^x` at 0) are included.
Fractional powers and logarithms are those of `x - x0 > 0`: such a series is the
expansion on the right of `x0` (`sqrt(-x)` at 0 is `im*sqrt(x)`, and `acos(1 + x)`
is `-im*sqrt(2)*sqrt(x) + …`, the value `N` gives above the branch cut).

At `±Inf`, an expansion may need an exponential or oscillating factor: terms
`exp(g)*(series)`, `cos(h)*(series)` and `sin(h)*(series)` with `g`, `h` growing at the
point. These come from `exp`, `sin` and `cos` of growing arguments and from the asymptotic
expansions of `erf`, `erfc`, `erfi`, `erfcx`, `expinti`, `expint`, `sinint`, `cosint`,
`besselj`, `bessely`, `besseli`, `besselk`, `airyai` and `airybi` (at `+Inf`) and `gamma`
(Stirling's series; `loggamma`'s is an ordinary series). Such a term is not dropped even when it is exponentially small:
`erf(x)` at `Inf` is `1 - exp(-x^2)*(1/(√π x) - …)`, not `1`. For `besseli` the term of
relative size `exp(-2x)` is left out, as is usual for real arguments.

Exponents may be irrational constants: `x + x^√2`, `1/(1 + x^π)` and `exp(x^(1/π))` expand
with terms such as `x^(2π)` and `x^(√2 + 1)`, ordered by their numeric values.

A power with an exponent that is not a rational number (`x^π`, `x^a`, `sin(x)^√2`) of a
base vanishing or infinite at the point, when it is a factor of the whole expression, is a
factor of the series: `x^a*(series)`, its terms and error counted in the series (`x^a*exp(x)`
is `x^a*(1 + x + x^2/2 + O(x^3))`), the base written `b = c t^k (b/(c t^k))` with
`t = x - x0`; at `±Inf` the factor is a power of `x` (of `-x` at `-Inf`). Series with the
same factor add, and a constant factor such as `x^π` is moved into the exponents when added
to a series without it; products, powers and `diff` combine the factors. A symbolic
exponent elsewhere (`x^a + x^2`, which term leads depending on `a`) is an error.

The result is a [`SymSeries`](@ref); `Sym(s)` gives the truncated expression.

```julia
series(sin(x), x)              # x - x^3/6 + x^5/120 + O(x^6)
series(1/sin(x), x, 0, 3)      # 1/x + x/6 + O(x^3)
series(sqrt(x + x^2), x, 0, 3) # sqrt(x) + x^(3/2)/2 - x^(5/2)/8 + O(x^3)
series(exp(x), x, 1, 3)        # ℯ + ℯ*(x - 1) + ℯ*(x - 1)^2/2 + O((x - 1)^3)
series(sqrt(x^2 + x), x, Inf, 3)   # x + 1/2 - 1/(8x) + 1/(16x^2) + O(1/x^3)
series(x^x, x, 0, 2)           # 1 + log(x)*x + O(x^2)
series(erf(x), x, Inf, 4)      # 1 + exp(-x^2)*(-1/(sqrt(π)*x) + 1/(2*sqrt(π)*x^3) + O(1/x^4))
series(x^a * exp(x), x, 0, 3)  # x^a*(1 + x + x^2/2 + O(x^3))
series(1 / (1 + x^π), x, 0, 7) # 1 - x^π + x^(2π) + O(x^7)
```

An expression with no series at the point (an essential singularity such as `exp(1/x)`
at 0, or a function without a derivative rule) gives an `ArgumentError`, as does a
reciprocal of a sum with such factors, as `1/(1 + exp(x))` at `Inf`.
"""
function series(ex, x::Sym, x0 = 0, n::Real = 6)
    _isnamed(x) || throw(ArgumentError("series: $x is not a symbol"))
    ex = _sinc_elementary(Sym(ex))
    # a symbolic ±∞ (what limit returns) is the point at infinity, not a finite point
    # Inf + t (series(1/(n+1), n, Sym(Inf)) gave 0; SymPy #16777)
    if x0 isa Sym && x0.kind == NUMBER
        v = _jlnum(x0.num)
        v isa Real && isinf(v) && (x0 = v)
    end
    x0 isa Sym && x0.kind != NUMBER && !_hassymbol(x0) && (try isinf(N(x0)) catch; false end) &&
        throw(ArgumentError("series: the point $x0 is infinite; use Inf or -Inf"))
    # powers b^e with e not a rational number (x^π, x^a) of a base vanishing or infinite at
    # the point: their t^(k e) (b ~ c t^k) in front of the series of the rest, whose base
    # tends to 1 (x^a (1 + x) is x^a (1 + x); SymPy #4047, #7168)
    pf = _power_prefactor(ex, x, x0)
    if pf !== nothing
        pre, rest = pf
        r = series(rest, x, x0, n)
        r isa SymSeries || (r = SymSeries(x, x0 isa Real && isinf(x0) ? x0 : Sym(x0), _psconst(Sym(r)), _Ex(rationalize(n))))
        return SymSeries(r.x, r.x0, r.ps, r.n, r.pcs, r.pre * pre)
    end
    sub = x0 isa Real && isinf(x0) ? (x0 > 0 ? inv(_T) : -inv(_T)) : Sym(x0) + _T
    ext = subs(ex, x => sub)
    order = _Ex(rationalize(n))
    atinf = x0 isa Real && isinf(x0)
    N = order
    # the terms below the order; exact (a polynomial, say) if nothing was cut off
    trim(p) = (keep = findall(<(order), p.e);
        _PS(p.e[keep], p.c[keep], p.o == _INF && length(keep) == length(p.e) ? _INF : order))
    for _ in 1:12
        r = try
            with(() -> atinf ? _aser(ext, N) : _asplain(_ser(ext, N)), _TPOS => atinf)
        catch err
            err isa _SeriesPrecision || rethrow()
            N += 2
            continue
        end
        o = minimum(p.o for p in vcat([r.plain], [q.p for q in r.pcs]))
        if o >= order
            plain, pcs = trim(r.plain), _Piece[]
            for q in r.pcs
                p = trim(q.p)
                if isempty(p.e) && isempty(q.G.e)    # cos(h) O(…) is in O(…): |cos h| ≤ 1
                    plain = _PS(plain.e, plain.c, min(plain.o, p.o))
                else
                    push!(pcs, _Piece(q.G, q.k, q.H, p))
                end
            end
            return SymSeries(x, atinf ? x0 : Sym(x0), plain, order, pcs)
        end
        N += max(order - o, 1 // 1) + 1
    end
    throw(ArgumentError("series: could not determine the series of $ex to order $n (too much cancellation)"))
end

# (t^α, the rest) for the factors b^e of ex with e not a rational number, b → 0 or ∞ at x0:
# b ~ c t^k (its series' leading term, k ≠ 0), b^e = t^(k e) c^e (b/(c t^k))^e (t = x - x0 > 0,
# or 1/x, -1/x at ±∞); nothing when ex has none
function _power_prefactor(ex::Sym, x::Sym, x0)
    pairs = ex.kind == MUL ? collect(zip(ex.args, ex.exps)) : ex.kind == POW ? [(ex.args[1], ex.args[2])] : Tuple{Sym,Sym}[]
    isempty(pairs) && return nothing
    atinf = x0 isa Real && isinf(x0)
    t = atinf ? (x0 > 0 ? inv(x) : -inv(x)) : x - Sym(x0)
    α, rest, found = ZERO, ex.kind == MUL ? Sym(ex.num) : ONE, false
    for (b, e) in pairs
        if has(b, x) && !has(e, x) && !(e.kind == NUMBER && isexact(e.num)) && !(b === E)
            lead = _series_lead(b, x, x0)
            if lead !== nothing && !iszero(lead[1])
                k, c = lead
                α += _sym(k) * e
                rest *= c^e * (b / (c * t^_sym(k)))^e
                found = true
                continue
            end
        end
        rest *= b^e
    end
    found || return nothing
    # (t^α as x^(-α) at +∞, (-x)^(-α) at -∞, where x and -x are positive)
    atinf ? (x0 > 0 ? x : -x)^(-α) : t^α, rest
end

# the leading term c t^k (k rational) of b's series at x0, or nothing (no plain leading term)
function _series_lead(b::Sym, x::Sym, x0)
    for n in (1, 2, 4, 8)
        sb = try series(b, x, x0, n) catch err; err isa ArgumentError || rethrow(); return nothing end
        (isempty(sb.pcs) && isone(sb.pre)) || return nothing
        isempty(sb.ps.e) && continue
        has(sb.ps.c[1], _LOGT) && return nothing
        return sb.ps.e[1], sb.ps.c[1]
    end
    nothing
end

# a series' plain terms, for code that reads them: an error with a power factor (x^a)
_plain_ps(s::SymSeries) = isone(s.pre) ? s.ps : throw(ArgumentError("series: a power factor $(s.pre)"))

# t and log(t) in terms of x
function _tvar(s::SymSeries)
    s.x0 isa Real && isinf(s.x0) && return s.x0 > 0 ? inv(s.x) : -inv(s.x)
    s.x - s.x0
end
# log(t) in terms of x: log(x - x0), or -log(x) at ∞ and -log(-x) at -∞
_logt(s::SymSeries) = s.x0 isa Real && isinf(s.x0) ? -log(s.x0 > 0 ? s.x : -s.x) : log(_tvar(s))
_coef(s::SymSeries, c::Sym) = subs(c, _LOGT => _logt(s))
# t^e in x: at ±∞ an irrational power as x^(-e) ((-x)^(-e)), x positive there ((1/x)^(-√2)
# stays as it is)
_tpow(s::SymSeries, e::_Ex) = s.x0 isa Real && isinf(s.x0) && !_isrational(e) ?
    (s.x0 > 0 ? s.x : -s.x)^(-_sym(e)) : _tvar(s)^_sym(e)
_term(s::SymSeries, e::_Ex, c::Sym) = _coef(s, c) * _tpow(s, e)

_psexpr(s::SymSeries, p::_PS) = _add_many(_term(s, e, c) for (e, c) in zip(p.e, p.c))
# the factor exp(g), exp(g) cos(h) or exp(g) sin(h) of a piece, in x
function _pfactor(s::SymSeries, q::_Piece)
    f = isempty(q.G.e) ? ONE : exp(_psexpr(s, q.G))
    q.k === :one && return f
    h = _psexpr(s, q.H)
    f * (q.k === :cos ? cos(h) : sin(h))
end
Sym(s::SymSeries) = s.pre * _add_many(vcat([_psexpr(s, s.ps)], [_pfactor(s, q) * _psexpr(s, q.p) for q in s.pcs]))
(s::SymSeries)(args...; kws...) = Sym(s)(args...; kws...)

# subs into a series: a parameter's value into the coefficients; x => c x^m (m > 0) about
# 0, the series of the composition, its order m times as high (cos(x) with x => x²;
# SymPy #19120). Other substitutions do not keep it a series about the point: an error
# (Sym(s) is the truncated expression).
function subs(s::SymSeries, ps::Pair...)
    for p in ps
        k, v = Sym(p.first), Sym(p.second)
        if _eq(k, s.x)
            s = _series_scale(s, v)
        elseif !has(v, s.x)
            sc(q::_PS) = _PS(q.e, Sym[subs(c, k => v) for c in q.c], q.o)
            s = _foldpre(SymSeries(s.x, s.x0, sc(s.ps), s.n, _Piece[_Piece(sc(q.G), q.k, sc(q.H), sc(q.p)) for q in s.pcs], subs(s.pre, k => v)))
        else
            throw(ArgumentError("subs: $k => $v does not keep the series a series in $(s.x); use subs(Sym(s), …) for the truncated expression"))
        end
    end
    s
end
# a series' coefficients simplified, its order kept (simplify took Sym(s), the truncated
# expression, dropping O(xⁿ) unannounced; SymPy #14622)
function simplify(s::SymSeries; assume = nothing)
    sc(q::_PS) = _PS(q.e, Sym[simplify(c; assume) for c in q.c], q.o)
    SymSeries(s.x, s.x0, sc(s.ps), s.n, _Piece[_Piece(q.G, q.k, q.H, sc(q.p)) for q in s.pcs], s.pre)
end
numer_denom(s::SymSeries) = throw(ArgumentError(
    "numer_denom: a series has no numerator and denominator; use numer_denom(Sym(s)) for the truncated expression"))

function _series_scale(s::SymSeries, v::Sym)
    bad() = throw(ArgumentError("subs: a series about 0 takes $(s.x) => c*$(s.x)^m (m > 0); use subs(Sym(s), …) for the truncated expression"))
    (s.x0 isa Sym && iszero(s.x0) && isempty(s.pcs) && isone(s.pre)) || bad()
    any(c -> has(c, _LOGT), s.ps.c) && bad()
    m = v.kind == MUL ? sum((q for (b, q) in zip(v.args, v.exps) if _eq(b, s.x)); init = ZERO) :
        v.kind == POW ? (_eq(v.args[1], s.x) ? v.args[2] : ZERO) : (_eq(v, s.x) ? ONE : ZERO)
    (m.kind == NUMBER && isexact(m.num) && m.num > N0) || bad()
    c = cancel(v / s.x^m)
    has(c, s.x) && bad()
    mr = _ex_of(m)
    ps = _PS(s.ps.e .* mr, Sym[ci * c^_sym(e) for (e, ci) in zip(s.ps.e, s.ps.c)],
             s.ps.o == _INF ? _INF : s.ps.o * mr)
    SymSeries(s.x, s.x0, ps, s.n * mr)
end

# a term as printed: the coefficient, then the power of (x - x0), not multiplied out
# (multiplying would distribute the coefficient over x - x0)
function _termstring(s::SymSeries, e::_Ex, c::Sym)
    atinf = s.x0 isa Real && isinf(s.x0)
    (iszero(e) || atinf || s.x0 isa Sym && iszero(s.x0)) && return string(_term(s, e, c))
    p = _tvar(s)^_sym(e)
    ps = e == 1 ? "(" * string(_tvar(s)) * ")" : string(p)
    k = _coef(s, c)
    isone(k) && return ps
    _eq(k, -ONE) && return "-" * ps
    ks = string(k)
    kp = k.kind == ADD ? "(" * ks * ")" : ks
    # a negative power as a quotient: -2/(x - 1), not -2*1/(x - 1) (found in the review of
    # SymPy #3730)
    startswith(ps, "1/") && return kp * ps[2:end]
    kp * "*" * ps
end

function _psstring(s::SymSeries, p::_PS)
    parts = String[]
    for (k, (e, c)) in enumerate(zip(p.e, p.c))
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
    p.o != _INF && push!(parts, " + " * _ostring(s, p.o))
    join(parts)
end
function _ostring(s::SymSeries, o)                   # at ±∞, in powers of 1/x either way
    ot = s.x0 isa Real && isinf(s.x0) ? inv(s.x) : _tvar(s)
    p = s.x0 isa Real && isinf(s.x0) && !_isrational(o) ? s.x^(-_sym(o)) : ot^_sym(o)
    "O(" * (o == 1 && ot.kind == ADD ? "(" * string(ot) * ")" : string(p)) * ")"
end

function Base.show(io::IO, s::SymSeries)
    if !isone(s.pre)                                   # x^a*(1 + x + O(x^3))
        (s.ps.o == _INF && isempty(s.pcs)) && return print(io, Sym(s))   # (exact: x^a, x^a*(x + 1))
        print(io, s.pre, "*(")
        show(io, SymSeries(s.x, s.x0, s.ps, s.n, s.pcs, ONE))
        print(io, ")")
        return
    end
    # with pieces, a plain part that is only O(…) goes last
    onlyo = !isempty(s.pcs) && isempty(s.ps.e)
    str = isempty(s.pcs) || !onlyo ? _psstring(s, s.ps) : ""
    for q in s.pcs
        piece = length(q.p.e) == 1 && q.p.o == _INF ? string(_pfactor(s, q) * _psexpr(s, q.p)) :
                isempty(q.p.e) ? string(_pfactor(s, q)) * "*" * _ostring(s, q.p.o) :
                string(_pfactor(s, q)) * "*(" * _psstring(s, q.p) * ")"
        str = isempty(str) ? piece : startswith(piece, "-") ? str * " - " * piece[2:end] : str * " + " * piece
    end
    onlyo && s.ps.o != _INF && (str *= " + " * _ostring(s, s.ps.o))
    print(io, str)
end
Base.show(io::IO, ::MIME"text/plain", s::SymSeries) = show(io, s)

# arithmetic: series about the same point, or a series and a constant
_same(a::SymSeries, b::SymSeries) = (_eq(a.x, b.x) && (a.x0 isa Sym && b.x0 isa Sym ? _eq(a.x0, b.x0) : isequal(a.x0, b.x0))) ||
    throw(ArgumentError("series about different points or in different variables"))
_lift(s::SymSeries, c) = (c = Sym(c); has(c, s.x) && throw(ArgumentError("add a series, not an expression in $(s.x)")); _psconst(c))
_with(a::SymSeries, ps::_PS, n = a.n, pre = a.pre) = _foldpre(SymSeries(a.x, a.x0, ps, n, _Piece[], pre))
_with(a::SymSeries, r::_AS, n = a.n, pre = a.pre) = _foldpre(SymSeries(a.x, a.x0, r.plain, n, r.pcs, pre))
# a power factor t^q with q a rational number (x^a with a => 1/2) moved into the exponents;
# with irr, q an irrational constant too (x^π, kept as a factor until added to a series
# without it)
function _foldpre(s::SymSeries, irr::Bool = false)
    isone(s.pre) && return s
    s.pre.kind == POW || return s
    β = s.pre.args[2]
    (β.kind == NUMBER && isexact(β.num)) || irr || return s
    q = _ex_of(β)
    q === nothing && return s
    atinf = s.x0 isa Real && isinf(s.x0)
    q = q * (atinf ? -1 : 1)
    sh(p::_PS) = _PS(p.e .+ q, p.c, p.o == _INF ? _INF : p.o + q)
    SymSeries(s.x, s.x0, sh(s.ps), s.n + q, _Piece[_Piece(r.G, r.k, r.H, sh(r.p)) for r in s.pcs], ONE)
end
# series with power prefactors (x^a (…)): added when equal (or both folded into the
# exponents, x^π (…) + x), multiplied together
_samepre(a::SymSeries, b::SymSeries) = _eq(a.pre, b.pre) ||
    throw(ArgumentError("series: adding series with different power factors ($(a.pre), $(b.pre)) is not supported"))
_as(s::SymSeries) = _AS(s.ps, s.pcs)
function Base.:+(a::SymSeries, b::SymSeries)
    _same(a, b)
    _eq(a.pre, b.pre) || ((a, b) = (_foldpre(a, true), _foldpre(b, true)))
    _samepre(a, b)
    _with(a, _asadd(_as(a), _as(b)), min(a.n, b.n))
end
Base.:*(a::SymSeries, b::SymSeries) = (_same(a, b); _with(a, _asmul(_as(a), _as(b)), min(a.n, b.n), a.pre * b.pre))
Base.:-(a::SymSeries) = _with(a, _asscale(_as(a), -ONE))
Base.:-(a::SymSeries, b::SymSeries) = a + (-b)
Base.inv(a::SymSeries) = _with(a, _asinv(_as(a), a.n), a.n, inv(a.pre))
Base.:/(a::SymSeries, b::SymSeries) = a * inv(b)
Base.:^(a::SymSeries, k::Integer) = _with(a, _aspowk(_as(a), k, a.n), a.n, a.pre^k)
function Base.:^(a::SymSeries, r::Rational)
    isempty(a.pcs) || throw(ArgumentError("series: a fractional power of a series with exponential or oscillating factors is not supported"))
    _with(a, _pspow(a.ps, Sym(r), a.n), a.n, a.pre^Sym(r))
end
for op in (:+, :*, :-, :/)
    @eval Base.$op(a::SymSeries, c::Union{Number,Sym}) = $op(a, _with(a, _lift(a, c), a.n, ONE))
    @eval Base.$op(c::Union{Number,Sym}, a::SymSeries) = $op(_with(a, _lift(a, c), a.n, ONE), a)
end
Base.:(==)(a::SymSeries, b::SymSeries) = _eq(a.x, b.x) && isequal(a.x0, b.x0) && _pseq(a.ps, b.ps) && _eq(a.pre, b.pre) &&
    length(a.pcs) == length(b.pcs) &&
    all(any(q.k === r.k && _pseq(q.G, r.G) && _pseq(q.H, r.H) && _pseq(q.p, r.p) for r in b.pcs) for q in a.pcs)

"""
    diff(s::SymSeries, x)

The series of the derivative, term by term; its order is one less (at a finite point).
"""
# d/dt (c t^e) = (e c + dc/dlog t) t^(e-1)
_psderiv(p::_PS) = _ps(p.e .- 1, [_sym(e) * c + diff(c, _LOGT) for (e, c) in zip(p.e, p.c)], p.o - 1)

# the derivative of Σ e^G trig(H) p, D the derivative of a series: e^G (trig(H) (DG p + Dp) + trig'(H) DH p)
function _asdiff(r::_AS, D)
    _asbuild() do add
        add(_PS0, :one, _PS0, D(r.plain))
        for q in r.pcs
            add(q.G, q.k, q.H, _psadd(_psmul(D(q.G), q.p), D(q.p)))
            q.k === :one && continue
            hp = _psmul(D(q.H), q.p)
            q.k === :cos ? add(q.G, :sin, q.H, _psscale(hp, MINUS_ONE)) : add(q.G, :cos, q.H, hp)
        end
    end
end

# ∫ s dx term by term, at a finite point and with no log t in the coefficients:
# ∫ c t^e dt = c t^(e+1)/(e+1), the order one more (integrate of a SymSeries was an error
# in RischIntegrate; SymPy #4831)
function _integrate_series(s::SymSeries, x::Sym)
    _eq(x, s.x) || throw(ArgumentError("integrate: the series is in $(s.x), not $x"))
    (s.x0 isa Real && isinf(s.x0)) && throw(ArgumentError("integrate: a series at ±∞ is not integrated term by term"))
    isempty(s.pcs) || throw(ArgumentError("integrate: a series with exponential terms is not integrated term by term"))
    isone(s.pre) || throw(ArgumentError("integrate: a series with a power factor $(s.pre) is not integrated term by term"))
    p = s.ps
    (any(e -> e == -1, p.e) || any(c -> has(c, _LOGT), p.c)) &&
        throw(ArgumentError("integrate: a series with a 1/(x - x0) or log term is not integrated term by term"))
    SymSeries(s.x, s.x0, _ps(p.e .+ 1, [c / _sym(e + 1) for (e, c) in zip(p.e, p.c)], p.o == _INF ? _INF : p.o + 1), s.n + 1)
end

function Base.diff(s::SymSeries, x::Sym)
    if !isone(s.pre)
        # t^α g: (t^α g)' = t^α (g' + α g/t) in x - x0; in a parameter, t^α (∂g + ∂α log(t) g)
        (s.x0 isa Real && isinf(s.x0)) && throw(ArgumentError("diff: a series at ±∞ with a power factor $(s.pre)"))
        isempty(s.pcs) || throw(ArgumentError("diff: a series with a power factor and exponential terms"))
        α = s.pre.kind == POW ? s.pre.args[2] : s.pre.kind == MUL ? s.pre.exps[1] : ONE
        g = s.ps
        if _eq(x, s.x)
            d = _psadd(_psderiv(g), _psscale(_ps(g.e .- 1, g.c, g.o - 1), α))
            return _with(s, d, s.n - 1)
        end
        dα = diff(α, x)
        d = _ps(g.e, [diff(c, x) + dα * _LOGT * c for c in g.c], g.o)
        return _with(s, d)
    end
    _eq(x, s.x) || return _with(s, _asdiff(_as(s), p -> _ps(p.e, [diff(c, x) for c in p.c], p.o)))
    d = _asdiff(_as(s), _psderiv)
    s.x0 isa Real && isinf(s.x0) || return _with(s, d, s.n - 1)
    # dt/dx: t = 1/x gives -t^2, t = -1/x gives t^2
    _with(s, _asmul(_asplain(_PS([2 // 1], [s.x0 > 0 ? -ONE : ONE], _INF)), d), s.n + 1)
end

"""
    SymJ.residue(f, x, x0 = 0)

The residue of `f` at `x = x0`: the coefficient of `(x - x0)^-1` in its Laurent series
there (0 where `f` is analytic), as SymPy's `residue`. At `x0 = ±∞`, the residue at
infinity, `-Res_{t=0} f(1/t)/t²`. An error where `f` has a branch point (fractional powers
or logarithms in its expansion: no Laurent series) or an essential singularity that
[`series`](@ref) does not expand. Not exported.

```julia
SymJ.residue(1 / sin(x), x, 0)                        # 1
SymJ.residue(exp(x) / x^3, x, 0)                      # 1/2
SymJ.residue((x^3 + 5) / ((x^4 - 1) * (x + 1)), x, -1)   # -9/4
SymJ.residue(1 / (x^2 + 1), x, im)                    # -im/2
```
"""
function residue(f, x::Sym, x0 = 0)
    _isnamed(x) || throw(ArgumentError("residue: $x is not a symbol"))
    p = x0 isa Sym ? x0 : Sym(x0)
    if _isinfpt(p)
        t = _fresh_dummy(Sym(f) + x)
        return -residue(subs(Sym(f), x => 1 / t) / t^2, t, 0)
    end
    s = series(f, x, p, 2)
    s isa SymSeries || return ZERO                    # (a constant)
    for (e, c) in zip(s.ps.e, s.ps.c)
        (isinteger(e) && !has(c, _LOGT)) || throw(ArgumentError(
            "residue: $f has a branch point at $x = $p (its expansion has fractional powers or logarithms), so no residue"))
    end
    i = findfirst(==(-1), s.ps.e)
    i === nothing ? ZERO : s.ps.c[i]
end
