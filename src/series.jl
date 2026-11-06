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
const _TPOS = ScopedValue(false)     # is _T positive (a series at ±∞)

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
    if hasneg && name in (:atan, :acot) && length(args) == 1
        # an infinite argument: atan(a) = ±π/2 - atan(1/a), the sign of a's leading term
        # (SymPy #27109); acot(a) = atan(1/a)
        v0, c0, _ = _pslead(a)
        w = _psinv(a, N)
        at = _psgeom(w, k -> isodd(k) ? Sym((-1)^((k - 1) ÷ 2) // k) : ZERO, N)
        name === :acot && return at
        # the sign of c0 t^v0: t > 0 at ±∞ (x = ±1/t); at a point, only for an even or
        # fractional power (atan(1/x) at 0 is ±π/2 - x by side)
        (_TPOS[] || iseven(numerator(v0)) || denominator(v0) > 1) ||
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
    # inverse functions at their branch points ±1: through square roots, whose series
    # have half-integer powers (acosh(1 + z) = √(2z)(1 - z/12 + …); SymPy #25175)
    if length(args) == 1 && (_eq(a0, ONE) || _eq(a0, MINUS_ONE))
        w = args[1]
        r = if name === :acos
            _eq(a0, ONE) ? 2asin(sqrt((1 - w) / 2)) : PI - 2asin(sqrt((1 + w) / 2))
        elseif name === :asin
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

# Li_s(a) for a → 1, with μ = -log a → 0: Γ(1 - s) μ^(s-1) + Σ ζ(s - k) (-μ)^k/k! for s
# not a positive integer; for s = n ≥ 1 the terms k = n - 1 are instead
# (-μ)^(n-1)/(n-1)! (H_(n-1) - log μ)
function _pspolylog(s::Sym, a::_PS, N::Rational{Int})
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
function _psbessel(name::Symbol, ν::Sym, u::_PS, N::Rational{Int})
    function JI(ν, alt)
        s = _val(u) * Rational{Int}(_numerator(ν.num), _denominator(ν.num))   # of u^ν
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

const _PS0 = _PS(Rational{Int}[], Sym[], _INF)
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

function _asinv(a::_AS, N::Rational{Int})
    isempty(a.pcs) && return _asplain(_psinv(a.plain, N))
    (length(a.pcs) == 1 && _isexact0(a.plain) && a.pcs[1].k === :one) ||
        throw(ArgumentError("series: the reciprocal of a sum with exponential or oscillating factors is not supported"))
    q = a.pcs[1]
    _AS(_PS0, [_Piece(_psscale(q.G, MINUS_ONE), :one, _PS0, _psinv(q.p, N))])
end

function _aspowk(a::_AS, k::Integer, N::Rational{Int})
    k == 0 && return _asconst(ONE)
    b = k > 0 ? a : _asinv(a, N)
    out = b
    for _ in 2:abs(k)
        out = _asmul(out, b)
    end
    out
end

# e^L: the negative powers of L stay in the factor e^G, the rest is a series
function _asexp(L::_PS, N::Rational{Int})
    neg = findall(<(0), L.e)
    isempty(neg) && return _asplain(_psexp(L, N))
    L.o > 0 || throw(_SeriesPrecision())
    rest = setdiff(eachindex(L.e), neg)
    G = _PS(L.e[neg], L.c[neg], _INF)
    _AS(_PS0, [_Piece(G, :one, _PS0, _psexp(_PS(L.e[rest], L.c[rest], L.o), N))])
end

# cos(u), sin(u) for u = H + r0 + h: cos(H + R) = cos H cos R - sin H sin R and so on
function _astrig(name::Symbol, u::_PS, N::Rational{Int})
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
function _psloggamma(u::_PS, N::Rational{Int})
    w = _psinv(u, N)
    L = _psadd(_psmul(_psadd(u, _psconst(Sym(-1 // 2))), _pslog(u, N)), _psscale(u, MINUS_ONE))
    tail = _psgeom(w, k -> isodd(k) ? (m = (k + 1) ÷ 2; Sym(SymNumber(_bernoulli(2m))) / (2m * (2m - 1))) : ZERO, N)
    _psadd(_psadd(L, _psconst(log(2PI) / 2)), tail)
end

# the series in e^G … pieces of ex at ±∞ (t > 0), to order N where possible
function _aser(ex::Sym, N::Rational{Int})
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

function _aspow(b::Sym, e::Sym, N::Rational{Int})
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
function _asfn(name::Symbol, args::Vector{Sym}, N::Rational{Int})
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
    n::Rational{Int}         # the order asked for (bounds the terms of an exact series)
    pcs::Vector{_Piece}      # at ±∞: terms with a factor e^g, e^g cos h or e^g sin h
end
SymSeries(x, x0, ps::_PS, n) = SymSeries(x, x0, ps, n, _Piece[])

"""
    series(ex, x, x0 = 0, n = 6)

The series of `ex` in `x` about `x0`, with the error `O((x - x0)^n)`: the terms of order
below `n` are exact. At `x0 = Inf` (or `-Inf`) it is a series in `1/x`, with error
`O((1/x)^n)`. Laurent series (negative powers), Puiseux series (fractional powers) and
logarithmic terms (`log(x)` in the coefficients, as for `x^x` at 0) are included.

At `±Inf`, an expansion may need an exponential or oscillating factor: terms
`exp(g)*(series)`, `cos(h)*(series)` and `sin(h)*(series)` with `g`, `h` growing at the
point. These come from `exp`, `sin` and `cos` of growing arguments and from the asymptotic
expansions of `erf`, `erfc`, `erfi`, `erfcx`, `expinti`, `expint`, `sinint`, `cosint`,
`besselj`, `bessely`, `besseli`, `besselk`, `airyai` and `airybi` (at `+Inf`) and `gamma`
(Stirling's series; `loggamma`'s is an ordinary series). Such a term is not dropped even when it is exponentially small:
`erf(x)` at `Inf` is `1 - exp(-x^2)*(1/(√π x) - …)`, not `1`. For `besseli` the term of
relative size `exp(-2x)` is left out, as is usual for real arguments.

The result is a [`SymSeries`](@ref); `Sym(s)` gives the truncated expression.

```julia
series(sin(x), x)              # x - x^3/6 + x^5/120 + O(x^6)
series(1/sin(x), x, 0, 3)      # 1/x + x/6 + O(x^3)
series(sqrt(x + x^2), x, 0, 3) # sqrt(x) + x^(3/2)/2 - x^(5/2)/8 + O(x^3)
series(exp(x), x, 1, 3)        # ℯ + ℯ*(x - 1) + ℯ*(x - 1)^2/2 + O((x - 1)^3)
series(sqrt(x^2 + x), x, Inf, 3)   # x + 1/2 - 1/(8x) + 1/(16x^2) + O(1/x^3)
series(x^x, x, 0, 2)           # 1 + log(x)*x + O(x^2)
series(erf(x), x, Inf, 4)      # 1 + exp(-x^2)*(-1/(sqrt(π)*x) + 1/(2*sqrt(π)*x^3) + O(1/x^4))
```

An expression with no series at the point (an essential singularity such as `exp(1/x)`
at 0, or a function without a derivative rule) gives an `ArgumentError`, as does a
reciprocal of a sum with such factors, as `1/(1 + exp(x))` at `Inf`.
"""
function series(ex, x::Sym, x0 = 0, n::Real = 6)
    _isnamed(x) || throw(ArgumentError("series: $x is not a symbol"))
    ex = _sinc_elementary(Sym(ex))
    sub = x0 isa Real && isinf(x0) ? (x0 > 0 ? inv(_T) : -inv(_T)) : Sym(x0) + _T
    ext = subs(ex, x => sub)
    order = Rational{Int}(rationalize(n))
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

# t and log(t) in terms of x
function _tvar(s::SymSeries)
    s.x0 isa Real && isinf(s.x0) && return s.x0 > 0 ? inv(s.x) : -inv(s.x)
    s.x - s.x0
end
# log(t) in terms of x: log(x - x0), or -log(x) at ∞ and -log(-x) at -∞
_logt(s::SymSeries) = s.x0 isa Real && isinf(s.x0) ? -log(s.x0 > 0 ? s.x : -s.x) : log(_tvar(s))
_coef(s::SymSeries, c::Sym) = subs(c, _LOGT => _logt(s))
_term(s::SymSeries, e::Rational{Int}, c::Sym) = _coef(s, c) * _tvar(s)^Sym(SymNumber(e))

_psexpr(s::SymSeries, p::_PS) = _add_many(_term(s, e, c) for (e, c) in zip(p.e, p.c))
# the factor exp(g), exp(g) cos(h) or exp(g) sin(h) of a piece, in x
function _pfactor(s::SymSeries, q::_Piece)
    f = isempty(q.G.e) ? ONE : exp(_psexpr(s, q.G))
    q.k === :one && return f
    h = _psexpr(s, q.H)
    f * (q.k === :cos ? cos(h) : sin(h))
end
Sym(s::SymSeries) = _add_many(vcat([_psexpr(s, s.ps)], [_pfactor(s, q) * _psexpr(s, q.p) for q in s.pcs]))
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
    "O(" * (o == 1 && ot.kind == ADD ? "(" * string(ot) * ")" : string(ot^Sym(SymNumber(o)))) * ")"
end

function Base.show(io::IO, s::SymSeries)
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
_with(a::SymSeries, ps::_PS, n = a.n) = SymSeries(a.x, a.x0, ps, n)
_with(a::SymSeries, r::_AS, n = a.n) = SymSeries(a.x, a.x0, r.plain, n, r.pcs)
_as(s::SymSeries) = _AS(s.ps, s.pcs)
Base.:+(a::SymSeries, b::SymSeries) = (_same(a, b); _with(a, _asadd(_as(a), _as(b)), min(a.n, b.n)))
Base.:*(a::SymSeries, b::SymSeries) = (_same(a, b); _with(a, _asmul(_as(a), _as(b)), min(a.n, b.n)))
Base.:-(a::SymSeries) = _with(a, _asscale(_as(a), -ONE))
Base.:-(a::SymSeries, b::SymSeries) = a + (-b)
Base.inv(a::SymSeries) = _with(a, _asinv(_as(a), a.n))
Base.:/(a::SymSeries, b::SymSeries) = a * inv(b)
Base.:^(a::SymSeries, k::Integer) = _with(a, _aspowk(_as(a), k, a.n))
function Base.:^(a::SymSeries, r::Rational)
    isempty(a.pcs) || throw(ArgumentError("series: a fractional power of a series with exponential or oscillating factors is not supported"))
    _with(a, _pspow(a.ps, Sym(r), a.n))
end
for op in (:+, :*, :-, :/)
    @eval Base.$op(a::SymSeries, c::Union{Number,Sym}) = $op(a, _with(a, _lift(a, c)))
    @eval Base.$op(c::Union{Number,Sym}, a::SymSeries) = $op(_with(a, _lift(a, c)), a)
end
Base.:(==)(a::SymSeries, b::SymSeries) = _eq(a.x, b.x) && isequal(a.x0, b.x0) && _pseq(a.ps, b.ps) &&
    length(a.pcs) == length(b.pcs) &&
    all(any(q.k === r.k && _pseq(q.G, r.G) && _pseq(q.H, r.H) && _pseq(q.p, r.p) for r in b.pcs) for q in a.pcs)

"""
    diff(s::SymSeries, x)

The series of the derivative, term by term; its order is one less (at a finite point).
"""
# d/dt (c t^e) = (e c + dc/dlog t) t^(e-1)
_psderiv(p::_PS) = _ps(p.e .- 1, [Sym(SymNumber(e)) * c + diff(c, _LOGT) for (e, c) in zip(p.e, p.c)], p.o - 1)

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

function Base.diff(s::SymSeries, x::Sym)
    _eq(x, s.x) || return _with(s, _asdiff(_as(s), p -> _ps(p.e, [diff(c, x) for c in p.c], p.o)))
    d = _asdiff(_as(s), _psderiv)
    s.x0 isa Real && isinf(s.x0) || return _with(s, d, s.n - 1)
    # dt/dx: t = 1/x gives -t^2, t = -1/x gives t^2
    _with(s, _asmul(_asplain(_PS([2 // 1], [s.x0 > 0 ? -ONE : ONE], _INF)), d), s.n + 1)
end
