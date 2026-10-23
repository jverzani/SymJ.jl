# Factorization into irreducibles over the integers (equivalently, the rationals).
#
# After the squarefree decomposition, each squarefree primitive factor is split:
#
# Univariate: Zassenhaus.
#   1. a small prime p keeping f squarefree mod p (of a few candidates, the one with
#      fewest modular factors, which keeps recombination cheap);
#   2. factor mod p: distinct-degree factorization, then Cantor–Zassenhaus for the
#      equal-degree parts;
#   3. Hensel lifting of the modular factors to p^k beyond twice lc(f) times the
#      Landau–Mignotte bound on factor coefficients;
#   4. recombination: subsets of the lifted factors, smallest first, each candidate
#      lc(f)·Π mod p^k (symmetric) tested by exact division.
#
# Multivariate: Kronecker substitution x_i -> t^(w_i), mixed-radix weights above each
# variable's degree. It is a ring map, injective on the monomials of every possible
# factor, so the factors of f are the decoded products of subsets of the univariate
# factors of its image, again tested by exact division. Correct, but the image has
# degree ~ Π (deg_i + 1), so it suits small multivariate inputs; Wang's EEZ algorithm
# (evaluation and multivariate Hensel lifting) is the fast method.

"""
    factor(p::Poly; gcd_method = :dense) -> (unit, [f₁ => m₁, f₂ => m₂, ...])
    factor(ex::Sym; gcd_method = :dense)

Factorization into irreducibles over the rationals: `p == unit * Π fᵢ^mᵢ` with each
`fᵢ` irreducible, primitive, of positive leading coefficient. Factors are sorted by
multiplicity, then degree. Coefficients must be exact.

For an expression the result is the factored expression; rational expressions are
reduced with [`cancel`](@ref) and their numerator and denominator factored, and, as
in `cancel`, function applications and non-integer powers are treated as variables:

```julia
factor(x^4 - 1)                 # (x - 1)*(x + 1)*(x^2 + 1)
factor(x^2*y^2 - 1)             # (x*y - 1)*(x*y + 1)
```

`gcd_method = :zippel` uses Zippel's sparse gcd (see [`gcd`](@ref)) for the gcds of the
squarefree decomposition and, for expressions, of `cancel`.

Multivariate inputs are factored through Kronecker substitution, whose cost grows
with the product of the degrees in each variable; very large inputs throw an
`ArgumentError` rather than run for a long time.
"""
factor(p::Poly{SymNumber}; gcd_method::Symbol = :dense) = _with_widening(_factor_poly, p; gcd_method)
function _factor_poly(p::Poly{SymNumber}; gcd_method::Symbol = :dense)
    unit, sqf = _factor_squarefree_poly(p; gcd_method)
    out = Pair{Poly{SymNumber},Int}[]
    for (q, m) in sqf, g in _factor_sqfree(q)
        push!(out, g => m)
    end
    sort!(out; by = fm -> (fm.second, degree(fm.first), string(Sym(fm.first))))
    unit, out
end

factor(ex::Sym; gcd_method::Symbol = :dense) = _with_widening(_factor_sym, ex; gcd_method)
function _factor_sym(ex::Sym; gcd_method::Symbol = :dense)
    n, d = numer_denom(_cancel(ex; gcd_method))
    atoms = Set{Sym}()
    (_poly_atoms!(atoms, n) && _poly_atoms!(atoms, d)) || return ex
    gs = sort!(collect(atoms); lt = _lt)
    num, den = _factor_poly(Poly{SymNumber}(n, gs); gcd_method), _factor_poly(Poly{SymNumber}(d, gs); gcd_method)
    _sqf_sym(num...) / _sqf_sym(den...)
end
factor(x::Number; gcd_method::Symbol = :dense) = Sym(x)

# irreducible factors of a squarefree, primitive, positive polynomial
function _factor_sqfree(q::Poly{SymNumber})
    act = [i for i in eachindex(q.gens) if _deg_in_all(q, i) > 0]
    isempty(act) && return Poly{SymNumber}[]
    total = degree(q)
    total == 1 && return [q]
    length(act) == 1 && return _factor_univariate(q, only(act))
    try
        return _factor_eez(q, act)
    catch err
        err isa OverflowError || rethrow()       # packed exponents of L^(r-1) q
        return _factor_kronecker(q, act)
    end
end

_deg_in_all(p::Poly, i::Int) = maximum((_exponent(m, p.bits, length(p.gens), i) for m in p.monos); init = 0)

## univariate

function _factor_univariate(q::Poly{SymNumber}, i::Int)
    n = length(q.gens)
    f = zeros(BigInt, _deg_in_all(q, i) + 1)
    for (m, c) in zip(q.monos, q.coefs)
        f[_exponent(m, q.bits, n, i) + 1] = big(_numerator(c))
    end
    s = _shift(q.bits, n, i)
    map(_zassenhaus(f)) do g
        M = eltype(q.monos)
        monos = M[M(e - 1) << s for e in length(g):-1:1 if !iszero(g[e])]
        coefs = SymNumber[SymNumber(g[e]) for e in length(g):-1:1 if !iszero(g[e])]
        _positive(Poly{SymNumber}(q.gens, q.bits, monos, coefs))
    end
end

# irreducible factors (primitive, positive leading coefficient) of a squarefree
# primitive f ∈ Z[x], given as coefficients low to high
function _zassenhaus(f::Vector{BigInt})
    n = length(f) - 1
    n <= 1 && return [f]
    p, modfacs = _choose_prime(f)
    length(modfacs) == 1 && return [f]
    if length(modfacs) >= _VANHOEIJ_MIN_FACTORS          # many modular factors: lattices
        fs = _vanhoeij(f, modfacs, p)
        fs === nothing || return fs
    end
    # p^k > 2 |lc| B, with B = 2^n ‖f‖₂ bounding the coefficients of any factor
    B = big(2)^n * (isqrt(sum(c -> c^2, f)) + 1)
    bound = 2 * abs(f[end]) * B
    k, pk = 1, big(p)
    while pk <= bound
        k += 1
        pk *= p
    end
    lifted = _hensel_lift(_zmod(f, pk), modfacs, p, k)
    _recombine(f, lifted, pk)
end

const _FACTOR_PRIMES = filter(_isprime32, 3:2:3000)

# a prime with lc(f) ≢ 0 and f squarefree mod p; of the first few, the one giving the
# fewest modular factors. Returns p and the monic modular factors.
function _choose_prime(f::Vector{BigInt})
    best = nothing
    tried = 0
    for p in _FACTOR_PRIMES
        mod(f[end], p) == 0 && continue
        fp = _umonic(_utrim!(Int[mod(c, p) for c in f]), p)
        _udeg(_ugcd(fp, _uderiv(fp, p), p)) == 0 || continue
        dd = _ddf(fp, p)
        nf = sum(_udeg(g) ÷ d for (g, d) in dd)
        if best === nothing || nf < best[1]
            best = (nf, p, fp, dd)
        end
        tried += 1
        (nf == 1 || tried >= 5) && break
    end
    best === nothing && error("no suitable prime for factoring")
    _, p, fp, dd = best
    facs = Vector{Int}[]
    for (g, d) in dd
        append!(facs, _edf(g, d, p))
    end
    p, facs
end

## arithmetic in Z_p[x] beyond modgcd.jl

_uderiv(a::Vector{Int}, p::Int) = _utrim!(Int[mod((i - 1) * a[i], p) for i in 2:length(a)])

function _usub(a::Vector{Int}, b::Vector{Int}, p::Int)
    r = zeros(Int, max(length(a), length(b)))
    for i in eachindex(a)
        r[i] = a[i]
    end
    for i in eachindex(b)
        r[i] = _submod(r[i], b[i], p)
    end
    _utrim!(r)
end

_umulmod(a::Vector{Int}, b::Vector{Int}, m::Vector{Int}, p::Int) = _udivrem(_umul(a, b, p), m, p)[2]

function _upowmod(a::Vector{Int}, e::Integer, m::Vector{Int}, p::Int)
    r, b = Int[1], _udivrem(a, m, p)[2]
    while e > 0
        isodd(e) && (r = _umulmod(r, b, m, p))
        e >>= 1
        e > 0 && (b = _umulmod(b, b, m, p))
    end
    r
end

# extended Euclid: (g, s, t) with s*a + t*b = g monic
function _uxgcd(a::Vector{Int}, b::Vector{Int}, p::Int)
    r0, r1 = a, b
    s0, s1 = Int[1], Int[]
    t0, t1 = Int[], Int[1]
    while !isempty(r1)
        q, r = _udivrem(r0, r1, p)
        r0, r1 = r1, r
        s0, s1 = s1, _usub(s0, _umul(q, s1, p), p)
        t0, t1 = t1, _usub(t0, _umul(q, t1, p), p)
    end
    c = invmod(r0[end], p)
    _uscale(r0, c, p), _uscale(s0, c, p), _uscale(t0, c, p)
end

# distinct-degree factorization of a monic squarefree f: [(g, d)], g the product of
# the irreducible factors of degree d
function _ddf(f::Vector{Int}, p::Int)
    out = Tuple{Vector{Int},Int}[]
    w = Int[0, 1]                         # x
    d = 0
    while 2 * (d + 1) <= _udeg(f)
        d += 1
        w = _upowmod(w, p, f, p)
        g = _ugcd(_usub(w, Int[0, 1], p), f, p)
        if _udeg(g) > 0
            push!(out, (g, d))
            f = _udivrem(f, g, p)[1]
            w = _udivrem(w, f, p)[2]
        end
    end
    _udeg(f) > 0 && push!(out, (f, _udeg(f)))
    out
end

# a small deterministic generator for the random choices of Cantor–Zassenhaus
mutable struct _XorShift
    s::UInt64
end
function _next!(r::_XorShift)
    x = r.s
    x ⊻= x << 13; x ⊻= x >> 7; x ⊻= x << 17
    r.s = x
end

# equal-degree factorization (p odd): the monic irreducible factors of degree d of g
function _edf(g::Vector{Int}, d::Int, p::Int, rng = _XorShift(0x9e3779b97f4a7c15))
    _udeg(g) == d && return [g]
    e = (big(p)^d - 1) ÷ 2
    while true
        a = _utrim!(Int[Int(_next!(rng) % UInt64(p)) for _ in 1:_udeg(g)])
        _udeg(a) < 1 && continue
        h = _ugcd(a, g, p)
        if 0 < _udeg(h) < _udeg(g)
            return vcat(_edf(h, d, p, rng), _edf(_udivrem(g, h, p)[1], d, p, rng))
        end
        b = _usub(_upowmod(a, e, g, p), Int[1], p)
        h = _ugcd(b, g, p)
        if 0 < _udeg(h) < _udeg(g)
            return vcat(_edf(h, d, p, rng), _edf(_udivrem(g, h, p)[1], d, p, rng))
        end
    end
end

## integer polynomials (Vector{BigInt}, low to high)

_ztrim!(a::Vector{BigInt}) = (while !isempty(a) && iszero(a[end]); pop!(a); end; a)
_zmod(a::Vector{BigInt}, m::BigInt) = _ztrim!(BigInt[mod(c, m) for c in a])
_zsym(a::Vector{BigInt}, m::BigInt) = _ztrim!(BigInt[(r = mod(c, m); r > m >> 1 ? r - m : r) for c in a])

function _zmul(a::Vector{BigInt}, b::Vector{BigInt})
    (isempty(a) || isempty(b)) && return BigInt[]
    r = zeros(BigInt, length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        r[i+j-1] += a[i] * b[j]
    end
    _ztrim!(r)
end

function _zsub(a::Vector{BigInt}, b::Vector{BigInt})
    r = zeros(BigInt, max(length(a), length(b)))
    for i in eachindex(a)
        r[i] += a[i]
    end
    for i in eachindex(b)
        r[i] -= b[i]
    end
    _ztrim!(r)
end

# a / b if b divides a exactly over Z, else nothing
function _zdivexact(a::Vector{BigInt}, b::Vector{BigInt})
    r = copy(a)
    db = length(b) - 1
    length(r) - 1 < db && return isempty(r) ? BigInt[] : nothing
    q = zeros(BigInt, length(r) - db)
    for i in (length(r) - 1):-1:db
        c, rem_ = divrem(r[i+1], b[end])
        iszero(rem_) || return nothing
        q[i-db+1] = c
        iszero(c) && continue
        for j in 0:db
            r[i-db+j+1] -= c * b[j+1]
        end
    end
    all(iszero, r) ? _ztrim!(q) : nothing
end

function _zprimitive(a::Vector{BigInt})
    g = reduce(gcd, a; init = big(0))
    a = isone(g) ? a : a .÷ g
    a[end] < 0 ? -a : a
end

## Hensel lifting

# F ≡ g*h mod p, g monic, lc(h) = lc(F), gcd(g, h) = 1 mod p: lift to F ≡ G*H mod p^k
# (linear lifting, one power of p per step)
function _hensel2(F::Vector{BigInt}, g::Vector{BigInt}, h::Vector{BigInt}, p::Int, k::Int)
    M = big(p)^k
    gp = _utrim!(Int[mod(c, p) for c in g])
    hp = _utrim!(Int[mod(c, p) for c in h])
    _, s, t = _uxgcd(gp, hp, p)
    pj = big(p)
    for _ in 1:(k - 1)
        e = _zmod(_zsub(F, _zmul(g, h)), M)
        ep = _utrim!(Int[mod(c ÷ pj, p) for c in e])
        if !isempty(ep)
            q, τ = _udivrem(_umul(t, ep, p), gp, p)
            σ = _utrim!(_usub(_umul(s, ep, p), _usub(Int[], _umul(q, hp, p), p), p))
            g = _zadd_scaled(g, τ, pj)
            h = _zadd_scaled(h, σ, pj)
        end
        pj *= p
    end
    _zmod(g, M), _zmod(h, M)
end

function _zadd_scaled(a::Vector{BigInt}, b::Vector{Int}, c::BigInt)
    r = zeros(BigInt, max(length(a), length(b)))
    for i in eachindex(a)
        r[i] = a[i]
    end
    for i in eachindex(b)
        r[i] += c * b[i]
    end
    _ztrim!(r)
end

# F ≡ lc(F) Π facs mod p (facs monic, coprime): monic lifts of the facs mod p^k
function _hensel_lift(F::Vector{BigInt}, facs::Vector{Vector{Int}}, p::Int, k::Int)
    M = big(p)^k
    l = F[end]
    length(facs) == 1 && return [_zmod(F .* invmod(l, M), M)]
    mid = length(facs) ÷ 2
    g = foldl((a, b) -> _umul(a, b, p), facs[1:mid])
    h = _uscale(foldl((a, b) -> _umul(a, b, p), facs[(mid+1):end]), Int(mod(l, p)), p)
    hb = BigInt.(h)
    hb[end] = l                            # exact leading coefficient
    G, H = _hensel2(F, BigInt.(g), hb, p, k)
    vcat(_hensel_lift(G, facs[1:mid], p, k), _hensel_lift(H, facs[(mid+1):end], p, k))
end

## recombination

# subsets of 1:n of size s, in lexicographic order
function _subsets(n::Int, s::Int)
    out = Vector{Int}[]
    idx = collect(1:s)
    s > n && return out
    while true
        push!(out, copy(idx))
        i = s
        while i >= 1 && idx[i] == n - s + i
            i -= 1
        end
        i == 0 && return out
        idx[i] += 1
        for j in (i + 1):s
            idx[j] = idx[j-1] + 1
        end
    end
end

function _recombine(f::Vector{BigInt}, lifted::Vector{Vector{BigInt}}, M::BigInt)
    out = Vector{BigInt}[]
    rest = collect(eachindex(lifted))
    s = 1
    while 2s <= length(rest)
        found = false
        for S in _subsets(length(rest), s)
            cand = BigInt[f[end]]
            for i in S
                cand = _zmod(_zmul(cand, lifted[rest[i]]), M)
            end
            g = _zprimitive(_zsym(cand, M))
            q = _zdivexact(f, g)
            q === nothing && continue
            push!(out, g)
            f = q
            deleteat!(rest, S)
            found = true
            break
        end
        found || (s += 1)
    end
    length(f) > 1 && push!(out, _zprimitive(f))
    out
end

## multivariate: Kronecker substitution

const _KRONECKER_MAX_DEGREE = 4000

function _factor_kronecker(q::Poly{SymNumber}, act::Vector{Int})
    n = length(q.gens)
    degs = [_deg_in_all(q, i) for i in act]
    w = ones(Int, length(act))
    for j in 2:length(act)
        w[j] = w[j-1] * (degs[j-1] + 1)
    end
    D = sum(degs .* w)
    D <= _KRONECKER_MAX_DEGREE || throw(ArgumentError(
        "multivariate factorization by Kronecker substitution would need degree $D; too large"))
    # the image in Z[t], factored (it need not be squarefree)
    t = Sym(:t)
    img = Dict{Int,SymNumber}()
    for (m, c) in zip(q.monos, q.coefs)
        e = sum(_exponent(m, q.bits, n, act[j]) * w[j] for j in eachindex(act))
        img[e] = c
    end
    F = Poly{SymNumber}([t], 64, UInt64[e for e in sort!(collect(keys(img)); rev = true)],
                        SymNumber[img[e] for e in sort!(collect(keys(img)); rev = true)])
    _, ufacs = _factor_poly(F)
    us = Poly{SymNumber}[]
    for (u, m) in ufacs, _ in 1:m
        push!(us, u)
    end
    # recombine: products of subsets, decoded, tested by division
    out = Poly{SymNumber}[]
    rest = collect(eachindex(us))
    s = 1
    while 2s <= length(rest) && !isone(q)
        found = false
        for S in _subsets(length(rest), s)
            prodp = foldl((a, i) -> a * us[rest[i]], S; init = one(F))
            g = _kronecker_decode(prodp, q, act, degs, w)
            g === nothing && continue
            ok, quo = divides(q, g)
            ok || continue
            push!(out, g)
            q = quo
            deleteat!(rest, S)
            found = true
            break
        end
        found || (s += 1)
    end
    q.monos[1] == 0 || push!(out, _positive(q))
    out
end

# the multivariate polynomial whose Kronecker image is u, or nothing if a term of u
# does not decode within the degree bounds
function _kronecker_decode(u::Poly{SymNumber}, q::Poly{SymNumber}, act, degs, w)
    n = length(q.gens)
    M = eltype(q.monos)
    monos, coefs = M[], SymNumber[]
    for (m, c) in zip(u.monos, u.coefs)
        e = Int(m)
        code = zero(M)
        for j in length(act):-1:1
            dj, e = divrem(e, w[j])
            dj > degs[j] && return nothing
            code |= M(dj) << _shift(q.bits, n, act[j])
        end
        push!(monos, code)
        push!(coefs, c)
    end
    perm = sortperm(monos; rev = true)
    _positive(Poly{SymNumber}(q.gens, q.bits, monos[perm], coefs[perm]))
end
