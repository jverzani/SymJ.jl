# Modular GCD of integer polynomials (Brown's dense algorithm, in the form of
# Geddes–Czapor–Labahn, Algorithms 7.1 MGCD / 7.2 PGCD).
#
# Over Z (`_gcd_modular`): compute images of the gcd modulo word-size primes, normalized
# so their leading coefficient is gcd(lc(a), lc(b)); discard unlucky primes by comparing
# leading monomials; combine images by Chinese remaindering; accept a candidate only
# when it divides both inputs. The trial division makes the result certain.
#
# Over Z_p (`_pgcd`): eliminate the last active variable x_k by evaluating it at points
# α, recurse, and rebuild x_k by Newton interpolation. Leading coefficients in Z_p[x_k]
# are handled by scaling each image to g(α), with g = gcd of the leading coefficients.
# Unlucky points give images of too large a degree, detected by leading monomial.
#
# Mod-p polynomials use the same packed lex monomials as `Poly`. With x_k the last
# active variable, two facts keep everything linear-time: evaluating x_k keeps terms
# sorted (the cleared field is the least significant in use), and terms sharing the
# same monomial in the other variables are contiguous runs.

const _MODPRIME_START = 2^31 - 1      # primes below 2^31: products fit in an Int

struct _Ctx
    bits::Int
    n::Int
    p::Int
end

# mod-p polynomial: strictly decreasing packed monomials, coefficients in 1:p-1
struct _MP
    monos::Vector{UInt64}
    coefs::Vector{Int}
end

_mpone() = _MP(UInt64[0], Int[1])
_mpiszero(a::_MP) = isempty(a.monos)

@inline _mulmod(a::Int, b::Int, p::Int) = (a * b) % p
@inline _submod(a::Int, b::Int, p::Int) = (d = a - b; d < 0 ? d + p : d)
@inline _addmod(a::Int, b::Int, p::Int) = (s = a + b; s >= p ? s - p : s)

## primes

function _isprime32(n::Int)
    n < 2 && return false
    for q in (2, 3, 5, 7)
        n % q == 0 && return n == q
    end
    d, s = n - 1, 0
    while iseven(d)
        d >>= 1
        s += 1
    end
    for a in (2, 3, 5, 7)            # deterministic for n < 3_215_031_751
        x = powermod(a, d, n)
        (x == 1 || x == n - 1) && continue
        composite = true
        for _ in 1:(s - 1)
            x = _mulmod(x, x, n)
            if x == n - 1
                composite = false
                break
            end
        end
        composite && return false
    end
    true
end

function _prevprime(n::Int)
    n -= 1
    while !_isprime32(n)
        n -= 1
    end
    n
end

## dense univariate polynomials mod p: Vector{Int}, low degree first, no trailing zeros

_utrim!(a::Vector{Int}) = (while !isempty(a) && a[end] == 0; pop!(a); end; a)
_udeg(a::Vector{Int}) = length(a) - 1

function _ueval(a::Vector{Int}, α::Int, p::Int)
    v = 0
    @inbounds for i in length(a):-1:1
        v = _addmod(_mulmod(v, α, p), a[i], p)
    end
    v
end

_uscale(a::Vector{Int}, c::Int, p::Int) = c == 1 ? a : _utrim!(Int[_mulmod(x, c, p) for x in a])
_umonic(a::Vector{Int}, p::Int) = isempty(a) ? a : _uscale(a, invmod(a[end], p), p)

function _umul(a::Vector{Int}, b::Vector{Int}, p::Int)
    (isempty(a) || isempty(b)) && return Int[]
    r = zeros(Int, length(a) + length(b) - 1)
    @inbounds for i in eachindex(a), j in eachindex(b)
        r[i+j-1] = _addmod(r[i+j-1], _mulmod(a[i], b[j], p), p)
    end
    _utrim!(r)
end

# quotient and remainder
function _udivrem(a::Vector{Int}, b::Vector{Int}, p::Int)
    r = copy(a)
    db = _udeg(b)
    _udeg(r) < db && return Int[], r
    q = zeros(Int, _udeg(r) - db + 1)
    inv = invmod(b[end], p)
    @inbounds for i in _udeg(r):-1:db
        c = _mulmod(r[i+1], inv, p)
        q[i-db+1] = c
        c == 0 && continue
        for j in 0:db
            r[i-db+j+1] = _submod(r[i-db+j+1], _mulmod(c, b[j+1], p), p)
        end
    end
    _utrim!(q), _utrim!(r)
end

# a := a mod b, in place
function _urem!(a::Vector{Int}, b::Vector{Int}, p::Int)
    db = _udeg(b)
    inv = invmod(b[end], p)
    @inbounds for i in _udeg(a):-1:db
        c = _mulmod(a[i+1], inv, p)
        c == 0 && continue
        for j in 0:db
            a[i-db+j+1] = _submod(a[i-db+j+1], _mulmod(c, b[j+1], p), p)
        end
    end
    resize!(a, min(length(a), db))
    _utrim!(a)
end

function _ugcd(a::Vector{Int}, b::Vector{Int}, p::Int)
    isempty(a) && return _umonic(b, p)
    isempty(b) && return _umonic(a, p)
    (length(a) == 1 || length(b) == 1) && return Int[1]
    a, b = copy(a), copy(b)
    while !isempty(b)
        _urem!(a, b, p)
        a, b = b, a
    end
    _umonic(a, p)
end

# gcd of a list, stopping early at 1
function _ugcd_all(qs::Vector{Vector{Int}}, p::Int)
    g = Int[]
    for q in qs
        g = _ugcd(g, q, p)
        length(g) == 1 && break
    end
    g
end

## multivariate mod-p polynomials

function _mp_from(P::Poly{SymNumber}, p::Int)
    monos, coefs = UInt64[], Int[]
    for (m, c) in zip(P.monos, P.coefs)
        r = c.tag == T_SMALL ? mod(c.n, p) : Int(mod(numerator(c.big::Rational{BigInt}), p))
        r == 0 && continue
        push!(monos, m)
        push!(coefs, r)
    end
    _MP(monos, coefs)
end

_mp_scale(a::_MP, c::Int, p::Int) = c == 1 ? a : _MP(a.monos, Int[_mulmod(x, c, p) for x in a.coefs])
_mp_monic(a::_MP, p::Int) = _mpiszero(a) ? a : _mp_scale(a, invmod(a.coefs[1], p), p)

_active(a::_MP) = reduce(|, a.monos; init = UInt64(0))

function _lastactive(mask::UInt64, ctx::_Ctx)
    for i in ctx.n:-1:1
        (mask >> _shift(ctx.bits, ctx.n, i)) & _mask(ctx.bits) != 0 && return i
    end
    0
end

# runs of terms sharing the monomial in the other variables, as dense polys in x_k
function _mp_split(a::_MP, k::Int, ctx::_Ctx)
    s = _shift(ctx.bits, ctx.n, k)
    fmask = _mask(ctx.bits) << s
    uppers, polys = UInt64[], Vector{Int}[]
    for (m, c) in zip(a.monos, a.coefs)
        u = m & ~fmask
        e = Int((m & fmask) >> s)
        if isempty(uppers) || uppers[end] != u
            push!(uppers, u)
            push!(polys, zeros(Int, e + 1))   # first term of a run has the top degree
        end
        polys[end][e+1] = c
    end
    uppers, polys
end

function _mp_join(uppers::Vector{UInt64}, polys::Vector{Vector{Int}}, k::Int, ctx::_Ctx)
    s = _shift(ctx.bits, ctx.n, k)
    monos, coefs = UInt64[], Int[]
    for (u, q) in zip(uppers, polys), e in length(q):-1:1
        q[e] == 0 && continue
        push!(monos, u | (UInt64(e - 1) << s))
        push!(coefs, q[e])
    end
    _MP(monos, coefs)
end

# a with x_k = α (x_k the last active variable, so the result stays sorted)
function _mp_eval(a::_MP, k::Int, α::Int, ctx::_Ctx)
    p, s = ctx.p, _shift(ctx.bits, ctx.n, k)
    fmask = _mask(ctx.bits) << s
    pw = Int[1]
    monos, coefs = UInt64[], Int[]
    am, ac = a.monos, a.coefs
    i, n = 1, length(am)
    @inbounds while i <= n
        u = am[i] & ~fmask
        v = 0
        while i <= n && am[i] & ~fmask == u
            e = Int((am[i] & fmask) >> s)
            while length(pw) <= e
                push!(pw, _mulmod(pw[end], α, p))
            end
            v = _addmod(v, _mulmod(ac[i], pw[e+1], p), p)
            i += 1
        end
        if v != 0
            push!(monos, u)
            push!(coefs, v)
        end
    end
    _MP(monos, coefs)
end

## gcd over Z_p, monic in lex order

function _pgcd(a::_MP, b::_MP, ctx::_Ctx)
    p = ctx.p
    _mpiszero(a) && return _mp_monic(b, p)
    _mpiszero(b) && return _mp_monic(a, p)
    k = _lastactive(_active(a) | _active(b), ctx)
    k == 0 && return _mpone()
    ua, qa = _mp_split(a, k, ctx)
    ub, qb = _mp_split(b, k, ctx)

    # contents in Z_p[x_k] and primitive parts
    # univariate in x_k: Euclid directly
    if ua == [0] && ub == [0]
        return _mp_join(UInt64[0], [_ugcd(qa[1], qb[1], p)], k, ctx)
    end
    conta, contb = _ugcd_all(qa, p), _ugcd_all(qb, p)
    c = _ugcd(conta, contb, p)
    # a primitive part free of the other variables is a unit
    (ua == [0] || ub == [0]) && return _mp_join(UInt64[0], [c], k, ctx)
    if length(conta) > 1
        qa = [_udivrem(q, conta, p)[1] for q in qa]
        a = _mp_join(ua, qa, k, ctx)
    end
    if length(contb) > 1
        qb = [_udivrem(q, contb, p)[1] for q in qb]
        b = _mp_join(ub, qb, k, ctx)
    end

    lca, lcb = qa[1], qb[1]
    g = _ugcd(lca, lcb, p)
    bound = min(maximum(_udeg, qa), maximum(_udeg, qb)) + _udeg(g)

    # interpolant: sorted other-variable monomials Hu, each with a dense poly in x_k
    Hu, Hh = UInt64[], Vector{Int}[]
    q = Int[1]                         # product of (x_k - α) over the points used
    lmH = typemax(UInt64)
    npts = 0
    α = Int(hash((p, k)) % UInt(p))
    for _ in 1:p
        α = α + 1 == p ? 0 : α + 1
        (_ueval(lca, α, p) == 0 || _ueval(lcb, α, p) == 0) && continue
        C = _pgcd(_mp_eval(a, k, α, ctx), _mp_eval(b, k, α, ctx), ctx)
        lmC = C.monos[1]
        lmC > lmH && continue                       # unlucky point
        if lmC < lmH                                 # all earlier points were unlucky
            Hu, Hh = UInt64[], Vector{Int}[]
            q, lmH, npts = Int[1], lmC, 0
        end
        Hu, Hh = _newton_step(Hu, Hh, C, _ueval(g, α, p), q, α, p)
        q = _umul(q, Int[p - α, 1], p)
        npts += 1
        if npts > bound
            keep = findall(!isempty, Hh)
            Hu, Hh = Hu[keep], Hh[keep]
            cont = _ugcd_all(Hh, p)
            length(cont) > 1 && (Hh = [_udivrem(h, cont, p)[1] for h in Hh])
            length(c) > 1 && (Hh = [_umul(h, c, p) for h in Hh])
            return _mp_monic(_mp_join(Hu, Hh, k, ctx), p)
        end
    end
    error("modular gcd ran out of evaluation points mod $p")
end

# H += (s*C - H(α)) / q(α) * q, merging the sorted monomials of H and C
function _newton_step(Hu, Hh, C::_MP, s::Int, q::Vector{Int}, α::Int, p::Int)
    w = invmod(_ueval(q, α, p), p)
    nu, nh = UInt64[], Vector{Int}[]
    sizehint!(nu, length(Hu) + length(C.monos))
    sizehint!(nh, length(Hu) + length(C.monos))
    i, j = 1, 1
    @inbounds while i <= length(Hu) || j <= length(C.monos)
        if j > length(C.monos) || (i <= length(Hu) && Hu[i] > C.monos[j])
            u, h, v = Hu[i], Hh[i], 0
            i += 1
        elseif i > length(Hu) || C.monos[j] > Hu[i]
            u, h, v = C.monos[j], Int[], _mulmod(C.coefs[j], s, p)
            j += 1
        else
            u, h, v = Hu[i], Hh[i], _mulmod(C.coefs[j], s, p)
            i += 1; j += 1
        end
        δ = _mulmod(_submod(v, _ueval(h, α, p), p), w, p)
        if δ != 0
            length(h) < length(q) && append!(h, zeros(Int, length(q) - length(h)))
            for t in eachindex(q)
                h[t] = _addmod(h[t], _mulmod(δ, q[t], p), p)
            end
            _utrim!(h)
        end
        push!(nu, u)
        push!(nh, h)
    end
    nu, nh
end

## gcd over Z

function _int_content(P::Poly{SymNumber})
    g = big(0)
    for c in P.coefs
        g = gcd(g, big(_numerator(c)))
        isone(g) && break
    end
    g
end

function _primitive(P::Poly{SymNumber})
    c = _int_content(P)
    _positive(isone(c) ? P : _scale(P, inv(SymNumber(c))))
end

_symmetric(r::Integer, m::Integer) = r > m >> 1 ? r - m : r

function _gcd_modular(a::Poly{SymNumber}, b::Poly{SymNumber})
    iszero(a) && return _positive(b)
    iszero(b) && return _positive(a)
    gs, bits, n = a.gens, a.bits, length(a.gens)
    c = SymNumber(gcd(_int_content(a), _int_content(b)))
    constant = _constpoly(SymNumber, gs, bits, c)
    a, b = _primitive(a), _primitive(b)
    (a.monos[1] == 0 || b.monos[1] == 0) && return constant
    lca, lcb = big(_numerator(a.coefs[1])), big(_numerator(b.coefs[1]))
    γ = gcd(lca, lcb)

    Hm = UInt64[]                    # candidate (gcd scaled to lc γ), coefficients mod M
    Hc = BigInt[]
    M = big(1)
    lmH = typemax(UInt64)
    p = _MODPRIME_START + 1
    while true
        p = _prevprime(p)
        (lca % p == 0 || lcb % p == 0) && continue
        ctx = _Ctx(bits, n, p)
        C = _pgcd(_mp_from(a, p), _mp_from(b, p), ctx)
        C.monos[1] == 0 && return constant            # lucky images are never too small
        lmC = C.monos[1]
        lmC > lmH && continue                          # unlucky prime
        C = _mp_scale(C, Int(mod(γ, p)), p)
        stable = false
        if lmC < lmH                                   # first prime, or all earlier unlucky
            Hm, Hc, M, lmH = C.monos, BigInt[_symmetric(x, p) for x in C.coefs], big(p), lmC
        else                                           # Chinese remaindering
            Hm2, Hc2, M = _crt(Hm, Hc, M, C, p)
            stable = Hm2 == Hm && Hc2 == Hc
            Hm, Hc = Hm2, Hc2
        end
        # Trial division only when the reconstruction looks finished: unchanged by the
        # last prime, or every coefficient far below the modulus (a wrong reconstruction
        # has residues spread over the whole range).
        stable || maximum(abs, Hc) <= M >> 8 || continue
        keep = findall(!iszero, Hc)
        G = _primitive(Poly{SymNumber}(gs, bits, Hm[keep], SymNumber.(Hc[keep])))
        divides(a, G)[1] && divides(b, G)[1] && return _scale(G, c)
    end
end

function _crt(Hm::Vector{UInt64}, Hc::Vector{BigInt}, M::BigInt, C::_MP, p::Int)
    monos, coefs = UInt64[], BigInt[]
    Mp = M * p
    w = invmod(Int(mod(M, p)), p)
    i, j = 1, 1
    while i <= length(Hm) || j <= length(C.monos)
        if j > length(C.monos) || (i <= length(Hm) && Hm[i] > C.monos[j])
            m, h, r = Hm[i], Hc[i], 0
            i += 1
        elseif i > length(Hm) || C.monos[j] > Hm[i]
            m, h, r = C.monos[j], big(0), C.coefs[j]
            j += 1
        else
            m, h, r = Hm[i], Hc[i], C.coefs[j]
            i += 1; j += 1
        end
        x = h + M * _mulmod(_submod(r, Int(mod(h, p)), p), w, p)
        x = _symmetric(mod(x, Mp), Mp)
        push!(monos, m)
        push!(coefs, x)
    end
    monos, coefs, Mp
end
