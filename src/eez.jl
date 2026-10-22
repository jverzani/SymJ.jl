# Multivariate factorization by Wang's EEZ approach: evaluation, univariate
# factorization, and multivariate Hensel lifting (Geddes–Czapor–Labahn, Algorithms 6.2
# to 6.4), over ℚ.
#
#   1. x: the active variable of least degree; the content in x (a polynomial in the
#      others) is factored on its own.
#   2. An integer point α for the other variables y with lc_x(q)(α) ≠ 0 and q(x, α)
#      squarefree; of a few such points, the one with fewest univariate factors. A
#      point with one factor proves q irreducible.
#   3. y -> y + α moves the point to 0: Taylor coefficients in y_j are plain
#      coefficients, and evaluation at the point drops terms.
#   4. Leading coefficients: L = lc_x(q) is imposed on every factor of L^(r-1) q (the
#      factors' own leading coefficients divide it); primitive parts at the end.
#   5. Lifting one y_j at a time, solving multivariate Diophantine equations over ℚ;
#      an exact check of the product decides success.
#   6. If α split a true factor, the full lift fails: then two-factor lifts of subsets
#      of the univariate factors, smallest first, so each factor found is irreducible.
# Kronecker substitution (factor.jl) remains the fallback and the reference.

const _P = Poly{SymNumber}

# terms with exponent e in variable j, that exponent removed (order is preserved)
function _coeff_of(p::_P, j::Int, e::Int)
    n = length(p.gens)
    unit = UInt64(e) << _shift(p.bits, n, j)
    monos, coefs = UInt64[], SymNumber[]
    for (m, c) in zip(p.monos, p.coefs)
        _exponent(m, p.bits, n, j) == e || continue
        push!(monos, m - unit)
        push!(coefs, c)
    end
    _P(p.gens, p.bits, monos, coefs)
end
_eval0(p::_P, j::Int) = _coeff_of(p, j, 0)
_ypow(p::_P, j::Int, m::Int) = UInt64(m) << _shift(p.bits, length(p.gens), j)

# p with y_j replaced by y_j + α
function _shift_var(p::_P, j::Int, α::Integer)
    iszero(α) && return p
    n = length(p.gens)
    s = _shift(p.bits, n, j)
    t = _MonoTable{SymNumber}(4 * length(p.monos))
    for (m, c) in zip(p.monos, p.coefs)
        e = _exponent(m, p.bits, n, j)
        base = m - (UInt64(e) << s)
        for k in 0:e
            _addto!(t, base + (UInt64(k) << s), c * SymNumber(binomial(big(e), k) * big(α)^(e - k)))
        end
    end
    monos, coefs = _contents(t)
    _sorted_poly(p.gens, p.bits, monos, coefs)
end

_lcoeff_x(p::_P, xi::Int) = _coeff_of(p, xi, _deg_in_all(p, xi))

function _set_lc(U::_P, xi::Int, newlc::_P)
    d = _deg_in_all(U, xi)
    xd = _ypow(U, xi, d)
    U - _mul_term(_coeff_of(U, xi, d), xd, N1) + _mul_term(newlc, xd, N1)
end

# gcd of the coefficients of p as a polynomial in x (a polynomial in the others)
function _content_x(p::_P, xi::Int)
    g = zero(p)
    for e in unique(_exponent(m, p.bits, length(p.gens), xi) for m in p.monos)
        g = gcd(g, _coeff_of(p, xi, e))
        _isunit(g) && break
    end
    g
end

# the primitive part in x of a factor found over ℚ: integer, x-content 1, positive
function _pp_x(U::_P, xi::Int)
    U = _integer_scaled(U)
    c = _content_x(U, xi)
    _positive(_isunit(c) ? U : divexact(U, c))
end

## multivariate Diophantine equations: Σ σ_i Π_{l≠i} a_l = c with deg_x σ_i < deg_x a_i

struct _DiophantBase
    us::Vector{_P}          # the univariate factors (all y at 0)
    s::Vector{_P}           # s_i = (Π_{l≠i} u_l)^(-1) mod u_i
end
function _DiophantBase(us::Vector{_P})
    s = map(eachindex(us)) do i
        b = foldl(*, (us[l] for l in eachindex(us) if l != i); init = one(us[i]))
        _qinvmod(b, us[i])
    end
    _DiophantBase(us, s)
end

function _mdiophant(as::Vector{_P}, c::_P, vars::Vector{Int}, d::Int, base::_DiophantBase)
    if isempty(vars)                             # univariate in x: as == base.us
        return [rem(c * base.s[i], base.us[i]) for i in eachindex(as)]
    end
    yj = vars[end]
    A = foldl(*, as)
    bs = [divexact(A, a) for a in as]
    anew = [_eval0(a, yj) for a in as]
    σ = _mdiophant(anew, _eval0(c, yj), vars[1:(end-1)], d, base)
    e = c - foldl(+, (σ[i] * bs[i] for i in eachindex(as)))
    for m in 1:d
        iszero(e) && break
        cm = _coeff_of(e, yj, m)
        iszero(cm) && continue
        Δ = _mdiophant(anew, cm, vars[1:(end-1)], d, base)
        ym = _ypow(c, yj, m)
        for i in eachindex(σ)
            δ = _mul_term(Δ[i], ym, N1)
            σ[i] += δ
            e -= δ * bs[i]
        end
    end
    σ
end

## multivariate Hensel lifting with imposed leading coefficients

# us: univariate factors (in x) of Q(x, 0) scaled so that Π us = Q(x, 0); lcs: the
# leading coefficient (in x) each lifted factor must have. The factors, or nothing.
function _mhensel(Q::_P, us::Vector{_P}, lcs::Vector{_P}, xi::Int, ys::Vector{Int})
    k = length(ys)
    As = Vector{_P}(undef, k + 1)
    As[k+1] = Q
    for j in k:-1:1
        As[j] = _eval0(As[j+1], ys[j])
    end
    foldl(*, us) == As[1] || return nothing
    dmax = maximum(_deg_in_all(Q, y) for y in ys; init = 0)
    base = _DiophantBase(us)
    U = copy(us)
    for j in 1:k
        yj = ys[j]
        U1 = copy(U)
        for i in eachindex(U)                    # impose lc, later variables at 0
            lci = lcs[i]
            for l in (j + 1):k
                lci = _eval0(lci, ys[l])
            end
            U[i] = _set_lc(U[i], xi, lci)
        end
        e = As[j+1] - foldl(*, U)
        for m in 1:_deg_in_all(As[j+1], yj)
            iszero(e) && break
            c = _coeff_of(e, yj, m)
            iszero(c) && continue
            σ = _mdiophant(U1, c, ys[1:(j-1)], dmax, base)
            ym = _ypow(Q, yj, m)
            for i in eachindex(U)
                U[i] += _mul_term(σ[i], ym, N1)
            end
            e = As[j+1] - foldl(*, U)
        end
        iszero(e) || return nothing
    end
    U
end

# the irreducible factors of qs (shifted: the point is 0), from the univariate
# factors us of qs(x, 0)
function _eez_split(qs::_P, us::Vector{_P}, xi::Int, ys::Vector{Int})
    length(us) == 1 && return [_positive(qs)]
    L = _lcoeff_x(qs, xi)
    L0 = foldl(_eval0, ys; init = L)                 # a constant
    l0 = L0.coefs[1]
    scale(u) = _scale(u, l0 / _lcoeff_x(u, xi).coefs[1])   # leading coefficient L(0)
    r = length(us)
    # all factors at once
    U = _mhensel(foldl(*, fill(L, r - 1); init = qs), map(scale, us), fill(L, r), xi, ys)
    U === nothing || return [_pp_x(u, xi) for u in U]
    # α split a true factor: two-factor lifts of subsets, smallest first
    for s in 1:(r ÷ 2), S in _subsets(r, s)
        rest = setdiff(1:r, S)
        uA = scale(foldl(*, us[S]))
        uB = scale(foldl(*, us[rest]))
        U2 = _mhensel(L * qs, [uA, uB], [L, L], xi, ys)
        U2 === nothing && continue
        g = _pp_x(U2[1], xi)
        return vcat([g], _eez_split(divexact(qs, g), us[rest], xi, ys))
    end
    [_positive(qs)]
end

## driver

function _factor_eez(q::_P, act::Vector{Int})
    xi = act[argmin([_deg_in_all(q, i) for i in act])]
    ys = [i for i in act if i != xi]
    c = _content_x(q, xi)
    if !_isunit(c)                                   # a factor free of x
        cf = [f for (f, _) in factor(c)[2]]
        return vcat(cf, _factor_sqfree(divexact(q, c)))
    end
    L = _lcoeff_x(q, xi)
    n = length(q.gens)
    # evaluation points: the one with fewest univariate factors of a few
    rng = _XorShift(0x2545f4914f6cdd1d)
    best = nothing
    good = 0
    for t in 1:60
        B = 1 + t ÷ 6
        α = [Int(_next!(rng) % UInt64(2B + 1)) - B for _ in ys]
        val(m) = prod((big(α[j])^_exponent(m, q.bits, n, ys[j]) for j in eachindex(ys)); init = big(1))
        sum((big(_numerator(c)) * val(m) for (m, c) in zip(L.monos, L.coefs)); init = big(0)) == 0 && continue
        u = zeros(BigInt, _deg_in_all(q, xi) + 1)
        for (m, c) in zip(q.monos, q.coefs)
            u[_exponent(m, q.bits, n, xi) + 1] += big(_numerator(c)) * val(m)
        end
        up = Poly{SymNumber}([q.gens[xi]], 64, UInt64[e - 1 for e in length(u):-1:1 if !iszero(u[e])],
                             SymNumber[SymNumber(u[e]) for e in length(u):-1:1 if !iszero(u[e])])
        gcd(up, diff(up, q.gens[xi])).monos[1] == 0 || continue    # squarefree
        fs = _zassenhaus(_zprimitive(_ztrim!(u)))
        length(fs) == 1 && return [q]                              # irreducible
        (best === nothing || length(fs) < length(best[2])) && (best = (α, fs))
        (good += 1) >= 3 && break
    end
    best === nothing && return _factor_kronecker(q, act)
    α, fs = best
    qs = q
    for (j, y) in enumerate(ys)
        qs = _shift_var(qs, y, α[j])
    end
    s = _shift(q.bits, n, xi)
    us = map(fs) do f
        ks = [e for e in length(f):-1:1 if !iszero(f[e])]
        _P(q.gens, q.bits, UInt64[UInt64(e - 1) << s for e in ks], SymNumber[SymNumber(f[e]) for e in ks])
    end
    gs = _eez_split(qs, us, xi, ys)
    map(gs) do g
        for (j, y) in enumerate(ys)
            g = _shift_var(g, y, -α[j])
        end
        _positive(g)
    end
end
