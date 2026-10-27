# Summation: SymJ.summation(f, k, a, b), the sum of f over the integers k = a, ..., b
# (b may be Inf).
#
# A chain of methods, as in SymPy:
#   - linearity: sums split, factors free of k taken out; a term free of k is f*(b - a + 1);
#   - numeric bounds with a short range: the terms added exactly;
#   - an antidifference F (F(k+1) - F(k) = f(k)), the sum F(b+1) - F(a):
#       polynomials in k by Bernoulli polynomials, B_{p+1}(k)/(p+1) for k^p;
#       hypergeometric terms (f(k+1)/f(k) rational in k: polynomials, r^k, factorials,
#       binomials, gamma) by Gosper's algorithm, which finds one or proves there is none;
#   - rational functions without one: partial fractions, Σ 1/(k+α)^m as generalized
#     harmonic numbers for integer shifts, polygamma functions otherwise, zeta values
#     for infinite sums;
#   - harmonic numbers (generalized too) times a term with an antidifference: by
#     parts; times 1/(k + α)^m of their own order: (H(q, m)^2 + H(q, 2m))/2;
#   - an infinite sum is the limit of the antidifference (`limit`, so convergence and
#     the cases of a parameter come from it);
#   - otherwise the sum is left unevaluated, as Sum(f, k, a, b).

"""
    SymJ.summation(f, k, a, b; assume = nothing)

The sum of `f` over the integers `k = a, a+1, ..., b`, in closed form where one is found.
`b` may be `Inf` (an infinite series, whose value is found as a limit, so it may be
`Inf`, or by cases of a parameter); `a` and `b` may be expressions. Assumptions about the
other symbols are given with `assume`, as for [`limit`](@ref).

Polynomials in `k` are summed by Faulhaber's formula, hypergeometric terms (polynomials
times `r^k`, factorials, binomials, rational functions that telescope) by Gosper's
algorithm, and other rational functions through partial fractions, as harmonic numbers
(generalized, `SymJ.harmonic(n, m)`) and polygamma functions (exact zeta values for even
powers). Harmonic numbers times a term with an antidifference are summed by parts. A
sum with no closed form is returned unevaluated, `Sum(f, k, a, b)`. Not exported.

```julia
SymJ.summation(k^2, k, 1, n)            # n^3/3 + n^2/2 + n/6
SymJ.summation(r^k, k, 0, n)            # (r^(n + 1) - 1)/(r - 1)
SymJ.summation(1/(k*(k + 1)), k, 1, Inf)  # 1
SymJ.summation(1/k^2, k, 1, Inf)        # π^2/6
```
"""
function summation(f, k::Sym, a, b; assume = nothing)
    _isnamed(k) || throw(ArgumentError("summation: $k is not a symbol"))
    f = Sym(f)
    infb = b isa Real && isinf(b)
    infb && b < 0 && throw(ArgumentError("summation: the upper bound is -Inf"))
    a isa Real && isinf(a) && throw(ArgumentError("summation from -Inf is not supported"))
    a, bs = Sym(a), infb ? nothing : Sym(b)
    _summation(f, k, a, bs, assume)
end

const _DIRECT_TERMS = 1000

# b === nothing: an infinite sum
function _summation(f::Sym, k::Sym, a::Sym, b, assume)
    iszero(f) && return ZERO
    # numeric bounds, few terms: add them
    if b !== nothing && a.kind == NUMBER && b.kind == NUMBER && _isint(a.num) && _isint(b.num) &&
       _toint(b.num) - _toint(a.num) <= _DIRECT_TERMS
        return _add_many(subs(f, k => Sym(i)) for i in _toint(a.num):_toint(b.num))
    end
    # a polynomial in k as a whole (Σ(2k + 1) is (n + 1)^2, not n(n + 1) + n + 1)
    P = try Poly{Sym}(expand(f), [k]) catch; nothing end
    P === nothing || return _from_antidifference(_poly_antidifference(P, k), k, a, b, assume)
    # the binomial theorem, and series with known sums
    r = _known_sum(f, k, a, b, assume)
    r === nothing || return r
    # linearity
    if f.kind == ADD
        terms = [Sym(c) * t for (c, t) in zip(f.coefs, f.args)]
        iszero(f.num) || push!(terms, Sym(f.num))
        parts = [_summation(t, k, a, b, assume) for t in terms]
        any(p -> p isa Cases, parts) && return _sum_cases(parts)
        return _add_many(parts)
    end
    c, t = _split_kfree(f, k)
    if !isone(c)
        s = _summation(t, k, a, b, assume)
        return s isa Cases ? _map_cases(v -> c * v, s) : c * s
    end
    if !has(f, k)                                    # a constant: (b - a + 1) copies
        b === nothing || return expand(f * (b - a + 1))
        return _summation_unevaluated(f, k, a, b)
    end
    r = _harmonic_sum(f, k, a, b, assume)
    r === nothing || return r
    F = _antidifference(f, k)
    F === nothing || return _from_antidifference(F, k, a, b, assume)
    if b !== nothing                                 # a definite sum: Zeilberger
        # (a pole or an unsolvable step along the way leaves the other methods to try)
        z = try _zeilberger_sum(f, k, a, b) catch err; err isa InterruptException && rethrow(); nothing end
        z === nothing || return z
    end
    r = _rational_sum(f, k, a, b)
    r === nothing || return r
    _summation_unevaluated(f, k, a, b)
end

_summation_unevaluated(f, k, a, b) = apply_function(:Sum, Sym[f, k, a, b === nothing ? Sym(Inf) : b])

# an unevaluated sum whose bounds become numbers (by substitution) is added up
register_function!(:Sum; eval = args -> begin
    length(args) == 4 || return nothing
    f, k, a, b = args
    (a.kind == NUMBER && b.kind == NUMBER && _isint(a.num) && _isint(b.num) &&
     _toint(b.num) - _toint(a.num) <= _DIRECT_TERMS) || return nothing
    _add_many(subs(f, k => Sym(i)) for i in _toint(a.num):_toint(b.num))
end)

# f = c * t with c free of k
function _split_kfree(f::Sym, k::Sym)
    f.kind == MUL || return (ONE, f)
    c, t = Sym(f.num), ONE
    for (base, p) in zip(f.args, f.exps)
        u = base^p
        has(u, k) ? (t *= u) : (c *= u)
    end
    c, t
end

# F(b+1) - F(a), or the limit of F at infinity minus F(a)
function _from_antidifference(F::Sym, k::Sym, a::Sym, b, assume)
    Fa = subs(F, k => a)
    b === nothing || return _present(subs(F, k => b + 1) - Fa)
    L = limit(F, k => Inf; assume)
    L isa Cases && return _map_cases(v -> v isa Sym && !_infinite(v) ? simplify(v - Fa) : v, L)
    _infinite(L) ? L : simplify(L - Fa)
end
_infinite(v) = v isa Sym && v.kind == NUMBER && isfloat(v.num) && isinf(Float64(v.num))

# a closed form, factored when that is smaller (n(n + 1)(2n + 1)/6)
function _present(e::Sym)
    s = simplify(e)
    f = try factor(expand(s)) catch; s end
    _nodes(f) <= _nodes(s) ? f : s
end

## known sums: geometric series, the binomial theorem, hypergeometric series

function _known_sum(f::Sym, k::Sym, a::Sym, b, assume)
    if b === nothing
        # p(k) r^k: converges exactly when |r| < 1, to -F(a) (the antidifference q(k) r^k
        # tends to 0); c r^k: c r^a/(1 - r)
        pg = _poly_geometric(f, k)
        if pg !== nothing
            p, r = pg
            v = has(p, k) ? (F = _gosper(f, k); F === nothing ? nothing : -subs(F, k => a)) : p * r^a / (1 - r)
            v === nothing || return _geometric_infinite(v, p, r, assume)
        end
    end
    if b !== nothing && a.kind == NUMBER && _isint(a.num) && 0 <= _toint(a.num) <= 10
        bt = _binomial_theorem(f, k, b)               # from 0, less the first terms
        bt === nothing || return expand(bt - _add_many(subs(f, k => Sym(j)) for j in 0:(_toint(a.num) - 1)))
    end
    b === nothing && return _hypergeometric_series(f, k, a)
    nothing
end

# f = p(k) r^k with p a polynomial and r free of k: (p, r), or nothing
function _poly_geometric(f::Sym, k::Sym)
    has(f, k) || return nothing
    ρ = cancel(subs(f, k => k + 1) / f)
    n, d = numer_denom(ρ)
    (_ispolyin(n, [k]) && _ispolyin(d, [k])) || return nothing
    dn, dd = _xdeg(n, k), _xdeg(d, k)
    dn == dd || return nothing
    r = cancel(_coef(n, k, dn) / _coef(d, k, dd))
    (has(r, k) || iszero(r)) && return nothing
    p = cancel(f / r^k)
    (has(r, k) || !_ispolyin(p, [k]) || !_ispolyin(numer_denom(p)[2], [k]) || has(numer_denom(p)[2], k)) && return nothing
    (p, r)
end

# f = c * r^k with c, r free of k: (c, r), or nothing
function _geometric(f::Sym, k::Sym)
    has(f, k) || return nothing
    ρ = cancel(subs(f, k => k + 1) / f)
    has(ρ, k) && return nothing
    c = cancel(f / ρ^k)
    has(c, k) && return nothing
    (c, ρ)
end

# the value s of a series p(k) r^k that converges exactly when -1 < r < 1, by cases
function _geometric_infinite(s::Sym, p::Sym, r::Sym, assume)
    if !_hassymbol(r)
        v = abs(_evalc(r))
        v < 1 - 1e-12 && return simplify(s)
        return _divergent(p, r)
    end
    A = _assumptions(assume)
    A === nothing && (A = Assumptions())
    for p in _named(r)
        _ismarkedreal(A, p) || get(A.notin, p, Int8(0)) >= _RR || (A = Assumptions(A, Membership(p, ℝ, false)))
    end
    _by_cases(A) do B
        rr = isempty(B.subst) ? r : subs(r, B.subst...)
        inside = _decide(Relation(:gt, rr, -ONE), B) && _decide(Relation(:lt, rr, ONE), B)
        inside ? simplify(isempty(B.subst) ? s : subs(s, B.subst...)) : ArgumentError("the series diverges")
    end
end
# |r| ≥ 1: to ±∞ for r ≥ 1 and a numeric polynomial factor of one sign, else no sum
function _divergent(p::Sym, r::Sym)
    if r.kind == NUMBER && r.num >= N1 && !_hassymbol(p)
        return Sym(_evalc(p) |> real |> sign > 0 ? Inf : -Inf)
    end
    throw(ArgumentError("summation: the series diverges"))
end

# Σ_{k=0}^{n} binomial(n, k) x^k = (1 + x)^n (x free of k, possibly 1)
function _binomial_theorem(f::Sym, k::Sym, n::Sym)
    fs = f.kind == MUL ? [(b^p) for (b, p) in zip(f.args, f.exps)] : [f]
    i = findfirst(u -> u.kind == FN && u.name === :binomial && _eq(u.args[1], n) && _eq(u.args[2], k), fs)
    i === nothing && return nothing
    rest = cancel(f / fs[i])
    g = rest == ONE ? (ONE, ONE) : _geometric(rest, k)
    g === nothing && return nothing
    c, x = g
    expand(c) == c && !has(c, k) || return nothing
    c * (1 + x)^n
end

# Σ_{k=a}^{∞} f(k) for a hypergeometric term whose ratio f(k+1)/f(k) is
# z/(k + 1) (0F0: exp(z)), z (k + α)/(k + 1) (1F0: (1 - z)^-α), or z/((k + 1)(k + b))
# with b = 1/2 or 3/2 (0F1: cosh, sinh of 2sqrt(z), cos, sin for z < 0), shifted to start
# at a
function _hypergeometric_series(f::Sym, k::Sym, a::Sym)
    ρ = _hyper_ratio(f, k)
    ρ === nothing && return nothing
    ρ = cancel(subs(ρ, k => k + a))
    t0 = subs(f, k => a)
    n, d = numer_denom(ρ)
    (_ispolyin(n, [k]) && _ispolyin(d, [k])) || return nothing
    dn, dd = _xdeg(n, k), _xdeg(d, k)
    # the ratio as z * Π(k + αᵢ) / ((k + 1) Π(k + βⱼ))
    lcn, lcd = _coef(n, k, dn), _coef(d, k, dd)
    z = cancel(lcn / lcd)
    roots_n = dn == 0 ? Sym[] : [expand(-r) for r in solve(n ~ 0, k)]
    roots_d = dd == 0 ? Sym[] : [expand(-r) for r in solve(d ~ 0, k)]
    (length(roots_n) == dn && length(roots_d) == dd) || return nothing
    i1 = findfirst(==(ONE), roots_d)
    i1 === nothing && return nothing                 # the (k + 1) of k!
    deleteat!(roots_d, i1)
    if isempty(roots_n) && isempty(roots_d)
        return simplify(t0 * exp(z))                                     # 0F0
    elseif length(roots_n) == 1 && isempty(roots_d)
        return simplify(t0 * (1 - z)^(-roots_n[1]))                     # 1F0, |z| < 1
    elseif isempty(roots_n) && length(roots_d) == 1 && roots_d[1] in (Sym(1 // 2), Sym(3 // 2))
        # 0F1(; b; z): cosh, sinh of 2 sqrt(z), or cos, sin of 2w for z = -w^2
        half = roots_d[1] == Sym(1 // 2)
        u = _square_root(expand(-4z))                # z = -(u/2)^2
        if u !== nothing
            return simplify(t0 * (half ? cos(u) : sin(u) / u))
        end
        u = _square_root(expand(4z))
        w = u === nothing ? 2sqrt(z) : u
        return simplify(t0 * (half ? cosh(w) : sinh(w) / w))
    end
    nothing
end

_map_cases(g, cs::Cases) = Cases{Any}([c.first => (c.second isa Exception ? c.second : g(c.second)) for c in cs.cases])
_sum_cases(parts) = throw(ArgumentError("summation: a sum of terms whose sums depend on cases is not supported yet; sum the terms separately"))

## harmonic numbers: summation by parts

# Σ_{k=a}^{b} g(k) H(k + c, m), c an integer, g with an antidifference G (G(k+1) - G(k) = g):
#   G(b+1) H(b+c+1, m) - G(a) H(a+c, m) - Σ_{k=a}^{b} G(k+1)/(k+c+1)^m
# by ΔH(k + c, m) = 1/(k+c+1)^m; the last sum is free of H, and H(b+c+1, m) is written
# H(b+c, m) + 1/(b+c+1)^m when smaller. A rational g without one is split into partial
# fractions: the polynomial part by parts, and 1/(k+α)^m (the same m, α an integer) by
#   Σ_{j=1}^{q} H(j, m)/j^m = (H(q, m)^2 + H(q, 2m))/2,
# H(k+c, m) moved to H(k+α, m) by finitely many terms. Other powers are Euler sums,
# without a closed form. Finite sums only, and only when every part closes.
function _harmonic_sum(f::Sym, k::Sym, a::Sym, b, assume)
    b === nothing && return nothing
    hs = _collect_fn!(Sym[], f, :harmonic)
    length(hs) == 1 || return nothing
    h = only(hs)
    c = expand(h.args[1] - k)
    m = _horder(h)
    (c.kind == NUMBER && _isint(c.num) && m.kind == NUMBER && _isint(m.num) && m.num >= N1) || return nothing
    g = cancel(f / h)
    (has(g, k) || !iszero(g)) && isempty(_collect_fn!(Sym[], g, :harmonic)) || return nothing
    # H(a + c) defined: a + c not a negative integer
    ac = expand(a + c)
    ac.kind == NUMBER && ac.num < N0 && return nothing
    G = has(g, k) ? _antidifference(g, k) : g * k
    G === nothing && return _harmonic_sum_pf(g, h, k, a, b, assume)
    # a polynomial G is free up to a constant: G(-c) = 0 cancels a factor k + c + 1
    _ispolyin(G, [k]) && (G = expand(G - subs(G, k => -c)))
    rest = _summation(cancel(subs(G, k => k + 1) / (k + c + 1)^m), k, a, b, assume)
    (rest isa Cases || _count_fn(rest, :Sum) > 0) && return nothing    # (no closed form)
    Gb = subs(G, k => b + 1)
    tail = -subs(G, k => a) * harmonic(ac, m) - rest
    # with H(b+c+1, m), or with H(b+c, m) and its last term apart, whichever is smaller
    q = expand(b + c)
    r1 = _present(Gb) * harmonic(q + 1, m) + _present(tail)
    r2 = _present(Gb) * harmonic(q, m) + _present(_harmonic_lower(Gb / (q + 1)^m + tail, q))
    _nodes(r2) <= _nodes(r1) ? r2 : r1
end
# H(q + 1, p) as H(q, p) + 1/(q + 1)^p, for every order p
function _harmonic_lower(e::Sym, q::Sym)
    hs = filter(h -> _eq(expand(h.args[1] - q), ONE), _collect_fn!(Sym[], e, :harmonic))
    isempty(hs) ? e : subs(e, (h => harmonic(q, _horder(h)) + 1 / (q + 1)^_horder(h) for h in hs)...)
end

# g rational in k without an antidifference, times h = H(k + c, m)
function _harmonic_sum_pf(g::Sym, h::Sym, k::Sym, a::Sym, b::Sym, assume)
    nu, de = numer_denom(g)
    (_ispolyin(nu, [k]) && _ispolyin(de, [k]) && has(de, k)) || return nothing
    pf = try apart(g, k) catch; return nothing end
    terms = pf.kind == ADD ? [Sym(t) * u for (t, u) in zip(pf.coefs, pf.args)] : [pf]
    pf.kind == ADD && !iszero(pf.num) && push!(terms, Sym(pf.num))
    c = _toint(expand(h.args[1] - k).num)
    m = _horder(h)
    poly = ZERO
    out = ZERO
    for t in terms
        if _ispolyin(t, [k])
            poly += t
            continue
        end
        w, u = _split_kfree(t, k)
        (u.kind == POW && u.args[2] == -m && _xdeg(u.args[1], k) == 1) || return nothing
        lin = u.args[1]
        c1 = _coef(lin, k, 1)
        isone(c1) || return nothing
        α = expand(_coef(lin, k, 0))
        (α.kind == NUMBER && _isint(α.num)) || return nothing
        αi = _toint(α.num)
        # H(k + c, m) = H(k + α, m) + D(k), D finitely many terms
        D = c >= αi ? _add_many(1 / (k + j)^m for j in (αi + 1):c) : -_add_many(1 / (k + j)^m for j in (c + 1):αi)
        T(q) = (harmonic(q, m)^2 + harmonic(q, 2m)) / 2
        lo, hi = expand(a + α - 1), expand(b + α)
        lo.kind == NUMBER && lo.num < N0 && return nothing
        out += w * (T(hi) - T(lo))
        iszero(D) && continue
        r = _summation(cancel(D / (k + α)^m), k, a, b, assume)
        (r isa Cases || _count_fn(r, :Sum) > 0) && return nothing
        out += w * r
    end
    if !iszero(poly)
        r = _harmonic_sum(poly * h, k, a, b, assume)
        r === nothing && return nothing
        out += r
    end
    _harmonic_combine(expand(out))
end

## antidifferences

function _antidifference(f::Sym, k::Sym)
    P = try Poly{Sym}(expand(f), [k]) catch; nothing end
    P === nothing || return _poly_antidifference(P, k)
    _gosper(f, k)
end

# Σ c_p k^p: Σ c_p B_{p+1}(k)/(p+1), B the Bernoulli polynomials (B(x+1) - B(x) = (p+1)x^p)
function _poly_antidifference(P, k::Sym)
    out = ZERO
    for (ev, cf) in zip(exponent_vectors(P), coefficients(P))
        p = ev[1]
        out += cf * _bernoulli_poly(p + 1, k) / (p + 1)
    end
    expand(out)
end

const _BERNOULLI = Rational{BigInt}[1]
function _bernoulli(n::Int)
    while length(_BERNOULLI) <= n
        m = length(_BERNOULLI)                      # B_m from B_0, ..., B_{m-1}
        s = sum(binomial(big(m + 1), j) * _BERNOULLI[j+1] for j in 0:(m - 1))
        push!(_BERNOULLI, -s / (m + 1))
    end
    _BERNOULLI[n+1]
end
_bernoulli_poly(n::Int, x::Sym) = _add_many(Sym(SymNumber(binomial(big(n), j) * _bernoulli(j))) * x^(n - j) for j in 0:n)

## Gosper's algorithm

# f(k+1)/f(k) as a rational function of k, factor by factor, or nothing when f is not a
# hypergeometric term in k
function _hyper_ratio(f::Sym, k::Sym)
    has(f, k) || return ONE
    fs = f.kind == MUL ? collect(zip(f.args, f.exps)) : [(f, ONE)]
    f.kind == POW && (fs = [(f.args[1], f.args[2])])
    r = ONE
    for (b, p) in fs
        has(b, k) || has(p, k) || continue
        if !has(p, k)
            (p.kind == NUMBER && _isint(p.num)) || return nothing
            q = _factor_ratio(b, k)
            q === nothing && return nothing
            r *= q^p
        else                                        # b^(c k + d), b free of k
            has(b, k) && return nothing
            c = expand(subs(p, k => k + 1) - p)
            has(c, k) && return nothing
            r *= b^c
        end
    end
    cancel(r)
end

function _factor_ratio(b::Sym, k::Sym)
    if b.kind == FN && length(b.args) <= 2
        u = b.args[end]
        if b.name in (:factorial, :gamma) && length(b.args) == 1
            m = expand(subs(u, k => k + 1) - u)
            (m.kind == NUMBER && _isint(m.num)) || return nothing
            s = _toint(m.num)
            base = b.name === :gamma ? u - 1 : u
            s >= 0 && return _prod_shift(base, s)
            return 1 / _prod_shift(base + s, -s)
        elseif b.name === :binomial && length(b.args) == 2
            n, m = b.args
            if !has(n, k) && expand(subs(m, k => k + 1) - m) == ONE          # binomial(n, k + c)
                return (n - m) / (m + 1)
            elseif !has(m, k) && expand(subs(n, k => k + 1) - n) == ONE      # binomial(k + c, m)
                return (n + 1) / (n + 1 - m)
            end
            # binomial(u, v) = u!/(v! (u - v)!), u and v with integer shifts in k
            rs = [_factor_ratio(factorial(w), k) for w in (n, m, expand(n - m))]
            any(isnothing, rs) && return nothing
            return rs[1] / (rs[2] * rs[3])
        end
        return nothing
    end
    # a rational function of k
    q = cancel(subs(b, k => k + 1) / b)
    nd = numer_denom(q)
    (_ispolyin(nd[1], [k]) && _ispolyin(nd[2], [k])) ? q : nothing
end
# (u + 1)(u + 2)...(u + s)
_prod_shift(u::Sym, s::Int) = prod((u + i for i in 1:s); init = ONE)

# an antidifference of the hypergeometric term f, or nothing (Gosper's algorithm)
function _gosper(f::Sym, k::Sym)
    ρ = _hyper_ratio(f, k)
    ρ === nothing && return nothing
    gf = _gosper_form(ρ, k)
    gf === nothing && return nothing
    A, B, C = gf
    # a(k) x(k+1) - b(k-1) x(k) = c(k) for a polynomial x
    Bm = expand(subs(B, k => k - 1))
    x = _gosper_x(expand(A), Bm, C, k)
    x === nothing && return nothing
    cancel(Bm * x / C * f)
end

# the Gosper form of a ratio: ρ = A(k)/B(k) * C(k+1)/C(k) with gcd(A(k), B(k+h)) = 1
# for every integer h ≥ 0; nothing if ρ is not a rational function of k
function _gosper_form(ρ::Sym, k::Sym)
    A, B = numer_denom(ρ)
    (_ispolyin(A, [k]) && _ispolyin(B, [k])) || return nothing
    C = ONE
    h = Sym(Symbol("_gosper_h"))
    R = try resultant(A, expand(subs(B, k => k + h)), k) catch; return nothing end
    for j in _nonneg_int_roots(R, h)
        g = gcd(A, expand(subs(B, k => k + j)))
        has(g, k) || continue
        A = cancel(A / g)
        B = cancel(B / subs(g, k => k - j))
        for i in 1:j
            C *= subs(g, k => k - i)
        end
    end
    A, B, expand(C)
end

# the nonnegative integer roots of R(h)
function _nonneg_int_roots(R::Sym, h::Sym)
    has(R, h) || return Int[]
    rs = try solve(R ~ 0, h) catch; return Int[] end
    sort!(unique([_toint(r.num) for r in rs if r.kind == NUMBER && _isint(r.num) && r.num >= N0]))
end

# x(k), or with `extra` unknowns in c (Zeilberger's coefficients): (x, their values)
function _gosper_x(a::Sym, b::Sym, c::Sym, k::Sym, extra::Vector{Sym} = Sym[])
    da, db, dc = _xdeg(a, k), _xdeg(b, k), _xdeg(c, k)
    la, lb = _coef(a, k, da), _coef(b, k, db)
    cand = Int[]
    if da != db || !iszero(expand(la - lb))
        push!(cand, dc - max(da, db))
    else
        push!(cand, dc - da + 1)
        q = cancel((_coef(b, k, db - 1) - _coef(a, k, da - 1)) / la)
        q.kind == NUMBER && _isint(q.num) && push!(cand, _toint(q.num))
    end
    d = maximum(cand)
    d < 0 && return nothing
    xs = [Sym(Symbol("_gosper_x", i)) for i in 0:d]
    xk = _add_many(xs[i+1] * k^i for i in 0:d)
    eq = expand(a * subs(xk, k => k + 1) - b * xk - c)
    unknowns = vcat(xs, extra)
    zero_all(e) = expand(subs(e, (v => ZERO for v in unknowns)...))
    if iszero(eq)
        x = zero_all(xk)
        return isempty(extra) ? x : (x, Dict{Sym,Sym}(v => ZERO for v in extra))
    end
    P = try Poly{Sym}(eq, [k]) catch; return nothing end
    eqs = [cf ~ 0 for cf in coefficients(P)]
    sols = try solve(eqs, unknowns) catch; return nothing end
    isempty(sols) && return nothing
    s = first(sols)
    x = zero_all(subs(xk, s...))                                     # free coefficients 0
    isempty(extra) && return x
    x, Dict{Sym,Sym}(v => zero_all(get(s, v, v)) for v in extra)
end
_xdeg(p::Sym, k::Sym) = has(p, k) ? degree(Poly{Sym}(expand(p), [k]), k) : 0
function _coef(p::Sym, k::Sym, d::Int)
    d < 0 && return ZERO
    has(p, k) || return d == 0 ? p : ZERO
    P = Poly{Sym}(expand(p), [k])
    for (ev, cf) in zip(exponent_vectors(P), coefficients(P))
        ev[1] == d && return cf
    end
    ZERO
end

## Zeilberger's algorithm: definite sums of hypergeometric terms

# Σ_{k=a}^{b} F(n, k) with b depending on n, when Gosper finds no antidifference:
# creative telescoping finds c_0(n), ..., c_J(n) and G with
#     Σ_j c_j(n) F(n+j, k) = G(n, k+1) - G(n, k),
# so that, the boundary terms vanishing, S(n) = Σ_k F(n, k) satisfies the recurrence
# Σ_j c_j(n) S(n+j) = 0. A first-order recurrence gives S as a product, from its value
# at a starting n. The closed form is checked against direct sums for several n; the
# boundaries are the summand's own (binomial(n, k) vanishing for k > n), so a sum whose
# bounds cut the summand off elsewhere gives no closed form here.

const _ZEIL_MAX_ORDER = 3

function _zeilberger_sum(f::Sym, k::Sym, a::Sym, b::Sym)
    ns = [s for s in _named(b) if !_eq(s, k)]
    length(ns) == 1 || return nothing
    n = only(ns)
    has(f, n) || return nothing
    S = _zeilberger_sum(f, k, a, b, n)
    S === nothing || return S
    # from a small positive start: the sum from 0, less its first terms
    if a.kind == NUMBER && _isint(a.num) && 0 < _toint(a.num) <= 10
        S = _zeilberger_sum(f, k, ZERO, b, n)
        first = try _add_many(subs(f, k => Sym(j)) for j in 0:(_toint(a.num) - 1)) catch; nothing end
        S === nothing || first === nothing || return S - first
    end
    nothing
end
function _zeilberger_sum(f::Sym, k::Sym, a::Sym, b::Sym, n::Sym)
    for J in 1:_ZEIL_MAX_ORDER
        z = _zeilberger(f, k, n, J)
        z === nothing && continue
        rec, R = z
        value = m -> _direct(f, k, n, a, b, m)
        # the summand's own bounds: the recurrence is homogeneous. First order: a product;
        # higher: the hypergeometric solutions (Hyper) that span S
        sol = J == 1 ? _first_order_solve(rec, f, k, n, a, b) : _recurrence_solve(rec, n, value)
        if sol !== nothing && _zeil_verified(sol[1], f, k, n, a, b, sol[2])
            return _merge_powers_terms(sol[1])
        end
        # other bounds: an inhomogeneous recurrence, its right side the boundary terms (a
        # sum of hypergeometric terms), made homogeneous
        g = _zeil_rhs(rec, R, f, k, n, a, b)
        g === nothing && return nothing
        if !iszero(g)
            hom = _annihilate(rec, g, n)
            sol = hom === nothing ? nothing : _recurrence_solve(hom, n, value)
            sol !== nothing && _zeil_verified(sol[1], f, k, n, a, b, sol[2]) && return _merge_powers_terms(sol[1])
        end
        # else d'Alembertian solutions (harmonic numbers), checked with their sums
        # unevaluated, else interlaced ones (Σ (-1)^k binomial(n, k)^2)
        avoid = vcat(_named(f), _named(a), _named(b))
        sol = _dalembert_solve(rec, g, n, value, avoid)
        sol === nothing && (sol = _liouvillian_solve(rec, g, n, value, avoid))
        sol === nothing && return nothing
        S, n0, raw = sol
        return _zeil_verified(raw, f, k, n, a, b, n0) ? _merge_powers_terms(S) : nothing
    end
    nothing
end

# the right side of Σ_j c_j S(n+j) = g for S(n) = Σ_{k=a}^{b(n)} F(n, k): summing the
# telescoper over k = a..b, g = G(n, b+1) - G(n, a) + Σ_j c_j Σ_{k=b(n)+1}^{b(n+j)} F(n+j, k)
# (the terms of S(n+j) past b(n)). G(n, b+1) is R(n, b+1) F(n, b) F(n, b+1)/F(n, b), the
# ratio as a rational function, so a zero of F there meets a pole of R. For an upper
# bound s n + t (s a positive integer) and a lower one free of n; nothing otherwise.
function _zeil_rhs(rec, R::Sym, F::Sym, k::Sym, n::Sym, a::Sym, b::Sym)
    has(a, n) && return nothing
    sl = expand(subs(b, n => n + 1) - b)
    (sl.kind == NUMBER && _isint(sl.num) && sl.num > N0) || return nothing
    s = _toint(sl.num)
    ρk = _hyper_ratio(F, k)
    ρk === nothing && return nothing
    Fb = subs(F, k => b)
    Gb = cancel(subs(cancel(subs(R, k => k + 1) * ρk), k => b)) * Fb
    Ga = subs(R, k => a) * subs(F, k => a)
    extra = ZERO
    for j in 1:(length(rec) - 1)
        for m in 1:(s * j)
            extra += rec[j+1] * subs(F, n => n + j, k => b + m)
        end
    end
    g = expand(simplify_gamma(Gb - Ga + extra))
    iszero(g) ? ZERO : g
end

# the recurrence coefficients [c_0, ..., c_J] (polynomials in n), or nothing
function _zeilberger(F::Sym, k::Sym, n::Sym, J::Int)
    ρn = _hyper_ratio(F, n)                         # F(n+1, k)/F(n, k)
    ρn === nothing && return nothing
    # r_j = F(n+j, k)/F(n, k)
    r = Sym[ONE]
    for j in 1:J
        push!(r, cancel(r[end] * subs(ρn, n => n + j - 1)))
    end
    cs = [Sym(Symbol("_zeil_c", j)) for j in 0:(J - 1)]
    lin = cancel(_add_many(cs[j+1] * r[j+1] for j in 0:(J - 1)) + r[J+1])     # c_J = 1
    N, Q = numer_denom(lin)
    # t(k) = N(k) F(n, k)/Q(k); the ratio of F/Q, its Gosper form, and Gosper's
    # equation with c(k) N(k) on the right, N linear in the c_j
    ρk = _hyper_ratio(F, k)
    ρk === nothing && return nothing
    ρ = cancel(ρk * Q / subs(Q, k => k + 1))
    gf = _gosper_form(ρ, k)
    gf === nothing && return nothing
    A, B, C = gf
    Bm = expand(subs(B, k => k - 1))
    sol = _gosper_x(expand(A), Bm, expand(C * N), k, cs)
    sol === nothing && return nothing
    x, cvals = sol
    coefs = vcat([cvals[c] for c in cs], ONE)
    all(iszero, coefs[1:end-1]) && return nothing
    # clear denominators: polynomial coefficients
    den = reduce((p, c) -> (d = numer_denom(c)[2]; cancel(p * d / gcd(p, d))), coefs; init = ONE)
    # the certificate: Σ_j c_j F(n+j, k) = G(n, k+1) - G(n, k), G = R F, the antidifference
    # of t = N F/Q from Gosper's x (scaled with the coefficients)
    R = cancel(den * Bm * subs(x, (c => cvals[c] for c in cs)...) / (C * Q))
    [expand(cancel(c * den)) for c in coefs], R
end

# c_0(n) S(n) + c_1(n) S(n+1) = 0: S(n) = S(n0) Π_{j=n0}^{n-1} ρ(j), ρ = -c_0/c_1
function _first_order_solve(rec, f::Sym, k::Sym, n::Sym, a::Sym, b::Sym)
    c0, c1 = rec
    ρ = cancel(-c0 / c1)
    # start past the nonnegative integer roots of c_1 (and of ρ's numerator)
    bad = Int[]
    for p in (numer_denom(ρ)...,)
        has(p, n) || continue
        for r in (try solve(p ~ 0, n) catch; Sym[] end)
            r.kind == NUMBER && _isint(r.num) && push!(bad, _toint(r.num))
        end
    end
    n0 = max(0, (isempty(bad) ? -1 : maximum(bad)) + 1)
    S0 = _direct(f, k, n, a, b, n0)
    S0 === nothing && return nothing
    _present_product(S0 * _hyper_product(ρ, n, n0)), n0
end

_direct(f, k, n, a, b, m) = (aa = subs(a, n => m); bb = subs(b, n => m);
    (aa.kind == NUMBER && bb.kind == NUMBER) ? _add_many(subs(f, n => m, k => Sym(i)) for i in _toint(aa.num):_toint(bb.num)) : nothing)

# Π_{j=n0}^{n-1} ρ(j) for a rational ρ: C^(n - n0) Π Γ(n + αᵢ)/Γ(n0 + αᵢ) / Π Γ(n + βᵢ)/Γ(n0 + βᵢ),
# the Γ of integer and half-integer shifts as factorials
function _hyper_product(ρ::Sym, n::Sym, n0::Int)
    num, den = numer_denom(ρ)
    C = ONE
    gam = ONE
    for (p, sgn) in ((num, 1), (den, -1))
        has(p, n) || (C *= sgn == 1 ? p : 1 / p; continue)
        lc = _coef(p, n, _xdeg(p, n))
        C *= sgn == 1 ? lc : 1 / lc
        roots = solve(p ~ 0, n)
        length(roots) == _xdeg(p, n) || (roots = _roots_with_multiplicity(p, n))
        for r in roots
            α = -r                                   # the factor (n + α)
            g = _gamma_shift(n, α) / _gamma_shift(Sym(n0), α)
            gam *= sgn == 1 ? g : 1 / g
        end
    end
    C^(n - n0) * gam
end
function _roots_with_multiplicity(p::Sym, n::Sym)
    q = factor(p)
    fs = q.kind == MUL ? collect(zip(q.args, q.exps)) : q.kind == POW ? [(q.args[1], q.args[2])] : [(q, ONE)]
    out = Sym[]
    for (b, e) in fs
        has(b, n) || continue
        for r in solve(b ~ 0, n), _ in 1:_toint(e.num)
            push!(out, r)
        end
    end
    out
end
# Γ(m + α) as factorial(m + α - 1); for a half-integer α = j + 1/2 as
# (2(m + j))! sqrt(π)/(4^(m + j) (m + j)!)
function _gamma_shift(m::Sym, α::Sym)
    if α.kind == NUMBER && isexact(α.num) && _denominator(α.num) == 2
        j = expand(m + α - Sym(1 // 2))
        return factorial(2j) * sqrt(PI) / (Sym(4)^j * factorial(j))
    end
    factorial(expand(m + α - 1))
end

# factorial(a)/(factorial(b) factorial(a - b)) written binomial(a, b), and factorials
# whose arguments differ by an integer combined (n!/(n - 1)! = n)
function _present_product(e::Sym)
    s = _merge_powers(simplify(_factorial_ratios(_gauss_multiplication(e))))
    r = _binomial_form(s)
    _nodes(r) < _nodes(s) ? r : s
end
# Gauss's multiplication theorem in a product: factorials (x - 1)! = Γ(x) at
# x = z + j/m for every j = 1, …, m-1, all to the same power, as
#   Π_j Γ(z + j/m) = (2π)^((m-1)/2) m^(1/2 - m z) Γ(m z)/Γ(z)
# (Γ(N + 1/3) Γ(N + 2/3) = 2π 3^(-1/2-3N) (3N)!/N!)
function _gauss_multiplication(e::Sym)
    e.kind == MUL || return e
    facs = [(i, b.args[1], p) for (i, (b, p)) in enumerate(zip(e.args, e.exps))
            if b.kind == FN && b.name === :factorial && length(b.args) == 1]
    length(facs) >= 1 || return e
    for (i, u, p) in facs
        x = expand(u + 1)                          # Γ(x)
        c = x.kind == ADD ? Sym(x.num) : ZERO
        (c.kind == NUMBER && !_isint(c.num) && isexact(c.num)) || continue
        m = Int(denominator(Rational{BigInt}(c.num)))
        m <= 12 || continue
        base = expand(x - c)                       # x = base + c
        fl = floor(Rational{BigInt}(c.num))
        z = expand(base + Sym(SymNumber(fl)))      # x = z + j/m, 0 < j < m
        # Γ(z + j/m) for every j, each to the power p
        idx = Int[]
        for j in 1:(m - 1)
            xj = expand(z + Sym(j // m))
            k = findfirst(t -> _eq(expand(t[2] + 1), xj) && _eq(t[3], p), facs)
            k === nothing && break
            push!(idx, facs[k][1])
        end
        length(idx) == m - 1 || continue
        rest = Sym(e.num) * prod((b^q for (l, (b, q)) in enumerate(zip(e.args, e.exps)) if !(l in idx)); init = ONE)
        g = (2 * PI)^(Sym(m - 1) / 2) * Sym(m)^(Sym(1 // 2) - m * z) * factorial(expand(m * z - 1)) / factorial(expand(z - 1))
        return _gauss_multiplication(rest * g^p)
    end
    e
end

# powers of numbers with the same exponent multiplied ((-4)^n/4^n = (-1)^n), and a
# rational coefficient that is a power of such a base taken in (2^(n-1)/2 = 2^(n-2))
function _merge_powers(e::Sym)
    e.kind == MUL || return e
    c = Sym(e.num)
    groups = Dict{Sym,Sym}()                       # exponent q (less its constant) => base
    other = Sym[]
    for (b, p) in zip(e.args, e.exps)
        if b.kind == NUMBER && isexact(b.num) && _hassymbol(p)
            p0 = p.kind == ADD ? Sym(p.num) : ZERO
            q = expand(p - p0)
            c *= b^p0
            # an integer multiple taken into the base: 3^(3n) = 27^n
            if q.kind == MUL && _isint(q.num) && 1 < abs(q.num) <= 12
                t = abs(q.num)
                b, q = b^t, expand(q / t)
            end
            if haskey(groups, -q)                  # b^(-q) = (1/b)^q
                groups[-q] *= 1 / b
            else
                groups[q] = get(groups, q, ONE) * b
            end
        else
            push!(other, b^p)
        end
    end
    isempty(groups) && return e
    # a sum's rational content taken into the coefficient (2^(n-3)(8n - 16) = 2^n (n - 2))
    for (i, t) in enumerate(other)
        t.kind == ADD && all(isexact, t.coefs) && isexact(t.num) || continue
        g = reduce(gcd, Rational{BigInt}.(t.coefs); init = Rational{BigInt}(t.num))
        (iszero(g) || isone(g)) && continue
        c *= g
        other[i] = expand(t / g)
    end
    pows = Sym[]
    for (q, b) in groups
        j = 0
        if b.kind == NUMBER && abs(b.num) != N1 && c.kind == NUMBER
            j = something(findfirst(t -> t != 0 && c == b^t, -12:12), 13) - 13
        end
        c /= b^j
        # (-1)^(-n) = (-1)^n
        b == -ONE && q.kind == MUL && q.num < N0 && (q = -q)
        push!(pows, b^(q + j))
    end
    r = c * prod(vcat(pows, other))
    _nodes(r) <= _nodes(e) ? r : e
end

function _factorial_ratios(e::Sym)
    e.kind == MUL || return e
    fs = [(b, p) for (b, p) in zip(e.args, e.exps)]
    out = Sym(e.num)
    used = falses(length(fs))
    for i in eachindex(fs)
        used[i] && continue
        b, p = fs[i]
        if b.kind == FN && b.name === :factorial && p.kind == NUMBER
            for j in eachindex(fs)
                (j == i || used[j]) && continue
                c, q = fs[j]
                (c.kind == FN && c.name === :factorial && q.kind == NUMBER && q == -p) || continue
                d = expand(b.args[1] - c.args[1])
                (d.kind == NUMBER && _isint(d.num) && 0 < abs(_toint(d.num)) <= 20) || continue
                # u!/v! with u = v + d
                m = _toint(d.num)
                r = m > 0 ? prod((c.args[1] + i for i in 1:m); init = ONE) : 1 / prod((b.args[1] + i for i in 1:(-m)); init = ONE)
                out *= r^p
                used[i] = used[j] = true
                break
            end
        end
        used[i] || (out *= b^p; used[i] = true)
    end
    out
end

function _binomial_form(e::Sym)
    e.kind == MUL || return e
    facts = [(b.args[1], p) for (b, p) in zip(e.args, e.exps) if b.kind == FN && b.name === :factorial]
    for (u, p) in facts
        p == ONE || continue
        for (v, q) in facts
            (q.kind == NUMBER && q.num < N0) || continue
            w = expand(u - v)
            # u!/(v! w!) with w = u - v: v! w! in the denominator (or v!^2 when w == v)
            if _eq(w, v) ? q == Sym(-2) : any(t -> _eq(t[1], w) && t[2] == -ONE, facts) && q == -ONE
                r = _binomial_form(cancel(e * factorial(v) * factorial(w) / factorial(u)))
                return r * binomial(u, _nodes(w) < _nodes(v) ? w : v)
            end
        end
    end
    e
end

# the closed form against direct sums at several n (parameters given values)
function _zeil_verified(S::Sym, f::Sym, k::Sym, n::Sym, a::Sym, b::Sym, n0::Int)
    # parameters at two sets of integer values (factorials of them stay exact)
    others = [s for s in _named(f) if !_eq(s, k) && !_eq(s, n)]
    for t in (5, 11)
        vals = [s => Sym(t + 3i) for (i, s) in enumerate(others)]
        for m in n0:(n0 + 6)
            d = _direct(f, k, n, a, b, m)
            d === nothing && continue
            d = isempty(vals) ? d : subs(d, vals...)
            v = subs(S, n => m)
            v = isempty(vals) ? v : subs(v, vals...)
            x, y = (try (_evalc(v), _evalc(d)) catch; return false end)
            isapprox(x, y; rtol = 1e-10, atol = 1e-12) || return false
        end
        isempty(others) && break
    end
    true
end

## rational functions: partial fractions and polygamma

# Σ 1/(k + α)^m over k = a..b, for each partial fraction c/(k + α)^m; nothing when a
# partial fraction is not of that form (an irreducible quadratic denominator)
function _rational_sum(f::Sym, k::Sym, a::Sym, b)
    n, d = numer_denom(f)
    (_ispolyin(n, [k]) && _ispolyin(d, [k]) && has(d, k)) || return nothing
    pf = try apart(f, k) catch; return nothing end
    terms = pf.kind == ADD ? [Sym(c) * t for (c, t) in zip(pf.coefs, pf.args)] : [pf]
    pf.kind == ADD && !iszero(pf.num) && push!(terms, Sym(pf.num))
    # a polynomial part summed on its own
    poly = _add_many(t for t in terms if _ispolyin(t, [k]))
    if !iszero(poly)
        rest = _add_many(t for t in terms if !_ispolyin(t, [k]))
        rp = _summation(poly, k, a, b, nothing)
        rr = _rational_sum(rest, k, a, b)
        (rr === nothing || rp isa Cases) && return nothing
        return rr isa Sym && rp isa Sym ? _harmonic_combine(rp + rr) : nothing
    end
    parts = Tuple{Sym,Sym,Int}[]                                 # (c, α, m)
    for t in terms
        c, u = _split_kfree(t, k)
        (u.kind == POW && u.args[2].kind == NUMBER && _isint(u.args[2].num) && u.args[2].num < N0) || return nothing
        m = -_toint(u.args[2].num)
        lin = u.args[1]
        _xdeg(lin, k) == 1 || return nothing
        c1 = _coef(lin, k, 1)
        α = expand(_coef(lin, k, 0) / c1)
        push!(parts, (c / c1^m, α, m))
    end
    b === nothing && return _rational_infinite(parts, a)
    _harmonic_combine(_add_many(_hurwitz_finite(c, α, m, a, b) for (c, α, m) in parts))
end

# harmonic numbers whose arguments differ by integers written with the smallest one:
# H(n + 2) = H(n + 1) + 1/(n + 2)
function _harmonic_combine(e::Sym)
    hs = Sym[]
    _collect_fn!(hs, e, :harmonic)
    length(hs) <= 1 && return _present(e)
    groups = Dict{Tuple{Sym,Sym},Vector{Tuple{Int,Sym}}}()        # (base, order)
    for h in hs
        u = h.args[1]
        c = u.kind == ADD ? Sym(u.num) : ZERO
        (c.kind == NUMBER && _isint(c.num)) || continue
        push!(get!(groups, (expand(u - c), _horder(h)), Tuple{Int,Sym}[]), (_toint(c.num), h))
    end
    subsd = Pair{Sym,Sym}[]
    for ((base, m), list) in groups
        c0 = minimum(first, list)
        h0 = harmonic(base + c0, m)
        for (c, h) in list
            c == c0 && continue
            push!(subsd, h => h0 + _add_many(1 / (base + j)^m for j in (c0 + 1):c))
        end
    end
    _present(isempty(subsd) ? e : subs(e, subsd...))
end
function _collect_fn!(out, e::Sym, name::Symbol)
    e.kind == FN && e.name === name && !any(h -> _eq(h, e), out) && push!(out, e)
    foreach(a -> _collect_fn!(out, a, name), e.args)
    foreach(a -> _collect_fn!(out, a, name), e.exps)
    out
end

# Σ_{k=a}^{b} c/(k + α)^m: harmonic numbers (of order m) for integer shifts, digamma and
# polygamma else
function _hurwitz_finite(c::Sym, α::Sym, m::Int, a::Sym, b::Sym)
    lo, hi = a + α, b + α + 1
    α.kind == NUMBER && _isint(α.num) && return c * (harmonic(b + α, m) - harmonic(a + α - 1, m))
    if m == 1
        return c * (apply_function(:digamma, Sym[hi]) - apply_function(:digamma, Sym[lo]))
    end
    s = (-1)^m / Sym(SymNumber(factorial(big(m - 1))))
    c * s * (apply_function(:polygamma, Sym[Sym(m - 1), lo]) - apply_function(:polygamma, Sym[Sym(m - 1), hi]))
end

# Σ_{k=a}^{∞} Σ c/(k + α)^m: the m = 1 terms converge together only when their
# coefficients add to 0; ζ(m, a + α) for m ≥ 2, exact for an integer shift and even m
function _rational_infinite(parts, a::Sym)
    ones = [(c, α) for (c, α, m) in parts if m == 1]
    total = _add_many(c for (c, _) in ones)
    if !isempty(ones) && !iszero(expand(total))
        total.kind == NUMBER && return Sym(total.num > N0 ? Inf : -Inf)
        return nothing
    end
    out = ZERO
    for (c, α) in ones                      # Σ c_i (ψ(∞) - ψ(a + α_i)), the ψ(∞) cancelling
        lo = expand(a + α)
        out += lo.kind == NUMBER && _isint(lo.num) && lo.num >= N1 ? -c * harmonic(lo - 1) :
            -c * apply_function(:digamma, Sym[lo])
    end
    for (c, α, m) in parts
        m == 1 && continue
        out += c * _hurwitz_zeta(m, expand(a + α))
    end
    simplify(out)
end

# ζ(m, q) = Σ_{j≥0} 1/(j + q)^m; for an integer q ≥ 1, ζ(m) - Σ_{j<q} j^-m, with
# ζ(2n) = (-1)^(n+1) B_{2n} (2π)^(2n) / (2 (2n)!)
function _hurwitz_zeta(m::Int, q::Sym)
    if q.kind == NUMBER && _isint(q.num) && q.num >= N1
        z = iseven(m) ? Sym(SymNumber((-1)^(m ÷ 2 + 1) * _bernoulli(m) * big(2)^m // (2 * factorial(big(m))))) * PI^m :
            apply_function(:zeta, Sym[Sym(m)])
        return z - _add_many(Sym(SymNumber(1 // big(j)^m)) for j in 1:(_toint(q.num) - 1))
    end
    apply_function(:zeta, Sym[Sym(m), q])
end

## harmonic numbers

"""
    SymJ.harmonic(n)
    SymJ.harmonic(n, m)

The `n`-th harmonic number `1 + 1/2 + ... + 1/n`, and the generalized one of order `m`,
`1 + 1/2^m + ... + 1/n^m` (`harmonic(n, 1)` is `harmonic(n)`): exact for an integer `n`,
otherwise left as it is. Not exported.
"""
harmonic(n) = apply_function(:harmonic, Sym[Sym(n)])
harmonic(n, m) = apply_function(:harmonic, Sym[Sym(n), Sym(m)])
_harmonic_float(n, m = 1) = isinteger(n) && n >= 0 ? sum(1 / Float64(i)^m for i in 1:Int(n); init = 0.0) :
    throw(DomainError(n, "harmonic: an integer argument for numeric evaluation"))
register_function!(:harmonic; julia = _harmonic_float,
    eval = args -> begin
        length(args) in (1, 2) || return nothing
        n = args[1]
        m = length(args) == 2 ? args[2] : ONE
        length(args) == 2 && isone(m) && return harmonic(n)
        (m.kind == NUMBER && _isint(m.num) && m.num >= N1) || return nothing
        (n.kind == NUMBER && _isint(n.num) && n.num >= N0 && _toint(n.num) <= 10_000) || return nothing
        p = _toint(m.num)
        Sym(SymNumber(sum(1 // big(i)^p for i in 1:_toint(n.num); init = big(0) // 1)))
    end)
# the order of a harmonic number, 1 for harmonic(n)
_horder(h::Sym) = length(h.args) == 2 ? h.args[2] : ONE
