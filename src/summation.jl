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
#   - rational functions without one: partial fractions, Σ 1/(k+α)^m as polygamma
#     functions, exact harmonic and zeta values for integer shifts;
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
and polygamma functions (exact zeta values for even powers). A sum with no closed form
is returned unevaluated, `Sum(f, k, a, b)`. Not exported.

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
    F = _antidifference(f, k)
    F === nothing || return _from_antidifference(F, k, a, b, assume)
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
    A, B = numer_denom(ρ)
    (_ispolyin(A, [k]) && _ispolyin(B, [k])) || return nothing
    # Gosper form: ρ = A(k)/B(k) * C(k+1)/C(k) with gcd(A(k), B(k+h)) = 1 for h ≥ 0
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
    C = expand(C)
    # a(k) x(k+1) - b(k-1) x(k) = c(k) for a polynomial x
    Bm = expand(subs(B, k => k - 1))
    x = _gosper_x(expand(A), Bm, C, k)
    x === nothing && return nothing
    cancel(Bm * x / C * f)
end

# the nonnegative integer roots of R(h)
function _nonneg_int_roots(R::Sym, h::Sym)
    has(R, h) || return Int[]
    rs = try solve(R ~ 0, h) catch; return Int[] end
    sort!(unique([_toint(r.num) for r in rs if r.kind == NUMBER && _isint(r.num) && r.num >= N0]))
end

function _gosper_x(a::Sym, b::Sym, c::Sym, k::Sym)
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
    iszero(eq) && return expand(subs(xk, (v => ZERO for v in xs)...))
    P = try Poly{Sym}(eq, [k]) catch; return nothing end
    eqs = [cf ~ 0 for cf in coefficients(P)]
    sols = try solve(eqs, xs) catch; return nothing end
    isempty(sols) && return nothing
    s = first(sols)
    expand(subs(subs(xk, s...), (v => ZERO for v in xs)...))      # free coefficients 0
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

## rational functions: partial fractions and polygamma

# Σ 1/(k + α)^m over k = a..b, for each partial fraction c/(k + α)^m; nothing when a
# partial fraction is not of that form (an irreducible quadratic denominator)
function _rational_sum(f::Sym, k::Sym, a::Sym, b)
    n, d = numer_denom(f)
    (_ispolyin(n, [k]) && _ispolyin(d, [k]) && has(d, k)) || return nothing
    pf = try apart(f, k) catch; return nothing end
    terms = pf.kind == ADD ? [Sym(c) * t for (c, t) in zip(pf.coefs, pf.args)] : [pf]
    pf.kind == ADD && !iszero(pf.num) && return nothing          # a polynomial part
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
    groups = Dict{Sym,Vector{Tuple{Int,Sym}}}()
    for h in hs
        u = h.args[1]
        c = u.kind == ADD ? Sym(u.num) : ZERO
        (c.kind == NUMBER && _isint(c.num)) || continue
        push!(get!(groups, expand(u - c), Tuple{Int,Sym}[]), (_toint(c.num), h))
    end
    subsd = Pair{Sym,Sym}[]
    for (base, list) in groups
        c0 = minimum(first, list)
        h0 = harmonic(base + c0)
        for (c, h) in list
            c == c0 && continue
            push!(subsd, h => h0 + _add_many(1 / (base + j) for j in (c0 + 1):c))
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

# Σ_{k=a}^{b} c/(k + α)^m: harmonic numbers for m = 1 and integer shifts, polygamma else
function _hurwitz_finite(c::Sym, α::Sym, m::Int, a::Sym, b::Sym)
    lo, hi = a + α, b + α + 1
    if m == 1
        α.kind == NUMBER && _isint(α.num) && return c * (harmonic(b + α) - harmonic(a + α - 1))
        return c * (apply_function(:digamma, Sym[hi]) - apply_function(:digamma, Sym[lo]))
    end
    s = (-1)^m / Sym(SymNumber(factorial(big(m - 1))))
    c * s * (apply_function(:polygamma, Sym[Sym(m - 1), hi]) - apply_function(:polygamma, Sym[Sym(m - 1), lo]))
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

The `n`-th harmonic number `1 + 1/2 + ... + 1/n`: exact for an integer `n`, otherwise
left as it is (`harmonic(n)`). Not exported.
"""
harmonic(n) = apply_function(:harmonic, Sym[Sym(n)])
register_function!(:harmonic; julia = n -> (isinteger(n) && n >= 0 ? sum(1 / i for i in 1:Int(n); init = 0.0) :
        throw(DomainError(n, "harmonic: an integer argument for numeric evaluation"))),
    eval = args -> begin
        n = args[1]
        (n.kind == NUMBER && _isint(n.num) && n.num >= N0 && _toint(n.num) <= 10_000) || return nothing
        Sym(SymNumber(sum(1 // big(i) for i in 1:_toint(n.num); init = big(0) // 1)))
    end)
