# Linear recurrences with polynomial coefficients: their hypergeometric solutions, by
# Petkovšek's algorithm Hyper, and d'Alembertian ones by reduction of order.
#
# For Σ_{i=0}^{d} p_i(n) y(n+i) = 0, a hypergeometric solution has
#     y(n+1)/y(n) = Z A(n)/B(n) C(n+1)/C(n)
# with A a monic factor of p_0(n), B one of p_d(n - d + 1), and C a polynomial. For each
# pair (A, B), Z is a nonzero root of the leading coefficients of
#     P_i(n) = p_i(n) Π_{j<i} A(n+j) Π_{j≥i}^{d-1} B(n+j),
# and C a polynomial solution of Σ_i Z^i P_i(n) C(n+i) = 0, its degree bounded by the
# indicial polynomial of the recurrence written in differences. A and B range over the
# factors over the rationals (irreducible factors of higher degree, and so solutions whose
# A or B needs algebraic roots, are not split).

const _HYPER_MAX_PAIRS = 400
const _HYPER_MAX_DEGREE = 40

# the hypergeometric solutions, as (ρ, C): y(n) = C(n) Π ρ(j), ρ = Z A/B; nothing when the
# search is too large
function _hyper(p::Vector{Sym}, n::Sym; algebraic::Bool = false)
    while length(p) > 1 && iszero(p[end])
        pop!(p)
    end
    while length(p) > 1 && iszero(p[1])
        popfirst!(p)
    end
    d = length(p) - 1
    d >= 1 || return Tuple{Sym,Sym}[]
    As = _monic_divisors(p[1], n; algebraic)
    Bs = _monic_divisors(expand(subs(p[end], n => n - d + 1)), n; algebraic)
    (As === nothing || Bs === nothing || length(As) * length(Bs) > _HYPER_MAX_PAIRS) && return nothing
    Z = Sym(Symbol("_hyper_Z"))
    out = Tuple{Sym,Sym}[]
    for A in As, B in Bs
        P = Sym[]
        for i in 0:d
            t = p[i+1]
            for j in 0:(i - 1)
                t *= subs(A, n => n + j)
            end
            for j in i:(d - 1)
                t *= subs(B, n => n + j)
            end
            push!(P, expand(t))
        end
        degs = [iszero(q) ? -1 : _xdeg(q, n) for q in P]
        m = maximum(degs)
        zpoly = _add_many(_coef(P[i+1], n, m) * Z^i for i in 0:d if degs[i+1] == m)
        has(zpoly, Z) || continue
        zs = try solve(expand(zpoly) ~ 0, Z) catch; continue end
        for z in unique(zs)
            iszero(z) && continue
            Q = Sym[expand(z^i * P[i+1]) for i in 0:d]
            for C in _poly_solutions(Q, n)
                ρ = cancel(z * A / B)
                r = cancel(ρ * subs(C, n => n + 1) / C)
                any(((ρ2, C2),) -> iszero(cancel(r - ρ2 * subs(C2, n => n + 1) / C2)), out) && continue
                push!(out, (ρ, C))
            end
        end
    end
    out
end

# the monic divisors of a polynomial in n over the rationals (its content dropped); with
# `algebraic`, its irreducible factors of higher degree split into linear ones over their
# roots (radicals), when solve finds them all
function _monic_divisors(p::Sym, n::Sym; algebraic::Bool = false)
    has(p, n) || return [ONE]
    f = factor(p)
    fs = f.kind == MUL ? collect(zip(f.args, f.exps)) : f.kind == POW ? [(f.args[1], f.args[2])] : [(f, ONE)]
    facs = Tuple{Sym,Int}[]
    for (b, e) in fs
        has(b, n) || continue
        (e.kind == NUMBER && _isint(e.num)) || return nothing
        lc = _coef(b, n, _xdeg(b, n))
        if algebraic && _xdeg(b, n) >= 2
            rs = try solve(b ~ 0, n) catch; Sym[] end
            if length(rs) == _xdeg(b, n)
                append!(facs, ((expand(n - r), _toint(e.num)) for r in rs))
                continue
            end
        end
        push!(facs, (expand(b / lc), _toint(e.num)))
    end
    divs = [ONE]
    for (b, e) in facs
        divs = [d * b^j for d in divs for j in 0:e]
    end
    [expand(d) for d in divs]
end

# a basis of the polynomial solutions C of Σ_i Q_i(n) C(n+i) = 0
function _poly_solutions(Q::Vector{Sym}, n::Sym)
    d = length(Q) - 1
    # in differences: E^i = Σ_j binomial(i, j) Δ^j, so r_j = Σ_i binomial(i, j) Q_i
    r = [expand(_add_many(binomial(i, j) * Q[i+1] for i in j:d)) for j in 0:d]
    nz = [j for j in 0:d if !iszero(r[j+1])]
    isempty(nz) && return Sym[]
    b = maximum(_xdeg(r[j+1], n) - j for j in nz)
    N = Sym(Symbol("_hyper_N"))
    ind = _add_many(_coef(r[j+1], n, _xdeg(r[j+1], n)) * prod((N - t for t in 0:(j - 1)); init = ONE)
                    for j in nz if _xdeg(r[j+1], n) - j == b)
    roots = try solve(expand(ind) ~ 0, N) catch; return Sym[] end
    ints = [_toint(v.num) for v in roots if v.kind == NUMBER && _isint(v.num) && v.num >= N0]
    isempty(ints) && return Sym[]
    D = maximum(ints)
    D <= _HYPER_MAX_DEGREE || return Sym[]
    # L(n^k) for k = 0..D, their coefficients as the columns of a matrix
    cols = [expand(_add_many(Q[i+1] * (n + i)^k for i in 0:d)) for k in 0:D]
    rows = maximum(c -> iszero(c) ? 0 : _xdeg(c, n), cols) + 1
    M = Sym[_coef(cols[k+1], n, j) for j in 0:(rows - 1), k in 0:D]
    B = LinearAlgebra.nullspace(M)
    [expand(_add_many(B[k+1, c] * n^k for k in 0:D)) for c in 1:size(B, 2)]
end

# Σ_i p_i(n) S(n+i) = 0 with S(m) given by `value(m)`: S as a combination of the
# hypergeometric solutions, its coefficients from initial values, checked at further n;
# nothing when the solutions do not span S
# the solution through Hyper's solutions over the rationals, else over the algebraic
# roots of the coefficients' factors (A = n - sqrt(2) dividing n^2 - 2)
function _recurrence_solve(p::Vector{Sym}, n::Sym, value)
    r = _recurrence_solve(p, n, value, false)
    r === nothing || return r
    any(c -> _has_irreducible_nonlinear(c, n), (p[1], p[end])) || return nothing
    _recurrence_solve(p, n, value, true)
end
function _has_irreducible_nonlinear(c::Sym, n::Sym)
    has(c, n) || return false
    f = factor(c)
    fs = f.kind == MUL ? f.args : f.kind == POW ? [f.args[1]] : [f]
    any(b -> has(b, n) && _xdeg(b, n) >= 2, fs)
end
function _recurrence_solve(p::Vector{Sym}, n::Sym, value, algebraic::Bool)
    sols = _hyper(copy(p), n; algebraic)
    (sols === nothing || isempty(sols)) && return nothing
    # start past the integer zeros and poles of the ratios and of the leading coefficient
    bad = Int[]
    for (ρ, _) in sols, q in numer_denom(ρ)
        _int_roots!(bad, q, n)
    end
    _int_roots!(bad, p[end], n)
    n0 = max(0, (isempty(bad) ? -1 : maximum(bad)) + 1)
    hs = Sym[_present_product(C * _hyper_product(ρ, n, n0)) for (ρ, C) in sols]
    r = length(hs)
    cs = [Sym(Symbol("_hyper_c", j)) for j in 1:r]
    S = _add_many(cs[j] * hs[j] for j in 1:r)
    eqs = Relation[]
    for m in n0:(n0 + r + 12)
        length(eqs) == r + 2 && break
        v = value(m)
        # (the solutions' values simplified: Γ of algebraic shifts, factorial(√2 + 2)/factorial(√2 - 1))
        v === nothing || push!(eqs, _add_many(cs[j] * simplify(subs(hs[j], n => m)) for j in 1:r) ~ v)
    end
    length(eqs) == r + 2 || return nothing
    sol = try solve(eqs, cs) catch; return nothing end
    (sol isa AbstractVector && !isempty(sol)) || return nothing
    s = first(sol)
    all(c -> haskey(s, c), cs) || return nothing
    simplify(subs(S, s...)), n0
end
function _int_roots!(out, q::Sym, n::Sym)
    has(q, n) || return out
    for v in (try solve(q ~ 0, n) catch; Sym[] end)
        v.kind == NUMBER && _isint(v.num) && push!(out, _toint(v.num))
    end
    out
end

## d'Alembertian solutions: reduction of order

# Σ p_i y(n+i) = g through iterated sums of hypergeometric terms (Abramov–Petkovšek): for
# a hypergeometric solution h of the homogeneous recurrence (Hyper), y = h z and
# z(n+i) = z(n) + Σ_{j<i} u(n+j), u = Δz, give
#     Σ_j s_j u(n+j) = g/h,   s_j = Σ_{i>j} p_i h(n+i)/h(n),
# one order lower (Σ_i p_i h(n+i)/h(n) = 0 drops z(n)). At order 0, u = g/s_0. Each level
# adds the solution h Σ_{m=n0}^{n-1} u(m), and the particular solution comes up the same
# way; the right side g need not be hypergeometric. The sums are left unevaluated for
# fitting the constants (added up directly at numbers), then summed. Nothing when a level has no
# hypergeometric solution.

# (particular, basis) from n0 on, or an Int: the start the solution needs instead
function _dalembert(p::Vector{Sym}, g::Sym, n::Sym, n0::Int, idx::Vector{Sym})
    d = length(p) - 1
    if d == 0
        bad = _int_roots!(Int[], p[1], n)
        any(>=(n0), bad) && return maximum(bad) + 1
        return (iszero(g) ? ZERO : _dal_present(g / p[1]), Sym[])
    end
    isempty(idx) && return nothing
    sols = _hyper(copy(p), n)
    (sols === nothing || isempty(sols)) && return nothing
    ρ, C = sols[argmin([_nodes(ρ) + _nodes(C) for (ρ, C) in sols])]
    # past the zeros and poles of h's ratio and the zeros of the leading coefficient
    bad = Int[]
    for q in (numer_denom(ρ)..., C, p[end])
        _int_roots!(bad, q, n)
    end
    any(>=(n0), bad) && return maximum(bad) + 1
    h = _present_product(C * _hyper_product(ρ, n, n0))
    r = cancel(ρ * subs(C, n => n + 1) / C)                # h(n+1)/h(n)
    q = Sym[p[1]]
    R = ONE
    for i in 1:d
        R = cancel(R * subs(r, n => n + i - 1))
        push!(q, cancel(p[i+1] * R))
    end
    D = reduce((a, c) -> (b = numer_denom(c)[2]; cancel(a * b / gcd(a, b))), q; init = ONE)
    s = Sym[expand(cancel(D * _add_many(q[i+1] for i in (j + 1):d))) for j in 0:(d - 1)]
    g2 = iszero(g) ? ZERO : _dal_present(D * g / h)
    sub = _dalembert(s, g2, n, n0, idx[2:end])
    sub === nothing && return nothing
    sub isa Int && return sub
    up, us = sub
    k = idx[1]
    Σn(u) = _summation_unevaluated(subs(u, n => k), k, Sym(n0), n - 1)   # (closed after the fit)
    (iszero(up) ? ZERO : h * Σn(up), Sym[h; [h * Σn(u) for u in us]])
end

# a summand at its smallest: cancelled, or factored with powers merged (2^(-k-1)/(k+1))
function _dal_present(e::Sym)
    c = _present_product(cancel(_merge_powers_terms(expand(e))))   # (-1)^(k+1) + (-1)^k = 0
    f = try _merge_powers(factor(c)) catch; c end
    _nodes(f) < _nodes(c) ? f : c
end

# Σ p_i y(n+i) = g with y(m) given by `value(m)`, through its d'Alembertian solutions, the
# constants from the values: (S, n0, S with its sums unevaluated); nothing when the
# solutions are not all d'Alembertian.
# `avoid`: names the summation indices must not take
function _dalembert_solve(p::Vector{Sym}, g::Sym, n::Sym, value, avoid)
    used = Set{Symbol}(Symbol(string(s)) for s in avoid)
    idx = Sym[Sym(s) for s in (:k, :j, :i, :m, :l, :r) if !(s in used)]
    n0 = 0
    for _ in 1:4
        res = _dalembert(p, g, n, n0, idx)
        res === nothing && return nothing
        res isa Int && (n0 = res; continue)
        yp, B = res
        r = length(B)
        cs = [Sym(Symbol("_dal_c", j)) for j in 1:r]
        eqs = Relation[]
        for m in n0:(n0 + r + 12)
            length(eqs) == r + 2 && break
            v = value(m)
            v === nothing && continue
            push!(eqs, _add_many(cs[j] * simplify(subs(B[j], n => m)) for j in 1:r) ~ v - simplify(subs(yp, n => m)))
        end
        length(eqs) == r + 2 || return nothing
        sol = try solve(eqs, cs) catch; return nothing end
        (sol isa AbstractVector && !isempty(sol)) || return nothing
        c = first(sol)
        all(x -> haskey(c, x), cs) || return nothing
        raw = yp + _add_many(c[cs[j]] * B[j] for j in 1:r)
        return simplify(_dal_close(raw)), n0, raw
    end
    nothing
end

# the unevaluated sums in e summed (closed forms, or Sum still), innermost first
function _dal_close(e::Sym)
    isempty(e.args) && return e
    e = _rebuild_node(e, Sym[_dal_close(a) for a in e.args])
    e.kind == FN && e.name === :Sum && length(e.args) == 4 ? summation(e.args...) : e
end

## rsolve

"""
    SymJ.rsolve(eq, y(n), inits...)

The solution of a linear recurrence with polynomial coefficients in `n`, such as
`y(n+2) ~ y(n+1) + y(n)`, as a combination of hypergeometric terms (Petkovšek's
algorithm Hyper). `inits` are initial values, `y(0) => 0, y(1) => 1`; without them the
answer is in terms of `y(0), …, y(d-1)` for a recurrence of order `d`. The right side
may also be a sum of hypergeometric terms in `n` (polynomials, `2^n`, `n!`, `(-1)^n`).
Not exported.

When the hypergeometric solutions are too few, the order is reduced through one of them
(d'Alembertian solutions): the answer is then in iterated sums of hypergeometric terms,
such as `harmonic(n)` or an unevaluated `Sum`, and the right side may be any function of
`n`. An error when the solution is not d'Alembertian (or Hyper's search is too large).

```julia
@syms n y()
SymJ.rsolve(y(n+2) ~ y(n+1) + y(n), y(n), y(0) => 0, y(1) => 1)   # Binet's formula
SymJ.rsolve(y(n+1) ~ 2y(n) + 1, y(n))                              # 2^n (y(0) + 1) - 1
SymJ.rsolve(y(n+1) - y(n) ~ 1/(n+1), y(n), y(0) => 0)              # harmonic(n)
```
"""
function rsolve(eq, yn::Sym, inits::Pair...)
    (yn.kind == FN && length(yn.args) == 1 && _isnamed(yn.args[1])) ||
        throw(ArgumentError("rsolve: $yn is not a function of a symbol, such as y(n)"))
    yname, n = yn.name, yn.args[1]
    e = eq isa Relation ? expand(eq.lhs - eq.rhs) : throw(ArgumentError("rsolve: $eq is not an equation"))
    # the coefficient of each y(n + s), and the terms free of y
    coefs = Dict{Int,Sym}()
    g = ZERO
    terms = e.kind == ADD ? [Sym(c) * t for (c, t) in zip(e.coefs, e.args)] : [e]
    e.kind == ADD && !iszero(e.num) && push!(terms, Sym(e.num))
    for t in terms
        ys = _yfactors(t, yname)
        if isempty(ys)
            g -= t                                 # Σ p_i y(n+i) = g
            continue
        end
        length(ys) == 1 || throw(ArgumentError("rsolve: the recurrence is not linear in $yname"))
        yt = only(ys)
        s = expand(yt.args[1] - n)
        (s.kind == NUMBER && _isint(s.num)) || throw(ArgumentError("rsolve: $yt is not $yname(n + integer)"))
        c = cancel(t / yt)
        _yfactors(c, yname) == Sym[] || throw(ArgumentError("rsolve: the recurrence is not linear in $yname"))
        coefs[_toint(s.num)] = get(coefs, _toint(s.num), ZERO) + c
    end
    isempty(coefs) && throw(ArgumentError("rsolve: no $yname in the equation"))
    lo, hi = extrema(keys(coefs))
    # y(n + lo) ... as y(n) ...: shift n by -lo
    p = Sym[expand(subs(get(coefs, i, ZERO), n => n - lo)) for i in lo:hi]
    g = expand(subs(g, n => n - lo))
    # denominators cleared: polynomial coefficients
    den = reduce((a, c) -> (q = numer_denom(c)[2]; cancel(a * q / gcd(a, q))), vcat(p, g); init = ONE)
    p = [expand(cancel(c * den)) for c in p]
    g = expand(cancel(g * den))
    all(c -> _ispolyin(c, [n]), p) || throw(ArgumentError("rsolve: the coefficients are not polynomials in $n"))
    d = length(p) - 1
    d >= 1 || throw(ArgumentError("rsolve: $eq does not relate values at different n"))
    # the initial values: given, or y(0), ..., y(d-1) themselves
    init = Dict{Int,Sym}()
    for (lhs, v) in inits
        (lhs isa Sym && lhs.kind == FN && lhs.name === yname && length(lhs.args) == 1 &&
         lhs.args[1].kind == NUMBER && _isint(lhs.args[1].num)) || throw(ArgumentError("rsolve: $lhs is not $yname(integer)"))
        init[_toint(lhs.args[1].num)] = Sym(v)
    end
    if isempty(init)
        for i in 0:(d - 1)
            init[i] = apply_function(yname, Sym[Sym(i)])
        end
    end
    length(init) == d || throw(ArgumentError("rsolve: a recurrence of order $d takes $d initial values"))
    # values by the recurrence itself
    vals = Dict{Int,Sym}(i => v for (i, v) in init)
    m0 = minimum(keys(init))
    function value(m)
        haskey(vals, m) && return vals[m]
        m < m0 + d && return nothing
        prev = [value(m - d + i) for i in 0:(d - 1)]
        any(isnothing, prev) && return nothing
        pd = subs(p[end], n => m - d)
        iszero(pd) && return nothing
        vals[m] = expand((subs(g, n => m - d) - _add_many(subs(p[i+1], n => m - d) * prev[i+1] for i in 0:(d - 1))) / pd)
    end
    # hypergeometric solutions (of the recurrence made homogeneous), else d'Alembertian ones
    hom = _annihilate(p, g, n)
    sol = hom === nothing ? nothing : _recurrence_solve(hom, n, value)
    sol === nothing && (sol = _dalembert_solve(p, g, n, value, _named(e)))
    sol === nothing && throw(ArgumentError("rsolve: the solution is not d'Alembertian (iterated sums of hypergeometric terms)"))
    S, n0 = sol
    Sv = length(sol) == 3 ? sol[3] : S             # (sums unevaluated: added up at numbers)
    # checked against the recurrence's own values
    for m in n0:(n0 + d + 4)
        v = value(m)
        v === nothing && continue
        iszero(expand(simplify(subs(Sv, n => m) - v))) ||
            throw(ArgumentError("rsolve: the solution found does not satisfy the recurrence"))
    end
    _merge_powers_terms(S)                         # (y's own index: the shift was only of n)
end

# 8·2^(n - 3) as 2^n, term by term
function _merge_powers_terms(e::Sym)
    e.kind == ADD || return _merge_powers(e)
    r = Sym(e.num) + _add_many(_merge_powers(Sym(c) * t) for (c, t) in zip(e.coefs, e.args))
    _nodes(r) <= _nodes(e) ? r : e
end

# Σ p_i y(n+i) = g, g a sum of hypergeometric terms, as a homogeneous recurrence that
# y satisfies: g's terms grouped into classes (two terms are alike when their ratio is
# rational in n; a class sums to one hypergeometric term), and each class annihilated in
# turn by L ↦ ρ_den(n) L(n+1) - ρ_num(n) L(n) for its ratio ρ = ρ_num/ρ_den, which maps
# the other classes to rational multiples of themselves. One order higher per class;
# nothing when g is not such a sum.
function _annihilate(p::Vector{Sym}, g::Sym, n::Sym)
    hom = p
    classes = _hyper_classes(g, n)
    classes === nothing && return nothing
    while !isempty(classes)
        h = popfirst!(classes)
        ρ = _hyper_ratio(h, n)
        ρ === nothing && return nothing
        rn, rd = numer_denom(ρ)
        d = length(hom) - 1
        new = [ZERO for _ in 0:(d + 1)]
        for i in 0:d
            new[i+2] += expand(rd * subs(hom[i+1], n => n + 1))
            new[i+1] -= expand(rn * hom[i+1])
        end
        hom = new
        # the other classes under the same operator: c ↦ c (ρ_den ρ_c - ρ_num), ρ_c = c(n+1)/c(n)
        new_classes = Sym[]
        for c in classes
            ρc = _hyper_ratio(c, n)
            ρc === nothing && return nothing
            q = cancel(rd * ρc - rn)
            iszero(q) || push!(new_classes, q * c)
        end
        classes = new_classes
    end
    # common factors of the coefficients removed
    gc = reduce(gcd, filter(!iszero, hom))
    has(gc, n) && (hom = [cancel(c / gc) for c in hom])
    [expand(c) for c in hom]
end

# the terms of g grouped by ratio (rational in n between members of a class), each class
# summed; nothing when a term is not hypergeometric
function _hyper_classes(g::Sym, n::Sym)
    iszero(g) && return Sym[]
    g = expand(g)
    terms = g.kind == ADD ? [Sym(c) * t for (c, t) in zip(g.coefs, g.args)] : [g]
    g.kind == ADD && !iszero(g.num) && push!(terms, Sym(g.num))
    classes = Vector{Sym}[]
    for t in terms
        _hyper_ratio(t, n) === nothing && return nothing
        i = findfirst(cl -> _alike(cl[1], t, n), classes)
        i === nothing ? push!(classes, [t]) : push!(classes[i], t)
    end
    [simplify_gamma(_add_many(cl)) for cl in classes]
end
function _alike(u::Sym, v::Sym, n::Sym)
    q = simplify_gamma(cancel(u / v))
    num, den = numer_denom(q)
    _ispolyin(num, [n]) && _ispolyin(den, [n])
end

function _yfactors(t::Sym, yname::Symbol)
    out = Sym[]
    _collect_fn!(out, t, yname)
    out
end
