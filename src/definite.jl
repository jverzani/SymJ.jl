# Definite integrals from an antiderivative F of f (the RischIntegrate extension supplies
# F). ∫ₐᵇ f dx is split at the singular points of f inside (a, b) (zeros of denominators
# and of logarithms' arguments, poles of tan, sec, cot, csc); on each piece the ends are
# F's values, or its one-sided limits at a singular point or at ±∞. Splitting also
# corrects a jump of F at such a point. A piece whose limit is infinite makes the integral
# ±∞ when all such pieces agree in sign, and an ArgumentError otherwise, as does a limit
# that does not exist. Decisions on parameters (is the singular point inside? the sign of
# a in a limit) go through _decide, so the answer is a Cases when the assumptions leave
# them open.
#
# With principal_value = true, a singular point s inside is crossed symmetrically, the
# Cauchy principal value lim_{ε→0⁺} (∫ₐ^{s-ε} + ∫_{s+ε}^b) f: its two one-sided values
# F(s⁻) - F(s⁺) become the one limit of F(s - ε) - F(s + ε); over (-∞, ∞), F(∞) - F(-∞)
# becomes the limit of F(R) - F(-R). The single limit is taken only when the one-sided
# values are not both finite, so an integral that converges keeps its value. A singular
# end is not crossed and is treated as before.

# ∫ f dx over (x, lo, hi); antider(f, x, A) is an antiderivative under the assumptions A
function _definite_integral(f, xab::Tuple, assume, antider; principal_value::Bool = false)
    length(xab) == 3 || throw(ArgumentError("integrate: give the limits as (x, a, b)"))
    x, lo, hi = Sym(xab[1]), _endpoint(xab[2]), _endpoint(xab[3])
    _isnamed(x) || throw(ArgumentError("integrate: $x is not a symbol"))
    A0 = _assumptions(assume)
    A0 === nothing && (A0 = Assumptions())
    # the other symbols are real, as for limit, and x too: conj(x) is x
    ps = Set{Sym}()
    for e in (Sym(f), lo, hi)
        _hassymbol(e) && union!(ps, free_symbols(e))
    end
    delete!(ps, x)
    f, lo, hi = (_realsyms(e, [x; collect(ps)]) for e in (Sym(f), lo, hi))
    isempty(ps) || (A0 = _with_fact(A0, FactAnd(Any[p ∈ ℝ for p in _sort_by_name!(collect(ps))])))
    v = _by_cases(A0) do A
        sub = collect(A.subst)
        g, l, h = isempty(sub) ? (Sym(f), lo, hi) : (subs(Sym(f), sub...), subs(lo, sub...), subs(hi, sub...))
        # failing that, once more with f simplified (cos² + sin² in a denominator, a 0)
        try
            _definite(g, x, l, h, antider, A, principal_value)
        catch err
            err isa InterruptException && rethrow()
            s = try simplify(g) catch; rethrow(err) end
            _eq(s, g) && rethrow()
            _definite(s, x, l, h, antider, A, principal_value)
        end
    end
    _tidy_constant(v)
end

# a value without symbols simplified, where that makes it smaller (radicals from the
# ends, as 12√2/(160 - 112√2) - …, to 1 - √2/2)
function _tidy_constant(v)
    (v isa Sym && !_hassymbol(v) && v.kind != NUMBER) || return v
    s = try simplify(v) catch; return v end
    _nodes(s) < _nodes(v) ? s : v
end

# an iterated integral whose inner integral has cases: over (x, a, b), x is between a and
# b, so a case that needs x outside, or at one point (x = 0), drops out. A case whose
# conditions on x hold throughout is integrated with its other conditions (on parameters)
# assumed; the answer is a Cases when such conditions are left. Conditions on x that hold
# on part of the interval only are an error.
function _definite_integral(cs::Cases, xab::Tuple, assume, antider; principal_value::Bool = false)
    length(xab) == 3 || throw(ArgumentError("integrate: give the limits as (x, a, b)"))
    x, lo, hi = Sym(xab[1]), _endpoint(xab[2]), _endpoint(xab[3])
    A0 = _assumptions(assume)
    A0 === nothing && (A0 = Assumptions())
    Ax = A0
    up = try _less(lo, hi, A0) catch err; err isa _NeedDecision || rethrow(); nothing end
    if up !== nothing
        a, b = up ? (lo, hi) : (hi, lo)
        _isinfpt(a) || (Ax = something(_with_fact(Ax, Relation(:gt, x, a)), Ax))
        _isinfpt(b) || (Ax = something(_with_fact(Ax, Relation(:lt, x, b)), Ax))
    end
    # cases on x ⋈ c, c strictly inside: the interval split at c, the pieces added
    if up !== nothing
        a, b = up ? (lo, hi) : (hi, lo)
        cuts = Sym[]
        for (conds, _) in cs, c in conds
            p = _fact_point(c, x)
            (p === nothing || any(q -> _eq(q, p), cuts)) && continue
            inside = try _less(a, p, A0) && _less(p, b, A0) catch err; err isa _NeedDecision || rethrow(); false end
            inside && _insert_sorted!(cuts, p, A0)
        end
        if !isempty(cuts)
            pts = [a; cuts; b]
            parts = [_definite_integral(cs, (x, pts[i], pts[i + 1]), A0, antider; principal_value) for i in 1:(length(pts) - 1)]
            all(v -> v isa Sym, parts) || throw(ArgumentError("integrate: the inner integral's cases split the interval of $x"))
            v = sum(parts)
            return up ? v : -v
        end
    end
    out = Pair{Vector{Any},Any}[]
    for (conds, v) in cs
        onx = [c for c in conds if _fact_has(c, x)]
        rest = Any[c for c in conds if !_fact_has(c, x)]
        any(c -> c isa Relation && c.op === :eq, onx) && continue
        hold = [_fact_holds(c, Ax) for c in onx]
        any(h -> h === false, hold) && continue
        all(h -> h === true, hold) || throw(ArgumentError(
            "integrate: the inner integral's cases ($(_condsstring(onx))) split the interval of $x"))
        B = isempty(rest) ? A0 : _with_fact(A0, rest)
        B === nothing && continue
        v isa Exception && throw(v)
        push!(out, rest => _definite_integral(v, xab, B, antider; principal_value))
    end
    isempty(out) && throw(ArgumentError("integrate: no case of the inner integral holds for $x between $lo and $hi"))
    length(out) == 1 && isempty(out[1].first) && return out[1].second
    T = typejoin(map(c -> typeof(c.second), out)...)
    Cases{T}(out)
end
# c for a condition x ⋈ c (or c ⋈ x), c free of x; else nothing
function _fact_point(c, x::Sym)
    c isa _NotFact && return _fact_point(c.f, x)
    c isa Relation || return nothing
    _eq(c.lhs, x) && !has(c.rhs, x) && return c.rhs
    _eq(c.rhs, x) && !has(c.lhs, x) && return c.lhs
    nothing
end
_fact_has(c, x::Sym) = c isa _NotFact ? _fact_has(c.f, x) : c isa Membership ? has(c.x, x) : has(c.lhs, x) || has(c.rhs, x)
_fact_holds(c, A) = c isa _NotFact ? (h = _fact_holds(c.f, A); h === nothing ? nothing : !h) : _holds(c, A)

_endpoint(v) = v isa Real && isinf(v) ? Sym(Float64(v)) : Sym(v)
_isinfpt(p::Sym) = p.kind == NUMBER && isfloat(p.num) && isinf(Float64(p.num))
_infsign(p::Sym) = Float64(p.num) > 0 ? 1 : -1

function _definite(f::Sym, x::Sym, lo::Sym, hi::Sym, antider, A::Assumptions, pv::Bool = false)
    _eq(lo, hi) && return ZERO
    has(f, x) || return _ordered_length(f, lo, hi, A)
    S = _singular_points(f, x, lo, hi)
    if isempty(S) && !_isinfpt(lo) && !_isinfpt(hi)
        F = _antiderivative(antider, f, x, A; numeric = true)
        F isa Sym || return _numeric_difference(F, f, lo, hi)
        return _sum_values([(1, _endvalue(F, x, hi, :-, A)), (-1, _endvalue(F, x, lo, :+, A))], A)
    end
    # in increasing order: ∫ₐᵇ = -∫ᵇₐ
    if _less(hi, lo, A)
        v = _definite(f, x, hi, lo, antider, A, pv)
        return v isa Sym ? -v : v
    end
    inside = Sym[]
    sing = Set{Sym}()                        # the singular points among the ends too
    for s in S
        (_eq(s, lo) || _eq(s, hi)) && (push!(sing, s); continue)
        # a parametric point is compared with a finite end, which reduces to facts about
        # the parameters (√(-a) > 0 is a < 0); with both ends infinite, its reality
        (!_hassymbol(s) || (_isinfpt(lo) && _isinfpt(hi))) && !_isreal_point(s, A) && continue
        (_isinfpt(lo) || _decide(Relation(:gt, s, lo), A)) && (_isinfpt(hi) || _decide(Relation(:lt, s, hi), A)) || continue
        any(t -> _eq(t, s), inside) || _insert_sorted!(inside, s, A)
    end
    union!(sing, inside)
    F = _antiderivative(antider, f, x, A; numeric = true)
    if !(F isa Sym)
        (_isinfpt(lo) || _isinfpt(hi)) && throw(ArgumentError(
            "integrate: the antiderivative of $f is known only numerically, which does not reach ±∞"))
        isempty(sing) || throw(ArgumentError(
            "integrate: the antiderivative of $f is known only numerically, and $f is singular on the interval"))
        return _numeric_difference(F, f, lo, hi)
    end
    terms = Tuple{Int,Sym}[]
    for s in inside                                  # F(s⁻) - F(s⁺)
        l, r = () -> _endvalue(F, x, s, :-, A, true), () -> _endvalue(F, x, s, :+, A, true)
        append!(terms, pv ? _pv_terms(l, r, () -> _pv_limit(F, x, s, A)) : [(1, l()), (-1, r())])
    end
    h, l = () -> _endvalue(F, x, hi, :-, A, hi in sing), () -> _endvalue(F, x, lo, :+, A, lo in sing)
    append!(terms, pv && _isinfpt(lo) && _isinfpt(hi) ? _pv_terms(h, l, () -> _pv_limit(F, x, nothing, A)) :
                   [(1, h()), (-1, l())])
    _sum_values(terms, A)
end

# the terms for a principal value: the one-sided values left() - right() when both are
# finite, else the single symmetric limit both()
function _pv_terms(left, right, both)
    lr = try
        (left(), right())
    catch err
        err isa ArgumentError || rethrow()
        nothing
    end
    lr !== nothing && !any(_isinfpt, lr) && return [(1, lr[1]), (-1, lr[2])]
    [(1, both())]
end

# lim_{ε→0⁺} F(s - ε) - F(s + ε), or for s = nothing lim_{R→∞} F(R) - F(-R)
function _pv_limit(F::Sym, x::Sym, s, A)
    t = _PV_VAR
    D, p, dir = s === nothing ? (subs(F, x => t) - subs(F, x => -t), Sym(Inf), :-) :
                                (subs(F, x => s - t) - subs(F, x => s + t), ZERO, :+)
    D = simplify(D; assume = A)       # the limit code may not see a hidden zero
    try
        _onesided(D, t, p, dir, A)
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("integrate: the principal value does not exist " *
                            (s === nothing ? "at ±∞" : "at $x = $s")))
    end
end
const _PV_VAR = Sym(:_definite_pv)

# ∫ c dx for c free of x
function _ordered_length(c::Sym, lo::Sym, hi::Sym, A)
    (_isinfpt(lo) || _isinfpt(hi)) || return c * (hi - lo)
    iszero(c) && return ZERO
    s = (_isinfpt(hi) ? _infsign(hi) : 0) - (_isinfpt(lo) ? _infsign(lo) : 0)
    s == 0 && throw(ArgumentError("integrate: the limits are both $(lo)"))
    _decide(Relation(:gt, c, ZERO), A) ? Sym(s * Inf) :
        _decide(Relation(:lt, c, ZERO), A) ? Sym(-s * Inf) : ZERO
end

# is the singular point s real (a numeric s by its value, else as decided)?
function _isreal_point(s::Sym, A)
    if !_hassymbol(s)
        v = try N(s) catch; return false end
        return v isa Real
    end
    r = _isreal(s, A)
    r === nothing ? _decide(Membership(s, NumberSet(_RR), false), A) : r
end

# p < q, the points finite or ±∞
function _less(p::Sym, q::Sym, A)
    _isinfpt(p) && return _infsign(p) < 0 && !(_isinfpt(q) && _infsign(q) < 0)
    _isinfpt(q) && return _infsign(q) > 0
    _decide(Relation(:lt, p, q), A)
end

function _insert_sorted!(v::Vector{Sym}, s::Sym, A)
    i = findfirst(t -> _decide(Relation(:lt, s, t), A), v)
    i === nothing ? push!(v, s) : insert!(v, i, s)
end

# simplified under the assumptions: limits must see that a/√(-a) + √(-a) is 0. With
# numeric = true, an antiderivative known only numerically (a callable, as RischIntegrate's
# root sums and ₂F₁ terms) is returned as it is.
function _antiderivative(antider, f::Sym, x::Sym, A; numeric::Bool = false)
    F = antider(f, x, A)
    numeric && !(F isa Sym) && applicable(F, 0.0) && return F
    (F isa Sym && !_hasfn(F, :Integral)) || throw(ArgumentError("integrate: no antiderivative found for $f"))
    _hassymbol(F - x) ? simplify(F; assume = A) : F
end

# F(hi) - F(lo), a float, for an antiderivative F known only numerically and continuous on
# [lo, hi] (no singular points there); the ends must be numbers
function _numeric_difference(F, f::Sym, lo::Sym, hi::Sym)
    (_hassymbol(lo) || _hassymbol(hi)) && throw(ArgumentError(
        "integrate: the antiderivative of $f is known only numerically; the limits must be numbers"))
    Sym(Float64(real(F(N(hi)) - F(N(lo)))))
end

# F at p from the side dir: its value, or its one-sided limit where that is needed
function _endvalue(F::Sym, x::Sym, p::Sym, dir::Symbol, A, singular::Bool = false)
    if !singular && !_isinfpt(p)
        v = try subs(F, x => p) catch err; err isa DivideError || rethrow(); nothing end
        # (a piecewise F at a parametric end, as ∫|x| at a/b, its case decided if A can)
        v !== nothing && _hasfn(v, :piecewise) && (v = refine(v; assume = A))
        v !== nothing && _finite_value(v) && return v
    end
    _onesided(F, x, p, dir, A)
end

function _finite_value(v::Sym)
    _hassymbol(v) && return true
    w = try N(v) catch; return false end
    isfinite(w)
end

function _onesided(F::Sym, x::Sym, p::Sym, dir::Symbol, A)
    x0 = _isinfpt(p) ? Float64(p.num) : p
    L = _limit_or(F, x, x0, dir, A)
    if L === nothing
        Gc = _neglogs(F, x, p, dir)          # log(u) with u < 0 near p: log(-u) + iπ
        if Gc !== nothing
            vs = [k => _limit_or(g, x, x0, dir, A) for (k, g) in Gc]
            if all(v -> v.second isa Sym, vs)
                if all(v -> !_isinfpt(v.second), vs)
                    L = _tidy_constant(sum(k * v for (k, v) in vs; init = ZERO))
                elseif length(vs) == 1 && isone(vs[1].first)
                    L = vs[1].second
                end
            end
        end
    end
    L === nothing && throw(ArgumentError("integrate: the integral does not converge: the antiderivative $F has no limit at $x = $p"))
    L
end

function _limit_or(F, x, x0, dir, A)
    L = try
        limit(F, x => x0; dir = dir === :+ ? :+ : :-, assume = A)
    catch err
        err isa Union{GruntzError,ArgumentError,DomainError} || rethrow()
        return nothing
    end
    L isa Cases || return L
    for (conds, v) in L                     # the case the assumptions (or a decision) give
        all(c -> c isa _NotFact ? !_decide(c.f, A) : _decide(c, A), conds) || continue
        return v isa Exception ? nothing : v
    end
    nothing
end

# F near p on the side dir with log(u), u real and negative there, written log(-u) + iπ:
# its terms grouped as k => g, k a constant (1, iπ, (iπ)², e^(iπa)) and g real near p, so
# F = Σ k g; nothing when there is no such log or the side cannot be tested numerically
function _neglogs(F::Sym, x::Sym, p::Sym, dir::Symbol)
    q = _isinfpt(p) ? (_infsign(p) > 0 ? 1e8 : -1e8) :
        (_hassymbol(p) ? nothing : (v = N(p); v isa Real ? v + (dir === :+ ? 1e-8 : -1e-8) * max(1, abs(v)) : nothing))
    q === nothing && return nothing
    changed = Ref(false)
    G = _replace_walk(F) do t
        (t.kind == FN && t.name === :log && length(t.args) == 1) || return nothing
        u = t.args[1]
        w = try N(subs(u, x => q)) catch; return nothing end
        (w isa Real && w < 0) || return nothing
        changed[] = true
        log(-u) + _MARK_IPI
    end
    changed[] || return nothing
    # G's terms grouped by a constant factor: a term linear in the mark, t₀ + c·mark, is
    # t₀ (factor 1) and c (factor iπ); one with the mark in an exponent, t ∝ e^(a·mark), is
    # t at mark = 0 with factor e^(iπa)
    G = expand(G)
    groups = Dict{Sym,Sym}()
    put!(k, t) = (groups[k] = get(groups, k, ZERO) + t)
    terms = G.kind == ADD ? [Sym(c) * t for (c, t) in zip(G.coefs, G.args)] : [G]
    G.kind == ADD && put!(ONE, Sym(G.num))
    for t in terms
        has(t, _MARK_IPI) || (put!(ONE, t); continue)
        cs = _mark_poly(t)
        if cs !== nothing
            for (j, c) in enumerate(cs)
                put!((IM * PI)^(j - 1), c)
            end
        else
            t0, c = subs(t, _MARK_IPI => ZERO), diff(t, _MARK_IPI)
            a = cancel(c / t)
            (has(a, x) || has(a, _MARK_IPI)) && return nothing
            put!(exp(a * IM * PI), t0)
        end
    end
    [k => g for (k, g) in groups if !iszero(g)]
end
const _MARK_IPI = Sym(:_definite_ipi)

# t as a polynomial in the mark (log(u)^2 gives its square): the coefficients c₀, c₁, …,
# or nothing
function _mark_poly(t::Sym)
    cs, d, f = Sym[], t, 1
    for j in 0:4
        push!(cs, subs(d, _MARK_IPI => ZERO) / f)
        d = diff(d, _MARK_IPI)
        f *= j + 1
        iszero(d) && return iszero(expand(t - sum(c * _MARK_IPI^(i - 1) for (i, c) in enumerate(cs)))) ? cs : nothing
    end
    nothing
end

# Σ cᵢ vᵢ, the vᵢ finite or ±∞
function _sum_values(terms, A::Assumptions)
    fin, infs = ZERO, Int[]
    for (c, v) in terms
        if _isinfpt(v)
            push!(infs, c * _infsign(v))
        else
            fin += c * v
        end
    end
    if !isempty(infs)
        all(==(infs[1]), infs) || throw(ArgumentError("integrate: the integral does not converge (∞ - ∞)"))
        return Sym(infs[1] * Inf)
    end
    fin = _neglog_consts(fin, A)
    _hassymbol(fin) && return fin
    s = simplify(fin)
    _nodes(s) <= _nodes(fin) ? s : fin
end

# log(u) for u negative (a number, or known from the assumptions) written log(-u) + iπ
# (the principal value), so the imaginary parts of the two ends cancel
_neglog_consts(v::Sym, A::Assumptions) = _replace_walk(v) do t
    (t.kind == FN && t.name === :log && length(t.args) == 1) || return nothing
    u = t.args[1]
    neg = u.kind == NUMBER ? u.num < N0 : _hassymbol(u) && _signs(u, A) == _NEG
    neg ? log(expand(-u)) + IM * PI : nothing
end

## singular points of the integrand

# the real points in [lo, hi] where f is singular: zeros of the bases of negative powers,
# of logarithms' arguments, of cos (tan, sec) and sin (cot, csc), of sinh (coth, csch)
function _singular_points(f::Sym, x::Sym, lo::Sym, hi::Sym)
    zs = Sym[]
    _zero_sets!(zs, f, x)
    pts = Sym[]
    for z in unique(zs)
        append!(pts, _zeros_in(z, x, lo, hi))
    end
    unique(pts)
end

function _zero_sets!(out, e::Sym, x::Sym)
    has(e, x) || return out
    k = e.kind
    if k == MUL || k == POW
        bs, ps = k == MUL ? (e.args, e.exps) : (e.args[1:1], e.args[2:2])
        for (b, p) in zip(bs, ps)
            # a negative power, or a symbolic one (x^n, singular at 0 for n < 0)
            has(b, x) && ((p.kind == NUMBER && p.num < N0) || (p.kind != NUMBER && !has(p, x))) && push!(out, b)
            _zero_sets!(out, b, x)
            _zero_sets!(out, p, x)
        end
        return out
    elseif k == FN && length(e.args) == 1
        u = e.args[1]
        e.name === :log && push!(out, u)
        e.name in (:tan, :sec) && push!(out, cos(u))
        e.name in (:cot, :csc) && push!(out, sin(u))
        e.name in (:coth, :csch) && push!(out, sinh(u))
    end
    foreach(a -> _zero_sets!(out, a, x), e.args)
    out
end

# the zeros of z in x (all of them in [lo, hi] for sin and cos of a linear argument)
function _zeros_in(z::Sym, x::Sym, lo::Sym, hi::Sym)
    if z.kind == FN && z.name in (:sin, :cos) && length(z.args) == 1
        u = z.args[1]
        k = expand(diff(u, x))
        (has(k, x) || iszero(k)) && throw(ArgumentError("integrate: cannot locate the singular points of the integrand (at the zeros of $z)"))
        (_hassymbol(lo) || _hassymbol(hi) || _hassymbol(subs(u, x => 0))) &&
            throw(ArgumentError("integrate: cannot locate the singular points of the integrand (at the zeros of $z) for symbolic limits"))
        (_isinfpt(lo) || _isinfpt(hi)) &&
            throw(ArgumentError("integrate: the integrand is singular at infinitely many points (the zeros of $z)"))
        off = z.name === :cos ? PI / 2 : ZERO               # u = off + nπ
        ua, ub = N(subs(u, x => lo)), N(subs(u, x => hi))
        n1, n2 = floor(Int, (min(ua, ub) - N(off)) / π) - 1, ceil(Int, (max(ua, ub) - N(off)) / π) + 1
        return Sym[expand(x - (u - (off + n * PI)) / k) for n in n1:n2]
    end
    rs = try solve(z ~ 0, x) catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("integrate: cannot locate the singular points of the integrand (the zeros of $z)"))
    end
    rs isa AbstractVector || throw(ArgumentError("integrate: cannot locate the singular points of the integrand (the zeros of $z)"))
    P = _trig_period(z, x)
    P === nothing && return Sym[r for r in rs]
    # periodic: solve gives one period's zeros; repeat the real ones across [lo, hi]
    vs = [(r, v) for r in rs for v in (try N(r) catch; nothing end,) if v isa Real]
    isempty(vs) && return Sym[r for r in rs]
    (_isinfpt(lo) || _isinfpt(hi)) &&
        throw(ArgumentError("integrate: the integrand is singular at infinitely many points (the zeros of $z)"))
    (_hassymbol(lo) || _hassymbol(hi)) &&
        throw(ArgumentError("integrate: cannot locate the singular points of the integrand (at the zeros of $z) for symbolic limits"))
    a, b = minmax(N(lo), N(hi))
    p = N(P)
    Sym[expand(r + n * P) for (r, v) in vs for n in (floor(Int, (a - v) / p) - 1):(ceil(Int, (b - v) / p) + 1)]
end

# 2π/k when x enters z only through sin and cos of arguments k x + c (one k), else nothing
function _trig_period(z::Sym, x::Sym)
    ks = Sym[]
    _trig_args!(ks, z, x) || return nothing
    (isempty(ks) || any(k -> k != ks[1], ks)) && return nothing
    2PI / (ks[1].num < N0 ? -ks[1] : ks[1])
end

function _trig_args!(ks, e::Sym, x::Sym)
    has(e, x) || return true
    if e.kind == FN && e.name in (:sin, :cos) && length(e.args) == 1
        k = expand(diff(e.args[1], x))
        (k.kind == NUMBER && !iszero(k)) || return false
        push!(ks, k)
        return true
    end
    e.kind == SYMBOL && return false
    all(a -> _trig_args!(ks, a, x), e.args) && (e.kind != MUL || all(q -> _trig_args!(ks, q, x), e.exps))
end

# the sign oracle for an integrator: the sign of a parameter expression under A, as 1 or
# -1, or nothing when A does not decide it
function _sign_oracle(A::Assumptions)
    function (p)
        m = _signs(Sym(p), A)
        m == _POS ? 1 : m == _NEG ? -1 : nothing
    end
end
