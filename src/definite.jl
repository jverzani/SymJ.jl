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
    # the other symbols are real, as for limit
    ps = Set{Sym}()
    for e in (Sym(f), lo, hi)
        _hassymbol(e) && union!(ps, free_symbols(e))
    end
    delete!(ps, x)
    isempty(ps) || (A0 = _with_fact(A0, FactAnd(Any[p ∈ ℝ for p in _sort_by_name!(collect(ps))])))
    _by_cases(A0) do A
        sub = collect(A.subst)
        g, l, h = isempty(sub) ? (Sym(f), lo, hi) : (subs(Sym(f), sub...), subs(lo, sub...), subs(hi, sub...))
        _definite(g, x, l, h, antider, A, principal_value)
    end
end

_endpoint(v) = v isa Real && isinf(v) ? Sym(Float64(v)) : Sym(v)
_isinfpt(p::Sym) = p.kind == NUMBER && isfloat(p.num) && isinf(Float64(p.num))
_infsign(p::Sym) = Float64(p.num) > 0 ? 1 : -1

function _definite(f::Sym, x::Sym, lo::Sym, hi::Sym, antider, A::Assumptions, pv::Bool = false)
    _eq(lo, hi) && return ZERO
    has(f, x) || return _ordered_length(f, lo, hi, A)
    S = _singular_points(f, x, lo, hi)
    if isempty(S) && !_isinfpt(lo) && !_isinfpt(hi)
        F = _antiderivative(antider, f, x, A)
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
    F = _antiderivative(antider, f, x, A)
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

# simplified under the assumptions: limits must see that a/√(-a) + √(-a) is 0
function _antiderivative(antider, f::Sym, x::Sym, A)
    F = antider(f, x, A)
    (F isa Sym && !_hasfn(F, :Integral)) || throw(ArgumentError("integrate: no antiderivative found for $f"))
    _hassymbol(F - x) ? simplify(F; assume = A) : F
end

# F at p from the side dir: its value, or its one-sided limit where that is needed
function _endvalue(F::Sym, x::Sym, p::Sym, dir::Symbol, A, singular::Bool = false)
    if !singular && !_isinfpt(p)
        v = try subs(F, x => p) catch err; err isa DivideError || rethrow(); nothing end
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
    if L === nothing && _isinfpt(p)
        G = _special_at_inf(F, x, x0, A)     # erf(u) → ±1 where u → ±∞, …
        G === nothing || (L = _limit_or(G, x, x0, dir, A))
    end
    if L === nothing
        Gc = _neglogs(F, x, p, dir)          # log(u) with u < 0 near p: log(-u) + iπ
        if Gc !== nothing
            L = _limit_or(Gc[1], x, x0, dir, A)
            L === nothing || _isinfpt(L) || (L = L + Gc[2])
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

# the values at ±∞ of special functions the limit code does not know: f(u) with u → ±∞
const _AT_INF = Dict{Symbol,Tuple{Any,Any}}(       # name => (value at +∞, value at -∞)
    :erf => (ONE, -ONE), :erfc => (ZERO, Sym(2)), :sinint => (PI / 2, -PI / 2),
    :cosint => (ZERO, nothing), :erfi => (Sym(Inf), Sym(-Inf)))
function _special_at_inf(F::Sym, x::Sym, x0, A)
    changed = Ref(false)
    G = _replace_walk(F) do t
        (t.kind == FN && length(t.args) == 1 && haskey(_AT_INF, t.name)) || return nothing
        Lu = _limit_or(t.args[1], x, x0, x0 > 0 ? :- : :+, A)
        (Lu !== nothing && _isinfpt(Lu)) || return nothing
        v = _AT_INF[t.name][_infsign(Lu) > 0 ? 1 : 2]
        v === nothing && return nothing
        changed[] = true
        v
    end
    changed[] ? G : nothing
end

# F as G + c near p on the side dir, G with log(u) written log(-u) where u, real, is
# negative there, and c the iπ those add; nothing when there is no such log or the side
# cannot be tested numerically
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
    c = expand(diff(G, _MARK_IPI))           # G is linear in the mark
    (expand(subs(G, _MARK_IPI => ZERO)), c * IM * PI)
end
const _MARK_IPI = Sym(:_definite_ipi)

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
    Sym[r for r in rs]
end

# the sign oracle for an integrator: the sign of a parameter expression under A, as 1 or
# -1, or nothing when A does not decide it
function _sign_oracle(A::Assumptions)
    function (p)
        m = _signs(Sym(p), A)
        m == _POS ? 1 : m == _NEG ? -1 : nothing
    end
end
