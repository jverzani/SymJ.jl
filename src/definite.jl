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
    # terms with an unknown function of x are held, ∫ₐᵇ f(x) dx, the rest integrated
    rest, held = _split_opaque(Sym(f), x)
    if !isempty(held)
        H = _add_many(c * _mkfn_node(:Integral, Sym[g, x, Sym(xab[2]), Sym(xab[3])]) for (c, g) in held)
        iszero(rest) && return H
        v = _definite_integral(rest, xab, assume, antider; principal_value)
        return v isa Cases ? _map_cases(w -> w + H, v) : v + H
    end
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
        # (a case where the integrand is not defined, b = 0 in exp(-(x - a)²/b): an error
        # for that case, not for the whole; SymPy #19672)
        g, l, h = isempty(sub) ? (Sym(f), lo, hi) : try
            (subs(Sym(f), sub...), subs(lo, sub...), subs(hi, sub...))
        catch err
            _isdivzero(err) || rethrow()
            throw(ArgumentError("integrate: the integrand is not defined where $(join(("$k = $w" for (k, w) in sub), ", "))"))
        end
        # the value refined under the case's assumptions (sin(nπ) = 0 for an integer n;
        # found in the review of SymPy #24000s)
        r = _definite_degenerate(g, x, l, h, A) do g, l, h
            # failing that, once more with f simplified (cos² + sin² in a denominator, a 0);
            # failing both, a known definite integral (no antiderivative needed)
            try
                try
                    _definite(g, x, l, h, antider, A, principal_value)
                catch err
                    err isa InterruptException && rethrow()
                    s = try simplify(g) catch; rethrow(err) end
                    _eq(s, g) && rethrow()
                    _definite(s, x, l, h, antider, A, principal_value)
                end
            catch err
                err isa InterruptException && rethrow()
                v = _known_definite(g, x, l, h)
                v === nothing && (v = _bose_fermi_definite(g, x, l, h, A))
                # (a Gaussian in real form once refined: |exp(w)|² with w complex; SymPy #26102)
                if v === nothing
                    Ax = _with_fact(A, x ∈ ℝ)
                    gr = Ax === nothing ? g : _refine(g, Ax)
                    v = _gaussian_definite(gr, x, l, h, A)
                end
                v === nothing && (v = _definite_fallback(g, x, l, h, antider, A, principal_value, 1))
                v === nothing ? rethrow() : v
            end
        end
        r isa Sym && _hassymbol(r) ? _refine(r, A) : r
    end
    _tidy_constant(v)
end

# The value run(g, l, h) of the integral, with a case for each parameter factor of its
# denominator that the integrand's own denominator does not have: where it vanishes the
# formula does not hold (∫₀¹ e^{ct} dt = (e^c - 1)/c, but 1 at c = 0; SymPy #27280), and
# the integral is computed again with that value substituted
function _definite_degenerate(run, g::Sym, x::Sym, l::Sym, h::Sym, A::Assumptions)
    v = run(g, l, h)
    (v isa Sym && _hassymbol(v)) || return v
    # (denominators inside functions and radicals too: √(1/M² + 1) is not defined at
    # M = 0, where a simpler form of 4/((M + M√(1/M² + 1))² - 1) has it)
    function nested_dens!(dens, e)
        _replace_walk(e) do t
            if t.kind == POW && t.args[2].kind == NUMBER && t.args[2].num < N0 && _hassymbol(t.args[1])
                push!(dens, t.args[1])
            elseif t.kind == MUL
                for (b, q) in zip(t.args, t.exps)
                    q.kind == NUMBER && q.num < N0 && _hassymbol(b) && push!(dens, b)
                end
            end
            nothing
        end
        dens
    end
    # the integrand's own, nested ones included: exp(-(x - a)²/b) is not defined at b = 0,
    # so that is no case to integrate again (it divided by zero; SymPy #19672)
    owndens = nested_dens!(Sym[numer_denom(g)[2]], g)
    own = unique(Sym[f for d in owndens for f in _param_factors(d, (x,))])
    dens = nested_dens!(Sym[numer_denom(cancel(v))[2]], v)
    D = Sym[]
    for d in dens, f in _param_factors(d, (x,))
        any(o -> cancel(o / f).kind == NUMBER, own) || any(o -> _eq(o, f), D) || push!(D, f)
    end
    isempty(D) && return v
    ps = _degenerate_subst(D, A, (x,))
    ps === nothing && return v
    _alternatives(ps) do p
        _definite_degenerate(run, subs(g, p), x, subs(l, p), subs(h, p), A)
    end
end

## definite integrals known without an antiderivative

# ∫₀^{π/2} sinᵃx cosᵇx dx = B((a+1)/2, (b+1)/2)/2 = Γ(p)Γ(q)/(2Γ(p+q)), for rational
# a, b > -1 (tan, cot, sec, csc as powers of sin and cos), times a constant; else nothing
function _known_definite(f::Sym, x::Sym, lo::Sym, hi::Sym)
    (iszero(lo) && _eq(hi, PI / 2)) || return nothing
    c, a, b = ONE, N0, N0
    fs = f.kind == MUL ? [(Sym(f.num), ONE); [(u, q) for (u, q) in zip(f.args, f.exps)]] : [(f, ONE)]
    for (u, q) in fs
        if !has(u, x) && !has(q, x)
            c *= u^q
            continue
        end
        if u.kind == POW && !has(u.args[2], x)
            u, q = u.args[1], expand(u.args[2] * q)
        end
        (q.kind == NUMBER && isexact(q.num)) || return nothing
        (u.kind == FN && length(u.args) == 1 && _eq(u.args[1], x)) || return nothing
        k = q.num
        u.name === :sin ? (a += k) : u.name === :cos ? (b += k) :
        u.name === :tan ? (a += k; b -= k) : u.name === :cot ? (a -= k; b += k) :
        u.name === :csc ? (a -= k) : u.name === :sec ? (b -= k) : return nothing
    end
    (a > SymNumber(-1) && b > SymNumber(-1)) || return nothing
    p, q = Sym((a + N1) / SymNumber(2)), Sym((b + N1) / SymNumber(2))
    gs = (_gamma_exact(p), _gamma_exact(q), _gamma_exact(p + q))
    c * gs[1] * gs[2] / (2gs[3])
end

# ∫₀^∞ c xᵖ/(e^(ax) ∓ 1) dx, a > 0 (decided): Γ(p+1) ζ(p+1)/a^(p+1) for -1 (Bose–Einstein,
# p > 0), (1 - 2^(-p)) Γ(p+1) ζ(p+1)/a^(p+1) for +1 (Fermi–Dirac, p > -1; log(2)/a at
# p = 0). SymPy #24431, where the residue criterion rightly finds no antiderivative.
function _bose_fermi_definite(f::Sym, x::Sym, lo::Sym, hi::Sym, A::Assumptions)
    (iszero(lo) && _isinfpt(hi) && _infsign(hi) > 0) || return nothing
    fs = f.kind == MUL ? [(u, q) for (u, q) in zip(f.args, f.exps)] : f.kind == POW ? [(f.args[1], f.args[2])] : [(f, ONE)]
    c = f.kind == MUL ? Sym(f.num) : ONE
    p, a, σ = ZERO, nothing, nothing
    for (u, q) in fs
        if !has(u, x) && !has(q, x)
            c *= u^q
        elseif _eq(u, x) && q.kind == NUMBER && isexact(q.num)
            p += q
        elseif a === nothing && q.kind == NUMBER && _isN1(-q.num) && u.kind == ADD && length(u.args) == 1 &&
               _isN1(u.coefs[1]) && abs(u.num) == N1 && u.args[1].kind == POW && u.args[1].args[1] === E
            w = u.args[1].args[2]
            a = cancel(w / x)
            has(a, x) && return nothing
            σ = u.num
        else
            return nothing
        end
    end
    a === nothing && return nothing
    _decide(Relation(:gt, a, ZERO), A) || return nothing
    s1 = p + 1
    if σ < N0
        p.num > N0 || return nothing
        return c * apply_function(:gamma, Sym[s1]) * apply_function(:zeta, Sym[s1]) / a^s1
    end
    p.num > SymNumber(-1) || return nothing
    iszero(p) && return c * log(Sym(2)) / a
    c * (1 - Sym(2)^(-p)) * apply_function(:gamma, Sym[s1]) * apply_function(:zeta, Sym[s1]) / a^s1
end

# Gaussian integrals, for a > 0 (decided by cases): ∫ℝ e^(-ax² + bx + c) dx =
# √(π/a) e^(b²/4a + c), and with b = 0, ∫₀^∞ xⁿ e^(-ax²) dx = Γ((n+1)/2)/(2a^((n+1)/2))
# (over ℝ twice that for even n, 0 for odd n). SymPy #28485.
function _gaussian_definite(f::Sym, x::Sym, lo::Sym, hi::Sym, A::Assumptions)
    whole = _isinfpt(lo) && _infsign(lo) < 0 && _isinfpt(hi) && _infsign(hi) > 0
    half = iszero(lo) && _isinfpt(hi) && _infsign(hi) > 0
    (whole || half) || return nothing
    fs = f.kind == MUL ? [(u, q) for (u, q) in zip(f.args, f.exps)] : f.kind == POW ? [(f.args[1], f.args[2])] : [(f, ONE)]
    c = f.kind == MUL ? Sym(f.num) : ONE
    n, Q = 0, nothing
    for (u, q) in fs
        if !has(u, x) && !has(q, x)
            c *= u^q
        elseif _eq(u, x) && q.kind == NUMBER && _isint(q.num) && q.num >= N0
            n += _toint(q.num)
        elseif u === E && Q === nothing && _ispoly_in(q, x) && degree(Poly{Sym}(expand(q), [x]), x) == 2
            Q = expand(q)
        else
            return nothing
        end
    end
    Q === nothing && return nothing
    a, b, c0 = -coeff(Q, x, 2), coeff(Q, x, 1), coeff(Q, x, 0)
    c *= exp(c0)
    if !_decide(Relation(:gt, a, ZERO), A)          # a ≤ 0: the integrand does not decay
        (_isreal(b, A) === true && _isreal(c, A) === true) ||
            throw(ArgumentError("integrate: $f diverges for $a ≤ 0"))
        sg = _signs(c, A)
        sg == _POS && return Sym(Inf)
        sg == _NEG && return Sym(-Inf)
        throw(ArgumentError("integrate: $f diverges for $a ≤ 0"))
    end
    if iszero(b)
        g = _gamma_exact(Sym((n + 1) // 2)) / (2 * a^Sym((n + 1) // 2))
        half && return c * g
        return isodd(n) ? ZERO : 2c * g
    end
    (whole && n == 0) || return nothing
    c * sqrt(PI / a) * exp(b^2 / (4a))
end

# Γ(a) for a rational a > 0: exact at integers and half-integers, else gamma(a)
function _gamma_exact(a::Sym)
    n = a.num
    if _isint(n)
        return Sym(factorial(big(_toint(n) - 1)))
    elseif _isint(n + n)
        m = _toint(n - SymNumber(1 // 2))
        return Sym(factorial(big(2m)) // (big(4)^m * factorial(big(m)))) * sqrt(PI)
    end
    apply_function(:gamma, Sym[a])
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
    # piecewise conditions free of x (b ≤ 1/2 in ∫ piecewise((a ≤ 1/2, b ≤ 1/2) => 1) da):
    # decided by cases on the parameters, then the pieces in x (SymPy #22023)
    if _hasfn(f, :piecewise)
        f = _replace_walk(f) do t
            (t.kind == FN && haskey(_COND_RELS, t.name) && !has(t, x) && _hassymbol(t)) || return nothing
            _decide(_relation_of(t), A) ? ONE : ZERO
        end
    end
    has(f, x) || return _ordered_length(f, lo, hi, A)
    v = _definite_by_pieces(f, x, lo, hi, antider, A, pv)
    v === nothing || return v
    v = _definite_param_pieces(f, x, lo, hi, antider, A, pv)
    v === nothing || return v
    # over (-a, a), finite: the odd terms, continuous there, integrate to 0 (x³cos(x/2)
    # √(4 - x²), which has no elementary antiderivative; SymPy #23246)
    if !_isinfpt(lo) && !_isinfpt(hi) && iszero(expand(lo + hi)) && !_eq(lo, hi)
        g = _drop_odd_terms(f, x, lo, hi)
        g === nothing || return iszero(g) ? ZERO : _definite(g, x, lo, hi, antider, A, pv)
    end
    # over the whole line: an odd f continuous there whose integral converges absolutely,
    # |f| ≤ C/x² far out (e^(-x²) sin x, sin x/(1 + x²); SymPy #21091)
    if _isinfpt(lo) && _isinfpt(hi) && _infsign(lo) < 0 && _infsign(hi) > 0 &&
       iszero(expand(subs(f, x => -x) + f)) && isempty(_singular_points(f, x, lo, hi)) && _abs_convergent_at_inf(f, x, A)
        return ZERO
    end
    S = _singular_points(f, x, lo, hi)
    F = _antiderivative(antider, f, x, A; numeric = true)
    # an antiderivative known only numerically that has an expression valid on the
    # interval (RischIntegrate's ₂F₁ terms, where their radicand is positive): that, so
    # the ends may be infinite (SymPy #24879)
    if !(F isa Sym) && _SYMBOLIC_ON[] !== nothing
        G = _SYMBOLIC_ON[](F, x, lo, hi)
        G isa Sym && (F = G)
    end
    # F may jump where f is smooth: 2atan(√(A²-1) tan(t/2)/(A-1))/√(A²-1) for 1/(A - cos t)
    # at t = π (#26283); split there too (a point where F is continuous changes nothing)
    F isa Sym && (S = _same_points(vcat(S, _jump_points(F, x, lo, hi))))
    if isempty(S) && !_isinfpt(lo) && !_isinfpt(hi)
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
        # a point k acos(u) + m: real for u² ≤ 1 (cos x = -C: SymPy #23613)
        ac = _hassymbol(s) ? _acos_affine(s) : nothing
        ac === nothing || _decide(Relation(:le, ac[1]^2, ONE), A) || continue
        (_isinfpt(lo) || _decide_pt(:gt, s, lo, A)) && (_isinfpt(hi) || _decide_pt(:lt, s, hi, A)) || continue
        any(t -> _eq(t, s), inside) || _insert_sorted!(inside, s, A)
    end
    union!(sing, inside)
    if !(F isa Sym)
        (_isinfpt(lo) || _isinfpt(hi)) && throw(ArgumentError(
            "integrate: the antiderivative of $f is known only numerically, which does not reach ±∞"))
        isempty(sing) || throw(ArgumentError(
            "integrate: the antiderivative of $f is known only numerically, and $f is singular on the interval"))
        return _numeric_difference(F, f, lo, hi)
    end
    # a simple pole inside, (x - s) f tending to a nonzero number: no integral, only a
    # principal value (F's two sides there can look equal: atan(±i) for 1/(C + cos x),
    # |C| < 1, came out 0; SymPy #23613)
    if !pv
        for s in inside
            _simple_pole(f, x, s, A) && throw(ArgumentError(
                "integrate: the integral does not converge: $f has a pole at $x = $s inside the interval"))
        end
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
    _decide_pt(:lt, p, q, A)
end

# s = k acos(u) + m, k and m numbers, u with parameters: (u, k, m), or nothing
function _acos_affine(s::Sym)
    as = Sym[]
    _replace_walk(s) do t
        t.kind == FN && t.name === :acos && length(t.args) == 1 && _hassymbol(t.args[1]) &&
            !any(b -> _eq(b, t), as) && push!(as, t)
        nothing
    end
    length(as) == 1 || return nothing
    a = only(as)
    k = expand(subs(s, a => ONE) - subs(s, a => ZERO))
    m = expand(subs(s, a => ZERO))
    (k.kind == NUMBER && !iszero(k.num) && !_hassymbol(m) && iszero(expand(s - k * a - m))) || return nothing
    a.args[1], k, m
end

# p op q (op :lt or :gt) for p or q a point k acos(u) + m and the other a number: as a
# condition on u, acos taking its values in [0, π] (acos(u) > 0 is u < 1), rather than one
# on acos(u); else as it is
function _decide_pt(op::Symbol, p::Sym, q::Sym, A)
    rel = Relation(op, p, q)
    d = expand(p - q)
    ac = _acos_affine(d)
    ac === nothing && return _decide(rel, A)
    u, k, m = ac
    t = try Float64(N(-m / k)) catch; return _decide(rel, A) end       # d op 0 ⇔ acos(u) op′ t
    gt = (op === :gt) == (k.num > N0)
    if gt                                                              # acos(u) > t
        t < 0 && return true
        t >= π && return false
        return _decide(Relation(:lt, u, cos(-m / k)), A)
    else                                                               # acos(u) < t
        t <= 0 && return false
        t > π && return true
        return _decide(Relation(:gt, u, cos(-m / k)), A)
    end
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

# set by the RischIntegrate extension: (F, x, lo, hi) -> an expression for the numeric
# antiderivative F valid on (lo, hi), or nothing
const _SYMBOLIC_ON = Ref{Any}(nothing)

# F(hi) - F(lo), a float, for an antiderivative F known only numerically and continuous on
# [lo, hi] (no singular points there); the ends must be numbers
function _numeric_difference(F, f::Sym, lo::Sym, hi::Sym)
    (_hassymbol(lo) || _hassymbol(hi)) && throw(ArgumentError(
        "integrate: the antiderivative of $f is known only numerically; the limits must be numbers"))
    Sym(Float64(real(F(N(hi)) - F(N(lo)))))
end

# F at p from the side dir: its value, or its one-sided limit where that is needed
function _endvalue(F::Sym, x::Sym, p::Sym, dir::Symbol, A, singular::Bool = false)
    G = _hyp_endpoint(F, x, p, dir, A)            # ₂F₁ terms by their values at p
    G === nothing || (F = G)
    if !singular && !_isinfpt(p)
        v = try subs(F, x => p) catch err; _isdivzero(err) || rethrow(); nothing end
        # (a piecewise F at a parametric end, as ∫|x| at a/b, its case decided if A can;
        # 0^x with x > 0 is 0)
        v !== nothing && (_hasfn(v, :piecewise) || _has_zero_power(v)) && (v = refine(v; assume = A))
        v !== nothing && !_has_zero_power(v) && _finite_value(v) && return v
    end
    _onesided(F, x, p, dir, A)
end

_has_zero_power(e::Sym) = (e.kind == POW && iszero(e.args[1])) ||
    (e.kind == MUL && any(iszero, e.args)) || any(_has_zero_power, e.args) ||
    (e.kind == MUL && any(_has_zero_power, e.exps))

# F with each ₂F₁(a, b; c; u) whose argument tends to u₀ at p (from the side dir)
# replaced by its value there: by continuity for |u₀| < 1, and by Gauss's sum
# Γ(c)Γ(c-a-b)/(Γ(c-a)Γ(c-b)) for u₀ = 1 when c - a - b > 0 (Abel's theorem). The limit
# code does not know ₂F₁, which made ∫₀¹ z^(x-1)(1-z)^(y-1) dz "divergent" (SymPy
# #26477, #26068). Only valid at p; nothing when F has no such term to replace.
function _hyp_endpoint(F::Sym, x::Sym, p::Sym, dir::Symbol, A)
    _hasfn(F, :_₂F₁) || return nothing
    x0 = _isinfpt(p) ? Float64(p.num) : p
    changed = Ref(false)
    G = _replace_walk(F) do t
        (t.kind == FN && t.name === :_₂F₁ && length(t.args) == 4 && has(t.args[4], x)) || return nothing
        a, b, c, u = t.args
        any(q -> has(q, x), (a, b, c)) && return nothing
        u0 = _limit_or(u, x, x0, dir, A)
        u0 isa Sym && !_isinfpt(u0) || return nothing
        if iszero(u0)
            changed[] = true
            return ONE
        elseif iszero(expand(u0 - 1))
            d = expand(c - a - b)
            _decide(Relation(:gt, d, ZERO), A) || return nothing
            changed[] = true
            Γ(w) = apply_function(:gamma, Sym[w])
            return Γ(c) * Γ(d) / (Γ(expand(c - a)) * Γ(expand(c - b)))
        end
        w = _hassymbol(u0) ? nothing : (try N(u0) catch; nothing end)
        (w isa Number && abs(w) < 1) || return nothing
        changed[] = true
        apply_function(:_₂F₁, Sym[a, b, c, u0])
    end
    changed[] ? G : nothing
end

function _finite_value(v::Sym)
    _hassymbol(v) && return true
    w = try N(v) catch; return false end
    isfinite(w)
end

function _onesided(F::Sym, x::Sym, p::Sym, dir::Symbol, A)
    x0 = _isinfpt(p) ? Float64(p.num) : p
    fail = Ref{Any}(nothing)
    L = _limit_or(F, x, x0, dir, A, fail)
    if L === nothing
        Gc = _neglogs(F, x, p, dir)          # log(u) with u < 0 near p: log(-u) + iπ
        if Gc !== nothing
            vs = [k => _limit_or(g, x, x0, dir, A) for (k, g) in Gc]
            if all(v -> v.second isa Sym, vs)
                if all(v -> !_isinfpt(v.second), vs)
                    L = _tidy_constant(sum(k * v for (k, v) in vs; init = ZERO))
                elseif count(v -> _isinfpt(v.second), vs) == 1 &&
                       (iv = findfirst(v -> _isinfpt(v.second), vs); isone(vs[iv].first))
                    # the real part infinite, the constants iπ beside it finite: ∫₁³
                    # dx/(x - 3), log(x - 3) = log(3 - x) + iπ, is -∞ (SymPy #23337)
                    L = vs[iv].second
                end
            end
        end
    end
    if L === nothing
        # a limit the limit code could not take (an unsupported function) is no proof
        # of divergence
        # (a limit the code could not take, an unsupported function or an undecided
        # sign, is no proof of divergence: only a limit that does not exist is)
        e = fail[]
        e isa GruntzError && e.kind !== :nolimit && throw(ArgumentError(
            "integrate: cannot find the limit of the antiderivative $F at $x = $p ($(sprint(showerror, e)))"))
        throw(ArgumentError("integrate: the integral does not converge: the antiderivative $F has no limit at $x = $p"))
    end
    L
end

# the limit, or nothing (the error, if one was the cause, left in fail[])
function _limit_or(F, x, x0, dir, A, fail = Ref{Any}(nothing); depth::Int = 0)
    L = try
        limit(F, x => x0; dir = dir === :+ ? :+ : :-, assume = A)
    catch err
        err isa Union{GruntzError,ArgumentError,DomainError} || rethrow()
        # a sign the limit could not split on: decided here, a case for each sign of a
        # parameter whose sign is open (e^(-ax) sin x at ∞ with a ≠ 0; SymPy #21054)
        if err isa GruntzError && err.kind === :fail && depth < 2 && A isa Assumptions
            ps = [p for p in free_symbols(Sym(F)) if !_eq(p, x) && !_eq(p, _MARK_IPI) &&
                  (m = _signs(p, A); m === nothing || (m != _POS && m != _NEG))]
            if !isempty(ps)
                p = first(ps)
                pos = _decide(Relation(:gt, p, ZERO), A)
                A2 = _with_fact(A, Relation(pos ? :gt : :le, p, ZERO))
                A2 === nothing || return _limit_or(F, x, x0, dir, A2, fail; depth = depth + 1)
            end
        end
        fail[] = err
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

# f = rest + Σ cᵢ gᵢ, the gᵢ the terms with an unknown function of x (f(x), a derivative
# of one, a held Sum or Integral), cᵢ their factors free of x: the integrator does not
# take those (RischNotImplemented on SymFunction; SymPy #25457, #26116, #26249), so
# their integrals are held and the rest integrated
function _split_opaque(f::Sym, x::Sym)
    held = Tuple{Sym,Sym}[]
    _has_opaque(f, x) || return f, held
    terms = f.kind == ADD ? [Sym(c) * a for (c, a) in zip(f.coefs, f.args)] : [f]
    rest = f.kind == ADD ? Sym(f.num) : ZERO
    for t in terms
        if !_has_opaque(t, x)
            rest += t
        elseif t.kind == MUL
            c = Sym(t.num) * prod((b^p for (b, p) in zip(t.args, t.exps) if !has(b, x) && !has(p, x)); init = ONE)
            push!(held, (c, cancel(t / c)))
        else
            push!(held, (ONE, t))
        end
    end
    rest, held
end

_has_opaque(e::Sym, x::Sym) = has(e, x) &&
    ((e.kind == FN && (e.name in (:Derivative, :Integral, :Sum, :Product) ||
                       (!haskey(FUNCTION_RULES, e.name) && !(e.name in _NO_CHAIN)))) ||
     any(a -> _has_opaque(a, x), e.args) || (e.kind == MUL && any(a -> _has_opaque(a, x), e.exps)))

# the singular points of an antiderivative F, where it may jump; none when they cannot be
# located (F is then trusted to be continuous, as before)
function _jump_points(F::Sym, x::Sym, lo::Sym, hi::Sym)
    try
        _singular_points(F, x, lo, hi)
    catch err
        err isa ArgumentError || rethrow()
        Sym[]
    end
end

function _zero_sets!(out, e::Sym, x::Sym)
    has(e, x) || return out
    k = e.kind
    if k == MUL || k == POW
        bs, ps = k == MUL ? (e.args, e.exps) : (e.args[1:1], e.args[2:2])
        for (b, p) in zip(bs, ps)
            # a negative power, or a symbolic one (x^n, singular at 0 for n < 0)
            # (the zeros of |u| are those of u)
            bz = b.kind == FN && b.name === :abs && length(b.args) == 1 ? b.args[1] : b
            has(b, x) && ((p.kind == NUMBER && p.num < N0) || (p.kind != NUMBER && !has(p, x))) && push!(out, bz)
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
    # none, by bounds: 1 + t^2 sin(t)^2 ≥ 1 (#27139)
    _sign_on(z, x, _num_end(lo, -Inf), _num_end(hi, Inf)) === nothing || return Sym[]
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
    rs = try _solve_generic(z ~ 0, x) catch err
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

# a limit of integration as a Float64, or `inf` when it is symbolic
function _num_end(e::Sym, inf::Float64)
    _hassymbol(e) && return inf
    v = try N(e) catch; return inf end
    v isa Real ? Float64(v) : inf
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

# a piecewise integrand (abs, sign, max, min, heaviside, piecewise) over numeric limits:
# the sum over the intervals between the points inside where a condition changes, each
# with its pieces fixed. A periodic condition (sin(2x) ≥ 0) changes in every period, and
# an antiderivative found by pieces knows only one period's points, so |sin(2x)| over
# [0, 2π] would come out as 2, not 4. nothing when there is nothing to split.
function _definite_by_pieces(f::Sym, x::Sym, lo::Sym, hi::Sym, antider, A, pv::Bool)
    (_isinfpt(lo) || _isinfpt(hi) || _hassymbol(lo) || _hassymbol(hi)) && return nothing
    g = piecewise(f)
    _haspieces(g, x) || return nothing
    a, b = N(lo), N(hi)
    (a isa Real && b isa Real) || return nothing
    a, b = minmax(a, b)
    # conditions of one sign on [a, b] by bounds (|3 + cos(3x)^5/2 - …| ≥ 1, #27283)
    g = _settle_conditions(g, x, a, b)
    _haspieces(g, x) || return _definite(g, x, lo, hi, antider, A, pv)
    cuts = _piece_cuts(g, x, a, b)
    cuts === nothing && return nothing
    if isempty(cuts)                         # the pieces are fixed on the interval
        gi = _pieces_at(g, x, Sym(SymNumber(rationalize(BigInt, (a + b) / 2; tol = (b - a) / 8))))
        return gi === nothing ? nothing : _definite(gi, x, lo, hi, antider, A, pv)
    end
    pts = Sym[(N(lo) <= N(hi) ? lo : hi); cuts; (N(lo) <= N(hi) ? hi : lo)]
    terms = Tuple{Int,Sym}[]
    for i in 1:(length(pts) - 1)
        p, q = pts[i], pts[i + 1]
        u, w = Float64(N(p)), Float64(N(q))
        m = Sym(SymNumber(rationalize(BigInt, (u + w) / 2; tol = (w - u) / 8)))
        gi = _pieces_at(g, x, m)
        gi === nothing && return nothing
        v = _definite(gi, x, p, q, antider, A, pv)
        v isa Sym || return nothing
        push!(terms, (1, v))
    end
    v = _sum_values(terms, A)
    N(lo) <= N(hi) ? v : -v
end

# the points strictly inside (a, b) where a condition of g in x changes, sorted, with a
# periodic condition's points repeated over the interval; nothing when they can't be found
function _piece_cuts(g::Sym, x::Sym, a::Real, b::Real)
    out = Sym[]
    ok = Ref(true)
    function visit(n::Sym)
        ok[] || return
        if n.kind == FN && (haskey(_COND_RELS, n.name) || n.name === :heaviside) && has(n, x)
            u = haskey(_COND_RELS, n.name) ? n.args[1] - n.args[2] : n.args[1]
            rs = try _solve_generic(u ~ 0, x) catch; (ok[] = false; return) end
            rs isa AbstractVector || (ok[] = false; return)
            # and where u is not defined, as it may change sign there (atan(1/x) at 0,
            # where the midpoint was then taken; SymPy #20755)
            for z in unique(_zero_sets!(Sym[], u, x))
                zs = try _solve_generic(z ~ 0, x) catch; continue end
                zs isa AbstractVector && (rs = vcat(rs, zs))
            end
            P = _trig_period(u, x)
            p = P === nothing ? nothing : N(P)
            for r in rs
                _hassymbol(r) && (ok[] = false; return)
                v = try N(r) catch; continue end
                v isa Real || continue
                ns = p === nothing ? (0:0) : ((floor(Int, (a - v) / p) - 1):(ceil(Int, (b - v) / p) + 1))
                for k in ns
                    c = p === nothing ? r : expand(r + k * P)
                    w = v + (p === nothing ? 0 : k * p)
                    (a < w < b) && !any(q -> isapprox(N(q), w; atol = 1e-12), out) && push!(out, c)
                end
            end
        end
        foreach(visit, n.args)
        n.kind == MUL && foreach(visit, n.exps)
    end
    visit(g)
    ok[] || return nothing
    sort!(out; by = N)
end

# the sign oracle for an integrator: the sign of a parameter expression under A, as 1 or
# -1, or nothing when A does not decide it
function _sign_oracle(A::Assumptions)
    function (p)
        m = _signs(Sym(p), A)
        m == _POS ? 1 : m == _NEG ? -1 : nothing
    end
end

# drop_removable_singularities = true: a case p = q whose value is the limit of the
# general formula as p → q (both sides) is a removable singularity of that formula, and
# is merged into it (∫₀¹ e^{ct} dt: (e^c - 1)/c, its value 1 at c = 0 the limit), the
# condition p ≠ q dropped from the general case. Repeated until none is left; a single
# case left with no conditions is its value.
function _drop_removable(v)
    v isa Cases || return v
    cs = collect(v.cases)
    changed = true
    while changed
        changed = false
        for (i, (ci, vi)) in enumerate(cs)
            k = findfirst(c -> c isa Relation && c.op === :eq && _isnamed(c.lhs), ci)
            k === nothing && continue
            eq = ci[k]
            rest = ci[[j for j in eachindex(ci) if j != k]]
            # the general case: the same conditions with p ≠ q in place of p = q
            j = findfirst(eachindex(cs)) do j
                j == i && return false
                cj = cs[j].first
                length(cj) == length(ci) || return false
                m = findfirst(c -> c isa Relation && c.op === :ne && _eq(c.lhs, eq.lhs) && _eq(c.rhs, eq.rhs), cj)
                m !== nothing && all(c -> any(r -> isequal(r, c), rest), cj[[t for t in eachindex(cj) if t != m]])
            end
            j === nothing && continue
            g = cs[j].second
            (g isa Sym && vi isa Sym) || continue
            _removable_at(g, eq.lhs, eq.rhs, vi) || continue
            cj = cs[j].first
            cs[j] = Any[c for c in cj if !(c isa Relation && c.op === :ne && _eq(c.lhs, eq.lhs) && _eq(c.rhs, eq.rhs))] => g
            deleteat!(cs, i)
            changed = true
            break
        end
    end
    length(cs) == 1 && isempty(cs[1].first) && return cs[1].second
    typeof(v)(cs)
end

# is val the limit of g as p → q, from both sides? (g simplified when the limit fails:
# -i e^{ib}/b + i e^{-ib}/b is 2 sin(b)/b)
function _removable_at(g::Sym, p::Sym, q::Sym, val::Sym)
    lim(h, dir) = try
        limit(h, p => q; dir = dir)
    catch err
        err isa Union{GruntzError,ArgumentError,DomainError} || rethrow()
        nothing
    end
    for dir in (:+, :-)
        L = lim(g, dir)
        L === nothing && (L = lim(simplify(g), dir))
        L === nothing && return false
        L isa Sym || return false
        d = cancel(L - val)
        iszero(d) || (!_hassymbol(d) && abs(N(d)) < 1e-12) || return false
    end
    true
end


# a piecewise integrand whose breakpoints in x move with one parameter y (a convolution
# of boxes: 0, 1, y - 1, y): the values of y where two of them (or one and a finite end)
# meet cut its line into intervals, on each of which their order is fixed; the integral
# on each, under that interval's assumption, makes a piecewise function of y (SymPy
# #23761). nothing when the points are not linear in one parameter.
const _PARAM_CUTS_MAX = 12
function _definite_param_pieces(f::Sym, x::Sym, lo::Sym, hi::Sym, antider, A, pv::Bool)
    (_hassymbol(lo) || _hassymbol(hi)) && return nothing
    g = piecewise(f)
    _haspieces(g, x) || return nothing
    bps = _all_breakpoints(g, x)
    (bps === nothing || !any(_hassymbol, bps)) && return nothing
    ys = unique(s for r in bps for s in free_symbols(r))
    length(ys) == 1 || return nothing
    y = only(ys)
    function lin(r)
        c1 = try coeff(expand(r), y, 1) catch; return nothing end
        c0 = expand(r - c1 * y)
        (c1.kind == NUMBER && c0.kind == NUMBER) ? (c1.num, c0.num) : nothing
    end
    ls = map(lin, bps)
    any(isnothing, ls) && return nothing
    ends = [p for p in (lo, hi) if !_isinfpt(p)]
    crit = Float64[]
    addc(v) = (isfinite(v) && !any(c -> isapprox(c, v; atol = 1e-12), crit)) && push!(crit, v)
    for i in eachindex(ls)
        a1, b1 = ls[i]
        for j in 1:(i - 1)
            a2, b2 = ls[j]
            a1 == a2 || addc(Float64((b2 - b1) / (a1 - a2)))
        end
        for e in ends
            iszero(a1) || addc(Float64((e.num - b1) / a1))
        end
    end
    isempty(crit) && return nothing
    length(crit) <= _PARAM_CUTS_MAX || return nothing
    sort!(crit)
    cy = [Sym(SymNumber(rationalize(BigInt, c; tol = 1e-14))) for c in crit]
    pieces = Pair[]
    for k in 0:length(crit)
        y0 = k == 0 ? crit[1] - 1 : k == length(crit) ? crit[end] + 1 : (crit[k] + crit[k+1]) / 2
        Ak = A
        k > 0 && (Ak = _with_fact(Ak, Relation(:gt, y, cy[k])))
        Ak === nothing && continue
        k < length(crit) && (Ak = _with_fact(Ak, Relation(:lt, y, cy[k+1])))
        Ak === nothing && continue
        # the points inside (lo, hi) at y0, in order
        l0, h0 = Float64(N(lo)), Float64(N(hi))
        vals = [(Float64(a) * y0 + Float64(b), r) for ((a, b), r) in zip(ls, bps)]
        inside = sort!([(v, r) for (v, r) in vals if l0 < v < h0]; by = first)
        pts = Tuple{Float64,Sym}[(l0, lo); inside; (h0, hi)]
        total = ZERO
        ysub = y => Sym(SymNumber(rationalize(BigInt, y0; tol = 1e-14)))
        for i in 1:(length(pts) - 1)
            (u, p), (w, q) = pts[i], pts[i+1]
            isapprox(u, w; atol = 1e-12) && continue
            m = isinf(u) ? w - 1 : isinf(w) ? u + 1 : (u + w) / 2
            gi = _pieces_at(g, x, Sym(SymNumber(rationalize(BigInt, m; tol = 1e-9))); extra = (ysub,))
            gi === nothing && return nothing
            v = _definite(gi, x, p, q, antider, Ak, pv)
            v isa Sym || return nothing
            total += v
        end
        c = k == length(crit) ? true : (y ≦ cy[k+1])
        push!(pieces, c => simplify(total))
    end
    isempty(pieces) && return nothing
    piecewise(pieces...)
end

# every breakpoint of g in x, numeric or with parameters; nothing when one can't be found
function _all_breakpoints(e::Sym, x::Sym)
    out = Sym[]
    ok = Ref(true)
    function visit(n::Sym)
        ok[] || return
        if n.kind == FN && (haskey(_COND_RELS, n.name) || n.name === :heaviside) && has(n, x)
            u = haskey(_COND_RELS, n.name) ? n.args[1] - n.args[2] : n.args[1]
            rs = try _solve_generic(u ~ 0, x) catch; (ok[] = false; return) end
            rs isa AbstractVector || (ok[] = false; return)
            for r in rs
                (_hassymbol(r) || _isreal_const(r) === true) || continue
                any(q -> iszero(expand(q - r)), out) || push!(out, r)
            end
            n.name === :heaviside && return
        end
        foreach(visit, n.args)
        n.kind == MUL && foreach(visit, n.exps)
    end
    visit(e)
    ok[] && !isempty(out) ? out : nothing
end

# f without its odd terms in x (f(-x) = -f(x)) that have no singular point on [lo, hi];
# nothing when there are none
function _drop_odd_terms(f::Sym, x::Sym, lo::Sym, hi::Sym)
    odd(t) = iszero(expand(subs(t, x => -x) + t))
    ok(t) = isempty(_singular_points(t, x, lo, hi))
    odd(f) && ok(f) && return ZERO
    e = f.kind == ADD ? f : expand(f)
    e.kind == ADD || return nothing
    terms = [Sym(c) * t for (c, t) in zip(e.coefs, e.args)]
    keep = [!(has(t, x) && odd(t) && ok(t)) for t in terms]
    all(keep) && return nothing
    _add_many(t for (t, k) in zip(terms, keep) if k) + Sym(e.num)
end

# f = n/d with d(s) = 0, d'(s) ≠ 0 and n(s) ≠ 0 (decided under A)
function _simple_pole(f::Sym, x::Sym, s::Sym, A)
    n, d = numer_denom(f)
    has(d, x) || return false
    at(e) = try simplify(subs(e, x => s)) catch err; err isa InterruptException && rethrow(); nothing end
    nonzero(v) = v isa Sym && (_hassymbol(v) ? (try _decide(Relation(:ne, v, ZERO), A) catch; false end) :
                               (w = try N(v) catch; nothing end; w isa Number && isfinite(w) && abs(w) > 1e-12))
    d0 = at(d)
    d0 isa Sym && iszero(d0) || return false
    nonzero(at(diff(d, x))) && nonzero(at(n))
end

# points that are one point merged: equal, or equal in value (at sample values of their
# parameters): acos(-C) and 2atan((C + 1)/√(1 - C²)) (SymPy #23613)
function _same_points(S::Vector{Sym})
    out = Sym[]
    for p in S
        any(q -> _eq(p, q) || _same_value(p, q), out) || push!(out, p)
    end
    out
end
function _same_value(p::Sym, q::Sym)
    syms = unique(vcat(free_symbols(p), free_symbols(q)))
    # (samples where both are real: the points are real there, branches agree)
    agree = 0
    for k in 1:(isempty(syms) ? 1 : 6)
        vals = [s => Sym(_SAMPLE_VALUES[mod1(k + 3j, length(_SAMPLE_VALUES))]) for (j, s) in enumerate(syms)]
        a, b = try (ComplexF64(N(subs(p, vals...))), ComplexF64(N(subs(q, vals...)))) catch; continue end
        (isfinite(a) && isfinite(b)) || continue
        isempty(syms) || (abs(imag(a)) < 1e-12 && abs(imag(b)) < 1e-12) || continue
        abs(a - b) <= 1e-10 * max(1, abs(a)) || return false
        agree += 1
    end
    agree >= (isempty(syms) ? 1 : 2)
end

# integrals over (-∞, ∞) or (0, ∞) with no elementary antiderivative, by known values and
# parts: a sum term by term; c sin(a x²), c cos(a x²) (Fresnel: ∫ sin(ax²) dx over the
# line is sign(a) √(π/(2|a|))); and n(x)/x² with n(0) = 0 by parts, ∫ n'(x)/x dx with the
# boundary terms -n/x (∫ (cos(x²) - cos x)/x² dx = π - √(2π); SymPy #22487). nothing when
# none applies. depth bounds the parts.
function _definite_fallback(g::Sym, x::Sym, l::Sym, h::Sym, antider, A, pv::Bool, depth::Int)
    line = _isinfpt(l) && _infsign(l) < 0 && _isinfpt(h) && _infsign(h) > 0
    half = iszero(l) && _isinfpt(h) && _infsign(h) > 0
    # to +∞ from a finite point: an integrand whose limit there is not 0 (infinite, or a
    # nonzero number) diverges, ±∞ by its sign (eˣ e^(log² x) on (1, ∞); SymPy #20831)
    if !_isinfpt(l) && _isinfpt(h) && _infsign(h) > 0
        L = try limit(g, x => Inf; assume = A) catch err; err isa InterruptException && rethrow(); nothing end
        if L isa Sym && !_hassymbol(L)
            _isinfpt(L) && !has(L, IM) && return L
            (L.kind == NUMBER && !iszero(L.num) && isreal(L.num)) && return L.num > N0 ? Sym(Inf) : Sym(-Inf)
        end
    end
    # ∫₀^∞ c xᵖ logᵏ(x) e^(-a x^q) dx = c dᵏ/dsᵏ [Γ(s/q)/(q a^(s/q))] at s = p + 1, for
    # a > 0, q > 0, p > -1 (e^(-x³): Γ(4/3); t log²(t) e^(-t): Γ''(2); SymPy #20685, #20043)
    if half
        v = _gamma_moment(g, x, A)
        v === nothing || return v
    end
    (line || half) || return nothing
    one(t) = try
        _definite(t, x, l, h, antider, A, pv)
    catch err
        err isa InterruptException && rethrow()
        _definite_fallback(t, x, l, h, antider, A, pv, depth)
    end
    if g.kind == ADD
        iszero(g.num) || return nothing                  # a constant: divergent
        vs = [one(Sym(c) * t) for (c, t) in zip(g.coefs, g.args)]
        all(v -> v isa Sym, vs) || return nothing
        return _sum_values([(1, v) for v in vs], A)
    end
    c, t = _split_kfree(g, x)
    if t.kind == FN && t.name in (:sin, :cos) && length(t.args) == 1
        a = expand(t.args[1] / x^2)
        if a.kind == NUMBER && !iszero(a.num)
            v = sqrt(PI / (2 * abs(a)))
            t.name === :sin && a.num < N0 && (v = -v)
            return c * (line ? v : v / 2)
        end
    end
    depth >= 1 || return nothing
    n = expand(g * x^2)
    (has(n, x) && !has(numer_denom(cancel(n))[2], x)) || return nothing
    iszero(simplify(subs(n, x => ZERO))) || return nothing
    B = -n / x
    ends = [_limit_or(B, x, Float64(h.num), :-, A), line ? _limit_or(B, x, -Inf, :+, A) : _limit_or(B, x, ZERO, :+, A)]
    all(e -> e isa Sym && !_isinfpt(e), ends) || return nothing
    rest = _definite_fallback_or(expand(diff(n, x) / x), x, l, h, antider, A, pv, depth - 1)
    rest isa Sym || return nothing
    simplify(ends[1] - ends[2] + rest)
end
function _gamma_moment(g::Sym, x::Sym, A)
    fs = g.kind == MUL ? [(b, p) for (b, p) in zip(g.args, g.exps)] :
         g.kind == POW ? [(g.args[1], g.args[2])] : [(g, ONE)]
    c = g.kind == MUL ? Sym(g.num) : ONE
    p, k, a, q = ZERO, 0, nothing, nothing
    for (b, e) in fs
        if !has(b, x) && !has(e, x)
            c *= b^e
        elseif _eq(b, x) && !has(e, x)
            p += e
        elseif b.kind == FN && b.name === :log && length(b.args) == 1 && _eq(b.args[1], x) &&
               e.kind == NUMBER && isexact(e.num) && _isint(e.num) && e.num > N0
            k += _toint(e.num)
        elseif b === E && has(e, x) && a === nothing
            w = expand(e)
            # -a x^q: a monomial in x
            (w.kind == MUL || w.kind == POW || _eq(w, x)) || return nothing
            qq = w.kind == MUL ? sum((ee for (bb, ee) in zip(w.args, w.exps) if _eq(bb, x)); init = ZERO) :
                 w.kind == POW ? (_eq(w.args[1], x) ? w.args[2] : return nothing) : ONE
            aa = cancel(-w / x^qq)
            has(aa, x) && return nothing
            a, q = aa, qq
        else
            return nothing
        end
    end
    (a === nothing || has(p, x)) && return nothing
    (q.kind == NUMBER && q.num > N0) || return nothing
    _decide(Relation(:gt, a, ZERO), A) || return nothing
    pv = _decide(Relation(:gt, p + 1, ZERO), A)
    pv === true || return nothing
    s = Sym(Symbol("_gm_s"))
    F = apply_function(:gamma, Sym[s / q]) / (q * a^(s / q))
    for _ in 1:k
        F = diff(F, s)
    end
    simplify(c * subs(F, s => p + 1))
end

function _definite_fallback_or(t::Sym, x, l, h, antider, A, pv, depth)
    try
        _definite(t, x, l, h, antider, A, pv)
    catch err
        err isa InterruptException && rethrow()
        _definite_fallback(t, x, l, h, antider, A, pv, depth)
    end
end

# |f| ≤ g with each sin, cos factor of f bounded by 1, and x² g bounded at ±∞
function _abs_convergent_at_inf(f::Sym, x::Sym, A)
    trig(b, p) = b.kind == FN && b.name in (:sin, :cos) && length(b.args) == 1 && p.kind == NUMBER && p.num > N0
    g = f.kind == MUL ? Sym(f.num) * prod((b^p for (b, p) in zip(f.args, f.exps) if !trig(b, p)); init = ONE) : f
    (_hasfn(g, :sin) || _hasfn(g, :cos)) && return false
    for x0 in (Inf, -Inf)
        L = try limit(x^2 * g, x => x0; assume = A) catch err; err isa InterruptException && rethrow(); return false end
        (L isa Sym && !_hassymbol(L) && !_isinfpt(L) && L.kind == NUMBER) || return false
    end
    true
end
