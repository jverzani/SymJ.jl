# solveset: every solution of an equation in one unknown, as a set: finitely many points,
# families {f(n) : n ∈ ℤ} for periodic and multivalued inverses, and a ConditionSet for
# what is not solved. `solve` gives finite lists (one period, principal branches) and is
# left as it is; see notes/solveset_scope.md.

"""
    SymJ.ImageSet

A family of solutions `{expr : n ∈ ℤ}` from [`SymJ.solveset`](@ref): `expr` in the bound
integer variables `vars`, which satisfy the relations `conds` (`n ≥ 0`; none: all of ℤ).
Not exported.
"""
struct ImageSet
    expr::Sym
    vars::Vector{Sym}
    conds::Vector{Relation}
end

"""
    SymJ.ConditionSet

The solutions of `eq` in `domain` that [`SymJ.solveset`](@ref) did not find in closed
form: `{x ∈ domain : eq}`. Not exported.
"""
struct ConditionSet
    x::Sym
    eq::Relation
    domain::NumberSet
end

"""
    SymJ.SolutionSet

The set [`SymJ.solveset`](@ref) gives: the union of finitely many `points`, of the
families in `families` ([`SymJ.ImageSet`](@ref)), and of the unsolved parts in
`unsolved` ([`SymJ.ConditionSet`](@ref)). Supports `in` (for a family not linear in its
variable, the variable solved for, or Lambert's branch found numerically:
`lambertw(1, -1) in {lambertw(1, n) | n ∈ ℤ}`), `isempty`, `∩` with an interval set (a
family linear in its variable cut to its finitely many members), and `collect` when it is
finite. Not exported.
"""
struct SolutionSet
    x::Sym
    points::Vector{Sym}
    families::Vector{ImageSet}
    unsolved::Vector{ConditionSet}
end

struct _SSUnsolved <: Exception end

const _SS_NAMES = (:n, :m, :k, :j, :l, :p, :q)

"""
    SymJ.solveset(eq, x; domain = ℝ)

Every solution of the equation `eq` in `x`, as a set ([`SymJ.SolutionSet`](@ref)): the
real ones, or with `domain = ℂ` the complex ones. Where `solve` gives one period or a
principal branch, `solveset` gives the families: `sin(x) ~ 1//2` has
`{π/6 + 2πn : n ∈ ℤ} ∪ {5π/6 + 2πn : n ∈ ℤ}`, and over ℂ, `exp(x) ~ 2` has
`{log(2) + 2πin : n ∈ ℤ}`. Nested functions are inverted in turn: `sin(sin(x)) ~ 0` is
`sin(x) ∈ {πm}`, which over ℝ is `sin(x) = 0`. An equation it cannot solve gives a
[`SymJ.ConditionSet`](@ref) rather than an error, and an inequality gives `solve`'s set.
The bound variables are named `n`, `m`, `k`, … (skipping names in `eq`). Over ℝ,
parameters are taken real and the answer comes by cases on them: `sin(x) ~ a` has the
families for `-1 ≤ a ≤ 1` and `∅` otherwise. Several trigonometric functions of one
angle with a parameter go through one of them (`tan(x) + cot(x) ~ a` in `tan(x)`) or
`tan(x/2)`, the angle `π` decided apart (`sin(x) + cos(x) ~ a`: `-√2 ≤ a ≤ √2, a ≠ -1`,
and `a = -1`). An equation through Lambert's W has every
branch: over ℂ `{W_k(a) : k ∈ ℤ}` (shown `lambertw(a, n)`), over ℝ the real branches
`W₀` and, on `[-1/e, 0)`, `W₋₁`. With Latexify loaded, the sets display as LaTeX. Not
exported.

```julia
SymJ.solveset(sin(x) ~ 1//2, x)              # {π/6 + 2πn} ∪ {5π/6 + 2πn}
SymJ.solveset(exp(x) ~ 2, x)                 # {log(2)}
SymJ.solveset(exp(x) ~ 2, x; domain = ℂ)     # {log(2) + 2πin}
SymJ.solveset(sin(x^2) ~ 0, x)               # {-√(πn), √(πn) : n ≥ 0}
S = SymJ.solveset(cos(x) ~ 0, x)
S ∩ solve([x ≧ 0, x ≦ 2PI], [x])             # {π/2, 3π/2}
SymJ.solveset(sin(x) ~ a, x)                 # a < -1 => ∅, -1 ≤ a ≤ 1 => {asin(a) + 2πn} ∪ …, a > 1 => ∅
SymJ.solveset(x*exp(x) ~ 1, x; domain = ℂ)   # {lambertw(1, n) | n ∈ ℤ}
```
"""
function solveset(eq::Relation, x::Sym; domain = ℝ)
    _isnamed(x) || throw(ArgumentError("solveset: $x is not a symbol"))
    (domain == ℝ || domain == ℂ) || throw(ArgumentError("solveset: the domain is ℝ or ℂ, not $domain"))
    if eq.op !== :eq
        domain == ℝ || throw(ArgumentError("solveset: an inequality is solved over ℝ only"))
        return solve(eq, x)
    end
    used = Set(Symbol(string(s)) for s in union(free_symbols(eq.lhs - eq.rhs), free_parameters(eq.lhs - eq.rhs)))
    names = [s for s in _SS_NAMES if !(s in used)]
    k = Ref(0)
    fresh() = (k[] += 1; Sym(k[] <= length(names) ? names[k[]] : Symbol("n", k[])))
    e = eq.lhs - eq.rhs
    unsolved = SolutionSet(x, Sym[], ImageSet[], [ConditionSet(x, eq, domain)])
    ps = Sym[s for s in union(free_symbols(e), free_parameters(e)) if !_eq(s, x)]
    if domain == ℝ && !isempty(ps)
        # parameters over ℝ are real; the answer by cases on them (sin(x) = a: |a| ≤ 1)
        B = something(_ASSUMING[], Assumptions())
        for p in ps
            _isreal(p, B) === true || (B = _with_fact(B, p ∈ ℝ))
        end
        run(C) = with(_SS_ASSUME => C) do
            try _ss(e, x, true, fresh, Sym[], Relation[], 0) catch err; err isa _SSUnsolved || rethrow(); nothing end
        end
        conv(m) = m isa Vector ? _ss_set(x, m, names) : unsolved
        r = _by_cases(run, B)
        r isa Cases || return conv(r)
        r = _map_cases(conv, r)
        return length(ps) == 1 ? _ss_merge_cases(r, ps[1]) : r
    end
    r = try
        _ss(e, x, domain == ℝ, fresh, Sym[], Relation[], 0)
    catch err
        err isa _SSUnsolved || rethrow()
        return unsolved
    end
    r isa Cases ? _map_cases(m -> _ss_set(x, m, names), r) : _ss_set(x, r, names)
end

# cases on one parameter with the same set merged, their conditions as one interval
# (a = -1, -1 < a < 1, a = 1: -1 ≤ a ≤ 1), less the points that are cases of their own
# (a ≠ -1); conditions in p alone (a² - 2 < 0) read as sets of its values. r as it is
# when a condition is not one in p alone
function _ss_merge_cases(r::Cases, p::Sym)
    full = SymSet([SymInterval(nothing, nothing, false, false)])
    function condset(c)
        c isa Relation || return nothing
        f = c.lhs - c.rhs
        all(v -> _eq(v, p), free_symbols(f)) || return nothing
        if _eq(c.lhs, p) && !_hassymbol(c.rhs)
            b = c.rhs
            c.op === :eq && return SymSet([SymInterval(b, b, true, true)])
            c.op === :ne && return SymSet([SymInterval(nothing, b, false, false), SymInterval(b, nothing, false, false)])
            c.op in (:lt, :le) && return SymSet([SymInterval(nothing, b, false, c.op === :le)])
            c.op in (:gt, :ge) && return SymSet([SymInterval(b, nothing, c.op === :ge, false)])
        end
        S = try solve(c, p; assume = (p ∈ ℝ,)) catch err; err isa ArgumentError || rethrow(); nothing end
        S isa SymSet && return S
        if S isa AbstractVector && c.op === :eq && all(v -> !_hassymbol(v), S)
            return SymSet([SymInterval(v, v, true, true) for v in _sort_roots!(Sym[Sym(v) for v in S])])
        end
        nothing
    end
    sets = SymSet[]
    for (conds, _) in r.cases
        S = full
        for c in conds
            T = condset(c)
            T === nothing && return r
            S = S ∩ T
        end
        push!(sets, S)
    end
    keys = String[]
    groups = Dict{String,Vector{Int}}()
    for (i, (_, v)) in enumerate(r.cases)
        k = string(v)
        haskey(groups, k) || (push!(keys, k); groups[k] = Int[])
        push!(groups[k], i)
    end
    ispoint(I) = I.lo !== nothing && I.hi !== nothing && _eq(I.lo, I.hi)
    bounds(I) = ispoint(I) ? Any[Relation(:eq, p, I.lo)] :
        Any[(I.lo === nothing ? Relation[] : [Relation(I.lclosed ? :ge : :gt, p, I.lo)])...;
            (I.hi === nothing ? Relation[] : [Relation(I.rclosed ? :le : :lt, p, I.hi)])...]
    # the points that are cases by themselves
    lone = Sym[S.intervals[1].lo for S in sets if length(S.intervals) == 1 && ispoint(S.intervals[1])]
    out = Tuple{Vector{Any},Any,Float64}[]
    lowest(I) = I.lo === nothing ? -Inf : real_part(N(I.lo))
    for k in keys
        is = groups[k]
        S = reduce(union, sets[is])
        isempty(S.intervals) && continue
        Is = S.intervals
        # gaps of single points, each a case of its own: the hull, less those points
        holes = Sym[]
        ok = all(eachindex(Is)[1:end-1]) do j
            a, b = Is[j], Is[j + 1]
            g = a.hi !== nothing && b.lo !== nothing && _eq(a.hi, b.lo) && !a.rclosed && !b.lclosed &&
                any(q -> _eq(q, a.hi), lone)
            g && push!(holes, a.hi)
            g
        end
        if ok
            H = SymInterval(Is[1].lo, Is[end].hi, Is[1].lclosed, Is[end].rclosed)
            push!(out, ([bounds(H); Any[Relation(:ne, p, h) for h in holes]], r.cases[is[1]].second, lowest(H)))
        else                                   # (a < -1 and a > 1: kept apart)
            for I in Is                               # (a² > 2: a < -√2, a > √2)
                push!(out, (bounds(I), r.cases[is[1]].second, lowest(I)))
            end
        end
    end
    length(out) == 1 && isempty(out[1][1]) && return out[1][2]
    sort!(out; by = c -> c[3])                 # (in the order of the parameter's values)
    T = typejoin(map(c -> typeof(c[2]), out)...)
    Cases{T}(Pair{Vector{Any},T}[c[1] => c[2] for c in out])
end

# the assumptions of the case being solved (parameters over ℝ), or nothing
const _SS_ASSUME = ScopedValue{Union{Nothing,Assumptions}}(nothing)

# members: (expression, bound variables, conditions on them)
const _SSMember = Tuple{Sym,Vector{Sym},Vector{Relation}}

function _ss_set(x::Sym, ms::Vector{_SSMember}, names)
    pts, fams = Sym[], ImageSet[]
    for (e, vs, cs) in _ss_merge(ms)
        vs = unique(Sym[v for v in vs if has(e, v)])
        cs = Relation[c for c in cs if any(v -> has(c.lhs - c.rhs, v), vs)]
        # each family's variables renamed n, m, … in order
        if !isempty(vs)
            new = Sym[Sym(i <= length(names) ? names[i] : Symbol("n", i)) for i in eachindex(vs)]
            tmp = Sym[Sym(Symbol("_ssr", i)) for i in eachindex(vs)]
            sub1, sub2 = [v => t for (v, t) in zip(vs, tmp)], [t => w for (t, w) in zip(tmp, new)]
            e = subs(subs(e, sub1...), sub2...)
            cs = Relation[Relation(c.op, subs(subs(c.lhs, sub1...), sub2...), subs(subs(c.rhs, sub1...), sub2...)) for c in cs]
            vs = new
        end
        if isempty(vs)
            any(p -> _ss_same(p, e), pts) || push!(pts, e)
        else
            any(f -> isequal(f.expr, e) && f.vars == vs, fams) || push!(fams, ImageSet(e, vs, cs))
        end
    end
    SolutionSet(x, _sort_roots!(pts), fams, ConditionSet[])
end

_ss_same(p::Sym, q::Sym) = isequal(p, q) || (!_hassymbol(p - q) && (v = try N(p - q) catch; nothing end; v isa Number && abs(v) < 1e-12))

# families a + b n merged: those with the same step b (and the same other variables and
# conditions) whose offsets, reduced mod b, are equally spaced around it are one family
# with the smaller step ({π/2 + 2πn} ∪ {3π/2 + 2πn} is {π/2 + πn}); conditions in n alone
# written n ⋈ c
function _ss_merge(ms::Vector{_SSMember})
    ms = _SSMember[(e, vs, Relation[_ss_normcond(c, vs) for c in cs]) for (e, vs, cs) in ms]
    lin(e, n) = (a = subs(e, n => 0); b = cancel(subs(e, n => 1) - a);
        iszero(expand(e - a - b * n)) && !_hassymbol(b) ? (a, b) : nothing)
    out = _SSMember[]
    used = falses(length(ms))
    for i in eachindex(ms)
        used[i] && continue
        e, vs, cs = ms[i]
        if isempty(vs) || any(c -> has(c.lhs - c.rhs, vs[end]), cs) || (ab = lin(e, vs[end])) === nothing
            push!(out, ms[i]); continue
        end
        n, (a, b) = vs[end], ab
        bn = abs(N(b))
        grp = [(i, a)]
        for j in (i + 1):length(ms)
            used[j] && continue
            e2, vs2, cs2 = ms[j]
            (length(vs2) == length(vs) && all(isequal.(vs2[1:end-1], vs[1:end-1])) && isequal(cs2, cs)) || continue
            ab2 = lin(subs(e2, vs2[end] => n), n)
            ab2 === nothing && continue
            q = cancel(ab2[2] / b)
            (q == 1 || q == -1) || continue
            push!(grp, (j, ab2[1]))
        end
        # offsets mod b; the subsets equally spaced by b/r around a period, largest first,
        # each one family of step b/r
        offs = Sym[]
        for (_, o) in grp
            v = try N(cancel(o / b)) catch err; err isa ArgumentError || rethrow(); nothing end
            (v isa Number && !(v isa Complex && abs(imag(v)) > 1e-12)) || (offs = nothing; break)
            push!(offs, expand(o - floor(Int, real_part(v) + 1e-9) * b))
        end
        if offs === nothing || length(grp) == 1
            push!(out, ms[i]); continue
        end
        for (j, _) in grp; used[j] = true; end
        vals = [mod(real_part(N(cancel(o / b))), 1.0) for o in offs]
        left = collect(eachindex(offs))
        for r in length(offs):-1:2
            again = true
            while again
                again = false
                for s0 in left
                    idx = Int[]
                    for k in 0:(r - 1)
                        t = mod(vals[s0] + k / r, 1.0)
                        j = findfirst(q -> q in left && min(abs(vals[q] - t), 1 - abs(vals[q] - t)) < 1e-9, eachindex(offs))
                        j === nothing && break
                        push!(idx, j)
                    end
                    length(idx) == r || continue
                    m0 = idx[argmin(vals[idx])]
                    push!(out, (offs[m0] + b / r * n, vs, cs))
                    setdiff!(left, idx)
                    again = true
                    break
                end
            end
        end
        for q in left
            push!(out, (offs[q] + b * n, vs, cs))
        end
        continue
        push!(out, ms[i])
    end
    out
end

# a + b n ⋈ 0 (a, b numbers) as n ⋈' c
function _ss_normcond(c::Relation, vs::Vector{Sym})
    d = expand(c.lhs - c.rhs)
    ns = [v for v in vs if has(d, v)]
    length(ns) == 1 || return c
    n = ns[1]
    a, b = subs(d, n => 0), cancel(subs(d, n => 1) - subs(d, n => 0))
    (iszero(expand(d - a - b * n)) && !_hassymbol(a) && !_hassymbol(b) && !iszero(b)) || return c
    bv = N(b)
    bv isa Real || return c
    op = bv > 0 ? c.op : get(Dict(:lt => :gt, :le => :ge, :gt => :lt, :ge => :le), c.op, c.op)
    Relation(op, n, cancel(-a / b))
end

# the solutions of e = 0 in x; bvs, bcs the bound variables (from outer inverses) e has
function _ss(e::Sym, x::Sym, real::Bool, fresh, bvs::Vector{Sym}, bcs::Vector{Relation}, depth::Int)
    depth > 6 && throw(_SSUnsolved())
    e = cancel(e)
    has(e, x) || throw(_SSUnsolved())                 # (no x: none or every value)
    gens = unique(_ss_gens!(Sym[], e, x))
    # algebraic in x: solve's roots, all of them
    if all(g -> _eq(g, x) || _ss_algebraic(g, x), gens)
        return _ss_finite(e, x, real, bvs, bcs, depth)
    end
    # trigonometric functions of multiples of one angle linear in x: solve's roots over one
    # period, each repeated by the period
    # (not with bound variables in e: their values are enumerated or kept by inverting)
    # (nor with parameters: the roots are not numbers to reduce modulo the period)
    P = any(v -> has(e, v), bvs) || any(s -> !_eq(s, x), free_symbols(e)) ? nothing : _ss_period(e, x)
    if P !== nothing
        r = _ss_solve(e, x, real, depth)
        return _ss_map(r) do rs
            n = fresh()
            Pn = N(P)
            out = _SSMember[]
            seen = Sym[]
            for ρ in rs
                v = try N(ρ) catch; throw(_SSUnsolved()) end
                v isa Number || throw(_SSUnsolved())
                q = floor(Int, real_part(v) / Pn + 1e-9)
                ρ0 = expand(ρ - q * P)
                any(s -> _ss_same(s, ρ0), seen) && continue
                push!(seen, ρ0)
                push!(out, (ρ0 + P * n, [bvs; n], copy(bcs)))
            end
            out
        end
    end
    # with parameters, two or more trigonometric functions of one angle: through tan(u/2)
    if P === nothing && !any(v -> has(e, v), bvs) && any(s -> !_eq(s, x), free_symbols(e))
        r = _ss_halfangle(e, x, real, fresh, bvs, bcs, depth)
        r === nothing || return r
    end
    # one periodic or multivalued generator g(u): its values y, g inverted with every
    # branch, then u = v solved
    tgens = Sym[g for g in gens if !_eq(g, x) && !_ss_algebraic(g, x)]
    if length(tgens) == 1 && _ss_invertible(tgens[1], x)
        G = tgens[1]
        Y = Sym(Symbol("_ssy", objectid(G)))
        eY = subs(e, G => Y)
        has(eY, x) && return _ss_lambert(e, x, real, fresh, bvs, bcs, depth)
        ys = try
            with(_SOLVE_STRICT => true) do
                _solve_eq(eY ~ 0, Y)
            end
        catch err
            err isa ArgumentError || rethrow()
            throw(_SSUnsolved())
        end
        ys isa AbstractVector || throw(_SSUnsolved())
        u = _ss_arg(G)
        out = _SSMember[]
        for y in ys
            # (over ℝ with parameters, the values of G real: their radicands ≥ 0, by cases;
            # tan x + cot x = a has none for |a| < 2)
            if real && _SS_ASSUME[] !== nothing && !any(v -> has(Sym(y), v), bvs)
                any(R -> _ss_decide_ge(R, Sym[]) === false, unique(_ss_radicands!(Sym[], Sym(y)))) && continue
            end
            for (v, nv, nc) in _ss_merge(_ss_invert(G, Sym(y), real, fresh, bvs, bcs))
                if _eq(u, x)
                    push!(out, (v, nv, nc))
                else
                    r = _ss(u - v, x, real, fresh, nv, nc, depth + 1)
                    r isa Cases && throw(_SSUnsolved())
                    append!(out, r)
                end
            end
        end
        return out
    end
    # nothing periodic: solve's roots are all of them
    all(g -> !_ss_periodic(g, x), gens) && return _ss_finite(e, x, real, bvs, bcs, depth)
    throw(_SSUnsolved())
end

real_part(v) = v isa Complex ? real(v) : v

# x outside its exponentials (x eˣ = a, 2ˣ = 4x): solve's roots through Lambert's W, each
# from (αx + β) e^(γx + δ) = K, a linear change of variable, so every branch of W gives
# one: over ℂ the branches W_k, k ∈ ℤ (checked at a few k); over ℝ, W₀(z) for z ≥ -1/e
# and W₋₁(z) for -1/e ≤ z < 0, decided by cases for a parameter z
function _ss_lambert(e::Sym, x::Sym, real::Bool, fresh, bvs, bcs, depth)
    any(v -> has(e, v), bvs) && throw(_SSUnsolved())
    gens = unique(_ss_gens!(Sym[], e, x))
    all(g -> _eq(g, x) || (_ss_isexp(g, x) && _ispositive_const(g.args[1]) && _ispolyin(g.args[2], [x])), gens) ||
        throw(_SSUnsolved())
    params = any(s -> !_eq(s, x), free_symbols(e))
    # numbers only, over ℝ: solve's real roots are all of them (both real branches)
    if real && !params
        rs = _ss_solve(e, x, true, depth)
        all(r -> _has_node(r, t -> t.kind == FN && t.name === :lambertw), rs) || throw(_SSUnsolved())
        return _SSMember[(r, copy(bvs), copy(bcs)) for r in rs]
    end
    rs = _ss_solve(e, x, false, depth)
    rs isa Cases && throw(_SSUnsolved())
    isW(t) = t.kind == FN && t.name === :lambertw
    out = _SSMember[]
    for r in rs
        ws = Sym[]
        _has_node(r, t -> (isW(t) && (any(w -> _eq(w, t), ws) || push!(ws, t)); false))
        isempty(ws) && throw(_SSUnsolved())
        z = ws[1].args[1]
        all(w -> _eq(w.args[1], z), ws) || throw(_SSUnsolved())
        branch(k) = _replace_walk(t -> isW(t) ? (iszero(k) ? lambertw(z) : apply_function(:lambertw, Sym[z, Sym(k)])) : nothing, r)
        if real
            A = _SS_ASSUME[]
            A === nothing && throw(_SSUnsolved())
            _decide(Relation(:ge, expand(z + exp(-ONE)), ZERO), A) || continue
            push!(out, (branch(0), copy(bvs), copy(bcs)))
            _decide(Relation(:lt, z, ZERO), A) && push!(out, (branch(-1), copy(bvs), copy(bcs)))
        else
            n = fresh()
            f = branch(n)
            # (each branch a root, at a few k, where there are no parameters)
            if !params
                for k in -2:2
                    v = try _evalc(subs(e, x => subs(f, n => k))) catch; throw(_SSUnsolved()) end
                    abs(v) < 1e-8 * max(1, abs(_evalc(subs(f, n => k)))) || throw(_SSUnsolved())
                end
            end
            any(m -> isequal(m[1], f), out) || push!(out, (f, [bvs; n], copy(bcs)))
        end
    end
    out
end

_ss_map(f, r) = r isa Cases ? _map_cases(f, r) : f(r)

function _ss_solve(e::Sym, x::Sym, real::Bool, depth::Int)
    r = try
        with(_SOLVE_STRICT => true) do
            A = _SS_ASSUME[]
            real ? solve(e ~ 0, x; assume = A === nothing ? (x ∈ ℝ,) : _with_fact(A, x ∈ ℝ)) : solve(e ~ 0, x)
        end
    catch err
        err isa ArgumentError || rethrow()
        throw(_SSUnsolved())
    end
    r isa Cases && depth > 0 && throw(_SSUnsolved())
    r isa Cases && return _map_cases(v -> v isa AbstractVector ? Sym[Sym(q) for q in v] : throw(_SSUnsolved()), r)
    r isa AbstractVector || throw(_SSUnsolved())
    Sym[Sym(q) for q in r]
end

# solve's roots; over ℝ with bound variables (an inner equation x² = 2πn), the complex
# roots kept where real, by conditions on the bound variables (their even roots'
# radicands ≥ 0)
function _ss_finite(e::Sym, x::Sym, real::Bool, bvs, bcs, depth)
    if real && any(v -> has(e, v), bvs)
        rs = _ss_solve(e, x, false, depth)
        out = _SSMember[]
        for ρ in rs
            cs = copy(bcs)
            ok = true
            for R in unique(_ss_radicands!(Sym[], ρ))
                d = _ss_decide_ge(R, bvs)
                d === false && (ok = false; break)
                d isa Relation && push!(cs, d)
            end
            ok && push!(out, (ρ, copy(bvs), cs))
        end
        return out
    end
    _ss_map(_ss_solve(e, x, real, depth)) do rs
        _SSMember[(ρ, copy(bvs), copy(bcs)) for ρ in rs]
    end
end

# R ≥ 0: true or false for a constant, else the relation (R linear in one bound variable)
function _ss_decide_ge(R::Sym, bvs)
    vs = [v for v in bvs if has(R, v)]
    if isempty(vs)
        # (a parameter's: by cases on it, over ℝ)
        A = _SS_ASSUME[]
        _hassymbol(R) && A !== nothing && return _decide(Relation(:ge, expand(R), ZERO), A)
        _hassymbol(R) && throw(_SSUnsolved())
        v = N(R)
        (v isa Real || abs(imag(v)) < 1e-12) || return false
        return real_part(v) >= -1e-12
    end
    length(vs) == 1 || throw(_SSUnsolved())
    Relation(:ge, R, ZERO)
end

function _ss_radicands!(out, e::Sym)
    if e.kind == POW || e.kind == MUL
        bs, ps = e.kind == MUL ? (e.args, e.exps) : (e.args[1:1], e.args[2:2])
        for (b, p) in zip(bs, ps)
            p.kind == NUMBER && isexact(p.num) && iseven(_denominator(p.num)) && _hassymbol(b) && push!(out, b)
        end
    end
    foreach(a -> _ss_radicands!(out, a), e.args)
    out
end

# the maximal parts of e in x other than sums, products and integer powers
function _ss_gens!(out, e::Sym, x::Sym)
    has(e, x) || return out
    if _eq(e, x)
        push!(out, e)
    elseif e.kind == ADD
        foreach(a -> _ss_gens!(out, a, x), e.args)
    elseif e.kind == MUL
        for (b, p) in zip(e.args, e.exps)
            if p.kind == NUMBER && _isint(p.num)
                _ss_gens!(out, b, x)
            else
                push!(out, b^p)
            end
        end
    elseif e.kind == POW && e.args[2].kind == NUMBER && _isint(e.args[2].num)
        _ss_gens!(out, e.args[1], x)
    else
        push!(out, e)
    end
    out
end

# a radical of an algebraic expression in x
_ss_algebraic(g::Sym, x::Sym) = g.kind == POW && g.args[2].kind == NUMBER && isexact(g.args[2].num) &&
    !has(g.args[2], x) && all(h -> _eq(h, x) || _ss_algebraic(h, x), _ss_gens!(Sym[], g.args[1], x))

const _SS_TRIG = (:sin, :cos, :tan, :cot, :sec, :csc)
const _SS_HYP = (:sinh, :cosh, :tanh)

_ss_isexp(g::Sym, x::Sym) = g.kind == POW && !has(g.args[1], x) && has(g.args[2], x)
_ss_periodic(g::Sym, x::Sym) = _ss_isexp(g, x) || (g.kind == FN && g.name in (_SS_TRIG..., _SS_HYP..., :log))
# e rational in sin, cos, tan, cot, sec, csc of one angle u = k x + c (k a number), two or
# more of them, with parameters (the period path needs numbers): sin u = 2T/(1 + T²),
# cos u = (1 - T²)/(1 + T²), …, in T = tan(u/2), one generator inverted with its branches;
# and u = π + 2πn, where T is undefined, when e vanishes there (by cases on a parameter,
# solved again at its value: sin x + cos x = a at a = -1). nothing when not of this form
function _ss_halfangle(e::Sym, x::Sym, real::Bool, fresh, bvs::Vector{Sym}, bcs::Vector{Relation}, depth::Int)
    us, names = Sym[], Set{Symbol}()
    _replace_walk(e) do t
        if t.kind == FN && t.name in _SS_TRIG && length(t.args) == 1 && has(t.args[1], x)
            any(u -> _eq(u, t.args[1]), us) || push!(us, t.args[1])
            push!(names, t.name)
        end
        nothing
    end
    (length(us) == 1 && length(names) >= 2) || return nothing
    u = us[1]
    k = expand(diff(u, x))
    (has(k, x) || _hassymbol(k) || iszero(k)) && return nothing
    Z = Sym(:_ssT)
    # one of the functions alone first: tan and cot in tan u, sin and csc in sin u, cos and
    # sec in cos u, sin in even powers only in cos u (cos in sin u) (tan x + cot x = a is
    # quadratic in tan x, quartic in tan(x/2))
    rational(q) = !has(q, x) && !_has_node(q, t -> t.kind == POW && has(t.args[1], Z) && t.args[2].kind == NUMBER && !_isint(t.args[2].num)) &&
                  !_has_node(q, t -> t.kind == MUL && any(((b, w),) -> has(b, Z) && w.kind == NUMBER && !_isint(w.num), zip(t.args, t.exps)))
    for (G, rules) in ((tan(u), [tan(u) => Z, cot(u) => 1 / Z]), (sin(u), [sin(u) => Z, csc(u) => 1 / Z]),
                       (cos(u), [cos(u) => Z, sec(u) => 1 / Z]),
                       (cos(u), [cos(u) => Z, sec(u) => 1 / Z, sin(u) => sqrt(1 - Z^2), csc(u) => 1 / sqrt(1 - Z^2)]),
                       (sin(u), [sin(u) => Z, csc(u) => 1 / Z, cos(u) => sqrt(1 - Z^2), sec(u) => 1 / sqrt(1 - Z^2)]))
        eZ = cancel(expand(subs(e, rules...)))
        rational(eZ) || continue
        return _ss(subs(eZ, Z => G), x, real, fresh, bvs, bcs, depth + 1)
    end
    rules = [sin(u) => 2Z / (1 + Z^2), cos(u) => (1 - Z^2) / (1 + Z^2), tan(u) => 2Z / (1 - Z^2),
             cot(u) => (1 - Z^2) / (2Z), sec(u) => (1 + Z^2) / (1 - Z^2), csc(u) => (1 + Z^2) / (2Z)]
    eZ = cancel(subs(e, rules...))
    has(eZ, x) && return nothing
    ms = _ss(subs(eZ, Z => tan(u / 2)), x, real, fresh, bvs, bcs, depth + 1)
    # u = π: sin u = 0, cos u = -1, tan u = 0, sec u = -1 (cot, csc undefined there)
    (:cot in names || :csc in names) && return ms
    v = cancel(subs(e, sin(u) => 0, cos(u) => -1, tan(u) => 0, sec(u) => -1))
    join(a, b) = _ss_map(a) do a_; _ss_map(b) do b_; _SSMember[a_; b_] end end
    iszero(v) && return join(ms, _ss(cos(u / 2), x, real, fresh, bvs, bcs, depth + 1))
    qs = Sym[q for q in free_symbols(v) if !_eq(q, x)]
    A = _SS_ASSUME[]
    (isempty(qs) || length(qs) > 1 || A === nothing) && return ms
    q = qs[1]
    r0s = try _solve_generic(v ~ 0, q) catch err; err isa ArgumentError || rethrow(); Sym[] end
    for r0 in (r0s isa AbstractVector && !(r0s isa Cases) ? r0s : Sym[])
        has(r0, q) && continue
        _decide(Relation(:eq, q - r0, ZERO), A) === true && return _ss(subs(e, q => r0), x, real, fresh, bvs, bcs, depth + 1)
    end
    ms
end

_ss_invertible(g::Sym, x::Sym) = (_ss_isexp(g, x) && _ispositive_const(g.args[1])) ||
    (g.kind == FN && length(g.args) == 1 && g.name in (_SS_TRIG..., _SS_HYP..., :log))
_ss_arg(g::Sym) = g.kind == POW ? g.args[2] : g.args[1]

# the period of e in x when e is a function of sin, cos, tan, cot, sec, csc of arguments
# k x + c (x nowhere else), the k rational multiples of one another; nothing otherwise
function _ss_period(e::Sym, x::Sym)
    ps = Sym[]
    ok = Ref(true)
    rest = _replace_walk(e) do t
        if t.kind == FN && t.name in _SS_TRIG && length(t.args) == 1 && has(t.args[1], x)
            k = expand(diff(t.args[1], x))
            (has(k, x) || _hassymbol(k) || iszero(k)) && (ok[] = false; return t)
            push!(ps, (t.name in (:tan, :cot) ? PI : 2PI) / (N(k) < 0 ? -k : k))
            return Sym(Symbol("_ssp", length(ps)))
        end
        nothing
    end
    (ok[] && !isempty(ps) && !has(rest, x)) || return nothing
    # rational in the trigonometric functions (no x in their arguments but linearly):
    # finitely many solutions in a period; sin(sin x) has infinitely many over ℂ
    phs = Sym[Sym(Symbol("_ssp", i)) for i in eachindex(ps)]
    all(g -> any(h -> _eq(g, h), phs), unique(_ss_gens!(Sym[], rest, phs[1]))) || return nothing
    for h in phs[2:end]
        all(g -> any(q -> _eq(g, q), phs) || !has(g, h), unique(_ss_gens!(Sym[], rest, h))) || return nothing
    end
    # P = the least common multiple of the periods: P₁ times the lcm of the ratios
    rs = SymNumber[]
    for p in ps
        q = cancel(p / ps[1])
        (q.kind == NUMBER && isexact(q.num)) || return nothing
        push!(rs, q.num)
    end
    L = reduce(lcm, (_numerator(r) for r in rs); init = big(1))
    G = reduce(gcd, (_denominator(r) for r in rs); init = big(0))
    ps[1] * Sym(SymNumber(L // G))
end

# g(u) = y inverted: the values of u, with the bound variables and conditions they bring;
# over ℝ only real values (y in g's real range; with a bound variable in y, the values of
# it for which y is, enumerated when finitely many)
function _ss_invert(g::Sym, y::Sym, real::Bool, fresh, bvs::Vector{Sym}, bcs::Vector{Relation})
    out = _SSMember[]
    # y's admissible instances: (y with bound variables substituted, conditions)
    function within(lo, hi; lc = true, hc = true)
        vs = [v for v in bvs if has(y, v)]
        if isempty(vs)
            if _hassymbol(y)
                # parameters (real): the bounds decided, by cases
                A = _SS_ASSUME[]
                A === nothing && throw(_SSUnsolved())
                okl = lo === nothing || _decide(Relation(lc ? :ge : :gt, expand(y - lo), ZERO), A)
                okh = hi === nothing || _decide(Relation(hc ? :le : :lt, expand(y - hi), ZERO), A)
                return okl && okh ? [(y, copy(bcs))] : Tuple{Sym,Vector{Relation}}[]
            end
            w = N(y)
            (w isa Real || abs(imag(w)) < 1e-12) || return Tuple{Sym,Vector{Relation}}[]
            w = real_part(w)
            okl = lo === nothing || (lc ? w >= lo - 1e-12 : w > lo + 1e-12)
            okh = hi === nothing || (hc ? w <= hi + 1e-12 : w < hi - 1e-12)
            return okl && okh ? [(y, copy(bcs))] : Tuple{Sym,Vector{Relation}}[]
        end
        length(vs) == 1 || throw(_SSUnsolved())
        local nb = vs[1]
        a, b = subs(y, nb => 0), cancel(subs(y, nb => 1) - subs(y, nb => 0))
        iszero(expand(y - a - b * nb)) || throw(_SSUnsolved())
        (_hassymbol(a) || _hassymbol(b)) && throw(_SSUnsolved())
        an, bn = real_part(N(a)), real_part(N(b))
        iszero(bn) && throw(_SSUnsolved())
        lo === nothing && hi === nothing && return [(y, copy(bcs))]
        if lo !== nothing && hi !== nothing
            l, h = minmax((lo - an) / bn, (hi - an) / bn)
            ks = ceil(Int, l - 1e-9):floor(Int, h + 1e-9)
            res = Tuple{Sym,Vector{Relation}}[]
            for k in ks
                yk = subs(y, nb => k)
                wk = real_part(N(yk))
                ((lc ? wk >= lo - 1e-12 : wk > lo + 1e-12) && (hc ? wk <= hi + 1e-12 : wk < hi - 1e-12)) || continue
                # the variable fixed at k: its conditions must hold there, and are then dropped
                all(c -> !has(c.lhs - c.rhs, nb) || c(nb => k) === true, bcs) || continue
                push!(res, (yk, Relation[c for c in bcs if !has(c.lhs - c.rhs, nb)]))
            end
            return res
        end
        # one side: a condition on the variable
        c = lo !== nothing ? Relation(lc ? :ge : :gt, y, Sym(lo)) : Relation(hc ? :le : :lt, y, Sym(hi))
        [(y, [bcs; c])]
    end
    mk(v, cs, extra = Sym[]) = push!(out, (v, [bvs; extra], cs))
    if g.kind == POW                                       # b^u (b > 0, constant)
        b = g.args[1]
        lb = _eq(b, E) ? ONE : log(b)
        if real
            for (yk, cs) in within(0, nothing; lc = false)
                mk(log(yk) / lb, cs)
            end
        else
            iszero(y) && return out
            n = fresh()
            mk(log(y) / lb + 2PI * IM * n / lb, copy(bcs), [n])
        end
        return out
    end
    f = g.name
    if f in (:sec, :csc)
        iszero(y) && return out
        return _ss_invert(apply_function(f === :sec ? :cos : :sin, g.args), 1 / y, real, fresh, bvs, bcs)
    end
    if f === :log
        if real
            for (yk, cs) in within(nothing, nothing)
                mk(exp(yk), cs)
            end
        else
            w = N(y)
            (w isa Number) || throw(_SSUnsolved())
            im_ = w isa Complex ? imag(w) : 0.0
            (-π < im_ <= π + 1e-12) && mk(exp(y), copy(bcs))
        end
        return out
    end
    n = fresh()
    if real
        if f === :sin || f === :cos
            for (yk, cs) in within(-1, 1)
                a = f === :sin ? asin(yk) : acos(yk)
                mk(a + 2PI * n, cs, [n])
                (yk == 1 || yk == -1) && continue
                mk(f === :sin ? PI - a + 2PI * n : -a + 2PI * n, cs, [n])
            end
        elseif f === :tan || f === :cot
            for (yk, cs) in within(nothing, nothing)
                mk((f === :tan ? atan(yk) : acot(yk)) + PI * n, cs, [n])
            end
        elseif f === :sinh
            for (yk, cs) in within(nothing, nothing)
                mk(asinh(yk), cs)
            end
        elseif f === :cosh
            for (yk, cs) in within(1, nothing)
                mk(acosh(yk), cs)
                yk == 1 || mk(-acosh(yk), cs)
            end
        elseif f === :tanh
            for (yk, cs) in within(-1, 1; lc = false, hc = false)
                mk(atanh(yk), cs)
            end
        end
    else
        c = copy(bcs)
        if f === :sin
            mk(asin(y) + 2PI * n, c, [n]); (y == 1 || y == -1) || mk(PI - asin(y) + 2PI * n, c, [n])
        elseif f === :cos
            mk(acos(y) + 2PI * n, c, [n]); (y == 1 || y == -1) || mk(-acos(y) + 2PI * n, c, [n])
        elseif f === :tan
            (y == IM || y == -IM) || mk(atan(y) + PI * n, c, [n])
        elseif f === :cot
            (y == IM || y == -IM) || mk(acot(y) + PI * n, c, [n])
        elseif f === :sinh
            mk(asinh(y) + 2PI * IM * n, c, [n]); mk(PI * IM - asinh(y) + 2PI * IM * n, c, [n])
        elseif f === :cosh
            mk(acosh(y) + 2PI * IM * n, c, [n]); (y == 1 || y == -1) || mk(-acosh(y) + 2PI * IM * n, c, [n])
        elseif f === :tanh
            (y == 1 || y == -1) || mk(atanh(y) + PI * IM * n, c, [n])
        end
    end
    out
end

## the set

# a case's set at values of the parameters (cs(a => 2) of a Cases from solveset)
_ssrel(c::Relation, vals) = Relation(c.op, subs(c.lhs, vals...), subs(c.rhs, vals...))
_subsvals(f::ImageSet, vals) = ImageSet(subs(f.expr, vals...), f.vars, Relation[_ssrel(c, vals) for c in f.conds])
_subsvals(c::ConditionSet, vals) = ConditionSet(c.x, _ssrel(c.eq, vals), c.domain)
_subsvals(S::SolutionSet, vals) = SolutionSet(S.x, _sort_roots!(Sym[subs(p, vals...) for p in S.points]),
    ImageSet[_subsvals(f, vals) for f in S.families], ConditionSet[_subsvals(c, vals) for c in S.unsolved])

Base.isempty(S::SolutionSet) = isempty(S.points) && isempty(S.families) && isempty(S.unsolved)

function Base.in(v, S::SolutionSet)
    w = Sym(v)
    any(p -> _ss_same(p, w), S.points) && return true
    for f in S.families
        length(f.vars) == 1 || throw(ArgumentError("membership in a family of more than one variable is not supported"))
        n = only(f.vars)
        any(kr -> all(c -> c(n => kr) === true, f.conds), _ss_indices(f.expr, n, w)) && return true
    end
    isempty(S.unsolved) || throw(ArgumentError("membership in a ConditionSet is not decided"))
    false
end

# the integers k with expr(k) = w: from a + b n directly; else n solved for (√(πn) = w),
# or, n only as Lambert W's branch, W solved for and its branch found numerically
# (lambertw(1, n) ∋ W₋₁(1)); an error when none of these applies
function _ss_indices(expr::Sym, n::Sym, w::Sym)
    isint(k) = (k isa Complex && abs(imag(k)) < 1e-9 && (k = real(k)); k isa Real && isfinite(k) && abs(k - round(k)) < 1e-9)
    near(p, q) = (d = try N(p - q) catch; nothing end; d isa Number && abs(d) < 1e-9 * max(1, abs(N(q))))
    a, b = subs(expr, n => 0), cancel(subs(expr, n => 1) - subs(expr, n => 0))
    if iszero(expand(expr - a - b * n))
        k = N((w - a) / b)
        return isint(k) ? [round(Int, real(k))] : Int[]
    end
    Ls = Sym[]
    _replace_walk(expr) do t
        t.kind == FN && t.name === :lambertw && length(t.args) == 2 && has(t.args[2], n) && !any(l -> _eq(l, t), Ls) && push!(Ls, t)
        nothing
    end
    if length(Ls) == 1 && !has(Ls[1].args[1], n)
        L = Ls[1]
        U = Sym(:_ssW)
        h = subs(expr, L => U)
        has(h, n) && throw(ArgumentError("membership in $(expr) is not supported"))
        us = try _solve_generic(h ~ w, U) catch err; err isa ArgumentError || rethrow(); throw(ArgumentError("membership in $(expr) is not supported")) end
        out = Int[]
        for uv in us
            un = try N(uv) catch; continue end
            k0 = round(Int, imag(complex(un)) / (2π))
            for k in (k0 - 2):(k0 + 2)
                near(apply_function(:lambertw, Sym[L.args[1], Sym(k)]), uv) || continue
                for j in _ss_indices(L.args[2], n, Sym(k))
                    j in out || push!(out, j)
                end
            end
        end
        return out
    end
    ks = try _solve_generic(expr ~ w, n) catch err; err isa ArgumentError || rethrow(); throw(ArgumentError("membership in $(expr) is not supported")) end
    (ks isa AbstractVector && !(ks isa Cases)) || throw(ArgumentError("membership in $(expr) is not supported"))
    out = Int[]
    for kv in ks
        k = try N(kv) catch; continue end
        isint(k) || continue
        j = round(Int, real(k))
        near(subs(expr, n => j), w) && !(j in out) && push!(out, j)
    end
    out
end

# the members in an interval set: a family linear in its one variable, cut to the
# finitely many in a bounded set
function Base.intersect(S::SolutionSet, T::SymSet)
    isempty(S.unsolved) || throw(ArgumentError("the intersection of a ConditionSet with an interval is not supported"))
    pts = Sym[p for p in S.points if p in T]
    for f in S.families, I in T.intervals
        (I.lo === nothing || I.hi === nothing) && throw(ArgumentError("a family cut to an unbounded interval is not finite"))
        length(f.vars) == 1 || throw(ArgumentError("a family of more than one variable cut to an interval is not supported"))
        n = only(f.vars)
        a, b = subs(f.expr, n => 0), cancel(subs(f.expr, n => 1) - subs(f.expr, n => 0))
        (iszero(expand(f.expr - a - b * n)) && !_hassymbol(a) && !_hassymbol(b)) ||
            throw(ArgumentError("a family not linear in its variable cut to an interval is not supported"))
        an, bn = real_part(N(a)), real_part(N(b))
        l, h = minmax((_rv(I.lo) - an) / bn, (_rv(I.hi) - an) / bn)
        for k in (ceil(Int, l) - 1):(floor(Int, h) + 1)
            all(c -> c(n => k) === true, f.conds) || continue
            q = expand(subs(f.expr, n => k))
            q in SymSet([I]) && !any(p -> _ss_same(p, q), pts) && push!(pts, q)
        end
    end
    SolutionSet(S.x, _sort_roots!(pts), ImageSet[], ConditionSet[])
end
Base.intersect(T::SymSet, S::SolutionSet) = S ∩ T

function Base.collect(S::SolutionSet)
    (isempty(S.families) && isempty(S.unsolved)) || throw(ArgumentError("the set is not finite (or not known to be)"))
    copy(S.points)
end

function Base.show(io::IO, f::ImageSet)
    print(io, "{", f.expr, " | ", join(string.(f.vars), ", "), " ∈ ℤ")
    for c in f.conds
        print(io, ", ", _condstring(c))
    end
    print(io, "}")
end
Base.show(io::IO, c::ConditionSet) = print(io, "{", c.x, " ∈ ", c.domain == ℝ ? "ℝ" : "ℂ", " | ", c.eq, "}")
function Base.show(io::IO, S::SolutionSet)
    parts = String[]
    isempty(S.points) || push!(parts, "{" * join(string.(S.points), ", ") * "}")
    append!(parts, string.(S.families))
    append!(parts, string.(S.unsolved))
    print(io, isempty(parts) ? "∅" : join(parts, " ∪ "))
end
Base.show(io::IO, ::MIME"text/plain", S::SolutionSet) = show(io, S)
