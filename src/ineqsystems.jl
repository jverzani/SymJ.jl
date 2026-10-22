# Systems of inequalities.
#
# One unknown: each relation is solved on its own (equations contribute their real
# roots as points) and the sets are intersected.
#
# Several unknowns, linear: the equations are solved first and substituted; Fourier–
# Motzkin elimination (exact) then removes the unknowns in the given order, recording
# the bounds of each in terms of the ones still to come. The result is a cylindrical
# description: the last unknown in an interval, each earlier one between bounds in the
# later ones. Nonlinear systems in several unknowns (cylindrical algebraic
# decomposition) are not supported.

"""
    CylindricalSet

The solution of a system of linear inequalities in several unknowns, as nested bounds
(a cylindrical description): the last unknown ranges over an interval, and each
earlier one lies between bounds in the later ones; unknowns fixed by equations are
given last. Not exported. Supports `in` for a point given as a `Dict`.
"""
struct CylindricalSet
    feasible::Bool
    levels::Vector{Tuple{Sym,Vector{Tuple{Sym,Bool}},Vector{Tuple{Sym,Bool}}}}  # (x, lowers, uppers); (bound, strict)
    equalities::Vector{Pair{Sym,Sym}}
end

Base.isempty(C::CylindricalSet) = !C.feasible

# a·v + k > 0 (strict) or ≥ 0, coefficients over the unknowns
struct _LinCon
    a::Vector{SymNumber}
    k::SymNumber
    strict::Bool
end

function _solve_ineq_system(eqs::AbstractVector{<:Relation}, vars::Vector{Sym})
    length(vars) == 1 && return _solve_ineq_1(eqs, only(vars))
    any(e -> e.op === :ne, eqs) && throw(ArgumentError(
        "solve: ≠ in a system of several unknowns is not supported (its solution is not convex)"))
    # equations first: pivots in terms of the free unknowns, substituted
    eqns = [e for e in eqs if e.op === :eq]
    ineqs = [e for e in eqs if e.op !== :eq]
    equalities = Pair{Sym,Sym}[]
    free = vars
    if !isempty(eqns)
        sol = solve(eqns, vars)
        sol === nothing && return CylindricalSet(false, [], [])
        equalities = [v => sol[v] for v in vars if haskey(sol, v)]
        free = [v for v in vars if !haskey(sol, v)]
        ineqs = [Relation(r.op, subs(r.lhs, sol), subs(r.rhs, sol)) for r in ineqs]
    end
    cons = [_lincon(r, free) for r in ineqs]
    _fourier_motzkin(cons, free, equalities)
end

function _solve_ineq_1(eqs, x::Sym)
    S = SymSet([SymInterval(nothing, nothing, false, false)])
    for r in eqs
        T = if r.op === :eq
            pts = _real_roots(solve(r, x))
            _set_from(pts, (v, atpoint) -> atpoint && any(p -> p == v, pts))
        else
            solve(r, x)
        end
        S = S ∩ T
    end
    S
end

# the relation as a·v + k > 0 or ≥ 0, with numeric coefficients
function _lincon(r::Relation, vars::Vector{Sym})
    A, b = _linear_system([r], vars)                       # lhs - rhs = A v - b
    all(e -> e.kind == NUMBER, Iterators.flatten((A, b))) || throw(ArgumentError(
        "solve: inequalities with symbolic coefficients in several unknowns are not supported yet"))
    a, k = [e.num for e in A[1, :]], -b[1].num
    # lhs - rhs = a·v + k; Lt, Le: -(lhs - rhs) > 0, ≥ 0; Gt, Ge: lhs - rhs > 0, ≥ 0
    s = r.op in (:lt, :le) ? NM1 : N1
    _LinCon(a .* s, k * s, r.op in (:lt, :gt))
end

# positive scaling to a canonical form: first nonzero coefficient ±1
function _normalize(c::_LinCon)
    i = findfirst(!iszero, c.a)
    i === nothing && return c
    s = inv(abs(c.a[i]))
    _LinCon(c.a .* s, c.k * s, c.strict)
end

# same direction: keep the tighter (smaller k; strict on a tie); constants decided
function _prune(cons::Vector{_LinCon})
    best = Dict{Vector{SymNumber},_LinCon}()
    for c0 in cons
        c = _normalize(c0)
        if all(iszero, c.a)
            (c.strict ? c.k > N0 : c.k >= N0) || return nothing    # infeasible
            continue
        end
        d = get(best, c.a, nothing)
        if d === nothing || c.k < d.k || (c.k == d.k && c.strict && !d.strict)
            best[c.a] = c
        end
    end
    collect(values(best))
end

function _fourier_motzkin(cons::Vector{_LinCon}, vars::Vector{Sym}, equalities)
    levels = Tuple{Sym,Vector{Tuple{Sym,Bool}},Vector{Tuple{Sym,Bool}}}[]
    cs = _prune(cons)
    cs === nothing && return CylindricalSet(false, levels, equalities)
    for j in eachindex(vars)
        pos = [c for c in cs if c.a[j] > N0]
        neg = [c for c in cs if c.a[j] < N0]
        rest = [c for c in cs if iszero(c.a[j])]
        # bound: v_j > (or ≥) -(k + Σ_{i≠j} a_i v_i)/a_j for a_j > 0; the reverse for a_j < 0
        bound(c) = (-(Sym(c.k) + _dot(c.a, vars, j)) / Sym(c.a[j]), c.strict)
        lowers = _dominant([bound(c) for c in pos], true)
        uppers = _dominant([bound(c) for c in neg], false)
        push!(levels, (vars[j], lowers, uppers))
        for p in pos, n in neg
            sp, sn = inv(p.a[j]), inv(-n.a[j])
            push!(rest, _LinCon(p.a .* sp .+ n.a .* sn, p.k * sp + n.k * sn, p.strict || n.strict))
        end
        cs = _prune(rest)
        cs === nothing && return CylindricalSet(false, levels, equalities)
    end
    CylindricalSet(true, levels, equalities)
end

_dot(a, vars, skip) = _add_many(Sym(a[i]) * vars[i] for i in eachindex(vars) if i != skip && !iszero(a[i]))

# of bounds differing by a constant, keep the tightest (lower: largest)
function _dominant(bs::Vector{Tuple{Sym,Bool}}, lower::Bool)
    out = Tuple{Sym,Bool}[]
    for (e, strict) in bs
        i = findfirst(o -> (d = cancel(e - o[1]); d.kind == NUMBER), out)
        if i === nothing
            push!(out, (e, strict))
        else
            d = cancel(e - out[i][1]).num
            tighter = lower ? d > N0 : d < N0
            (tighter || (iszero(d) && strict)) && (out[i] = (e, strict))
        end
    end
    out
end

## membership and printing

function Base.in(pt::AbstractDict, C::CylindricalSet)
    C.feasible || return false
    d = Dict{Sym,Sym}(Sym(k) => Sym(v) for (k, v) in pt)
    ok(e, v, strict, lower) = (w = N(subs(e, d)); lower ? (strict ? v > w : v >= w) : (strict ? v < w : v <= w))
    for (x, lo, hi) in C.levels
        v = N(d[x])
        all(b -> ok(b[1], v, b[2], true), lo) && all(b -> ok(b[1], v, b[2], false), hi) || return false
    end
    all(eq -> abs(N(d[eq.first]) - N(subs(eq.second, d))) < 1e-12, C.equalities)
end

function _show_level(io, x, lo, hi)
    uniform(bs) = isempty(bs) || all(b -> b[2] == bs[1][2], bs)
    if uniform(lo) && uniform(hi)
        side(bs, lower) = isempty(bs) ? (lower ? "-∞" : "∞") :
            length(bs) == 1 ? string(bs[1][1]) : string(lower ? max(first.(bs)...) : min(first.(bs)...))
        print(io, x, " ∈ ", isempty(lo) || lo[1][2] ? "(" : "[", side(lo, true), ", ",
              side(hi, false), isempty(hi) || hi[1][2] ? ")" : "]")
    else
        parts = vcat([string(b[1], b[2] ? " < " : " ≤ ", x) for b in lo],
                     [string(x, b[2] ? " < " : " ≤ ", b[1]) for b in hi])
        join(io, parts, ", ")
    end
end

function Base.show(io::IO, C::CylindricalSet)
    C.feasible || return print(io, "∅")
    parts = String[]
    for (x, lo, hi) in reverse(C.levels)
        push!(parts, sprint(_show_level, x, lo, hi))
    end
    for (v, e) in C.equalities
        push!(parts, string(v, " = ", e))
    end
    join(io, parts, ", ")
end
Base.show(io::IO, ::MIME"text/plain", C::CylindricalSet) = show(io, C)
