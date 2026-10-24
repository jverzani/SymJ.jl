# Subsets of the real line as finite unions of disjoint intervals with exact endpoints:
# the solutions of inequalities in one variable. Not exported (SymJ.SymSet).
#
# Every set is built the same way: sorted breakpoints split the line into open pieces
# and the breakpoints themselves; a predicate decides each piece (on an open piece by
# a sample, since nothing changes inside it); consecutive included pieces merge. That
# serves solving (the predicate is the inequality) and ∪, ∩, setdiff (membership in
# the operands). Endpoints must evaluate to real numbers.

"""
    SymInterval

An interval of the real line with exact endpoints: `lo`, `hi` (`nothing` for -∞, ∞) and
whether each end is closed. A point is the closed interval `[a, a]`.
"""
struct SymInterval
    lo::Union{Nothing,Sym}
    hi::Union{Nothing,Sym}
    lclosed::Bool
    rclosed::Bool
end

"""
    SymSet

A finite union of disjoint intervals with exact endpoints, sorted and merged; the
solutions of an inequality ([`solve`](@ref)). Supports `in`, `∪`, `∩`, `setdiff`,
`isempty`; the intervals are `S.intervals`. Calling tests membership: `S(2)`, or
`S(x => 2)` for a set from `solve(ineq, x)` (which records its variable).
"""
struct SymSet
    intervals::Vector{SymInterval}
    var::Union{Nothing,Sym}           # the variable it was solved for, if known
end
SymSet(ints::Vector{SymInterval}) = SymSet(ints, nothing)

(S::SymSet)(v::Union{Sym,Real}) = v in S
function (S::SymSet)(p::Pair)
    S.var !== nothing && _eq(Sym(p.first), S.var) ||
        throw(ArgumentError("this set is not known to be in $(p.first); call it with a value, S(v)"))
    p.second in S
end

_rv(e::Sym) = (v = N(e); v isa Real ? Float64(v) : throw(ArgumentError("$e is not a real number")))

Base.isempty(S::SymSet) = isempty(S.intervals)
Base.:(==)(A::SymSet, B::SymSet) = length(A.intervals) == length(B.intervals) &&
    all(_sameint(a, b) for (a, b) in zip(A.intervals, B.intervals))
_sameint(a::SymInterval, b::SymInterval) = a.lclosed == b.lclosed && a.rclosed == b.rclosed &&
    _sameend(a.lo, b.lo) && _sameend(a.hi, b.hi)
_sameend(a, b) = (a === nothing && b === nothing) || (a !== nothing && b !== nothing && a == b)

## membership

function _in_interval(v, I::SymInterval)          # v: a Sym (exact) or a Float64
    function side(e, closed, below)               # v on the right side of endpoint e
        e === nothing && return true
        v isa Sym && v == e && return closed
        x, y = v isa Sym ? _rv(v) : v, _rv(e)
        x == y ? closed : (below ? x > y : x < y)
    end
    side(I.lo, I.lclosed, true) && side(I.hi, I.rclosed, false)
end
Base.in(v::Union{Sym,Real}, S::SymSet) = (w = v isa Real ? Float64(v) : v; any(I -> _in_interval(w, I), S.intervals))

## construction from breakpoints and a predicate

# pred(v, point::Bool): is the piece containing v included? (v exact at a breakpoint,
# a Float64 sample inside an open piece)
function _set_from(bps::Vector{Sym}, pred)
    vals = map(_rv, bps)
    p = sortperm(vals)
    bps, vals = bps[p], vals[p]
    keep = [i == 1 || vals[i] != vals[i-1] for i in eachindex(vals)]   # distinct values
    bps, vals = bps[keep], vals[keep]
    m = length(bps)
    # pieces in order: open (−∞, b1), {b1}, (b1, b2), ..., {bm}, (bm, ∞)
    samples = m == 0 ? [0.0] :
        vcat([vals[1] - 1.0], [(vals[i] + vals[i+1]) / 2 for i in 1:(m-1)], [vals[m] + 1.0])
    incl_open = [pred(s, false) for s in samples]
    incl_pt = [pred(b, true) for b in bps]
    out = SymInterval[]
    lo, lclosed = nothing, false                # the interval being built, if open
    building = false
    for i in 0:m
        # open piece between bps[i] and bps[i+1]
        a = i == 0 ? nothing : bps[i]
        if incl_open[i+1]
            building || (building = true; lo = a; lclosed = false)
        elseif building
            push!(out, SymInterval(lo, a, lclosed, i > 0 && incl_pt[i]))
            building = false
        end
        i == m && break
        # the breakpoint bps[i+1]
        b = bps[i+1]
        if incl_pt[i+1]
            building || (building = true; lo = b; lclosed = true)
        elseif building
            push!(out, SymInterval(lo, b, lclosed, false))
            building = false
        end
    end
    building && push!(out, SymInterval(lo, nothing, lclosed, false))
    SymSet(out)
end

_endpoints(S::SymSet) = Sym[e for I in S.intervals for e in (I.lo, I.hi) if e !== nothing]

Base.union(A::SymSet, B::SymSet) = _set_from(vcat(_endpoints(A), _endpoints(B)), (v, _) -> v in A || v in B)
Base.intersect(A::SymSet, B::SymSet) = _set_from(vcat(_endpoints(A), _endpoints(B)), (v, _) -> v in A && v in B)
Base.setdiff(A::SymSet, B::SymSet) = _set_from(vcat(_endpoints(A), _endpoints(B)), (v, _) -> v in A && !(v in B))

## printing: (-∞, -1) ∪ [2, 3] ∪ {5}

function Base.show(io::IO, I::SymInterval)
    if I.lo !== nothing && I.hi !== nothing && I.lo == I.hi
        print(io, "{", I.lo, "}")
        return
    end
    print(io, I.lclosed ? "[" : "(", I.lo === nothing ? "-∞" : I.lo, ", ",
          I.hi === nothing ? "∞" : I.hi, I.rclosed ? "]" : ")")
end
function Base.show(io::IO, S::SymSet)
    isempty(S) && return print(io, "∅")
    join(io, S.intervals, " ∪ ")
end
Base.show(io::IO, ::MIME"text/plain", S::SymSet) = show(io, S)

## solving inequalities

# the real values among roots (exact)
_real_roots(rs) = filter(r -> (v = try _evalc(r) catch; return false end; abs(imag(v)) <= 1e-10 * max(1, abs(v))), rs)

function _solve_inequality(r::Relation, x::Sym)
    n, d = numer_denom(cancel(r.lhs - r.rhs))
    for e in (n, d), s in free_symbols(e) ∪ free_parameters(e)
        _eq(s, x) || throw(ArgumentError("solve: inequalities with other symbols or parameters are not supported yet"))
    end
    atoms = Set{Sym}()
    (_poly_atoms!(atoms, n) && _poly_atoms!(atoms, d)) ||
        throw(ArgumentError("solve: an inequality with inexact coefficients is not supported yet"))
    all(a -> _eq(a, x) || !has(a, x), atoms) ||
        throw(ArgumentError("solve: only polynomial and rational inequalities are supported yet"))
    zs = has(n, x) ? _real_roots(_solve_eq(n ~ 0, x)) : Sym[]
    ps = has(d, x) ? _real_roots(_solve_eq(d ~ 0, x)) : Sym[]
    pole(v) = any(p -> p == v, ps)
    _set_from(vcat(zs, ps), (v, atpoint) -> begin
        atpoint && pole(v) && return false
        w = v isa Sym ? v : Sym(v)
        val = r(x => w)
        val isa Bool ? val : false
    end) |> S -> SymSet(S.intervals, x)
end
