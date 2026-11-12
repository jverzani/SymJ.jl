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
solutions of an inequality ([`solve`](@ref)). Supports `in`, `∪`, `∩`, `setdiff`, `⊆`,
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

# (a real number written with complex parts, a root through ω, has an imaginary part of
# rounding size only)
_rv(e::Sym) = (v = N(e); v isa Real ? Float64(v) :
    abs(imag(v)) <= 1e-12 * max(1, abs(v)) ? Float64(real(v)) : throw(ArgumentError("$e is not a real number")))

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
Base.issubset(A::SymSet, B::SymSet) = isempty(setdiff(A, B))
Base.issubset(A::Union{SymSet,SymInterval}, B::Union{SymSet,SymInterval}) = _asset(A) ⊆ _asset(B)
_asset(S::SymSet) = S
_asset(I::SymInterval) = SymSet([I])

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
    S = _inequality_with_infinity(r, x)
    S === nothing || return S
    e = r.lhs - r.rhs
    # a piecewise function of x: piece by piece, each on its condition's set less the
    # earlier ones' (piecewise(x > 10 => 1, x ≤ 10 => 0) > 0 is x > 10; SymPy #19718)
    if _has_fn(e, (:piecewise, :abs, :sign, :max, :min, :heaviside))
        S = _solve_inequality_pieces(r, x)
        S === nothing || return S
    end
    n, d = numer_denom(cancel(e))
    for e in (n, d), s in free_symbols(e) ∪ free_parameters(e)
        _eq(s, x) || throw(ArgumentError("solve: inequalities with other symbols or parameters are not supported yet"))
    end
    atoms = Set{Sym}()
    (_poly_atoms!(atoms, n) && _poly_atoms!(atoms, d)) ||
        throw(ArgumentError("solve: an inequality with inexact coefficients is not supported yet"))
    if !all(a -> _eq(a, x) || !has(a, x), atoms)
        S = _solve_inequality_continuous(r, x, atoms)
        S === nothing && throw(ArgumentError("solve: only polynomial and rational inequalities are supported yet"))
        return S
    end
    # complex coefficients: n/d = (P + iQ)/|d|² for real x, with P, Q real polynomials; an
    # order holds only where the value is real, Q = 0 (sampled as real, (x + i)(x + 2i) < 0
    # held on the whole line; it holds at x = 0 alone; SymPy #13105)
    if _complex_coefs(n, x) || _complex_coefs(d, x)
        nr, ni = _reim_poly(n, x)
        dr, di = _reim_poly(d, x)
        P, Q, D2 = expand(nr * dr + ni * di), expand(ni * dr - nr * di), expand(dr^2 + di^2)
        if !iszero(Q)
            pts = has(Q, x) ? _real_roots(_solve_eq(Q ~ 0, x)) : Sym[]
            keep = Sym[p for p in pts if !iszero(simplify(subs(D2, x => p))) &&
                       (val = Relation(r.op, P / D2, ZERO)(x => p); val isa Bool && val)]
            return SymSet(SymInterval[SymInterval(p, p, true, true) for p in sort(unique(keep); by = _rv)], x)
        end
        n, d = P, D2
    end
    zs = has(n, x) ? _real_roots(_solve_eq(n ~ 0, x)) : Sym[]
    ps = has(d, x) ? _real_roots(_solve_eq(d ~ 0, x)) : Sym[]
    # the zeros of the denominators cancel removed (x = 1 in (x^2 - 1)/(x - 1)) are
    # excluded too
    for b in unique(_zero_sets!(Sym[], e, x)), h in _real_roots(_solve_eq(b ~ 0, x))
        any(p -> p == h, ps) || push!(ps, h)
    end
    pole(v) = any(p -> p == v, ps)
    rc = Relation(r.op, n / d, ZERO)            # what is evaluated: free of what cancel removed
    _set_from(vcat(zs, ps), (v, atpoint) -> begin
        atpoint && pole(v) && return false
        w = v isa Sym ? v : Sym(v)
        val = rc(x => w)
        val isa Bool ? val : false
    end) |> S -> SymSet(S.intervals, x)
end

# a polynomial in x with a non-real coefficient; its real and imaginary parts for real x
_complex_coefs(p::Sym, x::Sym) =
    has(p, x) ? any(k -> _isreal(coeff(p, x, k), Assumptions()) === false, 0:degree(Poly{Sym}(p, [x]), x)) :
    _isreal(p, Assumptions()) === false
function _reim_poly(p::Sym, x::Sym)
    has(p, x) || return (real(p), imag(p))
    cs = [coeff(p, x, k) for k in 0:degree(Poly{Sym}(p, [x]), x)]
    (sum(real(c) * x^(k - 1) for (k, c) in enumerate(cs)), sum(imag(c) * x^(k - 1) for (k, c) in enumerate(cs)))
end

# f ⋈ ±∞ with f finite: f < ∞ wherever f is defined, f < -∞ nowhere, and so on (x < -Inf
# raised DivisionByZero; SymPy #16582)
function _inequality_with_infinity(r::Relation, x::Sym)
    _isinfpt(r.lhs) == _isinfpt(r.rhs) && return nothing
    inf, f, op = _isinfpt(r.rhs) ? (r.rhs, r.lhs, r.op) :
        (r.lhs, r.rhs, get(Dict(:lt => :gt, :le => :ge, :gt => :lt, :ge => :le), r.op, r.op))
    op in (:lt, :le, :gt, :ge) || return nothing
    (_has_node(f, _isinfpt) || !has(f, x)) && return nothing
    holds = (op in (:lt, :le)) == (Float64(inf.num) > 0)
    holds || return _set_from(Sym[], (v, atpoint) -> false) |> S -> SymSet(S.intervals, x)
    d = denominator(cancel(f))
    ps = has(d, x) ? _real_roots(_solve_eq(d ~ 0, x)) : Sym[]
    _set_from(ps, (v, atpoint) -> !(atpoint && any(p -> p == v, ps))) |> S -> SymSet(S.intervals, x)
end

function _solve_inequality_pieces(r::Relation, x::Sym)
    e = r.lhs - r.rhs
    _hasfn(e, :piecewise) || (e = piecewise(e))           # |x|, sign, max, min as pieces
    pw = try piecewise_fold(e, x) catch err; err isa ArgumentError || rethrow(); return nothing end
    (pw.kind == FN && pw.name === :piecewise) || return nothing
    reals = SymSet([SymInterval(nothing, nothing, false, false)])
    function cond_set(c::Sym)
        c == ONE && return reals
        c == ZERO && return SymSet(SymInterval[])
        if c.kind == FN && haskey(_COND_RELS, c.name)
            q = _relation_of(c)
            q.op === :eq && return SymSet(_set_from(_real_roots(_solve_eq(q.lhs - q.rhs ~ 0, x)), (v, atpoint) -> atpoint).intervals)
            q.op === :ne && return setdiff(reals, cond_set(_mkfn_node(:eq, Sym[q.lhs, q.rhs])))
            return _solve_inequality(Relation(q.op, q.lhs, q.rhs), x)
        elseif c.kind == FN && c.name in (:And, :Or)
            ss = map(cond_set, c.args)
            return reduce(c.name === :And ? intersect : union, ss)
        elseif c.kind == FN && c.name === :Not
            return setdiff(reals, cond_set(only(c.args)))
        end
        throw(ArgumentError("solve: the condition $c is not supported"))
    end
    covered = SymSet(SymInterval[])
    out = SymSet(SymInterval[])
    try
        for (c, v) in _pieces(pw)
            Sc = setdiff(cond_set(c), covered)
            covered = union(covered, Sc)
            isempty(Sc) && continue
            Sv = if has(v, x)
                _solve_inequality(Relation(r.op, v, ZERO), x)
            else
                d = _decide(Relation(r.op, v, ZERO))
                d === true ? reals : d === false ? SymSet(SymInterval[]) : return nothing
            end
            out = union(out, intersect(Sc, Sv))
        end
    catch err
        err isa ArgumentError || rethrow()
        return nothing
    end
    SymSet(out.intervals, x)
end

# an inequality in radicals, logs and exponentials of x (no periodic functions, so the
# equation has finitely many roots): its breakpoints the real roots of the equation, the
# ends of the real domain (zeros of radicands and of log's arguments) and the poles; each
# piece decided at a point, outside the real domain excluded (x/√(1 - x²) > 2; SymPy
# #21173). nothing when not of that kind or the roots are not found.
# (atan and acot too, of such arguments in turn, acot(u) as atan(1/u), with a jump where
# u = 0: |acot x| on (-1, 1), atan(log x) > 0; SymPy #20755, #19977)
function _solve_inequality_continuous(r::Relation, x::Sym, atoms)
    acot2atan(t) = _replace_walk(t) do q
        q.kind == FN && q.name === :acot && length(q.args) == 1 ? atan(1 / acot2atan(q.args[1])) : nothing
    end
    lhs, rhs = acot2atan(r.lhs), acot2atan(r.rhs)
    e = lhs - rhs
    function okexpr(u)
        n, d = numer_denom(u)
        as = Set{Sym}()
        (_poly_atoms!(as, n) && _poly_atoms!(as, d)) || return false
        all(ok, as)
    end
    ok(a) = _eq(a, x) || !has(a, x) ||
        (a.kind == POW && a.args[2].kind == NUMBER && isexact(a.args[2].num) && _ispolyin(a.args[1], [x])) ||
        (a.kind == FN && a.name in (:log, :exp, :cbrt) && length(a.args) == 1 && _ispolyin(a.args[1], [x])) ||
        (a.kind == FN && a.name in (:log, :exp, :cbrt, :atan) && length(a.args) == 1 && okexpr(a.args[1]))
    n0, d0 = numer_denom(cancel(e))
    as0 = Set{Sym}()
    (_poly_atoms!(as0, n0) && _poly_atoms!(as0, d0) && all(ok, as0)) || return nothing
    rts = try _solve_eq(Relation(:eq, lhs, rhs), x) catch err; err isa ArgumentError || rethrow(); return nothing end
    rts isa AbstractVector || return nothing
    pts = Sym[q for q in _real_roots(rts) if !_hassymbol(q)]
    # the ends of the real domain and the jumps: zeros of radicands, of log's arguments, and
    # of the denominators inside atan (1/x in acot x), at any depth
    us = Sym[]
    _replace_walk(e) do a
        if a.kind == POW && has(a.args[1], x)
            push!(us, a.args[1])
        elseif a.kind == FN && a.name === :log && has(a.args[1], x)
            push!(us, a.args[1])
        elseif a.kind == FN && a.name === :atan && has(a.args[1], x)
            d = numer_denom(a.args[1])[2]
            has(d, x) && push!(us, d)
        end
        nothing
    end
    for u in unique(us)
        zs = try _solve_eq(u ~ 0, x) catch err; err isa ArgumentError || rethrow(); return nothing end
        zs isa AbstractVector || return nothing
        append!(pts, [q for q in _real_roots(zs) if !_hassymbol(q)])
    end
    for b in unique(_zero_sets!(Sym[], e, x))
        zs = try _solve_eq(b ~ 0, x) catch err; err isa ArgumentError || rethrow(); continue end
        zs isa AbstractVector && append!(pts, [q for q in _real_roots(zs) if !_hassymbol(q)])
    end
    rel = _REL_INFO[r.op][2]
    S = _set_from(unique(pts), (v, atpoint) -> begin
        w = v isa Sym ? v : Sym(v)
        val = try N(subs(e, x => w)) catch err; err isa Union{DivisionByZero,DomainError,ArgumentError} || rethrow(); return false end
        val isa Number || return false
        val isa Complex && (abs(imag(val)) <= 1e-12 * max(1, abs(val)) ? (val = real(val)) : return false)
        isfinite(val) || return false
        atpoint && abs(val) <= 1e-12 && (val = 0.0)
        rel(val, 0)
    end)
    SymSet(S.intervals, x)
end
