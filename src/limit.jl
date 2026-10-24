# Limits: SymJ's `limit`, computed by GruntzLimit (Gruntz's algorithm) on SymJ's
# expressions, with SymJ's assumptions and answers.
#
# The answer is gruntz_limit's generic one where there is one (constants that may vanish
# taken to be nonzero); only when the limit depends on a sign it cannot decide is it
# given by cases, through gruntz_cases.
#
# The assumptions are handed to GruntzLimit as `signs` (c => :positive and the like):
# equalities are substituted; a parameter's range gives its bounds and excluded points;
# linear relations and ≠ facts are signs of expressions. What `signs` cannot say (a gap
# in a range, integrality) is applied to the answer: cases whose conditions contradict
# the assumptions are dropped, and conditions they imply removed. Answers are Syms
# (±∞ as Sym(Inf)); `gruntz_cases`'s cases become a `Cases`, its conditions relations,
# and a limit that does not exist or is undefined an error that calling the case throws.
# As in GruntzLimit, the other symbols are real numbers.

"""
    limit(f, x => x0; dir = :both, assume = nothing)

The limit of `f` as `x` tends to `x0`: a number, an expression, or `±Inf`. The limit is
two-sided unless `dir` is `:+` (from above) or `:-` (from below). Computed by Gruntz's
algorithm ([GruntzLimit](https://github.com/jverzani/GruntzLimit.jl)), which handles
exponentials and logarithms that cancel in ways numerical evaluation cannot follow.

The other symbols are real numbers, and the answer is the generic one: a constant that
may vanish is taken to be nonzero. Facts about them can be given with `assume` (or in
an [`assuming`](@ref) block, see [`Assumptions`](@ref)). When the limit depends on the
sign of something the facts leave open, the answer is a `Cases` value, one limit for each
case, which can be called with values (`cs(a => 2)`). The limit is a `Sym`, `Sym(Inf)` or
`Sym(-Inf)` for an infinite one. A limit that does not exist (`sin(x)` at `Inf`) throws
a `GruntzError`.

```julia
limit(sin(x) / x, x => 0)                    # 1
limit((1 + 1 / x)^x, x => Inf)               # ℯ
limit(1 / x, x => 0; dir = :+)               # Inf
limit(x^a, x => Inf)                         # a < 0 => 0, a == 0 => 1, a > 0 => Inf
limit(x^a, x => Inf; assume = @assume(a > 1))   # Inf
```
"""
function limit(f, xx0::Pair; dir = :both, assume = nothing)
    x, x0 = Sym(xx0.first), xx0.second
    _isnamed(x) || throw(ArgumentError("limit: $x is not a symbol"))
    dir = dir in (:+, +, "+") ? :+ : dir in (:-, -, "-") ? :- : dir in (:both, "both") ? :both :
        throw(ArgumentError("limit: dir must be :both, :+ or :-, not $(repr(dir))"))
    f = Sym(f)
    x0 = x0 isa Sym && x0.kind == NUMBER && isfloat(x0.num) && isinf(Float64(x0.num)) ? Float64(x0.num) : x0
    A = _assumptions(assume)
    signs = Pair{Sym,Symbol}[]
    if A !== nothing
        for p in _named(f)
            _eq(p, x) || get(A.notin, p, Int8(0)) < _RR ||
                throw(ArgumentError("limit: the other symbols are taken to be real, and $p is not"))
        end
        sub = [k => v for (k, v) in A.subst if !_eq(k, x)]
        if !isempty(sub)
            f = subs(f, sub...)
            x0 isa Real || (x0 = subs(Sym(x0), sub...))
        end
        signs = _gruntz_signs(A, x)
    end
    x0g = x0 isa Real ? x0 : Sym(x0)
    # the generic answer, as gruntz_limit gives it (a constant that may vanish taken to be
    # nonzero, where the expression is defined); by cases only when the limit depends on a
    # sign that is not decided
    L = try
        GruntzLimit.gruntz_limit(f, x, x0g; dir, signs)
    catch err
        (err isa GruntzError && occursin("gruntz_cases", err.msg)) || rethrow()
        nothing
    end
    L === nothing || return _limit_value(L)
    cs = GruntzLimit.gruntz_cases(f, x, x0g; dir, signs)
    cases = Pair{Vector{Any},Any}[]
    for (conds, v) in cs
        rels = Any[_relation(c) for c in conds]
        if A !== nothing && !isempty(rels)
            B = _build(_flatten!(Any[], (A.facts, rels)))
            B === nothing && continue                     # excluded by the assumptions
            rels = _simplify_conds(rels, A)
        end
        push!(cases, rels => _limit_value(v))
    end
    isempty(cases) && throw(ArgumentError("limit: no case is consistent with the assumptions"))
    # the cases left cover what the assumptions allow: one answer throughout is the answer
    all(c -> _sameanswer(c.second, cases[1].second), cases) && (cases = [Any[] => cases[1].second])
    if length(cases) == 1 && isempty(cases[1].first)
        v = cases[1].second
        v isa Exception && throw(v)
        return v
    end
    T = typejoin(map(c -> typeof(c.second), cases)...)
    Cases{T}([c.first => c.second for c in cases])
end

const _GRUNTZ_REL = Dict(:< => :lt, :<= => :le, :(==) => :eq, :!= => :ne, :>= => :ge, :> => :gt)
_relation(c) = Relation(_GRUNTZ_REL[c.rel], Sym(c.lhs), Sym(c.rhs))

function _limit_value(v)
    v isa Exception && return v
    v === :nolimit && return GruntzLimit.GruntzError("the limit does not exist", :nolimit)
    v === :undefined && return GruntzLimit.GruntzError("the expression is not defined there", :undefined)
    Sym(v)
end

# the assumptions as GruntzLimit's signs: c => :positive, :nonnegative or :nonzero
function _gruntz_signs(A::Assumptions, x::Sym)
    out = Pair{Sym,Symbol}[]
    for (p, S) in A.sets
        _eq(p, x) && continue
        isempty(S.intervals) && continue
        I1, I2 = first(S.intervals), last(S.intervals)
        I1.lo === nothing || push!(out, expand(p - I1.lo) => (I1.lclosed ? :nonnegative : :positive))
        I2.hi === nothing || push!(out, expand(I2.hi - p) => (I2.rclosed ? :nonnegative : :positive))
        for (I, J) in zip(S.intervals[1:end-1], S.intervals[2:end])    # a single point left out
            I.hi !== nothing && J.lo !== nothing && I.hi == J.lo && !I.rclosed && !J.lclosed &&
                push!(out, expand(p - I.hi) => :nonzero)
        end
    end
    for c in A.lin
        e = expand(_add_many(Sym(v) * t for (t, v) in c.a) + c.k)
        has(e, x) || push!(out, e => (c.strict ? :positive : :nonnegative))
    end
    for d in A.nonzero
        has(d, x) || push!(out, d => :nonzero)
    end
    for (p, vs) in A.excluded
        _eq(p, x) && continue
        for v in vs
            push!(out, expand(p - v) => :nonzero)
        end
    end
    out
end
