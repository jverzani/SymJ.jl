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

# |u^r| as |u|^r for u and r real (here: free of i)
_abs_of_powers(e::Sym) = _replace_walk(e) do t
    t.kind == FN && t.name === :abs && length(t.args) == 1 && t.args[1].kind == POW || return nothing
    b, r = t.args[1].args
    (has(b, IM) || has(r, IM)) && return nothing
    abs(_abs_of_powers(b))^r
end

# does e depend on x: not through a derivative at a point whose bound variable is x (f'(3)
# is Derivative(f(x), x, 3))
_depends(e::Sym, x::Sym) = !has(e, x) ? false : _eq(e, x) ? true :
    _isderivat(e) ? (has(e.args[3], x) || (!_eq(e.args[2], x) && has(e.args[1], x))) :
    any(a -> _depends(a, x), e.args) || (e.kind == MUL && any(a -> _depends(a, x), e.exps))

# an unknown function (SymFunction(:f)(x), or a derivative of one) depending on x, or nothing
function _unknown_fn(e::Sym, x::Sym)
    _depends(e, x) || return nothing
    if e.kind == FN && (e.name === :Derivative ||
                        (!haskey(FUNCTION_RULES, e.name) && !(e.name in _NO_CHAIN)))
        return e
    end
    for a in e.args
        u = _unknown_fn(a, x)
        u === nothing || return u
    end
    if e.kind == MUL
        for a in e.exps
            u = _unknown_fn(a, x)
            u === nothing || return u
        end
    end
    nothing
end

function _unknown_consts!(out, e::Sym, x::Sym)
    if e.kind == FN && !_depends(e, x) && (e.name === :Derivative ||
                                      (!haskey(FUNCTION_RULES, e.name) && !(e.name in _NO_CHAIN)))
        any(c -> _eq(c, e), out) || push!(out, e)
        return out
    end
    foreach(a -> _unknown_consts!(out, a, x), e.args)
    e.kind == MUL && foreach(a -> _unknown_consts!(out, a, x), e.exps)
    out
end

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

`f` may itself be a `Cases` value (as a definite integral with a variable end gives): the
cases whose conditions on `x` hold near `x0` are kept, each with its limit (∫₁ᴹ as `M → ∞`
takes the case `M > 1`).

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
    f = _sinc_elementary(Sym(f))
    # |u^r| = |u|^r, the symbols real: |(-n/(n+1))^n| = (n/(n+1))^n (SymPy #26513)
    f = _abs_of_powers(f)
    # an unknown function of x: its limit is not known (it need not be continuous) (#27143)
    u = _unknown_fn(f, x)
    u === nothing || throw(ArgumentError("limit: $u is an unknown function of $x; its limit is not known"))
    # one free of x, f(2), is a constant: a symbol for the limit, put back after
    cs = Sym[]
    _unknown_consts!(cs, f, x)
    if !isempty(cs)
        ks = [Sym(Symbol("_lim_c", i)) for i in eachindex(cs)]
        L = limit(subs(f, (c => k for (c, k) in zip(cs, ks))...), x => x0; dir, assume)
        return L isa Sym ? subs(L, (k => c for (c, k) in zip(cs, ks))...) : L
    end
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
        # an integer variable at ±∞: sin(πv), cos(πv) without v's integer-valued polynomial part
        if get(A.dom, x, _CR) == _ZR && x0 isa Real && isinf(x0)
            f = _integer_trig_shift(f, x, x0, A)
        end
    end
    x0g = x0 isa Real ? x0 : Sym(x0)
    # a piecewise function whose conditions are free of x (a ≠ 0) is a split on the
    # parameters: the limit of each branch, under its condition
    P = _param_piecewise(f, x)
    if P !== nothing
        return piecewise((c => limit(subs(f, P => v), x => x0; dir, assume) for (c, v) in _pieces(P))...)
    end
    # piecewise functions, heaviside and dirac: the piece on each side
    if _haspieces(f, x)
        if x0g isa Real && isinf(x0g)
            f = _pieces_near(f, x, x0g, x0g > 0 ? :- : :+, assume)
        elseif dir === :both
            fp, fm = _pieces_near(f, x, x0g, :+, assume), _pieces_near(f, x, x0g, :-, assume)
            if !_eq(fp, fm)
                Lp, Lm = limit(fp, x => x0; dir = :+, assume), limit(fm, x => x0; dir = :-, assume)
                _sameanswer(Lp, Lm) || throw(GruntzError("one-sided limits differ", :nolimit))
                return Lp
            end
            f = fp
        else
            f = _pieces_near(f, x, x0g, dir, assume)
        end
    end
    # the generic answer, as gruntz_limit gives it (a constant that may vanish taken to be
    # nonzero, where the expression is defined); by cases only when the limit depends on a
    # sign that is not decided
    L = try
        GruntzLimit.gruntz_limit(f, x, x0g; dir, signs)
    catch err
        (err isa GruntzError && err.kind === :sign) || rethrow()
        nothing
    end
    L === nothing || return _limit_value(L)
    cs = GruntzLimit.gruntz_cases(f, x, x0g; dir, signs)
    cases = Pair{Vector{Any},Any}[]
    for (conds, v) in cs
        rels = Any[_relation(c) for c in conds]
        if A !== nothing && !isempty(rels)
            B = try _build(_flatten!(Any[], (A.facts, rels))) catch err; err isa ArgumentError || rethrow(); A end
            B === nothing && continue                     # excluded by the assumptions
            any(r -> _holds(r, A) === false, rels) && continue   # (2^y > y for y > 0)
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

# For an integer n: sin(πv) = (-1)^p sin(π(v - p)) and cos(πv) = (-1)^p cos(π(v - p)) for
# any integer-valued p; with p the polynomial part of v's expansion at ±∞ (integer
# coefficients), v - p has a limit: sin(π√(4n² + 1)) = sin(π(√(4n² + 1) - 2n)) → 0
# (SymPy #28313)
function _integer_trig_shift(f::Sym, n::Sym, x0::Real, A)
    _replace_walk(f) do e
        (e.kind == FN && e.name in (:sin, :cos) && length(e.args) == 1) || return nothing
        v = expand(e.args[1] / PI)
        has(v, n) || return nothing
        s = try series(v, n, x0, 2) catch err; err isa ArgumentError || rethrow(); return nothing end
        p = ZERO
        for (ex, c) in zip(s.ps.e, s.ps.c)
            ex > 0 && continue                       # the vanishing terms (t = 1/n)
            (isinteger(ex) && c.kind == NUMBER && _isint(c.num)) || return nothing
            p += c * (x0 > 0 ? n : -n)^Int(-ex)
        end
        iszero(p) && return nothing
        # (refining only the sign: refine would take the integer multiple of π back out)
        _refine((-1)^p, A) * apply_function(e.name, Sym[PI * expand(v - p)])
    end
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

## an answer by cases: its limit, and substitution into it

# The cases a limit in x keeps are those whose conditions on x hold near x0 on the side
# taken (M > 1 as M → ∞; x < c from below c); each kept case's answer has its limit,
# its conditions on other symbols kept. A condition on x that holds on one side only, with
# dir = :both, gives the two one-sided limits, which must agree.
function limit(cs::Cases, xx0::Pair; dir = :both, assume = nothing)
    x, x0 = Sym(xx0.first), xx0.second
    x0s = x0 isa Real && isinf(x0) ? Sym(Float64(x0)) : Sym(x0)
    if dir in (:both, "both") && !_isinfpt(x0s)
        Lp, Lm = limit(cs, xx0; dir = :+, assume), limit(cs, xx0; dir = :-, assume)
        _sameanswer(Lp, Lm) || throw(GruntzError("one-sided limits differ", :nolimit))
        return Lp
    end
    side = _isinfpt(x0s) ? (Float64(x0s.num) > 0 ? :- : :+) : (dir in (:+, +, "+") ? :+ : :-)
    out = Pair{Vector{Any},Any}[]
    for (conds, v) in cs
        keep = true
        rest = Any[]
        for c in conds
            _fact_has(c, x) || (push!(rest, c); continue)
            h = _eventually(c, x, x0s, side, assume)
            h === nothing && throw(ArgumentError("limit: cannot tell whether the case $(_condsstring([c])) holds as $x → $x0"))
            h || (keep = false; break)
        end
        keep || continue
        v isa Exception && throw(v)
        push!(out, rest => (v isa Sym ? limit(v, xx0; dir = side, assume) : v))
    end
    isempty(out) && throw(ArgumentError("limit: no case holds as $x → $x0"))
    length(out) == 1 && isempty(out[1].first) && return out[1].second
    T = typejoin(map(c -> typeof(c.second), out)...)
    Cases{T}(out)
end

# does the fact about x hold for x near x0 on the side (:+ above, :- below)? A fact d ⋈ 0
# holds there by the sign of d's one-sided limit, or for d = x - c with c = x0, by the side
function _eventually(c, x::Sym, x0::Sym, side::Symbol, assume)
    c isa _NotFact && return (h = _eventually(c.f, x, x0, side, assume); h === nothing ? nothing : !h)
    c isa Relation || return nothing
    d = expand(c.lhs - c.rhs)
    L = try limit(d, x => (_isinfpt(x0) ? Float64(x0.num) : x0); dir = side, assume) catch; return nothing end
    L isa Sym || return nothing
    s = _isinfpt(L) ? (Float64(L.num) > 0 ? 1 : -1) : _sign_known(L, assume)
    if s == 0
        # d → 0: its sign near x0 from the sign of d' at x0 (d = x - c: the side)
        _isinfpt(x0) && return nothing
        s = _sign_known(subs(diff(d, x), x => x0), assume)
        (s === nothing || s == 0) && return nothing
        side === :- && (s = -s)
        c.op in (:eq,) && return false
        c.op in (:ne,) && return true
    end
    s === nothing && return nothing
    c.op === :eq ? false : c.op === :ne ? true :
    c.op in (:gt, :ge) ? s > 0 : s < 0
end

# the sign of a constant, or of an expression from the assumptions: -1, 0, 1 or nothing
function _sign_known(v::Sym, assume)
    _hassymbol(v) || return _csign(v)
    A = _assumptions(assume)
    A === nothing && return nothing
    m = _signs(v, A)
    m === _POS ? 1 : m === _NEG ? -1 : m === _ZER ? 0 : nothing
end

"""
    subs(cs::Cases, x => v, ...)

Substitution into each case: its answer and its conditions. A condition that becomes
decided (a statement about numbers) is dropped when it holds, and its case when it does
not; the answer of the one case left with no conditions is returned as it is.
"""
function subs(cs::Cases, ps::Pair...)
    out = Pair{Vector{Any},Any}[]
    for (conds, v) in cs
        keep = true
        rest = Any[]
        for c in conds
            c2 = _subsfact_case(c, ps)
            h = _decided(c2)
            h === nothing ? push!(rest, c2) : (h || (keep = false; break))
        end
        keep || continue
        push!(out, rest => _subsvals(v, ps))
    end
    isempty(out) && throw(ArgumentError("subs: no case holds for $(join(ps, ", "))"))
    if length(out) == 1 && isempty(out[1].first)
        v = out[1].second
        v isa Exception && throw(v)
        return v
    end
    T = typejoin(map(c -> typeof(c.second), out)...)
    Cases{T}(out)
end
_subsfact_case(c::_NotFact, ps) = _NotFact(_subsfact_case(c.f, ps))
_subsfact_case(c::Relation, ps) = Relation(c.op, subs(c.lhs, ps...), subs(c.rhs, ps...))
_subsfact_case(c::Membership, ps) = Membership(subs(c.x, ps...), c.set, c.negated)
_subsfact_case(c, ps) = c
# true or false for a fact about numbers; nothing otherwise
function _decided(c)
    c isa _NotFact && return (h = _decided(c.f); h === nothing ? nothing : !h)
    if c isa Relation
        d = c.lhs - c.rhs
        _hassymbol(d) || return try _holds_at(c, ()) catch; nothing end
        # with its symbols real, as a condition on an order is (b^2 + 1 > 0)
        A = try Assumptions(FactAnd(Any[p ∈ ℝ for p in free_symbols(d)])) catch; return nothing end
        return try _holds(c, A) catch; nothing end
    elseif c isa Membership
        _hassymbol(c.x) && return nothing
        r = c.x in c.set
        return r isa Bool ? r != c.negated : nothing
    end
    nothing
end
