# Answers by cases: an operation that takes assumptions asks for decisions (does this
# fact hold?) with `_decide`. When the assumptions do not settle one, the operation is
# run again for each way it can go, the case's condition added to the assumptions:
#   - a fact c ⋈ 0 with c known to be real: one case for each sign of c the assumptions
#     allow (c < 0, c = 0, c > 0);
#   - otherwise two cases, the fact and its negation (`!(a ≫ 0)`, which with symbols
#     complex by default includes a not real). The negation is recorded as a decided
#     answer rather than added as a fact.
# A fact is first reduced to a simpler one where that is equivalent (√a > 0 to a > 0,
# 1/a > 0 to a > 0, 2a - 1 > 0 to a > 1/2). Cases with the same answer are merged, as in
# GruntzLimit (c < 0 and c = 0 to c ≤ 0). With one case and no conditions, the answer is
# returned as it is.

"""
    Cases

The answer of an operation that depends on facts the assumptions do not decide: a
vector of `conditions => answer` pairs, the conditions a vector of facts that all hold
in that case. Not exported.

Calling it with values for the symbols in the conditions, `cs(a => 2)`, gives the
answer for those values; indexing with them, `cs[a => 2]`, gives their case.

```julia
cs = solve(x^2 ~ a, x; assume = (x ≫ 0,))
cs[1]          # [a ≫ 0] => [sqrt(a)]
cs(a => 4)     # [2]
```
"""
struct Cases{T} <: AbstractVector{Pair{Vector{Any},T}}
    cases::Vector{Pair{Vector{Any},T}}
end
Base.size(cs::Cases) = size(cs.cases)
Base.getindex(cs::Cases, i::Int) = cs.cases[i]
Base.IndexStyle(::Type{<:Cases}) = IndexLinear()

function Base.getindex(cs::Cases, vals::Pair...)
    for c in cs.cases
        all(f -> _holds_at(f, vals), c.first) && return c
    end
    throw(ArgumentError("no case holds for $(join(vals, ", "))"))
end
function (cs::Cases)(vals::Pair...)
    v = cs[vals...].second
    v isa Exception && throw(v)
    _subsvals(v, vals)
end
_subsvals(v::Sym, vals) = subs(v, vals...)
_subsvals(v::AbstractVector{Sym}, vals) = Sym[_subsvals(u, vals) for u in v]
_subsvals(v::AbstractVector, vals) = map(u -> _subsvals(u, vals), v)
_subsvals(v, vals) = v

# shown as GruntzLimit shows its LimitCases: conditions with Julia's operators
function Base.show(io::IO, ::MIME"text/plain", cs::Cases)
    println(io, "Cases with ", length(cs.cases), " cases:")
    for (i, (conds, v)) in enumerate(cs.cases)
        print(io, "  ")
        isempty(conds) ? print(io, "always") : print(io, _condsstring(conds))
        print(io, "  =>  ")
        v isa Exception ? print(io, sprint(showerror, v)) :
            v isa AbstractVector{Sym} ? print(io, "[", join(v, ", "), "]") : show(io, v)
        i < length(cs.cases) && println(io)
    end
end

# the negation of a fact, as a condition (not added to the assumptions)
struct _NotFact
    f::Any
end
Base.show(io::IO, n::_NotFact) = print(io, "!(", n.f, ")")
Base.:(==)(a::_NotFact, b::_NotFact) = a.f == b.f
Base.hash(n::_NotFact, h::UInt) = hash(n.f, hash(:not, h))
const _ASCII_REL = Dict(:eq => "==", :ne => "!=", :lt => "<", :le => "<=", :gt => ">", :ge => ">=")
_condstring(c::Relation) = string(c.lhs, " ", _ASCII_REL[c.op], " ", c.rhs)
_condstring(c::_NotFact) = "!(" * _condstring(c.f) * ")"
_condstring(c) = string(c)
# a lower and an upper bound of the same expression chained: -3 < a < 3
function _condsstring(conds)
    parts, used = String[], falses(length(conds))
    for (i, c) in enumerate(conds)
        used[i] && continue
        used[i] = true
        if c isa Relation && c.op in (:gt, :ge)
            j = findfirst(k -> !used[k] && conds[k] isa Relation && conds[k].op in (:lt, :le) &&
                               _eq(conds[k].lhs, c.lhs), eachindex(conds))
            if j !== nothing
                used[j] = true
                lo, hi = c.op === :gt ? "<" : "<=", conds[j].op === :lt ? "<" : "<="
                push!(parts, string(c.rhs, " ", lo, " ", c.lhs, " ", hi, " ", conds[j].rhs))
                continue
            end
        end
        push!(parts, _condstring(c))
    end
    join(parts, " && ")
end

# does the condition hold for these values?
function _holds_at(f, vals)
    f isa _NotFact && return !_holds_at(f.f, vals)
    if f isa Membership
        v = subs(f.x, vals...)
        r = v in f.set
        r isa Bool || throw(ArgumentError("give values for the symbols in $f"))
        return r != f.negated
    end
    r = f::Relation
    d = subs(r.lhs - r.rhs, vals...)
    _hassymbol(d) && throw(ArgumentError("give values for the symbols in $r"))
    if r.op in (:lt, :le, :gt, :ge) && _isreal_const(d) !== true
        return false                         # an order fact needs a real number
    end
    v = _decide(Relation(r.op, d, ZERO))
    v isa Bool ? v : false
end

## decisions

struct _NeedDecision <: Exception
    f::Any                               # a Relation c ⋈ 0, or a Membership
end

# decided answers for the facts the assumptions cannot add (negations)
const _DECIDED = ScopedValue{Dict{Any,Bool}}(Dict{Any,Bool}())

# does the fact hold? Throws _NeedDecision when the assumptions do not say
function _decide(f, A::Assumptions)
    d = _DECIDED[]
    haskey(d, f) && return d[f]
    v = _holds(f, A)
    v === nothing || return v
    g = _reduce_fact(f)                       # an equivalent fact, perhaps decided
    haskey(d, g) && return d[g]
    g == f || (v = _holds(g, A)) === nothing || return v
    throw(_NeedDecision(f))
end

# true, false or nothing
function _holds(f, A::Assumptions)
    if f isa Membership
        q = f.set.rank == _ZR ? _isinteger(f.x, A) : f.set.rank == _RR ? _isreal(f.x, A) :
            f.set.rank == _CR ? true : _isrational(f.x, A)
        return q === nothing ? nothing : q != f.negated
    end
    r = f::Relation
    c = expand(r.lhs - r.rhs)
    if r.op in (:eq, :ne)
        nz = _isnonzero(c, A)
        return nz === nothing ? nothing : (r.op === :ne) == nz
    end
    m = _signs(c, A)
    if m === nothing
        return _isreal(c, A) === false ? false : nothing
    end
    ok = r.op === :lt ? _NEG : r.op === :le ? (_NEG | _ZER) : r.op === :gt ? _POS : (_POS | _ZER)
    m & ~ok == 0 ? true : m & ok == 0 ? false : nothing
end

function _isrational(e::Sym, A::Assumptions)
    _isinteger(e, A) === true && return true
    e.kind == NUMBER && return isexact(e.num)
    _isnamed(e) && return get(A.dom, e, _CR) <= _QR ? true : get(A.notin, e, Int8(0)) >= _QR ? false : nothing
    nothing
end

## reducing a fact to an equivalent simpler one

const _FALSE = Relation(:eq, ONE, ZERO)

# f ⇔ the result; f itself when nothing simpler is known
function _reduce_fact(f)
    if f isa Membership && f.set.rank == _RR && !f.negated
        # b^q, q not an integer, |q| < 1 (and k b^q, k real): real exactly when b ≥ 0
        c = f.x
        c.kind == MUL && length(c.args) == 1 && (c = c.args[1]^c.exps[1])
        if c.kind == POW && c.args[2].kind == NUMBER && isexact(c.args[2].num) &&
           !isone(_denominator(c.args[2].num)) && abs(c.args[2].num) < N1
            return _solved_for(Relation(c.args[2].num > N0 ? :ge : :gt, c.args[1], ZERO))
        end
        return f
    end
    f isa Relation || return f
    r = Relation(f.op, expand(f.lhs - f.rhs), ZERO)
    for _ in 1:8
        g = _reduce_once(r)
        g === nothing && break
        r = g
    end
    _solved_for(r)
end

function _reduce_once(r::Relation)
    c = r.lhs
    flip(op) = op === :lt ? :gt : op === :le ? :ge : op === :gt ? :lt : op === :ge ? :le : op
    if c.kind == MUL && length(c.args) == 1               # k * b^q
        k = c.num
        op = k < N0 ? flip(r.op) : r.op
        return Relation(op, c.args[1]^c.exps[1], ZERO)
    end
    if c.kind == POW && c.args[2].kind == NUMBER && isexact(c.args[2].num) && r.op in (:lt, :le, :gt, :ge)
        b, q = c.args[1], c.args[2].num
        if isone(_denominator(q)) && isodd(_numerator(q))
            # b^n, n odd: the sign of b (never zero for n < 0)
            return Relation(q < N0 && r.op === :le ? :lt : q < N0 && r.op === :ge ? :gt : r.op, b, ZERO)
        elseif !isone(_denominator(q)) && abs(q) < N1
            # b^q, q not an integer, |q| < 1: real and positive exactly when b > 0, zero
            # when b = 0 (q > 0), never negative
            r.op === :gt && return Relation(:gt, b, ZERO)
            r.op === :ge && return Relation(q > N0 ? :ge : :gt, b, ZERO)
            r.op === :le && q > N0 && return Relation(:eq, b, ZERO)
            return _FALSE
        end
    end
    nothing
end

# c ⋈ 0 with c linear in one symbol: the symbol against a value (2a - 1 > 0 to a > 1/2)
function _solved_for(r::Relation)
    vs = _named(r.lhs)
    length(vs) == 1 || return r
    p = only(vs)
    a, k = _linform(r.lhs)
    (length(a) == 1 && haskey(a, p)) || return r
    c = a[p]
    op = c < N0 && r.op in (:lt, :le, :gt, :ge) ?
        Dict(:lt => :gt, :le => :ge, :gt => :lt, :ge => :le)[r.op] : r.op
    Relation(op, p, -k / Sym(c))
end

## running an operation by cases

# A with the fact g, nothing if they contradict; A itself when g cannot be added (a
# membership of an expression, a relation not supported), its answer then only decided
function _with_fact(A::Assumptions, g)
    g isa Membership && !_isnamed(g.x) && return A
    try
        _build(_flatten!(Any[], (A.facts, g)))
    catch err
        err isa ArgumentError || rethrow()
        A
    end
end

const _MAX_CASES = 64

# run(A) gives the answer, asking for decisions with `_decide`; the answer, or a Cases
function _by_cases(run, A::Assumptions)
    count = Ref(0)
    leaves = _case_tree(run, A, Dict{Any,Bool}(), count)
    length(leaves) == 1 && isempty(leaves[1].first) &&
        return (v = leaves[1].second; v isa Exception ? throw(v) : v)
    T = typejoin(map(l -> typeof(l.second), leaves)...)
    Cases{T}([_simplify_conds(l.first, A) => l.second for l in leaves])
end

# a case's conditions without those the others imply (a ≥ 0 & a > 0 is a > 0)
function _simplify_conds(conds, A::Assumptions)
    out = Vector{Any}(conds)
    i = 1
    while i <= length(out)
        c = out[i]
        if !(c isa _NotFact)
            rest = [out[j] for j in eachindex(out) if j != i && !(out[j] isa _NotFact)]
            B = try _build(_flatten!(Any[], (A.facts, rest))) catch; nothing end
            if B !== nothing && _holds(c, B) === true
                deleteat!(out, i)
                continue
            end
        end
        i += 1
    end
    out
end

function _case_tree(run, A::Assumptions, decided::Dict{Any,Bool}, count::Ref{Int})
    need = try
        v = with(() -> run(A), _DECIDED => decided)
        return [Any[] => v]
    catch err
        err isa _NeedDecision || (err isa ArgumentError ? (return [Any[] => err]) : rethrow())
        err
    end
    count[] += 1
    count[] > _MAX_CASES && return [Any[] => ArgumentError("too many cases")]
    f = need.f
    g = _reduce_fact(f)
    branches = Pair{Any,Any}[]                       # condition => (assumptions, decided)
    c = g isa Relation ? expand(g.lhs - g.rhs) : nothing
    if g isa Relation && g.op in (:lt, :le, :gt, :ge) && _isreal(c, A) === true
        m = _signs(c, A)
        for (bit, op, s) in ((_NEG, :lt, -1), (_ZER, :eq, 0), (_POS, :gt, 1))
            m & bit == 0 && continue
            h = _solved_for(Relation(op, c, ZERO))
            B = _with_fact(A, h)
            B === nothing && continue
            d = copy(decided)
            d[f] = d[g] = _REL_INFO[g.op][2](s, 0)   # f ⇔ g, and g at this sign
            push!(branches, h => (B, d))
        end
    elseif g isa Relation && g.op in (:eq, :ne)
        for op in (:eq, :ne)
            h = _solved_for(Relation(op, c, ZERO))
            B = _with_fact(A, h)
            B === nothing && continue
            d = copy(decided); d[f] = d[g] = op === g.op
            push!(branches, h => (B, d))
        end
    else
        B = _with_fact(A, g)
        d1 = copy(decided); d1[f] = d1[g] = true
        B === nothing || push!(branches, g => (B, d1))
        d2 = copy(decided); d2[f] = d2[g] = false
        push!(branches, _NotFact(g) => (A, d2))
    end
    kids = [cond => _case_tree(run, B, d, count) for (cond, (B, d)) in branches]
    _merge_cases(kids)
end

# children with the same answers merged: all the same, no condition; two neighbouring
# signs, one condition (c < 0 and c = 0 give c ≤ 0)
function _merge_cases(kids)
    length(kids) == 1 && return [vcat(Any[kids[1].first], l.first) => l.second for l in kids[1].second]
    same(a, b) = length(a) == length(b) && all(_samecase(p, q) for (p, q) in zip(a, b))
    all(k -> same(k.second, kids[1].second), kids) && return kids[1].second
    if length(kids) == 3 && all(k -> k.first isa Relation, kids)
        lt, eq, gt = kids
        if same(lt.second, eq.second)
            kids = [Relation(:le, lt.first.lhs, lt.first.rhs) => lt.second, gt]
        elseif same(eq.second, gt.second)
            kids = [lt, Relation(:ge, gt.first.lhs, gt.first.rhs) => gt.second]
        elseif same(lt.second, gt.second)
            kids = [Relation(:ne, lt.first.lhs, lt.first.rhs) => lt.second, eq]
        end
    end
    out = Pair{Vector{Any},Any}[]
    for (cond, ls) in kids, l in ls
        push!(out, vcat(Any[cond], l.first) => l.second)
    end
    out
end
_samecase(p, q) = p.first == q.first && _sameanswer(p.second, q.second)
_sameanswer(a, b) = a isa Exception || b isa Exception ? (typeof(a) == typeof(b) && sprint(showerror, a) == sprint(showerror, b)) :
    isequal(a, b)

## with assumptions: the roots the assumptions allow for x, by cases where that depends
## on facts they do not decide

function _solve_assumed(eq::Relation, x::Sym, A::Assumptions)
    _isnamed(x) || throw(ArgumentError("solve: $x is not a symbol"))
    eq.op === :eq || return _solve_ineq_assumed(eq, x, A)
    _by_cases(A) do B
        sub = [k => w for (k, w) in B.subst if !_eq(k, x)]
        e = isempty(sub) ? eq : Relation(:eq, subs(eq.lhs, sub...), subs(eq.rhs, sub...))
        v = get(B.subst, x, nothing)
        if v !== nothing                          # x itself is given: does it solve eq?
            return _decide(Relation(:eq, expand(subs(e.lhs - e.rhs, x => v)), ZERO), B) ? Sym[v] : Sym[]
        end
        filter(r -> _allowed(r, x, B), _solve_eq(e, x))
    end
end

# may x take the value r?
function _allowed(r::Sym, x::Sym, A::Assumptions)
    S = get(A.sets, x, nothing)
    dom = get(A.dom, x, _CR)
    # (a set implies real values)
    (dom < _RR || (dom == _RR && S === nothing)) && !_decide(Membership(r, NumberSet(dom), false), A) && return false
    no = get(A.notin, x, Int8(0))
    no > 0 && !_decide(Membership(r, NumberSet(no), true), A) && return false
    S === nothing || _inset(r, S, A) || return false
    for v in get(A.excluded, x, Sym[])
        _decide(Relation(:ne, r, v), A) || return false
    end
    for c in A.lin                                # relations of x with other symbols
        any(t -> has(t, x), keys(c.a)) || continue
        val = expand(_add_many(Sym(a) * subs(t, x => r) for (t, a) in c.a) + c.k)
        _decide(Relation(c.strict ? :gt : :ge, val, ZERO), A) || return false
    end
    for d in A.nonzero
        has(d, x) || continue
        _decide(Relation(:ne, expand(subs(d, x => r)), ZERO), A) || return false
    end
    true
end

_inset(r::Sym, S::SymSet, A::Assumptions) = any(S.intervals) do I
    (I.lo === nothing || _decide(Relation(I.lclosed ? :ge : :gt, r, I.lo), A)) &&
        (I.hi === nothing || _decide(Relation(I.rclosed ? :le : :lt, r, I.hi), A))
end

# an inequality: its set within x's; parameters given values by equalities substituted
function _solve_ineq_assumed(r::Relation, x::Sym, A::Assumptions)
    get(A.dom, x, _CR) <= _QR && throw(ArgumentError("solve: inequalities in an integer or rational unknown are not supported"))
    sub = [k => w for (k, w) in A.subst if !_eq(k, x)]
    S = _solve_eq(isempty(sub) ? r : Relation(r.op, subs(r.lhs, sub...), subs(r.rhs, sub...)), x)
    T = get(A.sets, x, nothing)
    T === nothing ? S : SymSet((S ∩ T).intervals, x)
end
