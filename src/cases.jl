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
answer for those values; indexing with them, `cs[a => 2]`, gives their case. `subs`
substitutes into each case, keeping those whose conditions still can hold, and `limit`
takes the cases that hold near the point.

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
_subsvals(v::AbstractDict{Sym,Sym}, vals) = Dict{Sym,Sym}(k => subs(w, vals...) for (k, w) in v)
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
_condstring(c::FactAnd) = join(("(" * _condstring(g) * ")" for g in c.facts), " && ")
_condstring(c::FactOr) = join(("(" * _condstring(g) * ")" for g in c.facts), " || ")
_condstring(c) = string(c)
# a lower and an upper bound of the same expression chained: -3 < a < 3
function _condsstring(conds)
    parts, used = String[], falses(length(conds))
    for (i, c) in enumerate(conds)
        used[i] && continue
        used[i] = true
        if c isa Relation && c.op in (:gt, :ge, :lt, :le)       # in either order
            lower = c.op in (:gt, :ge)
            j = findfirst(k -> !used[k] && conds[k] isa Relation &&
                               conds[k].op in (lower ? (:lt, :le) : (:gt, :ge)) &&
                               _eq(conds[k].lhs, c.lhs), eachindex(conds))
            if j !== nothing
                used[j] = true
                l, u = lower ? (c, conds[j]) : (conds[j], c)
                lo, hi = l.op === :gt ? "<" : "<=", u.op === :lt ? "<" : "<="
                push!(parts, string(l.rhs, " ", lo, " ", l.lhs, " ", hi, " ", u.rhs))
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
    _hassymbol(e) || return _const_rational(e)
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

# c ⋈ 0 with c linear in a symbol: that symbol against the rest (2a - 1 > 0 to a > 1/2,
# a - b < 0 to a < b)
function _solved_for(r::Relation)
    vs = _named(r.lhs)
    isempty(vs) && return r
    a, k = _linform(r.lhs)
    i = findfirst(v -> haskey(a, v), vs)
    i === nothing && return r
    p = vs[i]
    c = a[p]
    rest = expand(r.lhs - Sym(c) * p)
    has(rest, p) && return r
    k = rest
    op = c < N0 && r.op in (:lt, :le, :gt, :ge) ?
        Dict(:lt => :gt, :le => :ge, :gt => :lt, :ge => :le)[r.op] : r.op
    Relation(op, p, -k / Sym(c))
end

## running an operation by cases

# the first factor of a product whose zero the assumptions leave open, or nothing
function _zero_factor(c::Sym, A::Assumptions)
    c.kind == MUL || return nothing
    count(_hassymbol, c.args) >= 2 || return nothing
    for (b, p) in zip(c.args, c.exps)
        _hassymbol(b) && p.kind == NUMBER && p.num > N0 || continue
        _isnonzero(b, A) === nothing && return b
    end
    nothing
end

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
    # the same answer in every case is the answer (a split on acos(A) > 0 for A > 1,
    # complex, gave 2π/√(A² - 1) in each; SymPy #26283)
    vs = [l.second for l in leaves]
    !any(v -> v isa Exception, vs) && all(v -> isequal(v, vs[1]), vs) && return vs[1]
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
        # an answer by cases of its own (one for each root of a degenerate factor): its
        # cases as leaves here
        v isa Cases && return [c.first => c.second for c in v.cases]
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
    elseif g isa Relation && g.op in (:eq, :ne) && (h = _zero_factor(c, A)) !== nothing
        # a product is zero when a factor is: split on one factor at a time (a = 0 or
        # not, then b), not on the product, whose zero is not a fact about the symbols
        for op in (:eq, :ne)
            B = _with_fact(A, Relation(op, h, ZERO))
            B === nothing || push!(branches, _solved_for(Relation(op, h, ZERO)) => (B, decided))
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
# equal answers, Cases compared case by case
_same_result(a, b) = a isa Cases && b isa Cases ?
    (length(a.cases) == length(b.cases) &&
     all(((p, q),) -> (isequal(p.first, q.first) || string(p.first) == string(q.first)) &&
                      _same_result(p.second, q.second), zip(a.cases, b.cases))) :
    _sameanswer(a, b)

_sameanswer(a, b) = a isa Exception || b isa Exception ? (typeof(a) == typeof(b) && sprint(showerror, a) == sprint(showerror, b)) :
    isequal(a, b)

## with assumptions: the roots the assumptions allow for x, by cases where that depends
## on facts they do not decide

function _solve_assumed(eq::Relation, x::Sym, A::Assumptions)
    _isnamed(x) || throw(ArgumentError("solve: $x is not a symbol"))
    eq.op === :eq || return _solve_ineq_assumed(eq, x, A)
    # floating point coefficients with parameters: the cases and solutions of the rationals
    # they stand for, the solutions' numbers floats again. Split as floats, the case
    # k0 = -k1 of (1.0k0 + 1.0k1)/(2.0(0.5k0 + 0.5k1)) left 0.0/0.0 (SymPy #26820).
    ex = eq.lhs - eq.rhs
    if !_isexactex(ex) && !_float_exponent(ex) && any(s -> !_eq(s, x), free_symbols(ex))
        q = try _exactfloats(ex) catch err; err isa ArgumentError || rethrow(); nothing end
        if q !== nothing && _isexactex(q)
            fl(v) = v isa AbstractVector{Sym} ? Sym[_floatify(r) for r in v] : v
            r = _solve_assumed(Relation(:eq, q, ZERO), x, A)
            return r isa Cases ? _map_cases(fl, r) : fl(r)
        end
    end
    _by_cases(A) do B
        sub = [k => w for (k, w) in B.subst if !_eq(k, x)]
        e = try
            isempty(sub) ? eq : Relation(:eq, subs(eq.lhs, sub...), subs(eq.rhs, sub...))
        catch err
            _isdivzero(err) || rethrow()
            return Sym[]                          # a denominator vanishes: not defined
        end
        v = get(B.subst, x, nothing)
        if v !== nothing                          # x itself is given: does it solve eq?
            return _decide(Relation(:eq, expand(subs(e.lhs - e.rhs, x => v)), ZERO), B) ? Sym[v] : Sym[]
        end
        rs = with(() -> _solve_case(e, x, B), _SOLVE_ASSUME => B)
        filter(r -> _allowed(r, x, B), _real_if_tiny(rs, x, B))
    end
end

## degenerate parameter values: coefficients that vanish
#
# The solutions of an equation with parameters hold where its coefficients do not
# vanish (a*x = b gives b/a for a ≠ 0). Each coefficient whose zero the assumptions
# leave open is decided by cases: the leading coefficients of a polynomial in x, until
# one is nonzero (a*x^2 + b*x + c: a ≠ 0, else b ≠ 0, else c ≠ 0); a denominator free of
# x (where it vanishes the equation is not defined); and otherwise the denominators of
# the solutions found (log(2)/a for exp(a*x) = 2). A case that sets a coefficient to zero
# is solved again with that value substituted.

# the factors of c's numerator that involve parameters but none of the unknowns, whose
# vanishing makes a case (each up to a constant factor)
# symbols that are not parameters: the variables of unknowns that are function calls
# (x for u(x): implicit differentiation divides by x without a case x = 0)
const _NOT_PARAMS = ScopedValue{Vector{Sym}}(Sym[])

function _param_factors(c::Sym, unknowns)
    out = Sym[]
    _hassymbol(c) || return out
    for part in (numer_denom(c)[1],)
        _hassymbol(part) || continue
        f = try factor(part) catch err; err isa ArgumentError || rethrow(); part end
        bs = f.kind == MUL ? [b for (b, p) in zip(f.args, f.exps) if p.kind == NUMBER && p.num > N0] :
             f.kind == POW && f.args[2].kind == NUMBER && f.args[2].num > N0 ? [f.args[1]] : [f]
        for b in bs
            _hassymbol(b) && !any(u -> has(b, u), unknowns) && !any(u -> has(b, u), _NOT_PARAMS[]) || continue
            atoms = Set{Sym}()
            _poly_atoms!(atoms, b)
            # a polynomial in the parameters, its coefficients constants such as √3
            # (a - √3 b; SymPy #24344)
            all(a -> _isnamed(a) || !_hassymbol(a), atoms) || continue
            b = _content_free(b)
            any(o -> (q = cancel(o / b); q.kind == NUMBER), out) || push!(out, b)
        end
    end
    out
end

# b without its numeric content and sign (2a^2 + 2b^2 and -a^2 - b^2 as a^2 + b^2)
function _content_free(b::Sym)
    b.kind == ADD && isexact(b.num) && all(isexact, b.coefs) || return b
    g = gcd(Rational{BigInt}[Rational{BigInt}(_numerator(c), _denominator(c)) for c in (iszero(b.num) ? b.coefs : [b.coefs; b.num])])
    b.coefs[1] < N0 && (g = -g)                   # the first term positive
    iszero(g) || isone(g) ? b : expand(b / Sym(SymNumber(g)))
end

# nothing when every factor of every expression in D is nonzero in this case; else the
# substitutions that make the first vanishing one zero: p => r for a factor linear in a
# parameter p (the assumptions substitute a = 0 or a = b themselves; this covers
# a*b = 1), one for each root otherwise (a^2 + b^2 = 0: a = ±i b)
function _degenerate_subst(D, B::Assumptions, unknowns)
    # the smallest factor first, solved for the parameter that gives the smallest value:
    # a large one substituted makes the system's coefficients large (a circuit's
    # determinant solved for one of its 13 parameters took hours; SymPy #24135)
    fs = Sym[]
    for g in D, f in _param_factors(g, unknowns)
        any(o -> _eq(o, f), fs) || push!(fs, f)
    end
    sort!(fs; by = f -> (_nodes(f), string(f)))
    for f in fs
        _decide(Relation(:ne, f, ZERO), B) && continue
        ps = _named(f)
        # f ≥ 0 in real parameters (a sum of squares, ζ²ω² + (Ω - ω₁)²): it vanishes where
        # its gradient does too; the real solutions of that system, each a substitution
        # (solved for one parameter, it had complex roots; SymPy #24518)
        sos = _nonneg_zero_set(f, ps, B)
        sos === nothing || return f => sos
        best = nothing
        for p in ps
            P = try Poly{Sym}(expand(f), [p]) catch; continue end
            degree(P, p) == 1 || continue
            c1, c0 = coeff(f, p, 1), coeff(f, p, 0)
            all(h -> _decide(Relation(:ne, h, ZERO), B), _param_factors(c1, unknowns)) || continue
            v = cancel(-c0 / c1)
            (best === nothing || (_nodes(v), string(p)) < (_nodes(best[2]), string(best[1]))) && (best = (p, v))
        end
        best === nothing || return f => [best[1] => best[2]]
        rs = isempty(ps) ? Sym[] : (try _solve_generic(f ~ 0, ps[1]) catch err; err isa ArgumentError || rethrow(); Sym[] end)
        isempty(rs) || return f => [ps[1] => r for r in rs]
        throw(ArgumentError("solve: the case $f = 0 is not supported"))
    end
    nothing
end

function _nonneg_zero_set(f::Sym, ps::Vector{Sym}, B::Assumptions)
    length(ps) >= 2 || return nothing
    all(p -> _isreal(p, B) === true, ps) || return nothing
    _ispoly_in(f, ps[1]) || return nothing
    all(p -> (P = try Poly{Sym}(expand(f), [p]) catch; nothing end; P !== nothing && degree(P, p) >= 2), ps) || return nothing
    (_nonneg(f, B) || _nonneg_quadratic(expand(f), ps, B, 2)) || return nothing
    eqs = Relation[Relation(:eq, f, ZERO); [Relation(:eq, diff(f, p), ZERO) for p in ps]]
    sols = try _solve_generic(eqs, ps) catch err; err isa ArgumentError || rethrow(); return nothing end
    out = Dict{Sym,Sym}[]
    for s in sols
        d = Dict{Sym,Sym}(p => v for (p, v) in s if !_eq(p, v))
        (isempty(d) || any(v -> has(v, IM), values(d))) && continue
        push!(out, d)
    end
    isempty(out) ? nothing : out
end

# a c p² + b p + c₀ ≥ 0 for real p when a > 0 and 4a c₀ - b² ≥ 0 (decided, or in turn
# such a quadratic in another parameter, to depth d)
function _nonneg_quadratic(f::Sym, ps::Vector{Sym}, B::Assumptions, d::Int)
    _nonneg(f, B) && return true
    d == 0 && return false
    for p in ps
        P = try Poly{Sym}(f, [p]) catch; continue end
        degree(P, p) == 2 || continue
        a, b, c = coeff(f, p, 2), coeff(f, p, 1), coeff(f, p, 0)
        _ispos(a, B) || continue
        _nonneg_quadratic(expand(4a * c - b^2), [q for q in ps if !_eq(q, p)], B, d - 1) && return true
    end
    false
end

# the answers under each alternative substitution for f = 0: one, or all of them when
# they are all empty (none has solutions); different answers depend on which root holds
function _alternatives(run, (f, ps))
    rs = try
        [run(p) for p in ps]
    catch err
        (err isa ArgumentError || _isdivzero(err)) && length(ps) > 1 || rethrow()
        throw(ArgumentError("solve: the case $f = 0 is not supported"))
    end
    length(rs) == 1 && return rs[1]
    all(r -> r isa AbstractVector && isempty(r), rs) && return rs[1]
    # the same answer under every alternative (the components ζ = 0 and ω = 0 of a sum of
    # squares ζ²ω² + (Ω - ω₁)²; SymPy #24518)
    all(r -> _same_result(r, rs[1]), rs) && return rs[1]
    # different answers for different roots (r² + rσ + σ² = 0 at r = σ(-1 ± i√3)/2): a
    # case for each, its condition the root
    conds = [p isa Pair ? Any[Relation(:eq, p.first, p.second)] :
             p isa AbstractDict ? Any[Relation(:eq, k, v) for (k, v) in p] : nothing for p in ps]
    any(isnothing, conds) && throw(ArgumentError("solve: the case $f = 0 is not supported"))
    Cases{Any}([c => r for (c, r) in zip(conds, rs)])
end

# the roots at which e stays defined in this case: an expression that must not vanish
# (a denominator, cos for tan) may vanish at a root for some parameter values
# (sin(x)/(cos(x) - a) at x = 0 when a = 1); each polynomial factor of its value there
# is decided, and a root where one vanishes is dropped
function _defined_roots(rs::AbstractVector, e::Sym, x::Sym, B::Assumptions)
    zs = unique(_zero_sets!(Sym[], e, x))
    (isempty(zs) || isempty(rs)) && return rs
    filter(rs) do r
        _hassymbol(r) || any(z -> any(u -> !_eq(u, x), _named(z)), zs) || return true
        for z in zs
            v = try subs(z, x => r) catch err; _isdivzero(err) && return false; rethrow() end
            _hassymbol(v) || continue
            # through a cubic's radicals: no polynomial factors in the parameters to decide
            # (and cancel over them took minutes; SymPy #25920)
            _has_high_root(v) && continue
            all(f -> _decide(Relation(:ne, f, ZERO), B), _param_factors(cancel(v), (x,))) || return false
        end
        true
    end
end

# |u| with u real (x real): the equation solved with |u| = u, keeping the roots where
# u ≥ 0, and with |u| = -u, keeping those where u ≤ 0 (|2x + 1| = x + 4; SymPy #28742)
function _solve_abs_split(e::Relation, x::Sym, B::Assumptions)
    _ismarkedreal(B, x) || return nothing
    g = e.lhs - e.rhs
    a = nothing
    _replace_walk(g) do t
        a === nothing && t.kind == FN && t.name === :abs && length(t.args) == 1 && has(t.args[1], x) &&
            _isreal(t.args[1], B) === true && (a = t)
        nothing
    end
    a === nothing && return nothing
    # a lone |u| = r is inverted directly (r ≥ 0 by cases)
    (_eq(e.lhs, a) && !has(e.rhs, x)) || (_eq(e.rhs, a) && !has(e.lhs, x)) && return nothing
    u = a.args[1]
    out = Sym[]
    for (w, op) in ((u, :ge), (-u, :le))
        rs = _solve_case(Relation(:eq, subs(e.lhs, a => w), subs(e.rhs, a => w)), x, B)
        rs isa AbstractVector || return nothing
        for r in rs
            _decide(Relation(op, expand(subs(u, x => r)), ZERO), B) || continue
            any(o -> iszero(expand(o - r)), out) || push!(out, r)
        end
    end
    _sort_roots!(out)
end

# Σ cᵢ uᵢ^(2kᵢ) = 0 with cᵢ > 0 and the uᵢ real (x real): every uᵢ = 0. The roots of one
# that can be solved, kept where the others vanish too ((x - 5)² + (sin x - sin 5)² = 0
# is x = 5; SymPy #24207). nothing when e is not such a sum.
function _solve_sum_of_squares(e::Relation, x::Sym, B::Assumptions)
    _ismarkedreal(B, x) || return nothing
    g = e.lhs - e.rhs
    g.kind == ADD && iszero(g.num) && length(g.args) >= 2 || return nothing
    us = Sym[]
    for (c, t) in zip(g.coefs, g.args)
        c > N0 || return nothing
        (t.kind == POW && t.args[2].kind == NUMBER && _isint(t.args[2].num) && iseven(_toint(t.args[2].num)) &&
         t.args[2].num > N0) || return nothing
        u = t.args[1]
        has(u, x) && _isreal(u, B) === true || return nothing
        push!(us, u)
    end
    # the most solvable first: polynomials, then the rest
    sort!(us; by = u -> (_ispoly_in(u, x) ? 0 : 1, _nodes(u)))
    for (i, u) in enumerate(us)
        rs = try _solve_case(Relation(:eq, u, ZERO), x, B) catch err; err isa ArgumentError || rethrow(); continue end
        rs isa AbstractVector || continue
        return [r for r in rs if all(j -> j == i || _vanishes(subs(us[j], x => r)), eachindex(us))]
    end
    nothing
end

# the one atom of n (a function call, power) that involves x, when there is one
function _single_x_atom(n::Sym, x::Sym)
    atoms = Set{Sym}()
    _poly_atoms!(atoms, n)
    bad = [a for a in atoms if has(a, x)]
    length(bad) == 1 && !_eq(bad[1], x) ? bad[1] : nothing
end

function _solve_case(e::Relation, x::Sym, B::Assumptions)
    if _has_fn(e.lhs - e.rhs, (:abs,))
        r = _solve_abs_split(e, x, B)
        r === nothing || return r
    end
    r = _solve_sum_of_squares(e, x, B)
    r === nothing || return r
    n, d = numer_denom(cancel(e.lhs - e.rhs))
    for f in _param_factors(d, (x,))
        _decide(Relation(:ne, f, ZERO), B) || return Sym[]
    end
    reduced = false
    if has(n, x) && _ispoly_in(n, x)
        while has(n, x)
            k = degree(Poly{Sym}(n, [x]), x)
            lc = coeff(n, x, k)
            all(f -> _decide(Relation(:ne, f, ZERO), B), _param_factors(lc, (x,))) && break
            n = expand(n - lc * x^k)              # lc vanishes in this case
            reduced = true
        end
    elseif has(n, x)
        # a polynomial in one function g of x (c cos(x) + c): its leading coefficients in g
        # likewise (c = 0: every x)
        g = _single_x_atom(n, x)
        if g !== nothing
            u = Sym(Symbol("_g", objectid(g)))
            m = expand(subs(n, g => u))
            if !has(m, x) && _ispoly_in(m, u)
                m0 = m
                while has(m, u)
                    k = degree(Poly{Sym}(m, [u]), u)
                    lc = coeff(m, u, k)
                    all(f -> _decide(Relation(:ne, f, ZERO), B), _param_factors(lc, (x, u))) && break
                    m = expand(m - lc * u^k)
                end
                if !has(m, u)
                    iszero(m) && throw(ArgumentError("solve: every value of $x satisfies $e"))
                    all(f -> _decide(Relation(:ne, f, ZERO), B), _param_factors(m, (x, u))) && return Sym[]
                    throw(ArgumentError("solve: every value of $x satisfies $e"))
                end
                _eq(m, m0) || (n = subs(m, u => g); reduced = true)
            end
        end
    end
    if !has(n, x)
        iszero(n) && throw(ArgumentError("solve: every value of $x satisfies $e"))
        all(f -> _decide(Relation(:ne, f, ZERO), B), _param_factors(n, (x,))) && return Sym[]
        throw(ArgumentError("solve: every value of $x satisfies $e"))
    end
    rs = reduced ? _in_domain(_solve_eq(n ~ 0, x), e.lhs - e.rhs, x) : _solve_eq(e, x)
    rs isa AbstractVector || return rs            # an inequality's set
    rs = _defined_roots(rs, e.lhs - e.rhs, x, B)
    reduced && return rs
    ps = _degenerate_subst([numer_denom(r)[2] for r in rs], B, (x,))
    ps === nothing && return rs
    _alternatives(ps) do p
        _solve_case(Relation(:eq, subs(e.lhs, p), subs(e.rhs, p)), x, B)
    end
end

# for a real x, float roots whose imaginary part is tiny beside their size (a multiple
# root split by rounding: 264 ± 6e-6im) are taken as real, with a warning
function _real_if_tiny(rs, x::Sym, A::Assumptions)
    rs isa Vector{Sym} || return rs
    get(A.dom, x, _CR) <= _RR || haskey(A.sets, x) || return rs
    out = Sym[]
    for r in rs
        w = (_hassymbol(r) || _isexactex(r)) ? nothing : try _evalc(r) catch; nothing end
        if w isa Complex && !iszero(imag(w)) && abs(imag(w)) <= cbrt(eps()) * abs(w)
            any(o -> !_hassymbol(o) && isapprox(N(o), real(w); rtol = cbrt(eps())), out) && continue
            @warn "solve: the root $r is taken as real ($(real(w))); its imaginary part is presumed rounding error"
            r = Sym(real(w))
        end
        push!(out, r)
    end
    out
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
