# Pattern matching: patterns with slots matched against expressions in canonical form
# (the design is notes/pattern_matching.md). A pattern is a Sym in which a slot is a
# private function node, _slot(x, flavor), and a head slot is _hslot(f, args...); it is
# built with SymJ's own arithmetic, so it is in the same canonical form as what it is
# matched against. It is compiled once to a tree of _PNode. Matching is backtracking in
# continuation-passing style: each matcher calls k(bindings) on success, and k returns
# true to stop the search.
#
# Sums and products are associative and commutative. Their pattern elements fall into
# literals (found by hash, no choice), other parts (each matching one subject element,
# tried only against elements with the right head), and slots, which take groups of the
# subject's elements: ~x and ~~~x one or more, ~!x and ~~x any number. Slots already
# bound take their value's elements; of the free ones, defaults go first (taking the
# number first, else nothing), then scaled slots (-~y), then slots with predicates, and
# the last takes what is left.

const _SL_ONE, _SL_DEF, _SL_STAR, _SL_PLUS, _SL_HEAD = 0x01, 0x02, 0x03, 0x04, 0x05

_slot(name::Symbol, fl::UInt8) = _mkfn_node(:_slot, Sym[Sym(name), Sym(Int(fl))])
_isslot(e::Sym) = e.kind == FN && e.name === :_slot
_hasslot(e::Sym) = (e.kind == FN && (e.name === :_slot || e.name === :_hslot)) ||
    any(_hasslot, e.args) || any(_hasslot, e.exps)

function _slotnames!(out::Set{Symbol}, e::Sym)
    e.kind == FN && (e.name === :_slot || e.name === :_hslot) && push!(out, e.args[1].name)
    foreach(a -> _slotnames!(out, a), e.args)
    foreach(a -> _slotnames!(out, a), e.exps)
    out
end

# ---- patterns from expressions

struct _PatBuild
    flavors::Dict{Symbol,UInt8}
    preds::Dict{Symbol,Any}
    order::Dict{Symbol,Int}      # where each slot first appears in the written pattern
end

# the argument of the outermost ~: (name, flavor, predicate expression or nothing)
function _slotspec(a)
    fl = _SL_ONE
    if a isa Expr && a.head === :call && length(a.args) == 2 && a.args[1] === :~
        a, fl = a.args[2], _SL_STAR
        if a isa Expr && a.head === :call && length(a.args) == 2 && a.args[1] === :~
            a, fl = a.args[2], _SL_PLUS
        end
    elseif a isa Expr && a.head === :call && length(a.args) == 2 && a.args[1] === :!
        a, fl = a.args[2], _SL_DEF
    end
    pred = nothing
    a isa Expr && a.head === :(::) && ((a, pred) = (a.args[1], a.args[2]))
    a isa Symbol || throw(ArgumentError("pattern: a slot is ~name (and ~x^2 is ~(x^2): write (~x)^2); got ~($a)"))
    a, fl, pred
end

function _pred_value(p)
    p isa Function && return p
    if p isa Symbol
        for M in (Main, Base, @__MODULE__)
            isdefined(M, p) && return getfield(M, p)
        end
    elseif p isa Expr && p.head === :.
        return Core.eval(Main, p)
    end
    throw(ArgumentError("pattern: unknown predicate $p"))
end

function _regslot!(P::_PatBuild, name::Symbol, fl::UInt8, pred)
    get(P.flavors, name, fl) == fl || throw(ArgumentError("pattern: ~$name is used as two kinds of slot"))
    P.flavors[name] = fl
    get!(P.order, name, length(P.order) + 1)
    pred === nothing || (P.preds[name] = _pred_value(pred))
    nothing
end

_to_pat(x::Sym, P) = x
_to_pat(x::Number, P) = Sym(x)
_to_pat(x::QuoteNode, P) = _to_pat(x.value, P)
_to_pat(x::Symbol, P) = x === :π || x === :pi ? PI : x === :ℯ ? E : x === :im ? Sym(im) : Sym(x)
_to_pat(x, P) = throw(ArgumentError("pattern: can't read $x"))

function _to_pat(ex::Expr, P::_PatBuild)
    ex.head === :call || throw(ArgumentError("pattern: can't read $ex"))
    f = ex.args[1]
    if f === :~ && length(ex.args) == 2
        name, fl, pred = _slotspec(ex.args[2])
        _regslot!(P, name, fl, pred)
        return _slot(name, fl)
    end
    args = Sym[_to_pat(a, P) for a in ex.args[2:end]]
    if f isa Expr && f.head === :call && length(f.args) == 2 && f.args[1] === :~
        name, fl, pred = _slotspec(f.args[2])
        fl == _SL_ONE || throw(ArgumentError("pattern: a head slot is written (~f)(…)"))
        _regslot!(P, name, _SL_HEAD, pred)
        return _mkfn_node(:_hslot, Sym[Sym(name); args])
    end
    f isa Function && return Sym(f(args...))
    f isa Symbol || throw(ArgumentError("pattern: can't read the function $f"))
    f === :// && length(args) == 2 && return args[1] / args[2]
    (f === :+ || f === :*) && length(args) == 1 && return args[1]
    fn = _pat_function(f)
    fn === nothing && return apply_function(f, args)
    try
        return Sym(fn(args...))
    catch err
        err isa MethodError || rethrow()
        return apply_function(f, args)
    end
end

function _pat_function(f::Symbol)
    M = @__MODULE__
    isdefined(M, f) || return nothing
    g = getfield(M, f)
    g isa Function && !(g isa Sym) ? g : nothing
end

# the function a head slot binds, for a function node's name
const _HEADS = Dict{Symbol,Any}()
_head_function(name::Symbol) = get!(_HEADS, name) do
    g = _pat_function(name)
    g === nothing || g isa AbstractSymbolic ? SymFunction(name) : g
end

# ---- compiled patterns

const _P_LIT, _P_SLOT, _P_HEAD, _P_FN, _P_POW, _P_ADD, _P_MUL = 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07

struct _PNode
    kind::UInt8
    sym::Sym                 # the pattern
    name::Symbol             # slot or function name
    flavor::UInt8
    pred::Any                # nothing, or the slot's predicate
    num::SymNumber           # a product's coefficient; a slot's scale in a sum (-~y: -1)
    kids::Vector{_PNode}     # arguments; [base, exponent]; a sum's or product's other parts
    lits::Vector{Sym}        # a sum's or product's parts without slots
    slots::Vector{_PNode}    # a sum's or product's slots
    divides::Bool            # a product's coefficient divides into its slots
end
_PN(kind, sym; name = NONAME, flavor = 0x00, pred = nothing, num = N1, kids = _PNode[],
    lits = Sym[], slots = _PNode[], divides = false) =
    _PNode(kind, sym, name, flavor, pred, num, kids, lits, slots, divides)

_zero_ok(sl::_PNode) = sl.flavor == _SL_DEF || sl.flavor == _SL_STAR

function _pcompile(p::Sym, P::_PatBuild)
    _hasslot(p) || return _PN(_P_LIT, p)
    if _isslot(p)
        name = p.args[1].name
        return _PN(_P_SLOT, p; name, flavor = UInt8(_toint(p.args[2].num)), pred = get(P.preds, name, nothing))
    end
    k = p.kind
    if k == FN && p.name === :_hslot
        name = p.args[1].name
        return _PN(_P_HEAD, p; name, pred = get(P.preds, name, nothing), kids = [_pcompile(a, P) for a in p.args[2:end]])
    elseif k == FN
        return _PN(_P_FN, p; name = p.name, kids = [_pcompile(a, P) for a in p.args])
    elseif k == POW
        return _PN(_P_POW, p; kids = [_pcompile(p.args[1], P), _pcompile(p.args[2], P)])
    end
    lits, others, slots = Sym[], _PNode[], _PNode[]
    if k == ADD
        p.num == N0 || push!(lits, Sym(p.num))
        for (c, a) in zip(p.coefs, p.args)
            if _isslot(a)
                s = _pcompile(a, P)
                push!(slots, _PN(_P_SLOT, a; name = s.name, flavor = s.flavor, pred = s.pred, num = c))
            else
                t = c == N1 ? a : Sym(c) * a
                _hasslot(t) ? push!(others, _pcompile(t, P)) : push!(lits, t)
            end
        end
        return _PN(_P_ADD, p; kids = _written_order(others, P), lits, slots = _slot_order(slots, P))
    end
    for (b, e) in zip(p.args, p.exps)
        if _isslot(b) && _eq(e, ONE)
            push!(slots, _pcompile(b, P))
        else
            f = _eq(e, ONE) ? b : b^e
            _hasslot(f) ? push!(others, _pcompile(f, P)) : push!(lits, f)
        end
    end
    _PN(_P_MUL, p; num = p.num, kids = _written_order(others, P), lits, slots = _slot_order(slots, P), divides = !isempty(slots))
end

# defaults first, then scaled slots, then slots with predicates, each as written: so the
# number in 3x goes to ~b in ~b*~x, whatever the slots' names
_slot_order(slots, P) = sort(slots; by = s -> (s.flavor == _SL_DEF ? 0 : s.num != N1 ? 1 : s.pred !== nothing ? 2 : 3, P.order[s.name]))
_written_order(kids, P) = sort(kids; by = q -> minimum((P.order[n] for n in _slotnames!(Set{Symbol}(), q.sym)); init = typemax(Int)))

"""
    SymJ.Pattern

A compiled pattern, made by [`SymJ.pattern`](@ref) and used with `match`, `eachmatch`
and `occursin`.
"""
struct Pattern
    ex::Sym
    node::_PNode
end
Base.show(io::IO, p::Pattern) = print(io, "pattern(", p.ex, ")")

"""
    SymJ.pattern(ex)

A pattern from the expression `ex`, usually quoted, in which `~x` marks a slot:

| Slot | In a function's arguments | In a sum or product |
|:--|:--|:--|
| `~x` | one argument | one or more terms (factors), bound to their sum (product) |
| `~!x` | one argument | any number; none gives 0 (1) |
| `~~x` | any number, bound to a vector | any number; none gives 0 (1) |
| `~~~x` | one or more, bound to a vector | one or more |

In an exponent, `~!n` is 1 when there is no power. A head slot `(~f)(~x)` matches any
function application and binds the function. `~x::pred` (with any flavor) accepts only
values for which `pred` returns `true`; a predicate's name is looked up in `Main`, then
`Base`, then SymJ, or a function can be interpolated (`:(~x::\$f)`).

The pattern is put in canonical form, as any expression is (`~x + ~x` is `2*~x`,
`exp(log(~x))` is `~x`), and means what that form means. A number divides into a slot:
`~x - ~y` matches `a - 2b` with `y => 2b`.

```julia
p = SymJ.pattern(:(~!a * sin(~x) + ~~r))
match(p, 3sin(t) + 1)          # Dict(:a => 3, :x => t, :r => 1)
eachmatch(SymJ.pattern(:(~x + ~y)), a + b + c)   # six matches
occursin(SymJ.pattern(:(sin(~x)^2)), 1 + sin(t)^2)   # true
```

Not exported.
"""
function pattern(ex)
    P = _PatBuild(Dict{Symbol,UInt8}(), Dict{Symbol,Any}(), Dict{Symbol,Int}())
    p = _to_pat(ex, P)
    gone = setdiff(keys(P.flavors), _slotnames!(Set{Symbol}(), p))
    isempty(gone) || throw(ArgumentError("pattern: the slot ~$(first(gone)) is lost in the pattern's canonical form, $p"))
    Pattern(p, _pcompile(p, P))
end

# ---- continuations: what is left to match once a part has matched, as a linked list
# of concrete goals rather than closures, so the search is fully typed

const _G_DONE, _G_MATCH, _G_ARGS, _G_OTHERS, _G_SLOTS = 0x01, 0x02, 0x03, 0x04, 0x05

struct _Goal
    kind::UInt8
    p::_PNode                    # _G_MATCH, _G_OTHERS, _G_SLOTS: the pattern
    s::Sym                       # _G_MATCH: the subject
    ps::Vector{_PNode}           # _G_ARGS: the patterns; _G_SLOTS: the slots in order
    as::Vector{Sym}              # _G_ARGS: the arguments; _G_OTHERS, _G_SLOTS: the elements
    i::Int
    j::Int
    used::BitVector
    f::Any                       # _G_DONE: called with the bindings, true to stop
    next::Union{_Goal,Nothing}
end

const _NOPN = _PN(_P_LIT, ZERO)
const _NOUSED = falses(0)
const _NOPS = _PNode[]
const _NOAS = Sym[]
_gdone(f) = _Goal(_G_DONE, _NOPN, ZERO, _NOPS, _NOAS, 0, 0, _NOUSED, f, nothing)
_gmatch(p, s, k) = _Goal(_G_MATCH, p, s, _NOPS, _NOAS, 0, 0, _NOUSED, nothing, k)
_gargs(ps, i, as, j, k) = _Goal(_G_ARGS, _NOPN, ZERO, ps, as, i, j, _NOUSED, nothing, k)
_gothers(p, i, s, elems, used, k) = _Goal(_G_OTHERS, p, s, _NOPS, elems, i, 0, used, nothing, k)
_gslots(p, order, i, s, elems, used, k) = _Goal(_G_SLOTS, p, s, order, elems, i, 0, used, nothing, k)

function _run(g::_Goal, B)::Bool
    K = g.kind
    K == _G_DONE && return g.f(B)::Bool
    k = g.next::_Goal
    K == _G_MATCH && return _pm(g.p, g.s, B, k)
    K == _G_ARGS && return _pm_args(g.ps, g.i, g.as, g.j, B, k)
    K == _G_OTHERS && return _pm_others(g.p, g.i, g.s, g.as, g.used, B, k)
    _pm_slots(g.p, g.ps, g.i, g.s, g.as, g.used, B, k)
end

# ---- bindings: a stack of (name, value), cut back on failure

struct _Binds
    names::Vector{Symbol}
    vals::Vector{Any}
end
_Binds() = _Binds(Symbol[], Any[])

function _getb(B::_Binds, n::Symbol)
    for i in eachindex(B.names)
        B.names[i] === n && return B.vals[i]
    end
    nothing
end

_same(a::Sym, b::Sym) = _eq(a, b)
_same(a::Vector{Sym}, b::Vector{Sym}) = length(a) == length(b) && all(i -> _eq(a[i], b[i]), eachindex(a))
_same(a, b) = isequal(a, b)

_predok(f, v) = try f(v) === true catch; false end

function _bind(B::_Binds, name::Symbol, v, pred, k::_Goal)
    old = _getb(B, name)
    old === nothing || return _same(old, v) && _run(k, B)
    pred === nothing || _predok(pred, v) || return false
    push!(B.names, name)
    push!(B.vals, v)
    r = _run(k, B)
    pop!(B.names)
    pop!(B.vals)
    r
end

_bdict(B::_Binds) = Dict{Symbol,Any}(zip(B.names, B.vals))

# ---- matching

function _pm(p::_PNode, s::Sym, B::_Binds, k::_Goal)::Bool
    K = p.kind
    K == _P_LIT && return _eq(p.sym, s) && _run(k, B)
    K == _P_SLOT && return _bind(B, p.name, s, p.pred, k)
    if K == _P_FN
        (s.kind == FN && s.name === p.name) || return false
        return _pm_args(p.kids, 1, s.args, 1, B, k)
    elseif K == _P_HEAD
        s.kind == FN || return false
        return _bind(B, p.name, _head_function(s.name), p.pred, _gargs(p.kids, 1, s.args, 1, k))
    elseif K == _P_POW
        return _pm_pow(p, s, B, k)
    end
    _pm_ac(p, s, B, k)
end

# a function's arguments, in order; sequence slots take runs, shortest first
function _pm_args(ps::Vector{_PNode}, i::Int, as::Vector{Sym}, j::Int, B::_Binds, k::_Goal)::Bool
    i > length(ps) && return j > length(as) && _run(k, B)
    p = ps[i]
    if p.kind == _P_SLOT && (p.flavor == _SL_STAR || p.flavor == _SL_PLUS)
        for n in (p.flavor == _SL_PLUS ? 1 : 0):(length(as) - j + 1)
            _bind(B, p.name, as[j:(j + n - 1)], p.pred, _gargs(ps, i + 1, as, j + n, k)) && return true
        end
        return false
    end
    j > length(as) && return false
    _pm(p, as[j], B, _gargs(ps, i + 1, as, j + 1, k))
end

# a power, then (for a default exponent ~!n) the whole subject to the power 1
function _pm_pow(p::_PNode, s::Sym, B::_Binds, k::_Goal)::Bool
    pb, pe = p.kids
    s.kind == POW && _pm(pb, s.args[1], B, _gmatch(pe, s.args[2], k)) && return true
    pe.kind == _P_SLOT && pe.flavor == _SL_DEF || return false
    _bind(B, pe.name, ONE, pe.pred, _gmatch(pb, s, k))
end

# the terms of a sum, its constant last; anything else is a one-term sum
function _add_parts(s::Sym)
    if s.kind == ADD
        out = Sym[c == N1 ? a : Sym(c) * a for (c, a) in zip(s.coefs, s.args)]
        s.num == N0 || push!(out, Sym(s.num))
        return out
    end
    s.kind == NUMBER && s.num == N0 ? Sym[] : Sym[s]
end

# the factors of a product and its coefficient
function _mul_parts(s::Sym)
    s.kind == MUL && return Sym[_eq(e, ONE) ? b : b^e for (b, e) in zip(s.args, s.exps)], s.num
    s.kind == NUMBER && return Sym[], s.num
    Sym[s], N1
end

# the subject's elements for a sum or product pattern, the number last; nothing when a
# product's coefficient can't match
function _ac_elems(p::_PNode, s::Sym)
    p.kind == _P_ADD && return _add_parts(s)
    out, d = _mul_parts(s)
    if p.divides
        q = d / p.num
        q == N1 || push!(out, Sym(q))
    else
        d == p.num || return nothing
    end
    out
end

function _pm_ac(p::_PNode, s::Sym, B::_Binds, k::_Goal)::Bool
    elems = _ac_elems(p, s)
    elems === nothing && return false
    used = falses(length(elems))
    for l in p.lits
        j = findfirst(i -> !used[i] && _eq(elems[i], l), eachindex(elems))
        j === nothing && return false
        used[j] = true
    end
    _pm_others(p, 1, s, elems, used, B, k)
end

function _headok(q::_PNode, e::Sym)
    q.kind == _P_FN && return e.kind == FN && e.name === q.name
    q.kind == _P_HEAD && return e.kind == FN
    q.kind == _P_POW && return e.kind == POW || (q.kids[2].kind == _P_SLOT && q.kids[2].flavor == _SL_DEF)
    true
end

# the parts that are not slots, each matching one element
function _pm_others(p::_PNode, i::Int, s::Sym, elems::Vector{Sym}, used::BitVector, B::_Binds, k::_Goal)::Bool
    if i > length(p.kids)
        order = sort(p.slots; by = sl -> _getb(B, sl.name) === nothing)   # bound slots first
        return _pm_slots(p, order, 1, s, elems, used, B, k)
    end
    q = p.kids[i]
    for j in eachindex(elems)
        (used[j] || !_headok(q, elems[j])) && continue
        used[j] = true
        r = _pm(q, elems[j], B, _gothers(p, i + 1, s, elems, used, k))
        used[j] = false
        r && return true
    end
    false
end

# a slot's value from a group of elements: their sum (divided by the slot's scale) or
# product. The group's terms (factors) are a subsequence of the subject's, so already in
# canonical order: the node is built from them directly.
function _group_value(p::_PNode, sl::_PNode, s::Sym, elems::Vector{Sym}, G::Vector{Int})
    issorted(G) || (G = sort(G))
    if p.kind == _P_ADD
        if s.kind == ADD
            n, c = length(s.args), N0
            idx = Int[]
            for j in G
                j > n ? (c = s.num) : push!(idx, j)
            end
            v = _mkadd(c, s.args[idx], s.coefs[idx])
        else
            v = isempty(G) ? ZERO : s
        end
        return sl.num == N1 ? v : _scale(v, N1 / sl.num)
    end
    if s.kind == MUL
        n, c = length(s.args), N1
        idx = Int[]
        for j in G
            j > n ? (c = elems[j].num) : push!(idx, j)
        end
        return _mkmul(c, s.args[idx], s.exps[idx])
    end
    prod((elems[j] for j in G); init = ONE)
end

_isneg(n::SymNumber) = try n < N0 catch; false end
_isneg(e::Sym) = (e.kind == NUMBER || e.kind == MUL) && _isneg(e.num)

function _pm_slots(p::_PNode, order::Vector{_PNode}, i::Int, s::Sym, elems::Vector{Sym}, used::BitVector, B::_Binds, k::_Goal)::Bool
    i > length(order) && return all(used) && _run(k, B)
    sl = order[i]
    old = _getb(B, sl.name)
    if old !== nothing                       # a bound slot takes its value's elements
        old isa Sym || return false
        parts = p.kind == _P_ADD ? _add_parts(sl.num == N1 ? old : Sym(sl.num) * old) : _mul_slot_parts(old)
        isempty(parts) && !_zero_ok(sl) && return false
        idx = Int[]
        for v in parts
            j = findfirst(jj -> !used[jj] && !(jj in idx) && _eq(elems[jj], v), eachindex(elems))
            j === nothing && return false
            push!(idx, j)
        end
        used[idx] .= true
        r = _pm_slots(p, order, i + 1, s, elems, used, B, k)
        used[idx] .= false
        return r
    end
    free = [j for j in eachindex(elems) if !used[j]]
    sort!(free; by = j -> elems[j].kind != NUMBER)    # the number first: ~a*~x against 2x is a => 2
    if i == length(order)                    # the last slot takes what is left
        isempty(free) && !_zero_ok(sl) && return false
        return _bind(B, sl.name, _group_value(p, sl, s, elems, free), sl.pred, k)
    end
    need = count(s -> !_zero_ok(s), @view order[(i + 1):end])
    maxn = length(free) - need
    maxn < 0 && return false
    sl.num != N1 && sort!(free; by = j -> _isneg(elems[j]) != _isneg(sl.num))   # -~y: negative terms first
    try_group(G) = begin
        used[G] .= true
        r = _bind(B, sl.name, _group_value(p, sl, s, elems, G), sl.pred, _gslots(p, order, i + 1, s, elems, used, k))
        used[G] .= false
        r
    end
    numj = sl.flavor == _SL_DEF ? findfirst(j -> elems[j].kind == NUMBER, free) : nothing
    if _zero_ok(sl)
        numj === nothing || maxn < 1 || try_group([free[numj]]) && return true
        try_group(Int[]) && return true
    end
    for m in 1:maxn
        _combinations(free, m) do G
            numj !== nothing && m == 1 && G[1] == free[numj] && return false
            try_group(G)
        end && return true
    end
    false
end

# a bound value's elements in a product: its factors and its coefficient
function _mul_slot_parts(v::Sym)
    fs, d = _mul_parts(v)
    d == N1 || push!(fs, Sym(d))
    fs
end

# f on each m-element subset of xs, in lexicographic order, until f returns true
function _combinations(f, xs::Vector{Int}, m::Int)
    n = length(xs)
    idx = collect(1:m)
    G = Vector{Int}(undef, m)
    while true
        for t in 1:m
            G[t] = xs[idx[t]]
        end
        f(copy(G)) && return true
        t = m
        while t >= 1 && idx[t] == n - m + t
            t -= 1
        end
        t == 0 && return false
        idx[t] += 1
        for u in (t + 1):m
            idx[u] = idx[u - 1] + 1
        end
    end
end

# ---- the interface

_matchsubject(ex) = ex isa Sym ? ex : Sym(ex)

"""
    match(p::SymJ.Pattern, ex)

The first match of the pattern `p` against `ex`, as a `Dict` from slot names to values,
or `nothing`. See [`SymJ.pattern`](@ref).
"""
function Base.match(p::Pattern, ex)
    out = Ref{Any}(nothing)
    _pm(p.node, _matchsubject(ex), _Binds(), _gdone(B -> (out[] = _bdict(B); true)))
    out[]
end

"""
    eachmatch(p::SymJ.Pattern, ex)

Every match of the pattern `p` against `ex`, as a vector of `Dict`s, in the order the
search finds them. See [`SymJ.pattern`](@ref).
"""
function Base.eachmatch(p::Pattern, ex)
    out = Dict{Symbol,Any}[]
    _pm(p.node, _matchsubject(ex), _Binds(), _gdone(B -> (d = _bdict(B); any(o -> isequal(o, d), out) || push!(out, d); false)))
    out
end

"""
    occursin(p::SymJ.Pattern, ex)

Whether the pattern `p` matches `ex` or a part of it: a term of a sum, a factor of a
product, an argument. See [`SymJ.pattern`](@ref).
"""
Base.occursin(p::Pattern, ex) = _occurs(p.node, _matchsubject(ex))

function _occurs(p::_PNode, s::Sym)
    _pm(p, s, _Binds(), _gdone(_ -> true)) && return true
    parts = s.kind == ADD ? _add_parts(s) : s.kind == MUL ? first(_mul_parts(s)) : s.args
    s.kind == ADD || s.kind == MUL ? any(t -> !_eq(t, s) && _occurs(p, t), parts) : any(t -> _occurs(p, t), parts)
end
