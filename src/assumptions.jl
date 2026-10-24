# Assumptions: facts about symbols, given for a call (`assume = …`) or a block
# (`assuming(…) do … end`), never stored on the symbols (see design.md).
#
# An `Assumptions` value holds the facts normalized:
#   - equalities solved for a symbol and substituted (`subst`), into the other facts
#     and into every query;
#   - for each symbol, the smallest number set known to contain it (`dom`), the
#     largest known not to (`notin`), and, for a real symbol, the set of its values
#     (`sets`, a SymSet) or values it does not take (`excluded`, for a complex one);
#   - relations among several symbols as linear constraints over atoms, a nonlinear
#     monomial counting as an atom (`a*b > 0` is linear in `a*b`).
# Queries are three-valued. The possible signs of a real expression are found by
# propagation through its structure (`_smask`), then, when that leaves more than one,
# by Fourier–Motzkin elimination over the linear constraints and the bounds the atoms
# get from their symbols' sets (`_fm_feasible`).

## number sets

"""
    NumberSet

One of the number sets `ℤ ⊂ ℚ ⊂ ℝ ⊂ ℂ`, also named `Integers`, `Rationals`, `Reals` and
`Complexes`, for facts such as `n ∈ ℤ` (or `n in Integers`) in [`Assumptions`](@ref). On
a number, `∈` is a `Bool` (`2 ∈ ℤ`); on a symbol it is a fact.
"""
struct NumberSet
    rank::Int8                        # 1 ℤ, 2 ℚ, 3 ℝ, 4 ℂ
end
const ℤ = NumberSet(1)
const ℚ = NumberSet(2)
const ℝ = NumberSet(3)
const ℂ = NumberSet(4)
# ASCII names (Mathematica's)
const Integers, Rationals, Reals, Complexes = ℤ, ℚ, ℝ, ℂ
const _SETNAMES = ("ℤ", "ℚ", "ℝ", "ℂ")
const _ZR, _QR, _RR, _CR = Int8(1), Int8(2), Int8(3), Int8(4)
Base.show(io::IO, S::NumberSet) = print(io, _SETNAMES[S.rank])

Base.in(x::Number, S::NumberSet) =
    S.rank == _CR || (isreal(x) && (S.rank == _RR || (isfinite(real(x)) &&
        !(x isa AbstractIrrational) && (S.rank == _QR || isinteger(real(x))))))

"""
    Membership

The fact `x ∈ S` (or `x ∉ S`) for a symbol `x` and a [`NumberSet`](@ref) `S`.
"""
struct Membership
    x::Sym
    set::NumberSet
    negated::Bool
end
Base.:!(m::Membership) = Membership(m.x, m.set, !m.negated)
Base.:(==)(a::Membership, b::Membership) = a.set == b.set && a.negated == b.negated && _eq(a.x, b.x)
Base.hash(m::Membership, h::UInt) = hash(m.x, hash(m.set.rank, hash(m.negated, h)))
Base.show(io::IO, m::Membership) = print(io, m.x, m.negated ? " ∉ " : " ∈ ", m.set)

function Base.in(x::Sym, S::NumberSet)
    x.kind == NUMBER && return S.rank >= _QR || isone(_denominator(x.num)) ||
        (isfloat(x.num) && isinteger(Float64(x.num)))
    _hassymbol(x) && return Membership(x, S, false)
    S.rank == _CR && return true
    r = _isreal_const(x)
    r === nothing && throw(ArgumentError("cannot decide whether $x ∈ $S"))
    S.rank == _RR && return r
    r || return false
    q = _const_rational(x)
    q === false && return false
    S.rank == _QR && q === true && return true
    v = _evalc(x)
    S.rank == _ZR && !isinteger(round(real(v); digits = 9)) && return false
    S.rank == _ZR && q === true && return true          # a rational at an integer
    throw(ArgumentError("cannot decide whether $x ∈ $S"))
end

# Is a constant rational? true or false where that is provable, else nothing:
#   - a product of powers of rationals: rational exactly when, over a coprime basis of
#     the bases (gcds, no factoring), each base's power is (√2·√3·√6 = 6, √2·√3 is not);
#   - not real: not rational;
#   - a sum: rational terms and one irrational term is irrational (1 + √2, (1 + √5)/2);
#     with several irrational terms, unknown;
#   - π or ℯ to a rational power times an algebraic number is transcendental (one
#     of them; π·ℯ is unknown). Functions such as log(2): unknown.
function _const_rational(e::Sym)::Union{Nothing,Bool}
    e = expand(e)
    e.kind == NUMBER && return isexact(e.num)
    _isreal_const(e) === false && return false             # not even real
    irrational = 0
    for t in (e.kind == ADD ? e.args : (e,))
        r = _monomial_rational(t)
        r === nothing && return nothing
        r || (irrational += 1)
    end
    irrational == 0 ? true : irrational == 1 ? false : nothing
end

function _monomial_rational(t::Sym)::Union{Nothing,Bool}
    fs = t.kind == MUL ? collect(zip(t.args, t.exps)) : t.kind == POW ? [(t.args[1], t.args[2])] :
         t.kind == CONSTANT ? [(t, ONE)] : return nothing
    exps = Dict{BigInt,Rational{BigInt}}()               # integer base => exponent
    transc = Symbol[]
    for (b, x) in fs
        (x.kind == NUMBER && isexact(x.num)) || return nothing
        q = Rational{BigInt}(_numerator(x.num), _denominator(x.num))
        if b.kind == CONSTANT && b.name in (:pi, :E)
            iszero(q) || push!(transc, b.name)
        elseif b.kind == NUMBER && isexact(b.num) && b.num > N0
            for (n, s) in ((_numerator(b.num), q), (_denominator(b.num), -q))
                isone(n) || (exps[BigInt(n)] = get(exps, BigInt(n), 0) + s)
            end
        else
            return nothing
        end
    end
    length(unique(transc)) > 1 && return nothing
    isempty(transc) || return false                     # π^q or ℯ^q times an algebraic number
    _coprime_basis!(exps)
    all(((n, q),) -> isinteger(q) || _perfect_root(n, Int(denominator(q))) !== nothing, exps)
end

# bases made pairwise coprime, exponents carried along (a^s b^t with g = gcd(a, b) is
# g^(s+t) (a/g)^s (b/g)^t)
function _coprime_basis!(exps::Dict{BigInt,Rational{BigInt}})
    while true
        ks = collect(keys(exps))
        pair = nothing
        for i in eachindex(ks), j in (i+1):length(ks)
            g = gcd(ks[i], ks[j])
            g > 1 && (pair = (ks[i], ks[j], g); break)
        end
        pair === nothing && break
        a, b, g = pair
        s, t = pop!(exps, a), pop!(exps, b)
        for (n, q) in ((g, s + t), (a ÷ g, s), (b ÷ g, t))
            isone(n) && continue
            exps[n] = get(exps, n, 0) + q
        end
        filter!(p -> !iszero(p.second), exps)
    end
    exps
end

## compound facts: & and | (Or only about one symbol), ! on any fact

const _NOTOP = Dict(:lt => :ge, :le => :gt, :gt => :le, :ge => :lt, :eq => :ne, :ne => :eq)
Base.:!(r::Relation) = Relation(_NOTOP[r.op], r.lhs, r.rhs)

struct FactAnd
    facts::Vector{Any}
end
struct FactOr
    facts::Vector{Any}
end
const _Fact = Union{Relation,Membership,FactAnd,FactOr}
_parts(f, T) = f isa T ? f.facts : Any[f]
Base.:&(a::_Fact, b::_Fact) = FactAnd(vcat(_parts(a, FactAnd), _parts(b, FactAnd)))
Base.:|(a::_Fact, b::_Fact) = FactOr(vcat(_parts(a, FactOr), _parts(b, FactOr)))
Base.:!(f::FactAnd) = FactOr(Any[!g for g in f.facts])
Base.:!(f::FactOr) = FactAnd(Any[!g for g in f.facts])
_showfact(io, f, sep) = join(io, (g isa Union{FactAnd,FactOr} ? "($g)" : string(g) for g in f.facts), sep)
Base.show(io::IO, f::FactAnd) = _showfact(io, f, " & ")
Base.show(io::IO, f::FactOr) = _showfact(io, f, " | ")

## the normalized facts

# Σ a[t] t + k > 0 (strict) or ≥ 0, t atoms, k a real constant
struct _ACon
    a::Dict{Sym,SymNumber}
    k::Sym
    strict::Bool
end

"""
    Assumptions(facts...)

Facts about symbols, normalized once so that operations taking them (`assume = A`, or
[`assuming`](@ref)) can reuse the work. A fact is:

- a relation: `x ≫ 0`, `a ≦ b + 1`, `x ~ 2`, `a ≶ 0` (≠). An order relation makes the
  symbols in it real.
- membership in a number set: `n ∈ ℤ`, `x ∈ ℝ`, `x ∉ ℤ` (`ℤ ⊂ ℚ ⊂ ℝ ⊂ ℂ`). Symbols are
  complex unless a fact says otherwise.
- a combination: `f & g`; `f | g` and `!f` for facts about one symbol
  (`!(0 ≦ x) | (x ≫ 1)`); `!` on a single relation or membership.

Several facts (or a vector of them) are all assumed. Relations among several symbols
are used as linear relations, a nonlinear monomial such as `a*b` taken as a quantity
of its own. Contradictory facts throw an `ArgumentError`. `Assumptions(A, facts...)`
adds facts to `A`.

```julia
A = Assumptions(x ≫ 0, n ∈ ℤ, a ≦ b)
```
"""
struct Assumptions
    facts::Vector{Any}                      # as given, flattened
    subst::Dict{Sym,Sym}                    # symbol => value, from equalities
    dom::Dict{Sym,Int8}                     # smallest set known to contain the symbol
    notin::Dict{Sym,Int8}                   # largest set known not to
    sets::Dict{Sym,SymSet}                  # a real symbol's values, when restricted
    excluded::Dict{Sym,Vector{Sym}}         # values a symbol does not take, when not real
    nonzero::Vector{Sym}                    # expressions in several symbols known ≠ 0
    lin::Vector{_ACon}                      # relations among several symbols
    cache::Dict{Sym,Union{Nothing,UInt8}}   # sign masks found
end

# with Julia's operators, as the facts are written in @assume
Base.show(io::IO, A::Assumptions) = (print(io, "Assumptions("); join(io, map(_condstring, A.facts), ", "); print(io, ")"))

Assumptions(A::Assumptions, facts...) = Assumptions(A.facts..., facts...)

function Assumptions(facts...)
    flat = _flatten!(Any[], facts)
    A = _build(flat)
    A === nothing && throw(ArgumentError("the assumptions are contradictory: " * join(string.(flat), ", ")))
    A
end

# the normalized facts, or nothing when they are contradictory
function _build(flat::Vector{Any})
    A = Assumptions(flat, Dict{Sym,Sym}(), Dict{Sym,Int8}(), Dict{Sym,Int8}(), Dict{Sym,SymSet}(),
                    Dict{Sym,Vector{Sym}}(), Sym[], _ACon[], Dict{Sym,Union{Nothing,UInt8}}())
    # membership facts first: equalities then avoid solving for those symbols
    mems = [f for f in flat if f isa Membership]
    for m in mems
        _isnamed(m.x) || throw(ArgumentError("membership facts are about symbols, not $(m.x)"))
    end
    # equalities: solved for a symbol where they can be, and substituted
    rels = Any[]
    for f in flat
        f isa Membership && continue
        if f isa Relation && f.op === :eq
            _add_equality!(A, f, Set(m.x for m in mems)) || push!(rels, f)
        else
            push!(rels, f)
        end
    end
    sub(e) = isempty(A.subst) ? e : expand(subs(e, A.subst...))
    for m in mems
        x = get(A.subst, m.x, nothing)
        if x !== nothing
            v = x in m.set
            v isa Bool || throw(ArgumentError("$(m.x) ∈ $(m.set) together with $(m.x) ~ $x is not supported"))
            v == m.negated && return nothing
            continue
        end
        if m.negated
            A.notin[m.x] = max(get(A.notin, m.x, Int8(0)), m.set.rank)
        else
            A.dom[m.x] = min(get(A.dom, m.x, _CR), m.set.rank)
        end
    end
    # the others, substituted; ≠ last, once it is known which symbols are real
    nes = Any[]
    for f in rels
        g = _subsfact(f, sub)
        if g isa Relation && g.op === :ne
            push!(nes, g)
        else
            _add_fact!(A, g) || return nothing
        end
    end
    for g in nes
        _add_ne!(A, g) || return nothing
    end
    _check!(A) || return nothing
    A
end

function _flatten!(out, fs)
    for f in fs
        if f isa Union{Tuple,AbstractVector}
            _flatten!(out, f)
        elseif f isa FactAnd
            _flatten!(out, f.facts)
        elseif f isa Bool
            f || push!(out, Relation(:eq, ONE, ZERO))     # false: contradicts
        elseif f isa Assumptions
            _flatten!(out, f.facts)
        elseif f isa _Fact
            push!(out, f)
        else
            throw(ArgumentError("not a fact: $f (use relations such as x ≫ 0, or x ∈ ℝ)"))
        end
    end
    out
end

_subsfact(f::Relation, sub) = Relation(f.op, sub(f.lhs), sub(f.rhs))
_subsfact(f::Membership, sub) = f
_subsfact(f::FactAnd, sub) = FactAnd(Any[_subsfact(g, sub) for g in f.facts])
_subsfact(f::FactOr, sub) = FactOr(Any[_subsfact(g, sub) for g in f.facts])

_named(e::Sym) = sort!(collect(free_symbols(e) ∪ free_parameters(e)); lt = _lt)
_ismarkedreal(A, x) = get(A.dom, x, _CR) <= _RR || haskey(A.sets, x)
_markreal!(A, x) = (A.dom[x] = min(get(A.dom, x, _CR), _RR))
_restrict!(A, x, S::SymSet) = (A.sets[x] = haskey(A.sets, x) ? A.sets[x] ∩ S : S; _markreal!(A, x))

# an equality solved for a symbol, added to subst (for a symbol with a membership fact,
# only a constant value, which the membership is checked against); false when it cannot
# be (nonlinear in several symbols: kept as a relation)
function _add_equality!(A::Assumptions, f::Relation, avoid)
    d = expand(isempty(A.subst) ? f.lhs - f.rhs : subs(f.lhs - f.rhs, A.subst...))
    vs = _named(d)
    isempty(vs) && return false                      # a constant: decided with the others
    # x ~ rhs is taken as giving x
    x = isempty(A.subst) ? f.lhs : expand(subs(f.lhs, A.subst...))
    if _isnamed(x) && !has(f.rhs, x)
        v = expand(isempty(A.subst) ? f.rhs : subs(f.rhs, A.subst...))
        if !(x in avoid) || !_hassymbol(v)
            _setsubst!(A, x, v)
            return true
        end
    end
    if length(vs) == 1
        x = only(vs)
        rs = try _solve_eq(d ~ 0, x) catch; return false end
        length(rs) == 1 || return false              # several values: a set
        _setsubst!(A, x, only(rs))
        return true
    end
    lf = _linform(d)
    lf === nothing && return false
    a, k = lf
    for x in vs
        c = get(a, x, nothing)
        c === nothing && continue
        rest = d - Sym(c) * x
        (has(rest, x) || (x in avoid && _hassymbol(rest))) && continue
        _setsubst!(A, x, expand(-rest / Sym(c)))
        return true
    end
    false
end
function _setsubst!(A::Assumptions, x::Sym, v::Sym)
    for (y, w) in A.subst
        A.subst[y] = expand(subs(w, x => v))
    end
    A.subst[x] = v
end

# add a fact (substituted); false if it is false
function _add_fact!(A::Assumptions, f)
    if f isa FactOr || f isa FactAnd
        vs = _factsyms(f)
        isempty(vs) && return _constfact(f) !== false        # (a substituted by a value)
        if length(vs) != 1
            f isa FactAnd && return all(g -> _add_fact!(A, g), f.facts)
            throw(ArgumentError("`|` of facts about different symbols is not supported: $f"))
        end
        x = only(vs)
        _restrict!(A, x, _factset(f, x))
        return true
    end
    f isa Membership && throw(ArgumentError("membership in a combination of facts is not supported: $f"))
    r = f::Relation
    d = expand(r.lhs - r.rhs)
    vs = _named(d)
    if isempty(vs)
        v = _decide(Relation(r.op, d, ZERO))
        return v isa Bool ? v : true
    end
    if r.op === :eq && length(vs) == 1              # several values: real roots as points
        x = only(vs)
        rs = _solve_eq(d ~ 0, x)
        isempty(_real_roots(rs)) && isempty(rs) && return false
        length(_real_roots(rs)) == length(rs) ||
            throw(ArgumentError("$r: values that are not all real are not supported"))
        _restrict!(A, x, _set_from(rs, (v, atpoint) -> atpoint))
        return true
    end
    foreach(x -> _markreal!(A, x), vs)
    if length(vs) == 1 && r.op !== :eq
        x = only(vs)
        S = try _solve_inequality(Relation(r.op, d, ZERO), x) catch; nothing end
        if S !== nothing
            _restrict!(A, x, S)
            return true
        end
    end
    lf = _linform(d)
    lf === nothing && throw(ArgumentError("$r is not supported as an assumption"))
    a, k = lf
    s = r.op in (:lt, :le) ? NM1 : N1
    sa = Dict{Sym,SymNumber}(t => c * s for (t, c) in a)
    if r.op === :eq
        push!(A.lin, _ACon(sa, k, false), _ACon(Dict(t => -c for (t, c) in sa), -k, false))
    else
        push!(A.lin, _ACon(sa, k * Sym(s), r.op in (:lt, :gt)))
    end
    true
end

# a combination of facts with no symbols: true, false, or nothing when undecided
function _constfact(f)
    if f isa Relation
        v = _decide(Relation(f.op, expand(f.lhs - f.rhs), ZERO))
        return v isa Bool ? v : nothing
    end
    f isa Membership && return nothing
    vs = map(_constfact, f.facts)
    f isa FactAnd ? (any(==(false), vs) ? false : all(==(true), vs) ? true : nothing) :
                    (any(==(true), vs) ? true : all(==(false), vs) ? false : nothing)
end

function _factsyms(f)
    f isa Relation && return _named(f.lhs - f.rhs)
    f isa Membership && return [f.x]
    unique!(reduce(vcat, (_factsyms(g) for g in f.facts); init = Sym[]))
end
function _factset(f, x::Sym)::SymSet
    f isa FactAnd && return reduce(∩, (_factset(g, x) for g in f.facts))
    f isa FactOr && return reduce(∪, (_factset(g, x) for g in f.facts))
    f isa Membership && throw(ArgumentError("membership in a combination of facts is not supported: $f"))
    r = f::Relation
    r.op === :eq && return _set_from(_real_roots(_solve_eq(r.lhs - r.rhs ~ 0, x)), (v, atpoint) -> atpoint)
    r.op === :ne && return setdiff(_REALS, _factset(!r, x))
    _solve_inequality(Relation(r.op, expand(r.lhs - r.rhs), ZERO), x)
end
const _REALS = SymSet([SymInterval(nothing, nothing, false, false)])

function _add_ne!(A::Assumptions, r::Relation)
    d = expand(r.lhs - r.rhs)
    vs = _named(d)
    if isempty(vs)
        v = _decide(Relation(:ne, d, ZERO))
        return v isa Bool ? v : true
    end
    if length(vs) == 1
        x = only(vs)
        if _ismarkedreal(A, x)
            _restrict!(A, x, _factset(Relation(:ne, d, ZERO), x))
        else
            for v in _solve_eq(d ~ 0, x)
                push!(get!(A.excluded, x, Sym[]), v)
            end
        end
        return true
    end
    push!(A.nonzero, d)
    true
end

# integers in their sets, empty sets, membership against non-membership, the linear
# relations feasible
function _check!(A::Assumptions)
    for (x, r) in A.notin
        get(A.dom, x, _CR) <= r && return false
    end
    for (x, S) in A.sets
        get(A.notin, x, Int8(0)) >= _RR && return false
        if get(A.dom, x, _CR) == _ZR
            S = _integer_hull(S)
            A.sets[x] = S
        end
        isempty(S) && return false
        for v in get(A.excluded, x, Sym[])
            v in S && (A.sets[x] = setdiff(S, _set_from([v], (u, atpoint) -> atpoint)))
        end
        isempty(A.sets[x]) && return false
    end
    for (x, vs) in A.excluded
        haskey(A.subst, x) && any(v -> _eq(v, A.subst[x]), vs) && return false
    end
    isempty(A.lin) || _fm_feasible(_with_bounds(A, copy(A.lin))) || return false
    true
end

# each interval narrowed to the integers in it ([ceil, floor]); empty ones dropped
function _integer_hull(S::SymSet)
    out = SymInterval[]
    for I in S.intervals
        lo = I.lo === nothing ? nothing : _intbound(I.lo, I.lclosed, true)
        hi = I.hi === nothing ? nothing : _intbound(I.hi, I.rclosed, false)
        lo !== nothing && hi !== nothing && lo > hi && continue
        push!(out, SymInterval(lo === nothing ? nothing : Sym(lo), hi === nothing ? nothing : Sym(hi),
                               lo !== nothing, hi !== nothing))
    end
    SymSet(out)
end
function _intbound(e::Sym, closed::Bool, lower::Bool)
    v = e.kind == NUMBER && isexact(e.num) ? Rational{BigInt}(_numerator(e.num), _denominator(e.num)) :
        _rv(e)
    n = lower ? ceil(BigInt, v) : floor(BigInt, v)
    n == v && !closed && (n += lower ? 1 : -1)
    n
end

## linear forms: Σ a[t] t + k, t atoms (named-containing terms without their numeric
## coefficient), k a constant; nothing when an atom's coefficient is not a number

function _linform(e::Sym)
    a = Dict{Sym,SymNumber}()
    k = ZERO
    terms = e.kind == ADD ? zip(e.coefs, e.args) : ((N1, e),)
    e.kind == ADD && (k = Sym(e.num))
    for (c, t) in terms
        if !_hassymbol(t)
            k += Sym(c) * t
            continue
        end
        if t.kind == MUL && !isone(t.num)
            c, t = c * t.num, t / Sym(t.num)
        end
        a[t] = get(a, t, N0) + c
    end
    filter!(p -> !iszero(p.second), a)
    a, k
end

## signs: masks of the possible signs of a real expression

const _NEG, _ZER, _POS = 0x01, 0x02, 0x04
const _ALLS = 0x07
const _SIGNBITS = ((_NEG, -1), (_ZER, 0), (_POS, 1))
_maskof(s::Int) = s < 0 ? _NEG : s == 0 ? _ZER : _POS

function _madd(m1::UInt8, m2::UInt8)
    out = 0x00
    for (b1, s1) in _SIGNBITS, (b2, s2) in _SIGNBITS
        (m1 & b1 == 0 || m2 & b2 == 0) && continue
        out |= s1 == 0 ? b2 : s2 == 0 ? b1 : s1 == s2 ? b1 : _ALLS
    end
    out
end
function _mmul(m1::UInt8, m2::UInt8)
    out = 0x00
    for (b1, s1) in _SIGNBITS, (b2, s2) in _SIGNBITS
        (m1 & b1 == 0 || m2 & b2 == 0) && continue
        out |= _maskof(s1 * s2)
    end
    out
end
_mneg(m::UInt8) = (m & _NEG == 0 ? 0x00 : _POS) | (m & _ZER) | (m & _POS == 0 ? 0x00 : _NEG)

# the sign of a constant: exact for numbers, numeric otherwise; nothing when it is not
# real or is too close to zero to tell (and not exactly zero)
function _csign(k::Sym)::Union{Nothing,Int}
    k.kind == NUMBER && return Int(sign(k.num))
    v = try _evalc(k) catch; return nothing end
    scale = max(1.0, abs(v))
    abs(imag(v)) > 1e-12 * scale && return nothing
    abs(real(v)) > 1e-12 * scale && return Int(sign(real(v)))
    iszero(expand(k)) ? 0 : nothing
end
_isreal_const(k::Sym) = (v = try _evalc(k) catch; return nothing end;
    abs(imag(v)) <= 1e-12 * max(1.0, abs(v)) ? true : false)

# the possible signs of e (a mask), or nothing when e is not known to be real
function _signs(e::Sym, A::Assumptions)::Union{Nothing,UInt8}
    isempty(A.subst) || (e = expand(subs(e, A.subst...)))
    haskey(A.cache, e) && return A.cache[e]
    m = _smask(e, A)
    # the structural answer is cached first: refining may ask about e again (a's sign
    # through the bounds of the atom a*b), and gets this sound, weaker answer
    A.cache[e] = m
    if m !== nothing && count_ones(m) > 1 && _hassymbol(e)
        m = _fm_refine(e, m, A)
    end
    if m !== nothing && count_ones(m) > 1 && e.kind == ADD     # a product of factors?
        f = try factor(e) catch; e end
        if f.kind == MUL
            mf = _smask(f, A)
            mf === nothing || (m &= mf)
        end
    end
    if m !== nothing && m & _ZER != 0 && count_ones(m) > 1 && _known_nonzero(e, A)
        m &= ~_ZER
    end
    A.cache[e] = m
    m
end

function _known_nonzero(e::Sym, A::Assumptions)
    any(d -> _eq(d, e) || _eq(d, -e), A.nonzero) && return true
    _isnamed(e) && return any(v -> iszero(v), get(A.excluded, e, Sym[]))
    false
end

function _setmask(S::SymSet)
    m = 0x00
    for I in S.intervals
        lo = I.lo === nothing ? -Inf : _rv(I.lo)
        hi = I.hi === nothing ? Inf : _rv(I.hi)
        (lo < 0) && (m |= _NEG)
        (hi > 0) && (m |= _POS)
        (lo < 0 < hi || (lo == 0 && I.lclosed) || (hi == 0 && I.rclosed)) && (m |= _ZER)
    end
    m
end

# structural: from the signs of the parts
function _smask(e::Sym, A::Assumptions)::Union{Nothing,UInt8}
    k = e.kind
    k == NUMBER && return _maskof(Int(sign(e.num)))
    if _isnamed(e)
        S = get(A.sets, e, nothing)
        S === nothing || return _setmask(S)
        _ismarkedreal(A, e) || return nothing
        return any(iszero, get(A.excluded, e, Sym[])) ? (_NEG | _POS) : _ALLS
    end
    if !_hassymbol(e)
        s = _csign(e)
        s === nothing && return _isreal_const(e) === true ? _ALLS : nothing
        return _maskof(s)
    end
    if k == ADD
        m = _maskof(Int(sign(e.num)))
        for (c, t) in zip(e.coefs, e.args)
            mt = _signs(t, A)
            mt === nothing && return nothing
            m = _madd(m, c < N0 ? _mneg(mt) : mt)
        end
        return m
    elseif k == MUL
        m = _maskof(Int(sign(e.num)))
        for (b, x) in zip(e.args, e.exps)
            mb = _powmask(b, x, A)
            mb === nothing && return nothing
            m = _mmul(m, mb)
        end
        return m
    elseif k == POW
        return _powmask(e.args[1], e.args[2], A)
    elseif k == FN
        return _fnmask(e, A)
    end
    nothing
end

function _powmask(b::Sym, x::Sym, A::Assumptions)
    if x.kind == NUMBER && isexact(x.num)
        mb = _signs(b, A)
        mb === nothing && return nothing
        q = x.num
        if isone(_denominator(q))
            n = _numerator(q)
            m = iseven(n) ? (_POS | (mb & _ZER)) : mb
            return n < 0 ? (m & ~_ZER == 0 ? nothing : m & ~_ZER) : m
        end
        mb & _NEG == 0 || return nothing             # a fractional power of a negative
        return q < N0 ? _POS : mb
    end
    if b.kind == CONSTANT && b.name === :E
        mx = _signs(x, A)
        return mx === nothing ? nothing : _POS
    end
    mb = _signs(b, A)
    (mb === nothing || mb != _POS) && return nothing
    _signs(x, A) === nothing ? nothing : _POS
end

function _fnmask(e::Sym, A::Assumptions)
    f = e.name
    length(e.args) == 1 || return nothing
    u = e.args[1]
    if f === :abs
        return _isnonzero(u, A) === true ? _POS : (_ZER | _POS)
    end
    mu = _signs(u, A)
    mu === nothing && return nothing
    f in (:sign, :atan, :asinh, :sinh, :tanh, :erf, :cbrt) && return mu
    f === :cosh && return _POS
    f in (:sin, :cos, :tan, :cot, :sec, :csc) && return _ALLS
    f === :floor && return (mu & _NEG) | (mu & (_ZER | _POS) == 0 ? 0x00 : _ZER) | (mu & _POS)
    f === :ceil && return (mu & _POS) | (mu & (_ZER | _NEG) == 0 ? 0x00 : _ZER) | (mu & _NEG)
    if f === :log
        mu == _POS || return nothing
        return _signs(u - 1, A)
    end
    nothing
end

## three-valued queries

# is e real? true, false, or nothing (unknown)
function _isreal(e::Sym, A::Assumptions)::Union{Nothing,Bool}
    _signs(e, A) === nothing || return true
    isempty(A.subst) || (e = expand(subs(e, A.subst...)))
    _isnamed(e) && return get(A.notin, e, Int8(0)) >= _RR ? false : nothing
    _hassymbol(e) || return _isreal_const(e)
    if e.kind == POW
        # b^q, q not an integer, |q| < 1: real exactly when b ≥ 0 (principal branch)
        b, q = e.args
        if q.kind == NUMBER && isexact(q.num) && !isone(_denominator(q.num)) && abs(q.num) < N1
            mb = _signs(b, A)
            mb !== nothing && mb & _NEG != 0 && mb & (_ZER | _POS) == 0 && return false
            _isreal(b, A) === false && return false
        end
    elseif e.kind == MUL || e.kind == ADD
        # real parts and one that is not real
        parts = e.kind == MUL ? [b^x for (b, x) in zip(e.args, e.exps)] : e.args
        rs = [_isreal(p, A) for p in parts]
        if count(==(false), rs) == 1 && all(!=(nothing), rs)
            e.kind == ADD && return false
            return all(i -> rs[i] == false || _isnonzero(parts[i], A) === true, eachindex(parts)) ? false : nothing
        end
    end
    nothing
end

# is e an integer?
function _isinteger(e::Sym, A::Assumptions)::Union{Nothing,Bool}
    isempty(A.subst) || (e = expand(subs(e, A.subst...)))
    k = e.kind
    k == NUMBER && return isexact(e.num) && isone(_denominator(e.num))
    _hassymbol(e) || return _const_rational(e) === false ? false : nothing
    if _isnamed(e)
        get(A.dom, e, _CR) == _ZR && return true
        get(A.notin, e, Int8(0)) >= _ZR && return false
        return nothing
    end
    if k == ADD
        isone(_denominator(e.num)) || return nothing
        all(c -> isexact(c) && isone(_denominator(c)), e.coefs) || return nothing
        return all(t -> _isinteger(t, A) === true, e.args) ? true : nothing
    elseif k == MUL
        isexact(e.num) && isone(_denominator(e.num)) || return nothing
        for (b, x) in zip(e.args, e.exps)
            (x.kind == NUMBER && isexact(x.num) && isone(_denominator(x.num)) && x.num > N0) || return nothing
            _isinteger(b, A) === true || return nothing
        end
        return true
    elseif k == POW
        x = e.args[2]
        (x.kind == NUMBER && isexact(x.num) && isone(_denominator(x.num)) && x.num > N0) || return nothing
        return _isinteger(e.args[1], A) === true ? true : nothing
    elseif k == FN && e.name in (:floor, :ceil) && length(e.args) == 1
        return _isreal(e.args[1], A) === true ? true : nothing
    end
    nothing
end

# is e nonzero?
function _isnonzero(e::Sym, A::Assumptions)::Union{Nothing,Bool}
    m = _signs(e, A)
    m === nothing || return m & _ZER == 0 ? true : m == _ZER ? false : nothing
    isempty(A.subst) || (e = expand(subs(e, A.subst...)))
    _known_nonzero(e, A) && return true
    if !_hassymbol(e)
        v = try _evalc(e) catch; return nothing end
        return abs(v) > 1e-12 ? true : iszero(expand(e)) ? false : nothing
    end
    if e.kind == MUL
        return all(b -> _isnonzero(b, A) === true, e.args) ? true : nothing
    elseif e.kind == POW
        b, x = e.args
        b.kind == CONSTANT && b.name === :E && return true
        return _isnonzero(b, A) === true ? true : nothing
    end
    nothing
end

## Fourier–Motzkin: is a set of constraints feasible? false only when certainly not

const _FM_MAX = 400

# the constraints with bounds on their atoms: from the atoms' signs, and from the sets
# of the symbols they are monotone functions of
function _with_bounds(A::Assumptions, cons::Vector{_ACon})
    atoms = Set{Sym}()
    for c in cons
        union!(atoms, keys(c.a))
    end
    for t in atoms
        m = _smask(t, A)
        if m !== nothing
            m & _NEG == 0 && push!(cons, _ACon(Dict(t => N1), ZERO, m & _ZER == 0))
            m & _POS == 0 && push!(cons, _ACon(Dict(t => NM1), ZERO, m & _ZER == 0))
        end
        _monotone_bounds!(cons, t, A)
    end
    cons
end

# t = f(x), f monotone, x a symbol in a set with hull (lo, hi): f(lo) ≤ t ≤ f(hi)
function _monotone_bounds!(cons, t::Sym, A::Assumptions)
    x, inc = _monotone_of(t, A)
    x === nothing && return
    S = A.sets[x]
    isempty(S.intervals) && return
    I1, I2 = first(S.intervals), last(S.intervals)
    for (e, closed, lower) in ((I1.lo, I1.lclosed, true), (I2.hi, I2.rclosed, false))
        e === nothing && continue
        v = try subs(t, x => e) catch; continue end
        _hassymbol(v) && continue
        _csign(v) === nothing && continue        # not a real value (log at 0, …)
        # t ≥ f(lo) for increasing f, t ≤ f(lo) for decreasing
        above = lower == inc
        push!(cons, above ? _ACon(Dict(t => N1), -v, !closed) : _ACon(Dict(t => NM1), v, !closed))
    end
end

# (x, increasing) when t is a monotone function of the symbol x on x's set
function _monotone_of(t::Sym, A::Assumptions)
    none = (nothing, false)
    if t.kind == POW
        b, q = t.args
        if b.kind == CONSTANT && b.name === :E
            return _isnamed(q) && haskey(A.sets, q) ? (q, true) : none
        end
        (_isnamed(b) && haskey(A.sets, b) && q.kind == NUMBER && isexact(q.num)) || return none
        m = _setmask(A.sets[b])
        n = q.num
        if isone(_denominator(n)) && isodd(_numerator(n)) && n > N0
            return (b, true)
        end
        m & _NEG == 0 || return none
        n < N0 && m & _ZER != 0 && return none
        return (b, n > N0)
    elseif t.kind == FN && length(t.args) == 1
        u = t.args[1]
        (_isnamed(u) && haskey(A.sets, u)) || return none
        t.name in (:atan, :asinh, :sinh, :tanh, :erf, :cbrt) && return (u, true)
        t.name === :log && _setmask(A.sets[u]) == _POS && return (u, true)
    elseif _isnamed(t) && haskey(A.sets, t)
        return (t, true)
    end
    none
end

# the possible signs of e narrowed by the linear relations
function _fm_refine(e::Sym, m::UInt8, A::Assumptions)
    lf = _linform(e)
    a, k = lf
    out = 0x00
    for (b, s) in _SIGNBITS
        m & b == 0 && continue
        cons = copy(A.lin)
        if s == 0
            push!(cons, _ACon(copy(a), k, false), _ACon(Dict(t => -c for (t, c) in a), -k, false))
        else
            sa = Dict(t => c * SymNumber(s) for (t, c) in a)
            push!(cons, _ACon(sa, k * s, true))
        end
        _fm_feasible(_with_bounds(A, cons)) && (out |= b)
    end
    out == 0x00 ? m : out             # all infeasible: the facts are contradictory; keep m
end

function _fm_feasible(cons::Vector{_ACon})
    cs = _fm_prune(cons)
    cs === nothing && return false
    while true
        atoms = Dict{Sym,Tuple{Int,Int}}()
        for c in cs, (t, v) in c.a
            p, n = get(atoms, t, (0, 0))
            atoms[t] = v > N0 ? (p + 1, n) : (p, n + 1)
        end
        isempty(atoms) && return true
        # the atom whose elimination makes the fewest new constraints
        t = argmin(t -> atoms[t][1] * atoms[t][2] - atoms[t][1] - atoms[t][2], collect(keys(atoms)))
        pos = [c for c in cs if get(c.a, t, N0) > N0]
        neg = [c for c in cs if get(c.a, t, N0) < N0]
        rest = [c for c in cs if !haskey(c.a, t)]
        length(rest) + length(pos) * length(neg) > _FM_MAX && return true    # give up
        for p in pos, n in neg
            sp, sn = inv(p.a[t]), inv(-n.a[t])
            a = Dict{Sym,SymNumber}()
            for (u, v) in p.a
                _eq(u, t) || (a[u] = get(a, u, N0) + v * sp)
            end
            for (u, v) in n.a
                _eq(u, t) || (a[u] = get(a, u, N0) + v * sn)
            end
            filter!(q -> !_eq(q.first, t) && !iszero(q.second), a)
            push!(rest, _ACon(a, p.k * Sym(sp) + n.k * Sym(sn), p.strict || n.strict))
        end
        cs = _fm_prune(rest)
        cs === nothing && return false
    end
end

# constants decided (nothing if one is false), duplicates in one direction reduced to
# the tighter
function _fm_prune(cons::Vector{_ACon})
    best = Dict{Vector{Pair{Sym,SymNumber}},_ACon}()
    out = _ACon[]
    for c in cons
        if isempty(c.a)
            s = _csign(c.k)
            s === nothing && continue
            (c.strict ? s > 0 : s >= 0) || return nothing
            continue
        end
        ks = sort!(collect(c.a); by = p -> p.first, lt = _lt)
        sc = inv(abs(ks[1].second))
        key = [t => v * sc for (t, v) in ks]
        cn = _ACon(Dict(key), c.k * Sym(sc), c.strict)
        d = get(best, key, nothing)
        if d === nothing
            best[key] = cn
            continue
        end
        s = _csign(expand(cn.k - d.k))
        if s === nothing
            push!(out, cn)
        elseif s < 0 || (s == 0 && cn.strict && !d.strict)
            best[key] = cn
        end
    end
    append!(out, values(best))
end

## the macro: facts written with Julia's comparison and Boolean operators

"""
    @assume(facts...)

[`Assumptions`](@ref) from facts written with Julia's own operators: `<`, `<=`, `>`,
`>=`, `==`, `!=` (and `≤`, `≥`, `≠`) make relations, chained comparisons such as
`1 < y <= 2` are split into their parts, and `&&`, `||` and `!` combine facts. Anything
else is used as it is, such as `n in Integers` or a fact made earlier.

```julia
A = @assume(x > 0 && 1 < y <= 2, n in Integers, a != 0)
solve(x^2 ~ 4, x; assume = @assume(x > 0))
assuming(@assume(x >= 1)) do … end
```
"""
macro assume(facts...)
    esc(Expr(:call, GlobalRef(SymJ, :Assumptions), map(_assume_fact, facts)...))
end

const _ASSUME_OPS = Dict(:< => :lt, :<= => :le, :≤ => :le, :> => :gt, :>= => :ge, :≥ => :ge,
                         :(==) => :eq, :!= => :ne, :≠ => :ne)
_assume_rel(op::Symbol, a, b) = Expr(:call, GlobalRef(SymJ, :_mkrel), QuoteNode(op), a, b)
_mkrel(op::Symbol, a, b) = Relation(op, Sym(a), Sym(b))

function _assume_fact(ex)
    ex isa Expr || return ex
    h, as = ex.head, ex.args
    h === :&& && return Expr(:call, GlobalRef(Base, :&), map(_assume_fact, as)...)
    h === :|| && return Expr(:call, GlobalRef(Base, :|), map(_assume_fact, as)...)
    if h === :comparison                         # a < b <= c: (a < b) & (b <= c)
        parts = [_assume_rel(_ASSUME_OPS[as[i]], as[i-1], as[i+1]) for i in 2:2:length(as)]
        return length(parts) == 1 ? parts[1] : Expr(:call, GlobalRef(Base, :&), parts...)
    end
    if h === :call
        f = as[1]
        f === :! && length(as) == 2 && return Expr(:call, GlobalRef(Base, :!), _assume_fact(as[2]))
        haskey(_ASSUME_OPS, f) && length(as) == 3 && return _assume_rel(_ASSUME_OPS[f], as[2], as[3])
    end
    ex
end

## scope: a block of assumptions, and the `assume` keyword

const _ASSUMING = ScopedValue{Union{Nothing,Assumptions}}(nothing)

"""
    assuming(f, facts...)
    assuming(facts...) do … end

Run `f` with the facts assumed by the operations that take assumptions, as if each
were given `assume = …`. Blocks nest: an inner block adds its facts to the outer
ones, and so does an `assume` keyword inside a block. The facts are normalized once,
when the block starts ([`Assumptions`](@ref)).

```julia
assuming(x ≫ 0) do
    solve(x^2 ~ 4, x)       # [2]
end
```
"""
function assuming(f, facts...)
    cur = _ASSUMING[]
    A = cur === nothing ? Assumptions(facts...) : Assumptions(cur, facts...)
    with(f, _ASSUMING => A)
end

"""
    @assuming facts... begin … end
    @assuming(facts..., begin … end)

Run the block with the facts assumed, the facts written as in [`@assume`](@ref);
`assuming(@assume(facts...)) do … end` in one step. The block's value is returned; as
with `do`, the block is a function body, so variables assigned in it are local to it.

```julia
@assuming x > 0 begin
    solve(x^2 ~ 4, x)       # [2]
end
```
"""
macro assuming(args...)
    (length(args) >= 1 && args[end] isa Expr && args[end].head === :block) ||
        throw(ArgumentError("@assuming facts... begin … end: the last argument must be a begin … end block"))
    facts, body = args[1:end-1], args[end]
    esc(Expr(:call, GlobalRef(SymJ, :assuming), Expr(:->, Expr(:tuple), body),
             Expr(:call, GlobalRef(SymJ, :Assumptions), map(_assume_fact, facts)...)))
end

# the assumptions for a call: the block's, with the keyword's added; nothing if none
function _assumptions(assume)
    cur = _ASSUMING[]
    (assume === nothing || (assume isa Union{Tuple,AbstractVector} && isempty(assume))) && return cur
    cur === nothing && return assume isa Assumptions ? assume : Assumptions(assume)
    Assumptions(cur, assume)
end
