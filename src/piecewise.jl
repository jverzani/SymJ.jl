# Piecewise functions, heaviside and dirac.
#
# `piecewise(c₁ => v₁, c₂ => v₂, …)` is an expression: the value of the first piece
# whose condition holds. It is a function node, `piecewise(c₁, v₁, c₂, v₂, …)`, its
# conditions themselves nodes: a relation, Lt(a, b) and its kin (Le, Gt, Ge, Eq, Ne),
# a conjunction And(c, d, …), or 1 (true) and 0 (false), to which a relation between
# constants evaluates. So a substitution picks the piece, and arithmetic, `diff`, `N`
# and `lambdify` take a piecewise as they take any function. Its derivative is taken
# piece by piece, a jump J at a breakpoint p adding J dirac(x - p), as for heaviside.
#
# As built, the pieces whose condition is false are dropped, a true condition ends the
# list, and a piece with the value of the true one after it is dropped. With no piece
# left the value is NaN, as for a value outside every condition.
#
# `heaviside(x)` (1/2 at 0, or `heaviside(x, h0)`) and `dirac(x)` (`dirac(x, k)`, its
# k-th derivative) are functions: diff(heaviside(x)) = dirac(x), diff(sign(x)) =
# 2 dirac(x). `piecewise(ex)` writes heaviside, sign, abs, max and min in ex as
# piecewise functions (for real arguments); `piecewise(cs)` turns cases into one.
#
# Unlike `Cases`, which is an answer for values of the parameters the assumptions do not
# settle, a piecewise function is an expression in its variable.

const _COND_OPS = Dict(:lt => :Lt, :le => :Le, :gt => :Gt, :ge => :Ge, :eq => :Eq, :ne => :Ne)
const _COND_RELS = Dict(v => k for (k, v) in _COND_OPS)
_iscond(c::Sym) = c.kind == FN && (haskey(_COND_RELS, c.name) || c.name === :And)

# a condition as a node: a relation, a Bool, a tuple or vector (all of them hold) or a
# node already
_cond(r::Relation) = apply_function(_COND_OPS[r.op], Sym[r.lhs, r.rhs])
_cond(b::Bool) = b ? ONE : ZERO
_cond(c::Sym) = (_iscond(c) || c == ONE || c == ZERO) ? c :
    throw(ArgumentError("piecewise: $c is not a condition"))
_cond(cs::Union{Tuple,AbstractVector}) = apply_function(:And, Sym[_cond(c) for c in cs])
_cond(c::FactAnd) = apply_function(:And, Sym[_cond(f) for f in c.facts])
_cond(c) = throw(ArgumentError("piecewise: $c is not a condition"))

# a condition as the conditions of pieces with one value: an Or as one piece for each
# of its parts, an And of Ors distributed over them (SymPy #22852)
_disjuncts(c::FactOr) = reduce(vcat, Any[_disjuncts(f) for f in c.facts])
function _disjuncts(c::Union{FactAnd,Tuple,AbstractVector})
    fs = c isa FactAnd ? c.facts : collect(c)
    any(f -> f isa FactOr, fs) || return Any[c]
    out = Any[Any[]]
    for f in fs
        ds = f isa FactOr ? _disjuncts(f) : Any[f]
        out = Any[vcat(o, Any[d]) for o in out for d in ds]
    end
    Any[Tuple(o) for o in out]
end
_disjuncts(c) = Any[c]

# a condition evaluated with floats (Lt(0.5, 1.0) is 1.0) as 1 or 0
_truth(c::Sym) = c.kind == NUMBER ? (iszero(c.num) ? ZERO : ONE) : c

# the node as a relation (for deciding it)
_relation_of(c::Sym) = Relation(_COND_RELS[c.name], c.args[1], c.args[2])

const _FLIP_OP = Dict(:lt => :gt, :le => :ge, :gt => :lt, :ge => :le, :eq => :eq, :ne => :ne)

for (op, name) in _COND_OPS
    f = _REL_INFO[op][2]
    register_function!(name; julia = (a, b) -> Float64(f(real(a), real(b))),
        eval = args -> begin
            length(args) == 2 || return nothing
            d = expand(args[1] - args[2])
            _hassymbol(d) || return (s = _csign(d); s === nothing ? nothing : f(s, 0) ? ONE : ZERO)
            # linear in one symbol: solved for it, 1 - x < 0 as x > 1
            vs = _named(d)
            length(vs) == 1 || return nothing
            v = only(vs)
            c = expand(subs(d, v => ONE) - subs(d, v => ZERO))
            (c.kind == NUMBER && !iszero(c.num) && iszero(expand(d - c * v - subs(d, v => ZERO)))) || return nothing
            rhs = expand(-subs(d, v => ZERO) / c)
            (_eq(args[1], v) && _eq(args[2], rhs)) && return nothing
            _mkfn_node(_COND_OPS[c.num > N0 ? op : _FLIP_OP[op]], Sym[v, rhs])
        end)
end

register_function!(:And; julia = (vs...) -> Float64(all(!iszero, vs)),
    eval = function (args)
        args = map(_truth, args)
        any(==(ZERO), args) && return ZERO
        rest = Sym[]
        for a in args
            a == ONE && continue
            a.kind == FN && a.name === :And ? append!(rest, a.args) : push!(rest, a)
        end
        isempty(rest) && return ONE
        length(rest) == 1 && return rest[1]
        length(rest) == length(args) && all(i -> _eq(rest[i], args[i]), eachindex(rest)) && return nothing
        _mkfn_node(:And, rest)
    end)

"""
    SymJ.piecewise(c₁ => v₁, c₂ => v₂, …)
    SymJ.piecewise(ex)
    SymJ.piecewise(cs::Cases)

A piecewise function: `v₁` where `c₁` holds, else `v₂` where `c₂` holds, and so on. A
condition is a relation (`x ≪ 0`), a tuple of them that all hold (`(0 ≦ x, x ≪ 1)`),
their combinations with `&` and `|` (`(x ≪ 0) | (x ≫ 1)`, made one piece for each
alternative), or `true`. Not exported.

With one expression, writes its `heaviside`, `sign`, `abs`, `max` and `min` as piecewise
functions (taking their arguments to be real); with cases, makes them one expression.

The derivative is taken piece by piece; where the function jumps by `J` at `p`, it has
the term `J*dirac(x - p)`, as `heaviside` does.

```julia
f = SymJ.piecewise(x ≪ 0 => -x, true => x^2)
f(x => -2), diff(f, x)          # 2, piecewise(x ≪ 0 => -1, true => 2x)
diff(SymJ.piecewise(x ≪ 0 => 0, true => 1), x)   # dirac(x)
```
"""
function piecewise(pairs::Pair...)
    isempty(pairs) && throw(ArgumentError("piecewise: no pieces"))
    args = Sym[]
    for (c, v) in pairs, d in _disjuncts(c)
        push!(args, _cond(d), Sym(v))
    end
    apply_function(:piecewise, args)
end

_pieces(e::Sym) = [(e.args[i], e.args[i+1]) for i in 1:2:length(e.args)]
_pwnode(ps) = _mkfn_node(:piecewise, Sym[x for p in ps for x in p])

register_function!(:piecewise;
    eval = function (args)
        (iseven(length(args)) && !isempty(args)) || return nothing
        ps = Tuple{Sym,Sym}[]
        for i in 1:2:length(args)
            c, v = _truth(args[i]), args[i+1]
            (c == ZERO || _infeasible_and(c)) && continue
            push!(ps, (c, v))
            c == ONE && break
        end
        isempty(ps) && return Sym(NaN)
        # a piece with the same value as the true one after it is not needed
        while length(ps) >= 2 && ps[end][1] == ONE && _eq(ps[end-1][2], ps[end][2])
            deleteat!(ps, length(ps) - 1)
        end
        # nor is a piece at a point, x = p, whose value there is that of the piece that
        # takes over at p (SymPy #27627)
        for i in length(ps)-1:-1:1
            _pw_point_redundant(ps, i) && deleteat!(ps, i)
        end
        ps[1][1] == ONE && return ps[1][2]
        length(ps) == length(args) ÷ 2 && all(i -> _eq(ps[i][1], args[2i-1]), eachindex(ps)) && return nothing
        _pwnode(ps)
    end)

# an And of relations that cannot all hold (x ≤ 0 and x > π; SymPy #25930)
function _infeasible_and(c::Sym)
    (c.kind == FN && c.name === :And && all(a -> a.kind == FN && haskey(_COND_RELS, a.name), c.args)) || return false
    B = try _build(_flatten!(Any[], (Any[], Any[_relation_of(a) for a in c.args]))) catch err
        err isa ArgumentError || rethrow()
        return false
    end
    B === nothing
end

# is piece i, with condition x = p, needless: the first later piece that holds at
# x = p has the same value there
function _pw_point_redundant(ps, i)
    c = ps[i][1]
    (c.kind == FN && c.name === :Eq) || return false
    l, r = c.args
    x, p = _isnamed(l) && !has(r, l) ? (l, r) : _isnamed(r) && !has(l, r) ? (r, l) : (nothing, nothing)
    x === nothing && return false
    at(e) = subs(e, x => p)
    try
        for j in i+1:length(ps)
            t = _truth(at(ps[j][1]))
            t == ZERO && continue
            t == ONE || return false              # undecided at p
            return _eq(at(ps[i][2]), at(ps[j][2]))
        end
    catch err
        err isa DivisionByZero || rethrow()       # sin(x)/x at 0: the piece is needed
    end
    false
end

# numeric value: the first piece whose condition holds, the others not evaluated
function _evalc_piecewise(e::Sym)
    for (c, v) in _pieces(e)
        iszero(real(_evalc(c))) || return _evalc(v)
    end
    complex(NaN)
end

## heaviside and dirac

_heaviside(x::Real, h0 = 0.5) = x < 0 ? zero(float(x)) : x > 0 ? one(float(x)) : float(h0)
_dirac(x::Real, k = 0) = iszero(x) ? Inf : zero(float(x))

"""
    SymJ.heaviside(x, h0 = 1/2)

The unit step: 0 for `x < 0`, 1 for `x > 0`, and `h0` at 0. Its derivative is
[`SymJ.dirac`](@ref). Not exported.
"""
heaviside(x) = apply_function(:heaviside, Sym[Sym(x)])
heaviside(x, h0) = apply_function(:heaviside, Sym[Sym(x), Sym(h0)])

"""
    SymJ.dirac(x, k = 0)

The Dirac delta (with `k`, its `k`-th derivative): 0 away from 0, the derivative of
[`SymJ.heaviside`](@ref), with `∫ f(x) dirac(x - a) dx = f(a) heaviside(x - a)`. Not
exported.
"""
dirac(x) = apply_function(:dirac, Sym[Sym(x)])
dirac(x, k::Integer) = apply_function(:dirac, Sym[Sym(x), Sym(k)])

register_function!(:heaviside; julia = _heaviside,
    eval = function (args)
        (1 <= length(args) <= 2) || return nothing
        h0 = length(args) == 2 ? args[2] : Sym(1 // 2)
        a = args[1]
        _hassymbol(a) && return nothing
        s = _csign(a)
        s === nothing ? nothing : s < 0 ? ZERO : s > 0 ? ONE : h0
    end,
    diff = (args, i) -> i == 1 ? dirac(args[1]) : ZERO)

register_function!(:dirac; julia = _dirac,
    eval = function (args)
        (1 <= length(args) <= 2) || return nothing
        length(args) == 2 && args[2] == ZERO && return dirac(args[1])
        a = args[1]
        _hassymbol(a) && return nothing
        s = _csign(a)
        s === nothing || s == 0 ? nothing : ZERO
    end,
    diff = (args, i) -> i == 1 ? dirac(args[1], (length(args) == 2 ? _toint(args[2].num) : 0) + 1) : ZERO)

# sign jumps by 2 at 0
let r = FUNCTION_RULES[:sign]
    register_function!(:sign; julia = r.julia, eval = r.eval, diff = (args, i) -> 2dirac(args[1]))
end

## rewriting as piecewise functions, and cases as one

function piecewise(ex::Sym)
    isempty(ex.args) && return ex
    e = _rebuild_node(ex, Sym[piecewise(a) for a in ex.args])
    e.kind == FN || return e
    a = isempty(e.args) ? ZERO : e.args[1]
    if e.name === :heaviside
        h0 = length(e.args) == 2 ? e.args[2] : Sym(1 // 2)
        return piecewise(a ≪ 0 => 0, (a ~ 0) => h0, true => 1)
    elseif e.name === :sign && length(e.args) == 1
        return piecewise(a ≪ 0 => -1, (a ~ 0) => 0, true => 1)
    elseif e.name === :abs && length(e.args) == 1
        return piecewise(a ≪ 0 => -a, true => a)
    elseif e.name in (:max, :min) && length(e.args) == 2
        b = e.args[2]
        return piecewise((e.name === :max ? a ≧ b : a ≦ b) => a, true => b)
    end
    e
end
piecewise(x::Number) = Sym(x)

function piecewise(cs::Cases)
    ps = Pair[]
    for (conds, v) in cs.cases
        v isa Exception && continue                       # undefined there: no piece
        v isa Sym || v isa Number || throw(ArgumentError("piecewise: the case answer $v is not an expression"))
        push!(ps, Any[_case_cond(c) for c in conds] => v)
    end
    isempty(ps) && return Sym(NaN)
    piecewise(ps...)
end
_case_cond(r::Relation) = r
_case_cond(c::_NotFact) = c.f isa Relation ? Relation(_NEGATE_OP[c.f.op], c.f.lhs, c.f.rhs) :
    throw(ArgumentError("piecewise: cannot write the condition $c"))
_case_cond(c) = throw(ArgumentError("piecewise: cannot write the condition $c"))
const _NEGATE_OP = Dict(:lt => :ge, :le => :gt, :gt => :le, :ge => :lt, :eq => :ne, :ne => :eq)

## refine: conditions the assumptions decide

function _refine_cond(e::Sym, A::Assumptions)
    haskey(_COND_RELS, e.name) || return nothing
    r = _relation_of(e)
    v = _holds(Relation(r.op, expand(r.lhs - r.rhs), ZERO), A)
    v === nothing ? nothing : _cond(v)
end

# an And of relations that contradict each other (with the assumptions): false
# (y ≥ 0 and y^3 < -8; SymPy #28323)
function _refine_and(e::Sym, A::Assumptions)
    all(a -> a.kind == FN && haskey(_COND_RELS, a.name), e.args) || return nothing
    rels = Any[_relation_of(a) for a in e.args]
    B = try _build(_flatten!(Any[], (A.facts, rels))) catch err
        err isa ArgumentError || rethrow()
        return nothing
    end
    B === nothing ? ZERO : nothing
end

# the sign of a heaviside's, sign's or abs's argument, where the assumptions decide it
function _refine_step(e::Sym, A::Assumptions)
    e.name in (:heaviside, :sign, :abs, :dirac) || return nothing
    u = e.args[1]
    # |u^r| = |u|^r for u and r real (|√x| = √|x|, |(-n/(n+1))^n|; SymPy #25779, #26513)
    if e.name === :abs && u.kind == POW && _isreal(u.args[1], A) === true && _isreal(u.args[2], A) === true &&
       !(u.args[2].kind == NUMBER && _isint(u.args[2].num))
        return abs(u.args[1])^u.args[2]
    end
    m = _signs(u, A)
    m === nothing && return nothing
    # dirac(u) = 0 where u ≠ 0 (SymPy #26691)
    e.name === :dirac && return m & _ZER == 0 ? ZERO : nothing
    if m == _POS
        return e.name === :abs ? u : ONE
    elseif m == _NEG
        return e.name === :heaviside ? ZERO : e.name === :sign ? -ONE : -u
    end
    nothing
end

## limits: the piece taken on one side of x0

# a piecewise node in e that involves x only in its values, not its conditions; or nothing
function _param_piecewise(e::Sym, x::Sym)
    if e.kind == FN && e.name === :piecewise && has(e, x) && all(((c, _),) -> !has(c, x), _pieces(e))
        return e
    end
    for a in (e.kind == MUL ? [e.args; e.exps] : e.args)
        p = _param_piecewise(a, x)
        p === nothing || return p
    end
    nothing
end

const _PIECEFNS = (:piecewise, :heaviside, :dirac)
_haspieces(e::Sym, x::Sym) = (e.kind == FN && e.name in _PIECEFNS && has(e, x)) ||
    any(a -> _haspieces(a, x), e.args) || (e.kind == MUL && any(a -> _haspieces(a, x), e.exps))

# f with each piecewise, heaviside and dirac in x replaced by its value as x → x0 from
# `side` (for x0 = ±∞, from below or above); an error when a condition is not decided
function _pieces_near(e::Sym, x::Sym, x0, side::Symbol, assume)
    _haspieces(e, x) || return e
    e = _rebuild_node(e, Sym[_pieces_near(a, x, x0, side, assume) for a in e.args])
    e.kind == FN && e.name in _PIECEFNS || return e
    if e.name === :piecewise
        for (c, v) in _pieces(e)
            _cond_near(c, x, x0, side, assume) && return v
        end
        return Sym(NaN)
    end
    s = _side_sign(e.args[1], x, x0, side, assume)
    if e.name === :heaviside
        return s < 0 ? ZERO : s > 0 ? ONE : length(e.args) == 2 ? e.args[2] : Sym(1 // 2)
    end
    s == 0 ? e : ZERO                                    # dirac
end

function _cond_near(c::Sym, x::Sym, x0, side, assume)
    c == ONE && return true
    c == ZERO && return false
    c.kind == FN && c.name === :And && return all(a -> _cond_near(a, x, x0, side, assume), c.args)
    r = _relation_of(c)
    _REL_INFO[r.op][2](_side_sign(r.lhs - r.rhs, x, x0, side, assume), 0)
end

# the sign (-1, 0, 1) of g for x near x0 on `side`: of its limit, or where that is 0, of
# g / t^j for t > 0 going to 0 (x - x0 or x0 - x, 1/|x| at ±∞)
function _side_sign(g::Sym, x::Sym, x0, side, assume)
    g = expand(g)
    iszero(g) && return 0
    inf = x0 isa Real && isinf(x0)
    t = inf ? (x0 > 0 ? 1 / x : -1 / x) : side === :+ ? x - x0 : x0 - x
    for j in 0:6
        L = limit(g / t^j, x => x0; dir = side, assume)
        L isa Sym || break
        s = _csign(L)
        # a limit with parameters: its sign by the assumptions (x heaviside(x) at L > 0;
        # SymPy #24042)
        if s === nothing
            A = _assumptions(assume)
            A === nothing && (A = _ASSUMING[])
            A === nothing && break
            s = _ispos(L, A) ? 1 : _ispos(-L, A) ? -1 : nothing
            s === nothing && break
        end
        s != 0 && return s
    end
    throw(ArgumentError("limit: cannot tell the sign of $g as $x → $x0 ($side)"))
end

## antiderivatives (with an integrator, `integ(f)`, from the RischIntegrate extension)

# ∫ e dx for e with piecewise, heaviside or dirac in x:
#   - a term g(x) dirac(m x + c, k), m a number: (-1)^k g⁽ᵏ⁾(a)/(m^k |m|) heaviside(x - a)
#     at a = -c/m;
#   - the rest, between consecutive breakpoints of its conditions: the integrand on each
#     interval (read at a point inside), its antiderivative, and constants that make the
#     whole continuous. The breakpoints are real numbers, or one expression in parameters
#     (heaviside(x - a)); others, which could not be ordered, are an error.
function _integrate_pieces(e::Sym, x::Sym, integ)
    # dirac inside a product with a sum: its terms apart (#27108)
    e.kind != ADD && _hasfn(e, :dirac) && (e = expand(e))
    terms = e.kind == ADD ? [Sym(c) * t for (c, t) in zip(e.coefs, e.args)] : [e]
    e.kind == ADD && !iszero(e.num) && push!(terms, Sym(e.num))
    deltas, rest = Sym[], Sym[]
    for t in terms
        r = _integrate_dirac(t, x)
        r === nothing ? push!(rest, t) : push!(deltas, r)
    end
    # the rest's abs, sign, max, min as pieces too; conditions of one sign everywhere
    # settled; terms without pieces integrated apart (|sin x| + x: the pieces' part periodic)
    g = _settle_conditions(piecewise(_add_many(rest)), x, -Inf, Inf)
    out = _add_many(deltas)
    if g.kind == ADD && _haspieces(g, x)
        plain = _add_many(Sym(c) * t for (c, t) in zip(g.coefs, g.args) if !_haspieces(t, x)) + Sym(g.num)
        g = _add_many(Sym(c) * t for (c, t) in zip(g.coefs, g.args) if _haspieces(t, x))
        iszero(plain) || (out = out + integ(plain))
    end
    _haspieces(g, x) || return iszero(g) ? out : out + integ(g)
    bps = _breakpoints(g, x)
    if bps === nothing
        # a condition quadratic in x, its roots in parameters: by the sign of the
        # discriminant, between the roots as heaviside steps (SymPy #27358); a sum, term
        # by term (each with one breakpoint in parameters, say)
        g2 = _quadratic_conditions(g, x)
        part(h) = _haspieces(h, x) ? _integrate_pieces(h, x, integ) : integ(h)
        # x only in one linear argument k x + c, c in parameters (|cos(x + y)|, #27234):
        # ∫ h(k x + c) dx = H(k x + c)/k, h free of the parameters
        lin = _linear_shift(g, x)
        if lin !== nothing
            k, c, h = lin
            return out + subs(part(h), x => k * x + c) / k
        end
        if g2 !== nothing
            return out + piecewise((c => part(expand(v)) for (c, v) in _pieces(g2))...)
        elseif g.kind == ADD
            return out + _add_many(part(Sym(c) * t) for (c, t) in zip(g.coefs, g.args)) +
                   (iszero(g.num) ? ZERO : integ(Sym(g.num)))
        end
        throw(ArgumentError("integrate: the points where the pieces of $g change could not be found and ordered"))
    end
    # a point inside each interval, the integrand there, its antiderivative
    pts = Sym[]
    if length(bps) == 1 && _hassymbol(bps[1])
        append!(pts, (bps[1] - 1, bps[1] + 1))
    else
        vals = [Float64(N(p)) for p in bps]
        for i in 0:length(bps)
            lo = i == 0 ? vals[1] - 2 : vals[i]
            hi = i == length(bps) ? vals[end] + 2 : vals[i+1]
            push!(pts, Sym(SymNumber(rationalize(BigInt, (lo + hi) / 2; tol = (hi - lo) / 8))))
        end
    end
    # (no piece holds on an interval: the integrand, and so the antiderivative, is not
    # defined there)
    Fs = Union{Sym,Nothing}[]
    numeric = !any(_hassymbol, bps)
    for (i, p) in enumerate(pts)
        gi = _pieces_at(g, x, p)
        F = gi === nothing ? nothing : integ(gi)
        if F isa Sym && numeric
            lo = i == 1 ? -Inf : Float64(N(bps[i-1]))
            hi = i == length(pts) ? Inf : Float64(N(bps[i]))
            F = _real_logs(F, x, p, lo, hi)
        end
        push!(Fs, F)
    end
    all(F -> F === nothing, Fs) &&
        throw(ArgumentError("integrate: no piece of $g holds anywhere"))
    last = findlast(F -> F !== nothing, Fs)
    # continuity at each breakpoint
    K = ZERO
    ps = Pair[]
    for i in 1:last
        if Fs[i] === nothing
            push!(ps, (x ≪ bps[i]) => Sym(NaN))
            continue
        end
        if i > 1 && Fs[i-1] !== nothing
            b = bps[i-1]
            # at a singular point (1/|sin x| at 0) the pieces are separate: nothing to match
            # (a limit that fails there, as of -log(x - 1) at 1 from below, counts as one)
            l = _side_value(Fs[i-1], x, b, :-)
            r = _side_value(Fs[i], x, b, :+)
            (l isa Sym && r isa Sym && _finite_value(l) && _finite_value(r)) && (K = simplify(K + l - r))
        end
        F = simplify(Fs[i] + K)
        push!(ps, (i == length(Fs) ? true : x ≪ bps[i]) => F)
    end
    # a periodic integrand (|sin x|): the pieces above hold for one period's breakpoints
    # only; repeated by the period, with the integral over a period per period passed
    if _periodic_conditions(g, x)
        P = _trig_period(g, x)
        R = P !== nothing && !any(_hassymbol, bps) && all(F -> F !== nothing, Fs) ?
            _periodic_antiderivative(ps, bps, x, P) : nothing
        R === nothing && throw(ArgumentError("integrate: the pieces of $g change at infinitely many points, and an antiderivative for all of them was not found"))
        return out + R
    end
    out + piecewise(ps...)
end

# on an interval (lo, hi) of real x where a linear u = c₁x + c₀ (numbers) is negative,
# log(u) as log(-u), a real antiderivative (they differ by the constant iπ): ∫ 1/|x| dx
# is -log(-x) for x < 0 (SymPy #22533)
function _real_logs(F::Sym, x::Sym, p::Sym, lo::Float64, hi::Float64)
    _hasfn(F, :log) || return F
    _replace_walk(F) do t
        (t.kind == FN && t.name === :log && length(t.args) == 1) || return nothing
        u = t.args[1]
        (has(u, x) && all(s -> _eq(s, x), free_symbols(u))) || return nothing
        c1 = try coeff(expand(u), x, 1) catch; return nothing end
        c0 = expand(u - c1 * x)
        (c1.kind == NUMBER && c0.kind == NUMBER && !iszero(c1.num)) || return nothing
        r = -Float64(c0.num) / Float64(c1.num)
        lo < r < hi && return nothing
        Float64(N(subs(u, x => p))) < 0 || return nothing
        log(-u)
    end
end

# (k, c, h) with g(x) = h(k x + c), k a nonzero number, c ≠ 0 in parameters and h free
# of them, k x + c the argument of a function in g; or nothing
function _linear_shift(g::Sym, x::Sym)
    args = Sym[]
    _fn_args!(args, g, x)
    for a in args
        k = expand(diff(a, x))
        (k.kind == NUMBER && !iszero(k)) || continue
        c = expand(a - k * x)
        _hassymbol(c) || continue
        h = expand(subs(g, x => (x - c) / k))
        any(v -> has(h, v), free_symbols(c)) && continue
        return k, c, h
    end
    nothing
end

function _fn_args!(out, e::Sym, x::Sym)
    has(e, x) || return out
    e.kind == FN && length(e.args) == 1 && !(e.name in _PIECEFNS) && push!(out, e.args[1])
    foreach(a -> _fn_args!(out, a, x), e.args)
    out
end

# does a condition of g change at infinitely many points: an argument periodic in x
function _periodic_conditions(g::Sym, x::Sym)
    has(g, x) || return false
    if g.kind == FN && (haskey(_COND_RELS, g.name) || g.name === :heaviside) && has(g, x)
        u = haskey(_COND_RELS, g.name) ? g.args[1] - g.args[2] : g.args[1]
        _trig_period(u, x) === nothing || return true
    end
    any(a -> _periodic_conditions(a, x), g.args) || (g.kind == MUL && any(a -> _periodic_conditions(a, x), g.exps))
end

# F(x) = I k + F₁(x - kP) for F₁ = piecewise(ps...), continuous on [c, c + P] with c the
# first breakpoint, k = floor((x - c)/P) and I = F₁(c + P) - F₁(c); each piece V(x - kP)
# as V(x) - k (V(x + P) - V(x)), a constant for a periodic integrand
function _periodic_antiderivative(ps, bps, x::Sym, P::Sym)
    c = bps[1]
    _isreal_const(c) === true || return nothing
    N(bps[end]) <= N(c + P) + 1e-12 || return nothing          # more than a period of points
    F1 = piecewise(ps...)
    l, r = _side_value(F1, x, c + P, :-), _side_value(F1, x, c, :+)
    (l isa Sym && r isa Sym) || return nothing
    # singular at the breakpoints (1/|sin x|): the pieces are not joined there, so no
    # constant per period either
    I = _finite_value(l) && _finite_value(r) ? simplify(l - r) : ZERO
    k = floor((x - c) / P)
    m = iszero(c) ? mod(x, P) : mod(x - c, P) + c
    out = Pair[]
    top = N(c + P) - 1e-12
    for (i, (cond, V)) in enumerate(ps)
        i == 1 && cond != ONE && continue                       # x < c: m ≥ c never
        D = simplify(subs(V, x => x + P) - V)
        has(D, x) && return nothing
        # c ≤ m < c + P: the piece reaching c + P is the last
        last = cond == ONE || (i <= length(bps) && N(bps[i]) >= top)
        push!(out, (last ? ONE : subs(_cond(cond), x => m)) => V - k * D)
        last && break
    end
    I * k + piecewise(out...)
end

# g = piecewise(u ⋚ 0 => v₁, true => v₀) with u quadratic in x, a numeric leading
# coefficient, as piecewise(D ≥ 0 => v₀ + (v₁ - v₀)(H(x - r₁) - H(x - r₂)), true => v₀)
# (the roles of v₀ and v₁ swapped when the condition holds outside the roots), for the
# discriminant D, conditions free of x and r₁ ≤ r₂ where D ≥ 0; else nothing
function _quadratic_conditions(g::Sym, x::Sym)
    (g.kind == FN && g.name === :piecewise && length(g.args) == 4 && g.args[3] == ONE) || return nothing
    c, v1, _, v0 = g.args
    (haskey(_COND_RELS, c.name) && !has(v1, x) && !has(v0, x)) || return nothing
    r = _relation_of(c)
    r.op in (:lt, :le, :gt, :ge) || return nothing
    u = expand(r.lhs - r.rhs)
    _ispolyin(u, [x]) || return nothing
    cs = [expand(subs(diff(u, x, k), x => 0) / factorial(k)) for k in 0:2]
    (iszero(expand(u - (cs[1] + cs[2] * x + cs[3] * x^2))) && cs[3].kind == NUMBER && !iszero(cs[3])) || return nothing
    C, B, A = cs
    D = expand(B^2 - 4A * C)
    has(D, x) && return nothing
    sq = sqrt(D)
    lo, hi = (-B - sq) / (2A), (-B + sq) / (2A)                  # sq ≥ 0
    r1, r2 = A.num > N0 ? (lo, hi) : (hi, lo)
    inside = (r.op in (:lt, :le)) == (A.num > N0)          # where u ⋚ 0 holds: between the roots
    a, b = inside ? (v0, v1) : (v1, v0)
    piecewise((D ≧ 0) => a + (b - a) * (heaviside(x - r1) - heaviside(x - r2)), true => a)
end

# F's one-sided limit at b, or (where the limit fails, as for complex exponentials) its
# value there when F is defined at b: an antiderivative is continuous where defined.
# Nothing when neither is known.
function _side_value(F::Sym, x::Sym, b::Sym, dir::Symbol)
    l = try limit(F, x => b; dir) catch err; err isa InterruptException && rethrow(); nothing end
    l isa Sym && return l
    v = try subs(F, x => b) catch err; _isdivzero(err) || rethrow(); return nothing end
    _defined_at(F, x, b) && _finite_value(v) ? v : nothing
end

# Σ J dirac(x - p) over the breakpoints p of the piecewise function e where its values
# from the left and right differ by J ≠ 0; zero when the breakpoints are not all known
# numbers or expressions in parameters, or repeat with a period
function _pw_jumps(e::Sym, x::Sym)
    any(i -> has(e.args[i], x), 1:2:length(e.args)) || return ZERO
    _periodic_conditions(e, x) && return ZERO
    bps = _breakpoints(e, x)
    bps === nothing && return ZERO
    out = ZERO
    for b in bps
        l, r = _side_value(e, x, b, :-), _side_value(e, x, b, :+)
        (l isa Sym && r isa Sym && _finite_value(l) && _finite_value(r)) || continue
        J = simplify(r - l)
        iszero(J) || (out = add(out, mul(J, dirac(x - b))))
    end
    out
end

function _integrate_dirac(t::Sym, x::Sym)
    fs = t.kind == MUL ? [(b, p) for (b, p) in zip(t.args, t.exps)] : [(t, ONE)]
    i = findfirst(((b, p),) -> b.kind == FN && b.name === :dirac && p == ONE, fs)
    i === nothing && return nothing
    d = fs[i][1]
    k = length(d.args) == 2 ? _toint(d.args[2].num) : 0
    u = expand(d.args[1])
    g = cancel(t / d)
    _haspieces(g, x) && return nothing
    m = expand(subs(u, x => ONE) - subs(u, x => ZERO))
    if !(m.kind == NUMBER && !iszero(m.num) && iszero(expand(u - m * x - subs(u, x => ZERO))))
        k == 0 ? (return _integrate_dirac_roots(g, u, x)) : return nothing
    end
    a = expand(-subs(u, x => ZERO) / m)
    gk = g
    for _ in 1:k
        gk = diff(gk, x)
    end
    (-1)^k * subs(gk, x => a) / (m^k * abs(m)) * heaviside(x - a)
end

# ∫ g(x) dirac(u(x)) dx = Σ g(r)/|u′(r)| heaviside(x - r) over the real roots r of u, all
# simple (u′(r) ≠ 0) and numeric, for a polynomial u (solve has all its roots; of
# sin(x) only those in one period); nothing otherwise
function _integrate_dirac_roots(g::Sym, u::Sym, x::Sym)
    _ispoly_in(u, x) || return nothing
    rs = try
        _solve_generic(u ~ 0, x)
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
    rs isa AbstractVector || return nothing
    du = diff(u, x)
    out = ZERO
    for r in rs
        _hassymbol(r) && return nothing
        v = try _evalc(r) catch; return nothing end
        abs(imag(v)) > 1e-12 * max(1, abs(v)) && continue        # not on the real line
        s = simplify(subs(du, x => r))
        _hassymbol(s) && return nothing
        abs(_evalc(s)) < 1e-12 && return nothing                   # a double root: no value
        out += subs(g, x => r) / abs(s) * heaviside(x - r)
    end
    out
end
_ispoly_in(u::Sym, x::Sym) = try
    Poly{Sym}(u, [x]); true
catch err
    err isa ArgumentError || rethrow()
    false
end

# the real numbers at which the conditions (and heaviside arguments) in e change, sorted
function _breakpoints(e::Sym, x::Sym)
    out = Sym[]
    ok = Ref(true)
    function visit(n::Sym)
        if n.kind == FN && (haskey(_COND_RELS, n.name) || n.name === :heaviside) && has(n, x)
            u = haskey(_COND_RELS, n.name) ? n.args[1] - n.args[2] : n.args[1]
            rs = try _solve_generic(u ~ 0, x) catch; (ok[] = false; Sym[]) end
            rs isa AbstractVector || (ok[] = false; return)
            for r in rs
                (_hassymbol(r) || _isreal_const(r) === true) && !any(q -> iszero(expand(q - r)), out) && push!(out, r)
            end
            n.name === :heaviside && return
        end
        foreach(visit, n.args)
        n.kind == MUL && foreach(visit, n.exps)
    end
    visit(e)
    (ok[] && !isempty(out)) || return nothing
    any(_hassymbol, out) && return length(out) == 1 ? out : nothing
    sort!(out; by = r -> real(_evalc(r)))
end

# e with each piecewise and heaviside in x replaced by its value at x = p (a point
# inside an interval, so no condition is on its boundary)
function _pieces_at(e::Sym, x::Sym, p::Sym; extra = ())
    _haspieces(e, x) || return e
    args = Sym[]
    for a in e.args
        r = _pieces_at(a, x, p; extra)
        r === nothing && return nothing
        push!(args, r)
    end
    e = _rebuild_node(e, args)
    e.kind == FN && e.name in _PIECEFNS || return e
    if e.name === :piecewise
        for (c, v) in _pieces(e)
            cv = _truth(subs(c, x => p, extra...))
            cv == ONE && return v
            cv == ZERO || return nothing
        end
        return nothing
    elseif e.name === :heaviside
        s = _csign(expand(subs(e.args[1], x => p, extra...)))
        return s === nothing ? nothing : s > 0 ? ONE : ZERO
    end
    nothing                                              # dirac in a product
end

## folding: an expression with piecewise functions as one piecewise function

"""
    SymJ.piecewise_fold(e, x; assume = nothing)
    SymJ.piecewise_fold(e; assume = nothing)

`e`, with piecewise functions (and `heaviside`, `sign`, `abs`, `max`, `min`) in `x`
inside sums, products and other functions, as one piecewise function of `x`: its value
between consecutive breakpoints (where a condition changes) and at each breakpoint,
neighbours with the same value merged. Without `x`, the one symbol in the conditions.
The breakpoints are ordered as numbers, or, with parameters (taken to be real), by the
assumptions (`assume`, or an [`assuming`](@ref) block); when they do not decide an
order the result is a [`SymJ.Cases`](@ref), a folded function for each (`a < b`,
`a = b`, `a > b`). Not exported.

```julia
@syms x a
SymJ.piecewise_fold(abs(x) + abs(x - 1), x)
# piecewise(x ≦ 0 => -2x + 1, x ≦ 1 => 1, true => 2x - 1)
SymJ.piecewise_fold(SymJ.piecewise(a ≫ 1 => a - 1//2, true => a^2 - a + 1//2) -
                    SymJ.piecewise(a ≫ 0 => 0, true => a^2))
# piecewise(a ≦ 0 => -a + 1/2, a ≦ 1 => a^2 - a + 1/2, true => a - 1/2)
```
"""
function piecewise_fold(e, x::Sym; assume = nothing)
    _isnamed(x) || throw(ArgumentError("piecewise_fold: $x is not a symbol"))
    A = _assumptions(assume)
    A === nothing && (A = Assumptions())
    g = piecewise(Sym(e))                    # heaviside, sign, abs, max, min as piecewise
    _hasfn(g, :dirac) && throw(ArgumentError("piecewise_fold: dirac has no value at its point"))
    # the parameters of the conditions are real (x ≪ a needs a real a): an order the
    # assumptions do not decide gives a case for each (a < b, a = b, a > b)
    ps = Set{Sym}()
    _pw_condsyms!(ps, g)
    delete!(ps, x)
    B = isempty(ps) ? A : _with_fact(A, FactAnd(Any[p ∈ ℝ for p in _sort_by_name!(collect(ps))]))
    _by_cases(C -> _piecewise_fold(g, x, C), B)
end

function _piecewise_fold(g::Sym, x::Sym, A::Assumptions)
    bps = _pw_breakpoints(g, x, A)
    isempty(bps) && return g
    # segments: (-∞, b₁), b₁, (b₁, b₂), …, bₖ, (bₖ, ∞), each with its value
    pts = Sym[]
    for i in 0:length(bps)
        p = i == 0 ? bps[1] - 1 : i == length(bps) ? bps[end] + 1 : (bps[i] + bps[i + 1]) / 2
        push!(pts, p)
        i < length(bps) && push!(pts, bps[i + 1])
    end
    vals = Sym[]
    for p in pts
        v = _pw_select(g, x, p, A)
        v === nothing && throw(ArgumentError("piecewise_fold: a condition of $g could not be decided at $x = $p"))
        push!(vals, simplify(v; assume = A))
    end
    # merge neighbours: an interval into the piece before it when the formulas agree, a
    # breakpoint when the formula before it gives its value there; a piece that starts
    # at a breakpoint takes the interval after it when that interval's formula gives the
    # breakpoint's value. A piece's condition is by its right end: x < b ending with an
    # interval, x ≤ b ending with the point b, true for the last.
    same(u, v) = iszero(simplify(u - v; assume = A))
    ps = Pair[]
    expr, ptval, ptb = nothing, nothing, nothing        # the piece being built
    for (j, v) in enumerate(vals)
        interval = isodd(j)
        b = interval ? nothing : bps[j ÷ 2]
        joins = j == 1 ? true :
            interval ? (expr !== nothing ? same(expr, v) : same(subs(v, x => ptb), subs(ptval, x => ptb))) :
            (expr !== nothing && same(subs(expr, x => b), subs(v, x => b)))
        if !joins
            # close the piece before segment j at the end of segment j - 1
            push!(ps, (isodd(j - 1) ? x ≪ bps[(j - 1 + 1) ÷ 2] : x ≦ bps[(j - 1) ÷ 2]) => (expr === nothing ? ptval : expr))
            expr, ptval, ptb = nothing, nothing, nothing
        end
        if interval
            expr === nothing && (expr = v)
        elseif expr === nothing
            ptval, ptb = subs(v, x => b), b             # (a point's value at the point)
        end
    end
    push!(ps, true => (expr === nothing ? ptval : expr))
    piecewise(ps...)
end
function piecewise_fold(e; assume = nothing)
    g = piecewise(Sym(e))
    vs = Set{Sym}()
    _pw_condsyms!(vs, g)
    length(vs) == 1 || throw(ArgumentError("piecewise_fold: give the variable; the conditions have $(isempty(vs) ? "no symbol" : join(_sort_by_name!(collect(vs)), ", "))"))
    piecewise_fold(g, only(vs); assume)
end

function _pw_condsyms!(vs, n::Sym)
    if n.kind == FN && n.name === :piecewise
        for (c, _) in _pieces(n)
            union!(vs, free_symbols(c))
        end
    end
    foreach(a -> _pw_condsyms!(vs, a), n.args)
    n.kind == MUL && foreach(a -> _pw_condsyms!(vs, a), n.exps)
    vs
end

# the real points where a condition in x changes, in increasing order (decided as
# numbers or under the assumptions)
function _pw_breakpoints(g::Sym, x::Sym, A::Assumptions)
    out = Sym[]
    function visit(n::Sym)
        if n.kind == FN && haskey(_COND_RELS, n.name) && has(n, x)
            rs = _solve_generic(n.args[1] - n.args[2] ~ 0, x)
            rs isa AbstractVector || throw(ArgumentError("piecewise_fold: the points where $(_relation_of(n)) changes were not found"))
            for r in rs
                r = _ode_tidy(r)
                (_hassymbol(r) || _isreal_const(r) === true) || continue
                any(q -> iszero(expand(q - r)), out) || push!(out, r)
            end
        end
        foreach(visit, n.args)
        n.kind == MUL && foreach(visit, n.exps)
    end
    visit(g)
    # insertion sort: the order of two points decided (by cases when not known); equal
    # points are one
    sorted = Sym[]
    for r in out
        k, dup = 1, false
        while k <= length(sorted)
            d = expand(r - sorted[k])
            _decide(Relation(:lt, d, ZERO), A) && break                  # r < sorted[k]
            _decide(Relation(:eq, d, ZERO), A) && (dup = true; break)
            k += 1
        end
        dup || insert!(sorted, k, r)
    end
    sorted
end

# a condition at x = p: true, false, or nothing (not decided)
function _pw_holds(c::Sym, x::Sym, p::Sym, A::Assumptions)
    c = _truth(subs(c, x => p))
    c == ONE && return true
    c == ZERO && return false
    if c.kind == FN && c.name === :And
        r = true
        for a in c.args
            t = _pw_holds(a, x, p, A)
            t === false && return false
            t === nothing && (r = nothing)
        end
        return r
    end
    (c.kind == FN && haskey(_COND_RELS, c.name)) || return nothing
    _decide(Relation(_COND_RELS[c.name], expand(c.args[1] - c.args[2]), ZERO), A)
end

# g with each piecewise in x replaced by its piece at x = p (conditions decided there)
function _pw_select(g::Sym, x::Sym, p::Sym, A::Assumptions)
    isempty(g.args) && return g
    args = Sym[]
    for a in g.args
        r = _pw_select(a, x, p, A)
        r === nothing && return nothing
        push!(args, r)
    end
    if g.kind == MUL
        exps = Sym[]
        for a in g.exps
            r = _pw_select(a, x, p, A)
            r === nothing && return nothing
            push!(exps, r)
        end
        return Sym(g.num) * prod((b^q for (b, q) in zip(args, exps)); init = ONE)
    end
    g = _rebuild_node(g, args)
    (g.kind == FN && g.name === :piecewise && has(g, x)) || return g
    for (c, v) in _pieces(g)
        t = _pw_holds(c, x, p, A)
        t === true && return v
        t === nothing && return nothing
    end
    Sym(NaN)
end

## printing: piecewise(x ≪ 0 => -x, true => x), conditions as relations

function _show_cond(io::IO, c::Sym)
    if c == ONE || c == ZERO
        print(io, c == ONE)
    elseif c.kind == FN && c.name === :And
        print(io, "(")
        join(io, (sprint(_show_cond, a) for a in c.args), ", ")
        print(io, ")")
    elseif c.kind == FN && haskey(_COND_RELS, c.name)
        show(io, _relation_of(c))
    else
        show(io, c)                 # not a condition (a pattern's slot; SymPy #25328)
    end
end
function _show_piecewise(io::IO, e::Sym)
    print(io, "piecewise(")
    for (i, (c, v)) in enumerate(_pieces(e))
        i > 1 && print(io, ", ")
        _show_cond(io, c)
        print(io, " => ")
        _show(io, v, 0)
    end
    print(io, ")")
end
