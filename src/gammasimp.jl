# Gamma and combinatorial simplification: factorials, binomials and gamma functions
# whose arguments differ by integers combined into one, with rational factors (SymPy's
# combsimp and gammasimp).
#
# Internally everything is a private Γ node (_Γ): factorial(u) = Γ(u + 1), binomial(u, v)
# = Γ(u + 1)/(Γ(v + 1) Γ(u - v + 1)). The arguments fall into classes modulo the integers;
# each Γ(g + m) of a class is written Γ(g + m₀)·(g + m₀)(g + m₀ + 1)…(g + m - 1), m₀ the
# smallest shift of the class, after which the expression is a rational function in one
# Γ per class and is cancelled. Two identities join classes: the reflection formula
# Γ(a)Γ(1 - a) = π/sin(πa) and the duplication formula Γ(a)Γ(a + 1/2) = 2^(1-2a)√π Γ(2a);
# each is applied as a substitution (Γ(1 - a) = π/(sin(πa)Γ(a))), and the candidate kept
# when smaller. The result is written back in factorials and binomials, or in gamma when
# the input had gamma.

"""
    SymJ.simplify_gamma(ex)

`ex` with its factorials, binomials and gamma functions combined: those whose arguments
differ by integers are written through one of them and cancelled, and the reflection
formula `gamma(a)gamma(1 - a) = π/sin(πa)` and duplication formula `gamma(a)gamma(a +
1/2) = 2^(1 - 2a)√π gamma(2a)` are applied where the result is smaller. The result is in
factorials and binomials, or in `gamma` when `ex` has `gamma` (SpecialFunctions'
`gamma`). As SymPy's `combsimp` and `gammasimp`. Not exported; used by
[`simplify`](@ref).

```julia
SymJ.simplify_gamma(factorial(n + 2) / factorial(n))             # (n + 1)(n + 2)
SymJ.simplify_gamma(binomial(n + 1, k) / binomial(n, k))         # (n + 1)/(n - k + 1)
SymJ.simplify_gamma(factorial(n) / (factorial(k) * factorial(n - k)))   # binomial(n, k)
SymJ.simplify_gamma(factorial(n + 1) - n * factorial(n))         # factorial(n)
SymJ.simplify_gamma(gamma(x) * gamma(1 - x))                     # π/sin(πx)
```
"""
function simplify_gamma(ex)
    e = Sym(ex)
    _hasgamma(e) || return e
    usegamma = !isempty(_collect_fn!(Sym[], e, :gamma))
    G = _gnormal(_to_Γ(e))
    best, score = e, _gscore(e)
    for c in (G, _gnormal(_greflect(G)), _gnormal(_gduplicate(G)))
        c = try cancel(c) catch; c end
        r, sc = _gabsorb(c, usegamma)
        sc < score && ((best, score) = (r, sc))
    end
    best
end

# fewer occurrences of gamma functions first, then fewer nodes
_gscore(e::Sym) = (_gcount(e), _nodes(e))
_gcount(e::Sym) = (e.kind == FN && e.name in (:factorial, :binomial, :gamma)) + sum(_gcount, e.args; init = 0) + sum(_gcount, e.exps; init = 0)

_hasgamma(e::Sym) = (e.kind == FN && (e.name in (:factorial, :binomial) || (e.name === :gamma && length(e.args) == 1))) ||
    any(_hasgamma, e.args) || any(_hasgamma, e.exps)

# Γ(u): a value for a positive integer u
_Γ(u::Sym) = u.kind == NUMBER && _isint(u.num) && u.num > N0 ? factorial(u - 1) : _mkfn_node(:_Γ, Sym[u])

function _to_Γ(e::Sym)
    isempty(e.args) && return e
    e.kind == MUL && return Sym(e.num) * prod((_to_Γ(b)^_to_Γ(p) for (b, p) in zip(e.args, e.exps)); init = ONE)
    args = Sym[_to_Γ(a) for a in e.args]
    if e.kind == FN
        e.name === :factorial && length(args) == 1 && return _Γ(args[1] + 1)
        e.name === :gamma && length(args) == 1 && return _Γ(args[1])
        if e.name === :binomial && length(args) == 2
            u, v = args
            return _Γ(u + 1) / (_Γ(v + 1) * _Γ(expand(u - v) + 1))
        end
    end
    _rebuild_node(e, args)
end

function _from_Γ(e::Sym, usegamma::Bool)
    isempty(e.args) && return e
    e.kind == MUL && return Sym(e.num) * prod((_from_Γ(b, usegamma)^_from_Γ(p, usegamma) for (b, p) in zip(e.args, e.exps)); init = ONE)
    args = Sym[_from_Γ(a, usegamma) for a in e.args]
    if e.kind == FN && e.name === :_Γ
        return usegamma ? apply_function(:gamma, args) : factorial(expand(args[1] - 1))
    end
    _rebuild_node(e, args)
end

# the Γ nodes of e
_gnodes(e::Sym) = _collect_fn!(Sym[], e, :_Γ)

# u = g + m with m the floor of u's rational constant: (g, m)
function _gsplit(u::Sym)
    c = u.kind == ADD ? u.num : u.kind == NUMBER ? u.num : N0
    isexact(c) || return (u, 0)
    q = Rational{BigInt}(_numerator(c), _denominator(c))
    m = Int(floor(q))
    expand(u - m), m
end

# Γ(y + s)/Γ(y), s an integer
_poch(y::Sym, s::Int) = s >= 0 ? prod((y + i for i in 0:(s - 1)); init = ONE) : 1 / prod((y + i for i in s:-1); init = ONE)

# each Γ written through the one of its class with the smallest shift, then cancelled
function _gnormal(e::Sym)
    gs = _gnodes(e)
    isempty(gs) && return e
    split = [_gsplit(g.args[1]) for g in gs]
    m0 = Dict{Sym,Int}()
    for (g, m) in split
        m0[g] = min(get(m0, g, m), m)
    end
    sub = Pair{Sym,Sym}[]
    for (node, (g, m)) in zip(gs, split)
        m == m0[g] && continue
        y = g + m0[g]
        push!(sub, node => _Γ(y) * _poch(y, m - m0[g]))
    end
    isempty(sub) ? e : subs(e, sub...)
end

# the classes g as (g, Γ node with the smallest shift)
function _gclasses(e::Sym)
    out = Dict{Sym,Tuple{Int,Sym}}()
    for node in _gnodes(e)
        g, m = _gsplit(node.args[1])
        (!haskey(out, g) || m < out[g][1]) && (out[g] = (m, node))
    end
    out
end

# Γ(b) = π/(sin(πa)Γ(a)) · Γ(b)/Γ(1 - a) for classes with a + b an integer; b the class
# with the larger key (1 - x, not x). sin(πa) = (-1)^m sin(πg) for a = g + m, and
# expanded when that is smaller (sin(π(x + 1/2)) = cos(πx)).
function _greflect(e::Sym)
    cl = collect(_gclasses(e))
    sub = Pair{Sym,Sym}[]
    done = Set{Sym}()
    for i in eachindex(cl), j in eachindex(cl)
        i == j && continue
        (gi, (mi, ni)), (gj, (mj, nj)) = cl[i], cl[j]
        (gi in done || gj in done) && continue
        _nodes(gj) > _nodes(gi) || (_nodes(gj) == _nodes(gi) && j > i) || continue
        a, b = ni.args[1], nj.args[1]
        t = expand(a + b)
        (t.kind == NUMBER && _isint(t.num)) || continue
        sn = sin(expand(PI * gi))
        sx = _expand_trig(sn)
        _nodes(sx) < _nodes(sn) && (sn = sx)
        push!(sub, nj => (isodd(mi) ? -PI : PI) / (sn * _Γ(a)) * _poch(1 - a, _toint(t.num) - 1))
        push!(done, gj)
    end
    isempty(sub) ? e : subs(e, sub...)
end

# Γ(a + 1/2) = 2^(1 - 2a)√π Γ(2a)/Γ(a) for classes a and a + 1/2
function _gduplicate(e::Sym)
    cl = _gclasses(e)
    sub = Pair{Sym,Sym}[]
    for (g, (m, node)) in cl
        h = expand(g + Sym(1 // 2))
        hg, hm = _gsplit(h)
        hm == 0 && haskey(cl, hg) || continue
        a = node.args[1]
        mb, nb = cl[hg]
        # Γ(b) = Γ(a + 1/2)·Γ(b)/Γ(a + 1/2), b = a + 1/2 + s
        s = _toint(expand(nb.args[1] - a - Sym(1 // 2)).num)
        push!(sub, nb => Sym(2)^(1 - 2a) * sqrt(PI) * _Γ(expand(2a)) / _Γ(a) * _poch(expand(a + Sym(1 // 2)), s))
    end
    isempty(sub) ? e : subs(e, sub...)
end

# a cancelled Γ form back in factorials (with binomials) or gamma, polynomial factors
# factored, and the score
function _gback(e::Sym, usegamma::Bool)
    r = _from_Γ(e, usegamma)
    r = something(try _gfactor(r) catch; nothing end, r)
    usegamma || (r = _gbinomials(r))
    n, d = numer_denom(r)
    if d.kind != NUMBER                              # (-n - 1)/(k - n - 1) as (n + 1)/(n + 1 - k)
        f = (-n) / (-d)
        (_nodes(f), _gnegs(-n) + _gnegs(-d)) < (_nodes(r), _gnegs(n) + _gnegs(d)) && (r = f)
    end
    r, _gscore(r)
end
_gnegs(p::Sym) = p.kind == ADD ? count(c -> c < N0, p.coefs) + (p.num < N0) : 0

# the classes' Γ(y) moved to Γ(y + s) = Γ(y)(y)ₛ where that absorbs factors ((n + 1)
# factorial(n) = factorial(n + 1)): every combination of shifts in -1:1 (-2:2 for one
# class), for up to four classes
function _gabsorb(e::Sym, usegamma::Bool)
    r, score = _gback(e, usegamma)
    nodes = [node for (_, (_, node)) in _gclasses(e)]
    (isempty(nodes) || length(nodes) > 4) && return r, score
    range = length(nodes) == 1 ? (-2:2) : (-1:1)
    for ss in Iterators.product((range for _ in nodes)...)
        all(iszero, ss) && continue
        sub = [node => _Γ(expand(node.args[1] + s)) / _poch(node.args[1], s) for (node, s) in zip(nodes, ss) if s != 0]
        c = try cancel(subs(e, sub...)) catch; continue end
        r2, sc = _gback(c, usegamma)
        sc < score && ((r, score) = (r2, sc))
    end
    r, score
end

# the polynomial factors of numerator and denominator factored, where they split
function _gfactor(e::Sym)
    n, d = numer_denom(e)
    split(p) = _ispolylike(p) && (f = factor(p); f.kind == MUL && (length(f.args) > 1 || f.exps[1] != ONE)) ? f : p
    fac(p) = p.kind == MUL ? Sym(p.num) * prod((split(b)^q for (b, q) in zip(p.args, p.exps)); init = ONE) : split(p)
    fac(n) / fac(d)
end

# factorial ratios as binomials, in each term
function _gbinomials(e::Sym)
    e.kind == ADD && return Sym(e.num) + _add_many(Sym(c) * _gbinomials(t) for (c, t) in zip(e.coefs, e.args))
    _binomial_form(e)
end
