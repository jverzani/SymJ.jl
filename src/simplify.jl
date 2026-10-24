# Simplification.
#
# `refine` applies the identities that need facts, decided by the assumptions core:
# sqrt(x^2) = |x|, |x| = x for x ≥ 0, (x^a)^b = x^(ab), log(exp(x)) = x, log(x^p) = p
# log(x), (-1)^(2n) = 1, sin(x + nπ) = (-1)^n sin(x), floor(n) = n. It is applied bottom
# up, each rule where its condition is known to hold, so the result is equal to the input
# for every value the assumptions allow.
#
# `simplify` refines, then searches: each specialized routine (simplify_rational,
# simplify_powers, simplify_trig, simplify_log, simplify_radicals, expansion,
# factoring) gives a candidate, and the one with the fewest
# nodes is kept, the input itself if none is smaller. Function arguments are simplified
# first. The routines are Julia functions on the canonical form; there is no general
# pattern matcher, as the canonical form already orders sums and products.

## refine

"""
    refine(ex; assume = nothing)

`ex` with the identities applied that hold under the assumptions (`assume`, or an
[`assuming`](@ref) block; see [`Assumptions`](@ref)):

- `sqrt(x^2) = abs(x)` for real `x`, `abs(x) = x` for `x ≥ 0` (`-x` for `x ≤ 0`), and
  `sign(x)` when the sign is known;
- `(x^a)^b = x^(a*b)` for `x ≥ 0`, `exp(u)^b = exp(b*u)` for real `u`, and
  `abs(x)^2 = x^2` for real `x`;
- `log(exp(u)) = u` for real `u`, `log(x^p) = p*log(x)` for `x > 0`;
- `(-1)^n = ±1` for an integer `n` of known parity, `sin`, `cos`, `tan` of `u + n*π` for
  an integer `n` (`sin(n*π) = 0`), and `floor(n) = ceil(n) = n`;
- equalities in the assumptions are substituted.

A relation is refined on both sides and is `true` or `false` when the assumptions decide
it. Symbols are complex unless the assumptions say otherwise, so without them only the
identities that hold for complex numbers apply.

```julia
refine(sqrt(x^2); assume = @assume(x > 0))      # x
refine(sqrt(x^2); assume = @assume(x in Reals)) # abs(x)
refine(sin(x + 2n*π); assume = @assume(n in Integers))   # sin(x)
refine(x > 0; assume = @assume(x > 1))          # true
```
"""
function refine(ex; assume = nothing)
    A = _assumptions(assume)
    _refine(Sym(ex), A === nothing ? Assumptions() : A)
end
function refine(r::Relation; assume = nothing)
    A = _assumptions(assume)
    A === nothing && (A = Assumptions())
    l, rr = _refine(r.lhs, A), _refine(r.rhs, A)
    v = _holds(Relation(r.op, expand(l - rr), ZERO), A)
    v === nothing ? Relation(r.op, l, rr) : v
end

function _refine(e::Sym, A::Assumptions)
    isempty(A.subst) || (e = subs(e, A.subst...))
    _refine_node(e, A)
end

function _refine_node(e::Sym, A::Assumptions)
    isempty(e.args) && return e
    n = e.kind == MUL ? Sym(e.num) * prod((_refine_node(b^p, A) for (b, p) in zip(e.args, e.exps)); init = ONE) :
        _rebuild_node(e, Sym[_refine_node(a, A) for a in e.args])
    r = _refine_rule(n, A)                        # a rule's result may itself refine
    r === nothing ? n : _refine_node(r, A)
end

# e with its arguments replaced (a product's are its factors' bases and exponents)
function _rebuild_node(e::Sym, args::Vector{Sym})
    k = e.kind
    k == ADD && return Sym(e.num) + _add_many(Sym(c) * a for (c, a) in zip(e.coefs, args))
    k == MUL && return Sym(e.num) * prod((b^p for (b, p) in zip(args, e.exps)); init = ONE)
    k == POW && return args[1]^args[2]
    k == FN && return apply_function(e.name, args)
    e
end

_nonneg(u, A) = (m = _signs(u, A); m !== nothing && m & _NEG == 0)
_nonpos(u, A) = (m = _signs(u, A); m !== nothing && m & _POS == 0)
_ispos(u, A) = _signs(u, A) == _POS

# one identity at the node, or nothing
function _refine_rule(e::Sym, A::Assumptions)
    k = e.kind
    if k == POW
        b, q = e.args
        if b.kind == NUMBER && b.num == NM1 && !(q.kind == NUMBER)
            _isinteger(expand(q / 2), A) === true && return ONE              # (-1)^(2n)
            _isinteger(expand((q - 1) / 2), A) === true && return -ONE
            return nothing
        end
        if b.kind == POW                                                     # (c^p)^q
            c, p = b.args
            _eq(c, E) && _isreal(p, A) === true && return exp(expand(p * q))   # exp(u)^q
            _nonneg(c, A) && _isreal(p, A) === true && return c^expand(p * q)
            if p.kind == NUMBER && isexact(p.num) && isone(_denominator(p.num)) && iseven(_numerator(p.num)) &&
               q.kind == NUMBER && _isreal(c, A) === true
                return abs(c)^expand(p * q)                                  # sqrt(x^2) = |x|
            end
        end
        # |u|^(2k) = u^(2k) for real u
        if b.kind == FN && b.name === :abs && q.kind == NUMBER && isexact(q.num) &&
           isone(_denominator(q.num)) && iseven(_numerator(q.num)) && _isreal(b.args[1], A) === true
            return b.args[1]^q
        end
        return nothing
    elseif k == FN && length(e.args) == 1
        u = e.args[1]
        f = e.name
        if f === :abs
            _nonneg(u, A) && return u
            _nonpos(u, A) && return -u
            if u.kind == MUL                         # factors of known sign come out
                out, rest = Sym(abs(u.num)), Sym[]
                for (b, p) in zip(u.args, u.exps)
                    t = b^p
                    _nonneg(t, A) ? (out *= t) : _nonpos(t, A) ? (out *= -t) : push!(rest, t)
                end
                (out != Sym(abs(u.num)) || !isone(abs(u.num))) && length(rest) < length(u.args) &&
                    return out * abs(prod(rest; init = ONE))
            end
            return nothing
        elseif f === :sign
            m = _signs(u, A)
            m == _POS && return ONE
            m == _NEG && return -ONE
            m == _ZER && return ZERO
            return nothing
        elseif f === :log
            if u.kind == POW
                b, p = u.args
                _eq(b, E) && _isreal(p, A) === true && return p              # log(exp(v)) = v
                _ispos(b, A) && _isreal(p, A) === true && return p * log(b)  # log(b^p), b > 0
            end
            return nothing
        elseif f in (:sin, :cos, :tan)
            return _refine_trig(f, u, A)
        elseif f in (:floor, :ceil)
            _isinteger(u, A) === true && return u
            if u.kind == ADD                          # integer terms come out
                n = [Sym(c) * t for (c, t) in zip(u.coefs, u.args) if _isinteger(Sym(c) * t, A) === true]
                isone(_denominator(u.num)) && push!(n, Sym(u.num))
                isempty(n) && return nothing
                m = _add_many(n)
                iszero(m) && return nothing
                return apply_function(f, Sym[expand(u - m)]) + m
            end
            return nothing
        end
    end
    nothing
end

# sin, cos, tan of u + mπ, m an integer: (-1)^m sin, (-1)^m cos, tan of the rest
function _refine_trig(f::Symbol, u::Sym, A::Assumptions)
    r = expand(u / PI)
    terms = r.kind == ADD ? [Sym(c) * t for (c, t) in zip(r.coefs, r.args)] : (r.kind == NUMBER ? Sym[] : Sym[r])
    c0 = r.kind == ADD ? Sym(r.num) : r.kind == NUMBER ? r : ZERO
    m = _add_many(t for t in terms if _isinteger(t, A) === true)
    iszero(m) && return nothing
    rest = expand(u - m * PI)
    s = f === :tan ? ONE : (-1)^m
    apply_function(f, Sym[rest]) * s
end

## simplify

"""
    simplify(ex; assume = nothing)

A simpler form of `ex`, by the number of nodes in the expression tree: `ex` is refined
(see [`refine`](@ref)), its function arguments simplified, and then several
transformations are tried and the smallest result kept (`ex` itself if none is smaller):

- rational functions over a common denominator, in lowest terms;
- powers: radicals of numbers combined (`sqrt(2)*sqrt(3) = sqrt(6)`), and of
  nonnegative quantities (`sqrt(x)*sqrt(y) = sqrt(x*y)` for `x, y ≥ 0`);
- trigonometric identities: a normal form modulo `sin^2 + cos^2 = 1`, multiple angles
  and sums of angles (`2sin(x)cos(x) = sin(2x)`, `sin(x)cos(y) + cos(x)sin(y) =
  sin(x + y)`);
- logarithms combined or expanded where their arguments are positive
  (`log(2) + log(3) = log(6)`);
- radicals denested (`sqrt(3 + 2sqrt(2)) = 1 + sqrt(2)`), denominators rationalized
  (`1/(1 + sqrt(2)) = sqrt(2) - 1`), radicands factored (`sqrt(x^2 + 2x + 1)`);
- expansion and factoring.

The transformations are also available on their own, unexported:
`SymJ.simplify_rational`, `SymJ.simplify_powers`, `SymJ.simplify_trig`,
`SymJ.simplify_log` and `SymJ.simplify_radicals`. Assumptions are given as for
`refine`.

```julia
simplify((x^2 - 1) / (x - 1))                    # x + 1
simplify((x + 1)^2 - x^2)                        # 2x + 1
simplify(sqrt(Sym(2)) * sqrt(Sym(3)))            # sqrt(6)
simplify(sqrt(x^2) / x; assume = @assume(x > 0)) # 1
simplify(1 - cos(x)^2)                           # sin(x)^2
```
"""
function simplify(ex; assume = nothing)
    A = _assumptions(assume)
    A === nothing && (A = Assumptions())
    _simplify(_refine(Sym(ex), A), A)
end
simplify(r::Relation; assume = nothing) =
    (v = refine(r; assume); v isa Bool ? v : Relation(v.op, simplify(v.lhs; assume), simplify(v.rhs; assume)))

function _simplify(e::Sym, A::Assumptions)
    _hassymbol(e) || e.kind == NUMBER || return _simplify_const(e, A)
    # function arguments first
    e = _map_fn_args(a -> _simplify(a, A), e)
    best, size = e, _nodes(e)
    powers(x) = simplify_powers(x; assume = A)
    logs(x) = simplify_log(x; assume = A)
    radicals(x) = simplify_radicals(x; assume = A)
    expansion(x) = _bigpower(x) ? nothing : expand(x)   # (x + y)^30 expanded never wins
    for g in (simplify_rational, powers, expansion, _factor_or_nothing, simplify_trig, logs, radicals,
              simplify_rational ∘ powers, expansion ∘ simplify_rational, expansion ∘ radicals)
        c = try g(e) catch err; err isa InterruptException && rethrow(); nothing end
        c isa Sym || continue
        c = _refine(c, A)
        n = _nodes(c)
        n < size && ((best, size) = (c, n))
    end
    best
end
# a constant: powers, logarithms and radicals of numbers
function _simplify_const(e::Sym, A::Assumptions)
    best = e
    for g in (simplify_powers, simplify_log, simplify_radicals)
        c = try g(e; assume = A) catch err; err isa InterruptException && rethrow(); e end
        c = try expand(c) catch; c end
        _nodes(c) < _nodes(best) && (best = c)
    end
    best
end

_bigpower(e::Sym) = (e.kind == POW && e.args[1].kind == ADD && e.args[2].kind == NUMBER && abs(e.args[2].num) > SymNumber(10)) ||
    (e.kind == MUL && any(i -> e.args[i].kind == ADD && e.exps[i].kind == NUMBER && abs(e.exps[i].num) > SymNumber(10), eachindex(e.args))) ||
    any(_bigpower, e.args)

function _factor_or_nothing(e::Sym)
    n, d = numer_denom(e)
    _ispolylike(n) && (_ispolylike(d) || d.kind == NUMBER) ? factor(n) / factor(d) : nothing
end
_ispolylike(e::Sym) = (atoms = Set{Sym}(); _poly_atoms!(atoms, e) && !isempty(atoms))

# the function applications in e with their arguments transformed
function _map_fn_args(g, e::Sym)
    isempty(e.args) && return e
    if e.kind == FN
        return apply_function(e.name, Sym[g(a) for a in e.args])
    elseif e.kind == MUL
        return Sym(e.num) * prod((_map_fn_args(g, b)^_map_fn_args(g, p) for (b, p) in zip(e.args, e.exps)); init = ONE)
    end
    _rebuild_node(e, Sym[_map_fn_args(g, a) for a in e.args])
end

# the size of an expression: its nodes, a product's exponents other than 1 and its
# coefficient other than ±1, a sum's coefficients other than ±1 and its constant term
function _nodes(e::Sym)
    k = e.kind
    k == NUMBER && return 1
    isempty(e.args) && return 1
    n = 1
    if k == ADD
        n += iszero(e.num) ? 0 : 1
        for (c, a) in zip(e.coefs, e.args)
            n += _nodes(a) + (abs(c) == N1 ? 0 : 1)
        end
    elseif k == MUL
        n += abs(e.num) == N1 ? 0 : 1
        for (b, p) in zip(e.args, e.exps)
            n += _nodes(b) + (p == ONE ? 0 : 1 + _nodes(p))
        end
    else
        for a in e.args
            n += _nodes(a)
        end
    end
    n
end

"""
    SymJ.simplify_rational(ex)

`ex` as one fraction in lowest terms: over a common denominator, with common factors of
the numerator and denominator cancelled ([`cancel`](@ref)). Not exported; used by
[`simplify`](@ref).
"""
simplify_rational(ex) = cancel(Sym(ex))

"""
    SymJ.simplify_powers(ex; assume = nothing)

`ex` with powers combined: radicals of positive numbers with the same exponent multiply
into one (`sqrt(2)*sqrt(3) = sqrt(6)`, `2^(1/3)*4^(1/3) = 2`), and so do powers of
quantities known to be nonnegative (`sqrt(x)*sqrt(y) = sqrt(x*y)` for `x, y ≥ 0`); nested
powers are taken apart as [`refine`](@ref) does. Not exported; used by
[`simplify`](@ref).
"""
function simplify_powers(ex; assume = nothing)
    A = assume isa Assumptions ? assume : (a = _assumptions(assume); a === nothing ? Assumptions() : a)
    _powers(_refine(Sym(ex), A), A)
end

function _powers(e::Sym, A::Assumptions)
    isempty(e.args) && return e
    if e.kind == MUL
        # factors grouped by a common non-integer exponent; within a group, positive
        # numbers and nonnegative quantities multiply under one power
        groups = Dict{Sym,Vector{Sym}}()
        other = Sym(e.num)
        for (b, p) in zip(e.args, e.exps)
            b2, p2 = _powers(b, A), _powers(p, A)
            intexp = p2.kind == NUMBER && isexact(p2.num) && isone(_denominator(p2.num))
            if !intexp && (b2.kind == NUMBER ? b2.num > N0 : _nonneg(b2, A))
                push!(get!(groups, p2, Sym[]), b2)
            else
                other *= b2^p2
            end
        end
        for (p, bs) in groups
            other *= prod(bs; init = ONE)^p
        end
        return other
    end
    _rebuild_node(e, Sym[_powers(a, A) for a in e.args])
end

## trigonometric simplification

"""
    SymJ.simplify_trig(ex)

`ex` with its trigonometric functions simplified. The arguments are grouped into base
angles (`2x` and `x` share the base `x`); `sin` and `cos` of multiples of a base are
written as polynomials in `S = sin(g)` and `C = cos(g)`, and `tan`, `cot`, `sec`, `csc`
through them. The resulting rational function is reduced modulo `S^2 + C^2 = 1` (once
eliminating `S^2`, once `C^2`) and cancelled, which is a normal form: an expression equal
to zero by trigonometric identities becomes `0`. The products are also linearized into
sines and cosines of multiples of the base (`2*sin(x)*cos(x)` is `sin(2x)`). The smallest
candidate is returned, with `sin/cos` written as `tan`. Not exported; used by
[`simplify`](@ref).

```julia
SymJ.simplify_trig(sin(x)^2 + cos(x)^2)          # 1
SymJ.simplify_trig(1 - cos(x)^2)                 # sin(x)^2
SymJ.simplify_trig(2*sin(x)*cos(x))              # sin(2x)
SymJ.simplify_trig(sin(2x) / sin(x))             # 2cos(x)
```
"""
function simplify_trig(ex)
    e = Sym(ex)
    bases = _trig_bases(e)
    isempty(bases) && return e
    length(bases) > 3 && return e                    # 2^k orders of elimination
    vars = [(Sym(Symbol("_trigS", i)), Sym(Symbol("_trigC", i))) for i in eachindex(bases)]
    p = _trig_poly(e, bases, vars)
    p === nothing && return e
    n, d = numer_denom(cancel(p))
    best, size = e, _nodes(e)
    back(q) = _trig_tan(subs(q, (v => f(g) for ((s, c), g) in zip(vars, bases) for (v, f) in ((s, sin), (c, cos)))...))
    for order in Iterators.product(((true, false) for _ in bases)...)
        ns = _trig_reduce(n, vars, order)
        ds = _trig_reduce(d, vars, order)
        for q in (cancel(ns / ds), ns / ds)
            r = try back(q) catch; continue end
            k = _nodes(r)
            k < size && ((best, size) = (r, k))
        end
    end
    # linearized: products of sines and cosines as sines and cosines of sums of multiples
    # of the bases (2 sin(x) cos(x) = sin(2x), sin(x) cos(y) + cos(x) sin(y) = sin(x + y))
    allv = Sym[v for sc in vars for v in sc]
    if _ispolyin(n, allv) && _ispolyin(d, allv)
        r = try _trig_tan(_trig_linear(n, vars, bases) / _trig_linear(d, vars, bases)) catch; nothing end
        r === nothing || (k = _nodes(r); k < size && ((best, size) = (r, k)))
    end
    best
end

const _TRIGFNS = (:sin, :cos, :tan, :cot, :sec, :csc)

# the base angles: arguments with their numeric content taken out, each with the gcd of
# the multiples it appears with (args x/2 and x: base x/2)
function _trig_bases(e::Sym)
    seen = Dict{Sym,Rational{BigInt}}()
    function walk(n::Sym)
        if n.kind == FN && n.name in _TRIGFNS && length(n.args) == 1
            q, h = _trig_split(n.args[1])
            seen[h] = haskey(seen, h) ? gcd(seen[h], q) : q
        end
        foreach(walk, n.args)
        foreach(walk, n.exps)
    end
    walk(e)
    Sym[expand(Sym(SymNumber(q)) * h) for (h, q) in seen]
end
# u = q * h, q the numeric content of u (2x + 2 is 2*(x + 1))
function _trig_split(u::Sym)
    q = if u.kind == MUL && isexact(u.num)
        Rational{BigInt}(_numerator(u.num), _denominator(u.num))
    elseif u.kind == ADD && isexact(u.num) && all(isexact, u.coefs)
        g = reduce(gcd, (Rational{BigInt}(_numerator(c), _denominator(c)) for c in u.coefs);
                   init = Rational{BigInt}(_numerator(u.num), _denominator(u.num)))
        iszero(g) ? one(Rational{BigInt}) : abs(g)
    else
        one(Rational{BigInt})
    end
    q, expand(u / Sym(SymNumber(q)))
end

# e with trigonometric functions of multiples of the bases as polynomials (and
# quotients) in S, C; nothing if an argument is not a multiple of a base
function _trig_poly(e::Sym, bases::Vector{Sym}, vars)
    ok = Ref(true)
    function tr(n::Sym)
        isempty(n.args) && return n
        if n.kind == FN && n.name in _TRIGFNS && length(n.args) == 1
            u = n.args[1]
            i = findfirst(g -> (r = expand(u / g); r.kind == NUMBER && isexact(r.num) && isone(_denominator(r.num))), bases)
            i === nothing && (ok[] = false; return n)
            k = Int(_numerator(expand(u / bases[i]).num))
            S, C = vars[i]
            s, c = _chebyshev(k, S, C)
            f = n.name
            return f === :sin ? s : f === :cos ? c : f === :tan ? s / c : f === :cot ? c / s : f === :sec ? 1 / c : 1 / s
        end
        n.kind == MUL && return Sym(n.num) * prod((tr(b)^tr(p) for (b, p) in zip(n.args, n.exps)); init = ONE)
        _rebuild_node(n, Sym[tr(a) for a in n.args])
    end
    r = tr(e)
    ok[] ? r : nothing
end

# (sin(k g), cos(k g)) as polynomials in S = sin(g), C = cos(g)
function _chebyshev(k::Int, S::Sym, C::Sym)
    k < 0 && ((s, c) = _chebyshev(-k, S, C); return (-s, c))
    k == 0 && return (ZERO, ONE)
    s0, c0, s1, c1 = ZERO, ONE, S, C
    for _ in 2:k
        s0, c0, s1, c1 = s1, c1, expand(2C * s1 - s0), expand(2C * c1 - c0)
    end
    s1, c1
end

# p reduced modulo S^2 + C^2 = 1 for each base: S^(2m+r) -> S^r (1 - C^2)^m when the
# order says S, else C^(2m+r) -> C^r (1 - S^2)^m
function _trig_reduce(p::Sym, vars, order)
    p = expand(p)
    for ((S, C), elimS) in zip(vars, order)
        has(p, S) || has(p, C) || continue
        P = try Poly{Sym}(p, [S, C]) catch; return p end
        a, b = elimS ? (S, C) : (C, S)
        ia = elimS ? 1 : 2
        acc = ZERO
        for (ev, cf) in zip(exponent_vectors(P), coefficients(P))
            ka, kb = ev[ia], ev[3 - ia]
            acc += cf * a^(ka % 2) * (1 - b^2)^(ka ÷ 2) * b^kb
        end
        p = expand(acc)
    end
    p
end

_ispolyin(p::Sym, vs::Vector{Sym}) = (try Poly{Sym}(expand(p), vs); true catch; false end)

# a polynomial in Sᵢ = sin(gᵢ), Cᵢ = cos(gᵢ) as a combination of sines and cosines of
# k·g = Σ kᵢ gᵢ: with Sᵢ = (zᵢ - 1/zᵢ)/(2i), Cᵢ = (zᵢ + 1/zᵢ)/2 it is a Laurent
# polynomial in zᵢ = exp(i gᵢ), whose terms z^k and z^-k pair into cos(k·g), sin(k·g)
function _trig_linear(p::Sym, vars, bases::Vector{Sym})
    zs = [Sym(Symbol("_trigz", i)) for i in eachindex(bases)]
    p = expand(p)
    D = Int[]
    for (S, C) in vars
        P = Poly{Sym}(p, [S, C])
        push!(D, maximum(sum, exponent_vectors(P); init = 0))
    end
    sub = [v => w for ((S, C), z) in zip(vars, zs) for (v, w) in ((S, (z^2 - 1) / (2IM * z)), (C, (z^2 + 1) / (2z)))]
    q = expand(subs(p, sub...) * prod(z^d for (z, d) in zip(zs, D)))
    Q = Poly{Sym}(q, zs)
    a = Dict{Vector{Int},Sym}()
    for (ev, cf) in zip(exponent_vectors(Q), coefficients(Q))
        a[collect(ev) .- D] = cf
    end
    out = get(a, zeros(Int, length(zs)), ZERO)
    for (k, ak) in a
        i = findfirst(!iszero, k)
        (i === nothing || k[i] < 0) && continue             # each pair once, from k > 0
        am = get(a, -k, ZERO)
        θ = expand(sum(kj * g for (kj, g) in zip(k, bases)))
        out += expand(ak + am) * cos(θ) + expand(IM * (ak - am)) * sin(θ)
    end
    expand(out)
end

# sin(u)^a * cos(u)^-a written tan(u)^a (cot for the reverse)
function _trig_tan(e::Sym)
    isempty(e.args) && return e
    if e.kind == MUL
        fs = [(_trig_tan(b), p) for (b, p) in zip(e.args, e.exps)]
        out = Sym(e.num)
        used = falses(length(fs))
        for (i, (b, p)) in enumerate(fs)
            used[i] && continue
            if b.kind == FN && b.name in (:sin, :cos)
                other = b.name === :sin ? :cos : :sin
                j = findfirst(j -> !used[j] && j != i && fs[j][1].kind == FN && fs[j][1].name === other &&
                                   _eq(fs[j][1].args[1], b.args[1]) && iszero(expand(fs[j][2] + p)), eachindex(fs))
                if j !== nothing
                    used[i] = used[j] = true
                    out *= apply_function(b.name === :sin ? :tan : :cot, b.args)^p
                    continue
                end
            end
            used[i] = true
            out *= b^p
        end
        return out
    end
    _rebuild_node(e, Sym[_trig_tan(a) for a in e.args])
end

## logarithms

"""
    SymJ.simplify_log(ex; assume = nothing)

`ex` with its logarithms combined or expanded, whichever is smaller:

- combined, `c₁ log(u₁) + c₂ log(u₂) = log(u₁^c₁ u₂^c₂)` for `uᵢ > 0` and real `cᵢ`
  (`log(2) + log(3) = log(6)`; positive numbers need no assumptions);
- expanded, `log(a*b) = log(a) + log(b)` for `a > 0` (a positive factor does not change
  the argument of the product, so `b` may be anything), `log(a^p) = p*log(a)` for `a > 0`
  and real `p`.

Not exported; used by [`simplify`](@ref).
"""
function simplify_log(ex; assume = nothing)
    A = assume isa Assumptions ? assume : (a = _assumptions(assume); a === nothing ? Assumptions() : a)
    e = _refine(Sym(ex), A)
    best = e
    for c in (_log_combine(e, A), _log_combine(_log_expand(e, A), A), _log_expand(e, A))
        _nodes(c) < _nodes(best) && (best = c)
    end
    best
end

_logpos(u::Sym, A) = u.kind == NUMBER ? u.num > N0 : _ispos(u, A)

function _log_combine(e::Sym, A::Assumptions)
    isempty(e.args) && return e
    if e.kind == ADD
        args = Sym[_log_combine(a, A) for a in e.args]
        prodarg, rest = ONE, ZERO
        n = 0
        for (c, t) in zip(e.coefs, args)
            k, l = t.kind == MUL && length(t.args) == 1 && _eq(t.exps[1], ONE) ? (t.num, t.args[1]) : (N1, t)
            if l.kind == FN && l.name === :log && length(l.args) == 1 && _logpos(l.args[1], A)
                prodarg *= l.args[1]^Sym(c * k)
                n += 1
            else
                rest += Sym(c) * t
            end
        end
        n >= 1 && return Sym(e.num) + rest + (n == 1 && isone(prodarg) ? ZERO : log(prodarg))
        return Sym(e.num) + rest
    end
    if e.kind == MUL && length(e.args) == 1 && _eq(e.exps[1], ONE)        # c log(u) = log(u^c)
        l = _log_combine(e.args[1], A)
        if l.kind == FN && l.name === :log && _logpos(l.args[1], A) && isexact(e.num) && isone(_denominator(e.num))
            return log(l.args[1]^Sym(e.num))
        end
    end
    e.kind == MUL && return Sym(e.num) * prod((_log_combine(b, A)^_log_combine(p, A) for (b, p) in zip(e.args, e.exps)); init = ONE)
    _rebuild_node(e, Sym[_log_combine(a, A) for a in e.args])
end

function _log_expand(e::Sym, A::Assumptions)
    isempty(e.args) && return e
    if e.kind == FN && e.name === :log && length(e.args) == 1
        u = _log_expand(e.args[1], A)
        if u.kind == MUL
            out, rest = ZERO, Sym(u.num)
            _logpos(Sym(u.num), A) && !isone(u.num) ? (out += log(Sym(u.num)); rest = ONE) : nothing
            for (b, p) in zip(u.args, u.exps)
                if _logpos(b, A) && _isreal(p, A) === true
                    out += p * log(b)
                else
                    rest *= b^p
                end
            end
            return out + (isone(rest) ? ZERO : log(rest))
        elseif u.kind == POW && _logpos(u.args[1], A) && _isreal(u.args[2], A) === true
            return u.args[2] * log(u.args[1])
        end
        return log(u)
    end
    e.kind == MUL && return Sym(e.num) * prod((_log_expand(b, A)^_log_expand(p, A) for (b, p) in zip(e.args, e.exps)); init = ONE)
    _rebuild_node(e, Sym[_log_expand(a, A) for a in e.args])
end

## radicals

"""
    SymJ.simplify_radicals(ex; assume = nothing)

`ex` with radicals simplified:

- denested: `sqrt(a + b*sqrt(c)) = sqrt(x) + sqrt(y)` when `a^2 - b^2 c` is the square
  of a rational (`sqrt(3 + 2sqrt(2)) = 1 + sqrt(2)`);
- denominators rationalized: `1/(a + b*sqrt(c))` multiplied by the conjugate
  (`1/(1 + sqrt(2)) = sqrt(2) - 1`);
- radicands factored, so that a power under a root comes out as `refine` allows
  (`sqrt(x^2 + 2x + 1) = abs(x + 1)` for real `x`, `x + 1` for `x ≥ -1`).

Not exported; used by [`simplify`](@ref).
"""
function simplify_radicals(ex; assume = nothing)
    A = assume isa Assumptions ? assume : (a = _assumptions(assume); a === nothing ? Assumptions() : a)
    _refine(_radicals(_refine(Sym(ex), A), A), A)
end

function _radicals(e::Sym, A::Assumptions)
    isempty(e.args) && return e
    if e.kind == MUL
        return Sym(e.num) * prod((_radicals(b^p, A) for (b, p) in zip(e.args, e.exps)); init = ONE)
    end
    if e.kind == POW
        b, q = _radicals(e.args[1], A), e.args[2]
        if q.kind == NUMBER && isexact(q.num)
            if q.num == SymNumber(1 // 2) && (d = _denest(b)) !== nothing
                return d
            end
            if isone(_denominator(q.num)) && q.num < N0 && (c = _conjugate(b)) !== nothing
                conj, norm = c                                  # 1/b = conj/norm
                return expand((conj / norm)^Sym(-q.num))
            end
            if !isone(_denominator(q.num)) && _hassymbol(b) && _ispolylike(b)
                f = try factor(b) catch; b end                  # a power under the root
                return f^q
            end
        end
        return b^q
    end
    _rebuild_node(e, Sym[_radicals(a, A) for a in e.args])
end

# b = p + r sqrt(c) with rational p, r, c (c not a square): (p, r, c), or nothing
function _quadratic_surd(b::Sym)
    b.kind == ADD && length(b.args) == 1 && isexact(b.num) || return nothing
    t, r = b.args[1], b.coefs[1]
    s = t.kind == MUL && length(t.args) == 1 ? (r *= t.num; t.args[1]^t.exps[1]) : t
    (s.kind == POW && s.args[1].kind == NUMBER && isexact(s.args[1].num) && s.args[2] == Sym(1 // 2)) || return nothing
    isexact(r) || return nothing
    (b.num, r, s.args[1].num)
end

# sqrt(p + r sqrt(c)) = sqrt(x) + sign(r) sqrt(y), x + y = p, 4xy = r^2 c, when
# p^2 - r^2 c is a rational square
function _denest(b::Sym)
    t = _quadratic_surd(b)
    t === nothing && return nothing
    p, r, c = t
    p > N0 || return nothing
    disc = Sym(p * p - r * r * c)
    disc.kind == NUMBER && disc.num >= N0 || return nothing
    s = sqrt(disc)
    s.kind == NUMBER || return nothing                         # not a rational square
    x, y = (Sym(p) + s) / 2, (Sym(p) - s) / 2
    (x.num >= N0 && y.num >= N0) || return nothing
    _powers(sqrt(x) + (r > N0 ? 1 : -1) * sqrt(y), Assumptions())   # sqrt(3/2) as sqrt(6)/2
end

# (p - r sqrt(c), p^2 - r^2 c) for b = p + r sqrt(c)
function _conjugate(b::Sym)
    t = _quadratic_surd(b)
    t === nothing && return nothing
    p, r, c = t
    norm = Sym(p * p - r * r * c)
    iszero(norm) && return nothing
    Sym(p) - Sym(r) * sqrt(Sym(c)), norm
end
