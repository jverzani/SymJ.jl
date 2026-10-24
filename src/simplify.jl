# Simplification.
#
# `refine` applies the identities that need facts, decided by the assumptions core:
# sqrt(x^2) = |x|, |x| = x for x ≥ 0, (x^a)^b = x^(ab), log(exp(x)) = x, log(x^p) = p
# log(x), (-1)^(2n) = 1, sin(x + nπ) = (-1)^n sin(x), floor(n) = n. It is applied bottom
# up, each rule where its condition is known to hold, so the result is equal to the input
# for every value the assumptions allow.
#
# `simplify` refines, then searches: each specialized routine (simplify_rational,
# simplify_powers, expansion, factoring) gives a candidate, and the one with the fewest
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
- expansion and factoring.

The transformations are also available on their own, unexported:
`SymJ.simplify_rational` and `SymJ.simplify_powers`. Assumptions are given as for
`refine`.

```julia
simplify((x^2 - 1) / (x - 1))                    # x + 1
simplify((x + 1)^2 - x^2)                        # 2x + 1
simplify(sqrt(Sym(2)) * sqrt(Sym(3)))            # sqrt(6)
simplify(sqrt(x^2) / x; assume = @assume(x > 0)) # 1
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
    expansion(x) = _bigpower(x) ? nothing : expand(x)   # (x + y)^30 expanded never wins
    for g in (simplify_rational, powers, expansion, _factor_or_nothing,
              simplify_rational ∘ powers, expansion ∘ simplify_rational)
        c = try g(e) catch err; err isa InterruptException && rethrow(); nothing end
        c isa Sym || continue
        c = _refine(c, A)
        n = _nodes(c)
        n < size && ((best, size) = (c, n))
    end
    best
end
function _simplify_const(e::Sym, A::Assumptions)
    c = try simplify_powers(e; assume = A) catch; e end
    _nodes(c) < _nodes(e) ? c : e
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
