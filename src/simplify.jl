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
# nodes is kept, the input itself if none is smaller. Function arguments (with exp(u)'s
# u and the radicand of a root) are simplified first. The routines are Julia functions
# on the canonical form; for expressions with trigonometric or hyperbolic functions, a
# few identities the trig normal form does not reach are rewrite rules (simprules.jl).

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
- `real`, `imag`, `conj`, `angle` and `abs` from the real and imaginary parts the
  assumptions give (`real(x) = x` and `conj(x) = x` for real `x`,
  `real(exp(im*x)) = cos(x)`, `abs(x + im*y) = sqrt(x^2 + y^2)` for real `x`, `y`);
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
    if k == FN && (haskey(_COND_RELS, e.name) || e.name in (:heaviside, :sign, :abs))
        r = haskey(_COND_RELS, e.name) ? _refine_cond(e, A) : _refine_step(e, A)
        r === nothing || return r
    end
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
            # |exp(w)| = exp(real(w)); otherwise sqrt(re^2 + im^2) from u's parts
            if u.kind == POW && _eq(u.args[1], E)
                return exp(_refine_node(real(u.args[2]), A))
            end
            if _isreal(u, A) !== true && (p = _reim(u, A)) !== nothing
                m2 = expand(p[1]^2 + p[2]^2)
                m2 = simplify_hyperbolic(simplify_trig(m2))
                return sqrt(m2)
            end
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
        elseif f in (:real, :imag, :conj, :angle)
            return _refine_parts(f, u, A)
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

# real, imag, conj, angle of u, from u's real and imaginary parts under the assumptions
function _refine_parts(f::Symbol, u::Sym, A::Assumptions)
    p = _reim(u, A)
    p === nothing && return nothing
    re, im = expand(p[1]), expand(p[2])
    f === :real && return re
    f === :imag && return im
    f === :conj && return iszero(im) ? u : re - IM * im
    # angle: of a real, 0 or π by its sign; else atan(im/re) by the quadrant
    if iszero(im)
        m = _signs(re, A)
        m == _POS && return ZERO
        m == _NEG && return PI
        return nothing
    end
    sr, si = _signs(re, A), _signs(im, A)
    (sr in (_POS, _NEG, _ZER) && si in (_POS, _NEG)) || return nothing
    sr == _ZER && return si == _POS ? PI / 2 : -PI / 2
    a = atan(im / re)
    sr == _POS ? a : si == _POS ? a + PI : a - PI
end

# (real part, imaginary part) of u under the assumptions, or nothing
function _reim(u::Sym, A::Assumptions)
    _isreal(u, A) === true && return (u, ZERO)
    p = _const_parts(u)
    p === nothing || return p
    _isnamed(u) && return nothing
    k = u.kind
    mulp((a, b), (c, d)) = (a * c - b * d, a * d + b * c)
    if k == ADD
        re, im = Sym(u.num), ZERO
        for (c, t) in zip(u.coefs, u.args)
            q = _reim(t, A)
            q === nothing && return nothing
            re += Sym(c) * q[1]
            im += Sym(c) * q[2]
        end
        return (re, im)
    elseif k == MUL
        acc = (Sym(u.num), ZERO)
        for (b, e) in zip(u.args, u.exps)
            q = _reim(b^e, A)
            q === nothing && return nothing
            acc = mulp(acc, q)
        end
        return acc
    elseif k == POW
        b, e = u.args
        if _eq(b, E)                                           # exp(a + ib)
            q = _reim(e, A)
            q === nothing && return nothing
            return (exp(q[1]) * cos(q[2]), exp(q[1]) * sin(q[2]))
        end
        if e.kind == NUMBER && isexact(e.num) && isone(_denominator(e.num)) && N0 < e.num <= SymNumber(8)
            q = _reim(b, A)
            q === nothing && return nothing
            acc = (ONE, ZERO)
            for _ in 1:Int(_numerator(e.num))
                acc = mulp(acc, q)
            end
            return acc
        end
        return nothing
    elseif k == FN && length(u.args) == 1
        f = u.name
        f in _REALFNS && return (u, ZERO)
        q = _reim(u.args[1], A)
        q === nothing && return nothing
        a, b = q
        f === :conj && return (a, -b)
        f === :sin && return (sin(a) * cosh(b), cos(a) * sinh(b))
        f === :cos && return (cos(a) * cosh(b), -sin(a) * sinh(b))
        f === :sinh && return (sinh(a) * cos(b), cosh(a) * sin(b))
        f === :cosh && return (cosh(a) * cos(b), sinh(a) * sin(b))
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
  sin(x + y)`), exponentials of imaginary arguments (`cos(x) + im*sin(x) = exp(im*x)`);
- hyperbolic identities, likewise, with exponentials (`cosh(x)^2 - sinh(x)^2 = 1`,
  `exp(x) + exp(-x) = 2cosh(x)`);
- logarithms combined or expanded where their arguments are positive
  (`log(2) + log(3) = log(6)`);
- radicals denested (`sqrt(3 + 2sqrt(2)) = 1 + sqrt(2)`), denominators rationalized
  (`1/(1 + sqrt(2)) = sqrt(2) - 1`), radicands factored (`sqrt(x^2 + 2x + 1)`);
- expansion and factoring.

The transformations are also available on their own, unexported:
`SymJ.simplify_rational`, `SymJ.simplify_powers`, `SymJ.simplify_trig`,
`SymJ.simplify_hyperbolic`, `SymJ.simplify_log` and `SymJ.simplify_radicals`.
Assumptions are given as for `refine`.

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
    function consider(g)
        c = try g(e) catch err; err isa InterruptException && rethrow(); nothing end
        c isa Sym || return
        c = _refine(c, A)
        n = _nodes(c)
        n < size && ((best, size) = (c, n))
    end
    for g in (simplify_rational, powers, expansion, _factor_or_nothing, simplify_trig, simplify_hyperbolic, logs, radicals, simplify_gamma,
              simplify_rational ∘ powers, expansion ∘ simplify_rational, expansion ∘ radicals)
        consider(g)
    end
    if _has_fn(e, _CIRCFNS)       # identities simplify_trig's normal form does not reach
        for g in (_simp_trig_rules, expand_trig, simplify_rational ∘ expand_trig)
            consider(g)
        end
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

# the function applications in e with their arguments transformed, counting exp(u)
# (stored as ℯ^u) and the radicand of b^(p/q) as arguments
function _map_fn_args(g, e::Sym)
    isempty(e.args) && return e
    if e.kind == FN
        return apply_function(e.name, Sym[g(a) for a in e.args])
    elseif e.kind == MUL
        return Sym(e.num) * prod((_map_pow_args(g, b, p) for (b, p) in zip(e.args, e.exps)); init = ONE)
    elseif e.kind == POW
        return _map_pow_args(g, e.args[1], e.args[2])
    end
    _rebuild_node(e, Sym[_map_fn_args(g, a) for a in e.args])
end
_map_pow_args(g, b::Sym, p::Sym) =
    b === E ? b^g(p) :
    p.kind == NUMBER && !isinteger(p) ? g(b)^p :
    _map_fn_args(g, b)^_map_fn_args(g, p)

# the size of an expression: its nodes, a product's exponents other than 1 and its
# coefficient other than ±1, a sum's coefficients other than ±1 and its constant term;
# a negative exponent counts one more, a division (as SymPy's count_ops counts it)
_division_cost(p::Sym) = p.kind == NUMBER && p.num < N0 ? 1 : 0

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
            n += _nodes(b) + (p == ONE ? 0 : 1 + _nodes(p)) + _division_cost(p)
        end
    elseif k == POW
        n += _nodes(e.args[1]) + _nodes(e.args[2]) + _division_cost(e.args[2])
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

## trigonometric and hyperbolic simplification
#
# One method for both, with a mode: the functions of a base angle g are polynomials in
# S and C (sin g, cos g; or sinh g, cosh g), reduced modulo the relation (S² + C² = 1;
# C² - S² = 1), and linearized through z = exp(i g) (or exp(g)). Exponentials join in:
# exp(i k g) = cos(kg) + i sin(kg), exp(k g) = cosh(kg) + sinh(kg); and the exponential
# form itself, z substituted back as exp(i g) or exp(g), is a candidate.

struct _Circ
    fns::NTuple{6,Symbol}          # sin cos tan cot sec csc, or the hyperbolic ones
    hyp::Bool
end
const _TRIGMODE = _Circ((:sin, :cos, :tan, :cot, :sec, :csc), false)
const _HYPMODE = _Circ((:sinh, :cosh, :tanh, :coth, :sech, :csch), true)
const _TRIGFNS = _TRIGMODE.fns

"""
    SymJ.simplify_trig(ex)

`ex` with its trigonometric functions simplified. The arguments are grouped into base
angles (`2x` and `x` share the base `x`); `sin` and `cos` of multiples of a base are
written as polynomials in `S = sin(g)` and `C = cos(g)`, and `tan`, `cot`, `sec`, `csc`
through them, as are exponentials of imaginary multiples (`exp(i*x) = cos(x) + i*sin(x)`).
The resulting rational function is reduced modulo `S^2 + C^2 = 1` (once eliminating
`S^2`, once `C^2`) and cancelled, which is a normal form: an expression equal to zero by
trigonometric identities becomes `0`. The products are also linearized into sines and
cosines of multiples of the base (`2*sin(x)*cos(x)` is `sin(2x)`), and written in
exponentials. The smallest candidate is returned, with `sin/cos` written as `tan`. Not
exported; used by [`simplify`](@ref).

```julia
SymJ.simplify_trig(sin(x)^2 + cos(x)^2)          # 1
SymJ.simplify_trig(1 - cos(x)^2)                 # sin(x)^2
SymJ.simplify_trig(2*sin(x)*cos(x))              # sin(2x)
SymJ.simplify_trig((exp(im*x) + exp(-im*x)) / 2) # cos(x)
```
"""
simplify_trig(ex) = _circ_simplify(Sym(ex), _TRIGMODE)

"""
    SymJ.simplify_hyperbolic(ex)

`ex` with its hyperbolic functions and exponentials simplified, as
[`SymJ.simplify_trig`](@ref) does for trigonometric ones: polynomials in `sinh` and
`cosh` of base arguments reduced modulo `cosh^2 - sinh^2 = 1`, linearized, or written in
exponentials, whichever is smallest (`cosh(x)^2 - sinh(x)^2 = 1`, `2sinh(x)cosh(x) =
sinh(2x)`, `cosh(x) + sinh(x) = exp(x)`, `exp(x) + exp(-x) = 2cosh(x)`). Not exported;
used by [`simplify`](@ref).
"""
simplify_hyperbolic(ex) = _circ_simplify(Sym(ex), _HYPMODE)

"""
    SymJ.expand_trig(ex)

`ex` with its trigonometric and hyperbolic functions of sums and integer multiples
expanded, bottom up, into functions of the terms: `sin(x + y) = sin(x)cos(y) +
cos(x)sin(y)`, `cos(2x) = 2cos(x)^2 - 1`, `tan(x + y) = (tan(x) + tan(y))/(1 -
tan(x)tan(y))`, and likewise `cot`, `sec`, `csc` and the hyperbolic functions. A
multiple is a polynomial in the sine and cosine of the base (Chebyshev), as is the
tangent's in the tangent; a rational multiple (`x/2`) is left. The reverse is
[`SymJ.simplify_trig`](@ref). Not exported.

```julia
SymJ.expand_trig(sin(2x))         # 2*sin(x)*cos(x)
SymJ.expand_trig(cos(3x))         # 4*cos(x)^3 - 3*cos(x)
SymJ.expand_trig(sinh(x + y))     # sinh(x)*cosh(y) + cosh(x)*sinh(y)
SymJ.expand_trig(tan(2x))         # 2*tan(x)/(1 - tan(x)^2)
```
"""
expand_trig(ex) = _expand_trig(Sym(ex))

# top down: a function is expanded on its argument as it is, then inside what remains
function _expand_trig(e::Sym)
    isempty(e.args) && return e
    if e.kind == FN && length(e.args) == 1
        for M in (_TRIGMODE, _HYPMODE)
            e.name in M.fns && return _expand_circ(e.name, e.args[1], M)
        end
    end
    e.kind == MUL && return Sym(e.num) * prod((_expand_trig(b)^_expand_trig(p) for (b, p) in zip(e.args, e.exps)); init = ONE)
    _rebuild_node(e, Sym[_expand_trig(a) for a in e.args])
end

# f(u) expanded: u a sum is split into its first term and the rest, an integer multiple
# k*h of h written through f(h)
function _expand_circ(f::Symbol, u::Sym, M::_Circ)
    s, c, t, ct, se, cs = M.fns
    σ = M.hyp ? 1 : -1                       # cos(a + b) = cos cos + σ sin sin, …
    f === se && return 1 / _expand_circ(c, u, M)
    f === cs && return 1 / _expand_circ(s, u, M)
    ex(g, v) = _expand_circ(g, v, M)
    if u.kind == ADD
        a = iszero(u.num) ? Sym(u.coefs[1]) * u.args[1] : Sym(u.num)
        b = u - a
        f === s && return expand(ex(s, a) * ex(c, b) + ex(c, a) * ex(s, b))
        f === c && return expand(ex(c, a) * ex(c, b) + σ * ex(s, a) * ex(s, b))
        f === t && return (ta = ex(t, a); tb = ex(t, b); (ta + tb) / (1 + σ * ta * tb))
        return (ca = ex(ct, a); cb = ex(ct, b); (ca * cb + σ) / (ca + cb))
    end
    if u.kind == MUL && isexact(u.num) && isone(_denominator(u.num)) && abs(u.num) > N1
        k = Int(_numerator(u.num))
        h = u / Sym(u.num)
        if f in (s, c)
            sn, cn = _chebyshev(k, ex(s, h), ex(c, h))
            return f === s ? sn : cn
        end
        # tan(k h) (cot) by the sum rule, k - 1 times, cancelled
        x1 = ex(f, h)
        r = x1
        for _ in 2:abs(k)
            r = f === t ? (r + x1) / (1 + σ * r * x1) : (r * x1 + σ) / (r + x1)
        end
        r = cancel(r)
        return k < 0 ? -r : r
    end
    apply_function(f, Sym[_expand_trig(u)])
end

function _circ_simplify(e::Sym, M::_Circ)
    bases = _circ_bases(e, M)
    isempty(bases) && return e
    length(bases) > 3 && return e                    # 2^k orders of elimination
    vars = [(Sym(Symbol("_circS", i)), Sym(Symbol("_circC", i))) for i in eachindex(bases)]
    p = _circ_poly(e, bases, vars, M)
    p === nothing && return e
    n, d = numer_denom(cancel(p))
    best, size = e, _nodes(e)
    consider(r) = (k = _nodes(r); k < size && ((best, size) = (r, k)))
    s, c = M.fns[1], M.fns[2]
    back(q) = _circ_tan(subs(q, (v => apply_function(f, Sym[g]) for ((S, C), g) in zip(vars, bases)
                                       for (v, f) in ((S, s), (C, c)))...), M)
    for order in Iterators.product(((true, false) for _ in bases)...)
        ns = _circ_reduce(n, vars, order, M)
        ds = _circ_reduce(d, vars, order, M)
        for q in (cancel(ns / ds), ns / ds)
            r = try back(q) catch; continue end
            consider(r)
        end
    end
    allv = Sym[v for sc in vars for v in sc]
    if _ispolyin(n, allv) && _ispolyin(d, allv)
        # linearized: products as functions of sums of multiples of the bases
        r = try _circ_tan(_circ_linear(n, vars, bases, M) / _circ_linear(d, vars, bases, M), M) catch; nothing end
        r === nothing || consider(r)
        # in exponentials: z = exp(i g) (exp(g)) substituted back
        r = try _circ_exp(n, vars, bases, M) / _circ_exp(d, vars, bases, M) catch; nothing end
        r === nothing || consider(cancel(r))
    end
    best
end

# exp(u) in the mode: (multiple, h) with u = multiple * i * h (trig) or multiple * h
# (hyperbolic); nothing if u is not of that kind
function _circ_expsplit(u::Sym, M::_Circ)
    if M.hyp
        has(u, IM) && return nothing
        return _trig_split(u)
    end
    v = expand(u / IM)
    has(v, IM) && return nothing
    _trig_split(v)
end

# the exponentials in e as (u, node): exp(u) nodes, and a product's factors ℯ^u
function _circ_exps!(out, n::Sym)
    if n.kind == POW && _eq(n.args[1], E)
        push!(out, n.args[2])
    elseif n.kind == MUL
        for (b, p) in zip(n.args, n.exps)
            _eq(b, E) && push!(out, p)
        end
    end
    out
end

# the base angles: arguments with their numeric content taken out, each with the gcd of
# the multiples it appears with (args x/2 and x: base x/2)
function _circ_bases(e::Sym, M::_Circ)
    seen = Dict{Sym,Rational{BigInt}}()
    nfn = Ref(0)
    add!(q, h) = (seen[h] = haskey(seen, h) ? gcd(seen[h], q) : q)
    function walk(n::Sym)
        if n.kind == FN && n.name in M.fns && length(n.args) == 1
            nfn[] += 1
            add!(_trig_split(n.args[1])...)
        end
        for u in _circ_exps!(Sym[], n)
            t = _circ_expsplit(u, M)
            t === nothing || add!(t...)
        end
        foreach(walk, n.args)
        foreach(walk, n.exps)
    end
    walk(e)
    # exponentials alone are left to the canonical form unless they could combine
    nfn[] == 0 && length(seen) == 0 && return Sym[]
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

# the multiple k of a base that u is, or nothing
_circ_multiple(u::Sym, g::Sym) = (r = expand(u / g); r.kind == NUMBER && isexact(r.num) && isone(_denominator(r.num)) ? Int(_numerator(r.num)) : nothing)

# e with the functions (and exponentials) of multiples of the bases as polynomials and
# quotients in S, C; nothing if an argument is not a multiple of a base
function _circ_poly(e::Sym, bases::Vector{Sym}, vars, M::_Circ)
    ok = Ref(true)
    function find(u)
        for (i, g) in enumerate(bases)
            k = _circ_multiple(u, g)
            k === nothing || return i, k
        end
        nothing
    end
    function expof(u)                       # exp(u), u = k i g (trig) or k g
        t = _circ_expsplit(u, M)
        t === nothing && return nothing
        f = find(expand(Sym(SymNumber(t[1])) * t[2]))
        f === nothing && return nothing
        i, k = f
        S, C = vars[i]
        s, c = _chebyshev(k, S, C)
        M.hyp ? c + s : c + IM * s
    end
    function tr(n::Sym)
        isempty(n.args) && return n
        if n.kind == POW && _eq(n.args[1], E)
            r = expof(n.args[2])
            r === nothing || return r
        end
        if n.kind == FN && n.name in M.fns && length(n.args) == 1
            f = find(n.args[1])
            f === nothing && (ok[] = false; return n)
            i, k = f
            S, C = vars[i]
            s, c = _chebyshev(k, S, C)
            j = findfirst(==(n.name), M.fns)
            return (s, c, s / c, c / s, 1 / c, 1 / s)[j]
        end
        if n.kind == MUL
            out = Sym(n.num)
            for (b, p) in zip(n.args, n.exps)
                r = _eq(b, E) ? expof(p) : nothing
                out *= r === nothing ? tr(b)^tr(p) : r
            end
            return out
        end
        _rebuild_node(n, Sym[tr(a) for a in n.args])
    end
    r = tr(e)
    ok[] ? r : nothing
end

# (f(k g) for f the sine, the cosine) as polynomials in S, C: the recurrence
# f(kg) = 2C f((k-1)g) - f((k-2)g) holds for sin, cos and sinh, cosh alike
function _chebyshev(k::Int, S::Sym, C::Sym)
    k < 0 && ((s, c) = _chebyshev(-k, S, C); return (-s, c))
    k == 0 && return (ZERO, ONE)
    s0, c0, s1, c1 = ZERO, ONE, S, C
    for _ in 2:k
        s0, c0, s1, c1 = s1, c1, expand(2C * s1 - s0), expand(2C * c1 - c0)
    end
    s1, c1
end

# p reduced modulo the relation for each base: eliminating S^2 (S^2 = 1 - C^2, or
# C^2 - 1 hyperbolic) when the order says S, else C^2 (1 - S^2, or 1 + S^2)
function _circ_reduce(p::Sym, vars, order, M::_Circ)
    p = expand(p)
    for ((S, C), elimS) in zip(vars, order)
        has(p, S) || has(p, C) || continue
        P = try Poly{Sym}(p, [S, C]) catch; return p end
        a, b = elimS ? (S, C) : (C, S)
        ia = elimS ? 1 : 2
        sq = M.hyp ? (elimS ? b^2 - 1 : 1 + b^2) : 1 - b^2      # a^2 in terms of b
        acc = ZERO
        for (ev, cf) in zip(exponent_vectors(P), coefficients(P))
            ka, kb = ev[ia], ev[3 - ia]
            acc += cf * a^(ka % 2) * sq^(ka ÷ 2) * b^kb
        end
        p = expand(acc)
    end
    p
end

_ispolyin(p::Sym, vs::Vector{Sym}) = (try Poly{Sym}(expand(p), vs); true catch; false end)

# S, C through z = exp(i g) (exp(g)): the Laurent polynomial in the zᵢ, times Π zᵢ^Dᵢ,
# and the Dᵢ
function _circ_laurent(p::Sym, vars, M::_Circ)
    zs = [Sym(Symbol("_circz", i)) for i in eachindex(vars)]
    p = expand(p)
    D = Int[]
    for (S, C) in vars
        P = Poly{Sym}(p, [S, C])
        push!(D, maximum(sum, exponent_vectors(P); init = 0))
    end
    den = M.hyp ? 2 : 2IM
    sub = [v => w for ((S, C), z) in zip(vars, zs) for (v, w) in ((S, (z^2 - 1) / (den * z)), (C, (z^2 + 1) / (2z)))]
    q = expand(subs(p, sub...) * prod(z^d for (z, d) in zip(zs, D)))
    Poly{Sym}(q, zs), D
end

# a polynomial in Sᵢ, Cᵢ as a combination of the sine and cosine of k·g = Σ kᵢ gᵢ: the
# terms z^k and z^-k pair into cos(k·g) and sin(k·g) (cosh, sinh)
function _circ_linear(p::Sym, vars, bases::Vector{Sym}, M::_Circ)
    Q, D = _circ_laurent(p, vars, M)
    a = Dict{Vector{Int},Sym}()
    for (ev, cf) in zip(exponent_vectors(Q), coefficients(Q))
        a[collect(ev) .- D] = cf
    end
    s, c = M.fns[1], M.fns[2]
    out = get(a, zeros(Int, length(vars)), ZERO)
    done = Set{Vector{Int}}()
    for k0 in keys(a)
        i = findfirst(!iszero, k0)
        i === nothing && continue
        k = k0[i] > 0 ? k0 : -k0                             # each pair once, by its k > 0
        k in done && continue
        push!(done, k)
        ak, am = get(a, k, ZERO), get(a, -k, ZERO)
        θ = expand(sum(kj * g for (kj, g) in zip(k, bases)))
        cs = expand(ak + am)
        sn = M.hyp ? expand(ak - am) : expand(IM * (ak - am))
        out += cs * apply_function(c, Sym[θ]) + sn * apply_function(s, Sym[θ])
    end
    expand(out)
end

# the polynomial in exponentials: z = exp(i g) (exp(g)) substituted back
function _circ_exp(p::Sym, vars, bases::Vector{Sym}, M::_Circ)
    Q, D = _circ_laurent(p, vars, M)
    out = ZERO
    for (ev, cf) in zip(exponent_vectors(Q), coefficients(Q))
        θ = expand(sum((e - d) * g for (e, d, g) in zip(ev, D, bases)))
        out += cf * exp(M.hyp ? θ : IM * θ)
    end
    expand(out)
end

# s(u)^a * c(u)^-a written as the tangent (cotangent for the reverse), and a negative
# power of the cosine (sine) left as a positive power of the secant (cosecant):
# 1/cosh(x)^2 is sech(x)^2. The sizes are the same; this is the usual presentation.
function _circ_recip(b::Sym, p::Sym, M::_Circ)
    if b.kind == FN && length(b.args) == 1 && b.name in (M.fns[1], M.fns[2]) &&
       p.kind == NUMBER && p.num < N0
        return apply_function(b.name === M.fns[2] ? M.fns[5] : M.fns[6], b.args)^Sym(-p.num)
    end
    b^p
end

function _circ_tan(e::Sym, M::_Circ)
    isempty(e.args) && return e
    s, c, t, ct = M.fns[1], M.fns[2], M.fns[3], M.fns[4]
    if e.kind == POW
        return _circ_recip(_circ_tan(e.args[1], M), e.args[2], M)
    end
    if e.kind == MUL
        fs = [(_circ_tan(b, M), p) for (b, p) in zip(e.args, e.exps)]
        out = Sym(e.num)
        used = falses(length(fs))
        for (i, (b, p)) in enumerate(fs)
            used[i] && continue
            if b.kind == FN && b.name in (s, c) && p.kind == NUMBER
                other = b.name === s ? c : s
                j = findfirst(j -> !used[j] && j != i && fs[j][1].kind == FN && fs[j][1].name === other &&
                                   _eq(fs[j][1].args[1], b.args[1]) && fs[j][2].kind == NUMBER &&
                                   sign(fs[j][2].num) == -sign(p.num), eachindex(fs))
                if j !== nothing
                    used[i] = used[j] = true
                    # sin^a cos^b with a, b of opposite signs: tan^m (cot^m) times the rest
                    a, bb = b.name === s ? (p.num, fs[j][2].num) : (fs[j][2].num, p.num)
                    m = min(abs(a), abs(bb))
                    u = b.args
                    out *= apply_function(a > N0 ? t : ct, u)^Sym(m)
                    ra, rb = a - (a > N0 ? m : -m), bb + (a > N0 ? m : -m)
                    iszero(ra) || (out *= _circ_recip(apply_function(s, u), Sym(ra), M))
                    iszero(rb) || (out *= _circ_recip(apply_function(c, u), Sym(rb), M))
                    continue
                end
            end
            used[i] = true
            out *= _circ_recip(b, p, M)
        end
        return out
    end
    _rebuild_node(e, Sym[_circ_tan(a, M) for a in e.args])
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

# force: every base taken as positive (but a number) and every exponent as real
function _log_expand(e::Sym, A::Assumptions, force::Bool = false)
    isempty(e.args) && return e
    pos(b) = b.kind == NUMBER ? b.num > N0 : force || _ispos(b, A)
    realp(p) = force || _isreal(p, A) === true
    if e.kind == FN && e.name === :log && length(e.args) == 1
        u = _log_expand(e.args[1], A, force)
        if u.kind == MUL
            out, rest = ZERO, Sym(u.num)
            _logpos(Sym(u.num), A) && !isone(u.num) ? (out += log(Sym(u.num)); rest = ONE) : nothing
            for (b, p) in zip(u.args, u.exps)
                if pos(b) && realp(p)
                    out += p * _log_expand(log(b), A, force)
                else
                    rest *= b^p
                end
            end
            return out + (isone(rest) ? ZERO : log(rest))
        elseif u.kind == POW && pos(u.args[1]) && realp(u.args[2])
            return u.args[2] * _log_expand(log(u.args[1]), A, force)
        end
        return log(u)
    end
    e.kind == MUL && return Sym(e.num) * prod((_log_expand(b, A, force)^_log_expand(p, A, force) for (b, p) in zip(e.args, e.exps)); init = ONE)
    _rebuild_node(e, Sym[_log_expand(a, A, force) for a in e.args])
end

"""
    SymJ.expand_log(ex; assume = nothing, force = false)

`ex` with its logarithms expanded: `log(a*b) = log(a) + log(b)` and `log(a^p) = p*log(a)`
where the assumptions decide that `a > 0` and `p` is real (a positive number always
comes out: `log(2x) = log(2) + log(x)`); factors that are not decided stay in one
logarithm. With `force = true`, every symbolic base is taken as positive and every
exponent as real, as SymPy's `force` does. The reverse, combining, is part of
[`SymJ.simplify_log`](@ref). Not exported.

```julia
SymJ.expand_log(log(x^2 * y); assume = @assume(x > 0, y > 0))   # 2log(x) + log(y)
SymJ.expand_log(log(x / y); force = true)                       # log(x) - log(y)
SymJ.expand_log(log(x * y))                                     # log(x*y): undecided
```
"""
function expand_log(ex; assume = nothing, force::Bool = false)
    A = assume isa Assumptions ? assume : (a = _assumptions(assume); a === nothing ? Assumptions() : a)
    _log_expand(Sym(ex), A, force)
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
