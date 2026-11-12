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
# a membership, or a combination of facts: decided under the assumptions, else the parts
# refined (refine((n = 0) | (n ≥ 1)) and refine(1/(1 + d) ∈ ℤ) were MethodErrors; SymPy
# #18389, #19145)
function refine(f::Union{FactAnd,FactOr,Membership}; assume = nothing)
    A = _assumptions(assume)
    A === nothing && (A = Assumptions())
    if f isa Membership
        v = _holds(f, A)
        return v === nothing ? f : v
    end
    parts = Any[g isa Relation || g isa Membership || g isa FactAnd || g isa FactOr ? refine(g; assume) : g for g in f.facts]
    if f isa FactOr
        any(p -> p === true, parts) && return true
        rest = Any[p for p in parts if p !== false]
        isempty(rest) && return false
        return length(rest) == 1 ? only(rest) : simplify(FactOr(rest); assume)
    else
        any(p -> p === false, parts) && return false
        rest = Any[p for p in parts if p !== true]
        isempty(rest) && return true
        return length(rest) == 1 ? only(rest) : simplify(FactAnd(rest); assume)
    end
end
function refine(r::Relation; assume = nothing)
    A = _assumptions(assume)
    A === nothing && (A = Assumptions())
    l, rr = _refine(r.lhs, A), _refine(r.rhs, A)
    # an infinite side: compared as it is, not through l - r (∞ - ∞ is NaN, so ∞ ≤ ∞ was
    # false and ∞ = ∞ became NaN = NaN; SymPy #19245)
    if _isinfpt(l) || _isinfpt(rr)
        # a real (so finite) side against ±∞: decided (x < ∞ for real x; SymPy #17113)
        if _isinfpt(l) != _isinfpt(rr)
            f, inf = _isinfpt(rr) ? (l, rr) : (rr, l)
            if !_has_node(f, _isinfpt) && _isreal(f, A) === true
                s = Float64(inf.num) > 0 ? 1 : -1
                fl, fr = _isinfpt(rr) ? (0.0, s * Inf) : (s * Inf, 0.0)   # f stands for 0
                return _REL_INFO[r.op][2](fl, fr)
            end
        end
        (_hassymbol(l) || _hassymbol(rr)) && return Relation(r.op, l, rr)
        a, b = _evalc(l), _evalc(rr)
        (isreal(a) && isreal(b)) || return r.op === :eq ? a == b : r.op === :ne ? a != b : Relation(r.op, l, rr)
        return _REL_INFO[r.op][2](real(a), real(b))
    end
    v = _holds(Relation(r.op, expand(l - rr), ZERO), A)
    v === nothing || return v
    # the relation reduced to a simpler equivalent one: √x < 0 is false, the principal
    # root never negative, with or without assumptions (SymPy #14815)
    g = _reduce_fact(Relation(r.op, l, rr))
    if g isa Relation && !(g.op === r.op && _eq(expand(g.lhs - g.rhs), expand(l - rr)))
        if !_hassymbol(g.lhs) && !_hassymbol(g.rhs)
            a, b = _evalc(g.lhs), _evalc(g.rhs)
            isreal(a) && isreal(b) && return _REL_INFO[g.op][2](real(a), real(b))
        else
            w = _holds(Relation(g.op, expand(g.lhs - g.rhs), ZERO), A)
            w === nothing || return w
        end
    end
    Relation(r.op, l, rr)
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

# exp(iπ w) with w = c₀ + Σ cⱼ tⱼ: a term with tⱼ an integer and cⱼ an even integer drops
# out, one with cⱼ odd is (-1)^tⱼ, and c₀ is reduced mod 2; nothing if nothing changes
function _refine_exp_ipi(u::Sym, A::Assumptions)
    w = expand(u / (IM * PI))
    c0, terms = if w.kind == ADD
        w.num, collect(zip(w.coefs, w.args))
    elseif w.kind == NUMBER
        w.num, Tuple{SymNumber,Sym}[]
    elseif w.kind == MUL
        N0, [(w.num, w / Sym(w.num))]
    else
        N0, [(N1, w)]
    end
    changed = false
    sign = ONE
    rest = ZERO
    for (c, t) in terms
        if isexact(c) && isone(_denominator(c)) && _isinteger(t, A) === true
            changed = true
            isodd(_numerator(c)) && (sign *= NM1^t)
        else
            rest += Sym(c) * t
        end
    end
    if isexact(c0) && (c0 >= SymNumber(2) || c0 < N0)
        changed = true
        c0 = c0 - SymNumber(2) * SymNumber(floor(_numerator(c0) // _denominator(c0) // 2))
    end
    # exp(iπ c₀) = cos(π c₀) + i sin(π c₀) where both are exact (c₀ = 1, 1/2, 1/3, …)
    k = ONE
    if isexact(c0) && !iszero(c0)
        cs, sn = cos(PI * Sym(c0)), sin(PI * Sym(c0))
        if !_has_fn(cs, (:cos,)) && !_has_fn(sn, (:sin,))
            changed = true
            k, c0 = expand(cs + IM * sn), N0
        end
    end
    changed || return nothing
    sign * k * exp(IM * PI * (Sym(c0) + rest))
end

# one identity at the node, or nothing
function _refine_rule(e::Sym, A::Assumptions)
    k = e.kind
    # ±∞ plus real terms is ±∞ (x + ∞ with x real; SymPy #27605); a real term is finite
    if k == ADD && isfloat(e.num) && isinf(Float64(e.num)) && all(t -> _isreal(t, A) === true, e.args)
        return Sym(e.num)
    end
    if k == FN && (haskey(_COND_RELS, e.name) || e.name in (:heaviside, :sign, :abs, :dirac))
        r = haskey(_COND_RELS, e.name) ? _refine_cond(e, A) : _refine_step(e, A)
        r === nothing || return r
    end
    k == FN && e.name === :atan2 && return _refine_atan2(e.args[1], e.args[2], A)
    k == FN && e.name === :And && return _refine_and(e, A)
    # angle(u) with Re u > 0: atan(Im u/Re u) (angle(a + iF) for a > 0, F real; SymPy
    # #15945)
    if k == FN && e.name === :angle && length(e.args) == 1 && _hassymbol(e.args[1])
        u = e.args[1]
        re, im_ = _refine(real(u), A), _refine(imag(u), A)
        if !_has_fn(re, (:real, :imag)) && !_has_fn(im_, (:real, :imag))
            _ispos(re, A) && return atan(cancel(im_ / re))
            # real parts of either sign: the two-argument arctangent (SymPy #11988)
            _isreal(re, A) === true && _isreal(im_, A) === true && return atan(im_, re)
        end
    end
    # mod(u + k m, m) with k m/m integer-valued: mod(u, m) (mod(n y + z, n) for an integer
    # y; SymPy #12592)
    if k == FN && e.name === :mod && length(e.args) == 2 && e.args[1].kind == ADD && _hassymbol(e.args[2])
        u, m = e.args
        drop = [_isinteger(cancel(Sym(c) * t / m), A) === true for (c, t) in zip(u.coefs, u.args)]
        if any(drop)
            return apply_function(:mod, Sym[Sym(u.num) + _add_many(Sym(c) * t for (c, t, d) in zip(u.coefs, u.args, drop) if !d), m])
        end
    end
    # mod(Σ cᵢ tᵢ + c₀, m) for integers m > 0, cᵢ and integer-valued tᵢ: the cᵢ, c₀ mod m
    # (mod(3n, 2) = mod(n, 2) for integer n; SymPy #15493)
    if k == FN && e.name === :mod && length(e.args) == 2
        u, m = e.args
        if m.kind == NUMBER && _isint(m.num) && m.num > N0 && _hassymbol(u)
            terms = u.kind == ADD ? collect(zip(u.coefs, u.args)) : u.kind == MUL ? [(u.num, u / Sym(u.num))] : [(N1, u)]
            c0 = u.kind == ADD ? u.num : N0
            if isexact(c0) && _isint(c0) && all(((c, t),) -> isexact(c) && _isint(c) && _isinteger(t, A) === true, terms)
                M = _toint(m.num)
                v = Sym(mod(_toint(c0), M)) + _add_many(Sym(mod(_toint(c), M)) * t for (c, t) in terms)
                _eq(v, u) || return apply_function(:mod, Sym[v, m])
            end
        end
    end
    # max and min without the arguments the assumptions show dominated (max(x, 0) = x
    # for x > 0; SymPy #16926)
    if k == FN && e.name in (:max, :min) && length(e.args) >= 2
        op = e.name === :max ? :ge : :le
        keep = trues(length(e.args))
        for i in eachindex(e.args), j in eachindex(e.args)
            (i == j || !keep[i] || !keep[j]) && continue
            d = expand(e.args[i] - e.args[j])
            (_isreal(e.args[i], A) === true && _isreal(e.args[j], A) === true) || continue
            _holds(Relation(op, d, ZERO), A) === true && (keep[j] = false)
        end
        all(keep) || return count(keep) == 1 ? e.args[findfirst(keep)] : apply_function(e.name, e.args[keep])
    end
    # floor(u) = n where n ≤ u < n + 1 by the assumptions, ceil alike, for small n
    # (floor(-1/x) = -1 for an integer x > 0; SymPy #21429)
    if k == FN && e.name in (:floor, :ceil) && length(e.args) == 1 && _hassymbol(e.args[1]) && _isreal(e.args[1], A) === true
        u = e.args[1]
        for n in -2:2
            N_ = Sym(n)
            ok = e.name === :floor ? (_holds(Relation(:ge, u, N_), A) === true && _holds(Relation(:lt, u, N_ + 1), A) === true) :
                                     (_holds(Relation(:gt, u, N_ - 1), A) === true && _holds(Relation(:le, u, N_), A) === true)
            ok && return N_
        end
    end
    if k == POW
        b, q = e.args
        if b.kind == NUMBER && b.num == NM1 && !(q.kind == NUMBER)
            _isinteger(expand(q / 2), A) === true && return ONE              # (-1)^(2n)
            _isinteger(expand((q - 1) / 2), A) === true && return -ONE
            return nothing
        end
        _eq(b, E) && return _refine_exp_ipi(q, A)                            # exp(2πi n) = 1
        # 0^q = 0 for Re q > 0 (z^x at z = 0, x > 0, in a definite integral; SymPy #26477)
        iszero(b) && _ispos(q, A) && return ZERO
        # (c y^k)^q = (c (-1)^k)^q (-y)^(kq) for y < 0 and c (-1)^k > 0: (-x^3)^(4/3) = x^4
        # for x < 0 (SymPy #25724)
        if b.kind == MUL && length(b.args) == 1 && b.exps[1].kind == NUMBER && _isint(b.exps[1].num) &&
           q.kind == NUMBER && isexact(q.num) && !_isint(q.num) && _signs(b.args[1], A) == _NEG
            kk = _toint(b.exps[1].num)
            c = Sym(b.num) * (iseven(kk) ? 1 : -1)
            # ((-y)^q itself, for kk = 1 and c = 1, is already in this form)
            c.kind == NUMBER && isexact(c.num) && c.num > N0 && !(kk == 1 && isone(c.num)) &&
                return c^q * (-b.args[1])^(kk * q)
        end
        # b^(k log a / log b) = a^k for a, b > 0 (SymPy #26006)
        if _ispos(b, A) && _has_fn(q, (:log,))
            t = expand(q * log(b))
            if t.kind == MUL
                for (f, p) in zip(t.args, t.exps)
                    (f.kind == FN && f.name === :log && length(f.args) == 1 && _isone(p) && !_eq(f.args[1], b) &&
                     _ispos(f.args[1], A)) || continue
                    kq = cancel(t / f)
                    has(kq, log(b)) && continue
                    return f.args[1]^kq
                end
            end
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
        b.kind == MUL && q.kind == NUMBER && isexact(q.num) && !isone(_denominator(q.num)) &&
            return _refine_powprod(b, q, A)
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
            # |b^q| = |b|^re(q) e^(-im(q) arg b) for a constant b off the positive axis, q's
            # parts known (|i^n| = 1 for real n; SymPy #14260)
            if u.kind == POW && !_hassymbol(u.args[1]) && !_eq(u.args[1], E)
                b, q = u.args
                ab = _angle(b)
                if ab !== nothing && !iszero(ab)
                    p = _reim(q, A)
                    p === nothing || return abs(b)^p[1] * exp(-p[2] * ab)
                end
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
        elseif f in keys(_INVERSE_OF) && u.kind == FN && u.name === _INVERSE_OF[f] && length(u.args) == 1
            return _refine_inverse(f, u.args[1], A)
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

# (∏ uᵢ)^q, q a fraction: (uv)^q = u^q v^q for u ≥ 0, so a real factor to an even power
# comes out as |c| (√(R²(b² + 1)) = |R|√(b² + 1)), and a positive number with it; other
# factors stay, so √(xy) is not split
function _refine_powprod(b::Sym, q::Sym, A::Assumptions)
    out, rest = ONE, Sym[]
    for (c, p) in zip(b.args, b.exps)
        if p.kind == NUMBER && isexact(p.num) && isone(_denominator(p.num)) && iseven(_numerator(p.num)) &&
           _isreal(c, A) === true
            out *= abs(c)^expand(p * q)
        else
            push!(rest, c^p)
        end
    end
    length(rest) == length(b.args) && return nothing
    sign(b.num) > 0 ? out * Sym(b.num)^q * prod(rest; init = ONE)^q : out * (Sym(b.num) * prod(rest; init = ONE))^q
end

# the radicals of e, under A, with their real factors to even powers brought out (_refine_powprod)
function _split_radicals(e::Sym, A::Assumptions)
    isempty(e.args) && return e
    e.kind == MUL &&
        return Sym(e.num) * prod((_split_pow(_split_radicals(b, A), p, A) for (b, p) in zip(e.args, e.exps)); init = ONE)
    n = _rebuild_node(e, Sym[_split_radicals(a, A) for a in e.args])
    n.kind == POW ? _split_pow(n.args[1], n.args[2], A) : n
end

# b^q, q a fraction, with b's real factors to even powers out: √(x²) = |x|, √(4x²c) = 2|x|√c
function _split_pow(b::Sym, q::Sym, A::Assumptions)
    (q.kind == NUMBER && isexact(q.num) && !isone(_denominator(q.num))) || return b^q
    if b.kind == POW
        c, p = b.args
        p.kind == NUMBER && isexact(p.num) && isone(_denominator(p.num)) && iseven(_numerator(p.num)) &&
            _isreal(c, A) === true && return abs(c)^expand(p * q)
    elseif b.kind == MUL
        r = _refine_powprod(b, q, A)
        r === nothing || return r
    end
    b^q
end

# atan(y, x) as atan(y/x) moved to the quadrant the signs of x and y give
function _refine_atan2(y::Sym, x::Sym, A::Assumptions)
    my, mx = _signs(y, A), _signs(x, A)
    (my === nothing || mx === nothing) && return nothing          # both must be real
    mx == _POS && return _atan2_value(y, x, 0, 1)
    if mx == _NEG
        my & _NEG == 0 && return _atan2_value(y, x, 1, -1)
        my == _NEG && return _atan2_value(y, x, -1, -1)
    elseif mx == _ZER
        my == _POS && return _atan2_value(y, x, 1, 0)
        my == _NEG && return _atan2_value(y, x, -1, 0)
        my == _ZER && return ZERO
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
        # 1/w^n = conj(w)^n/|w|^(2n): x²/(z + i zR) with real symbols (SymPy #26102)
        if e.kind == NUMBER && isexact(e.num) && isone(_denominator(e.num)) && -SymNumber(8) <= e.num < N0
            q = _reim(b^Sym(-e.num), A)
            q === nothing && return nothing
            iszero(q[2]) && return nothing
            m = expand(q[1]^2 + q[2]^2)
            return (q[1] / m, -q[2] / m)
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
        # acos, asin of a real number outside [-1, 1] (Julia's values, approached from
        # above the real axis): acos(2) = -i acosh(2), asin(2) = π/2 + i acosh(2)
        if f in (:acos, :asin) && iszero(b)
            if _holds(Relation(:ge, a, ONE), A) === true
                return f === :acos ? (ZERO, -acosh(a)) : (PI / 2, acosh(a))
            elseif _holds(Relation(:le, a, -ONE), A) === true
                return f === :acos ? (PI, -acosh(-a)) : (-PI / 2, acosh(-a))
            elseif _holds(Relation(:ge, a, -ONE), A) === true && _holds(Relation(:le, a, ONE), A) === true
                return (u, ZERO)
            end
        end
    end
    nothing
end

# f(g(v)) = v for g the function f inverts, where v is in f's range (SymPy #28376):
# atanh(tanh(v)), asinh(sinh(v)) for real v; asin(sin(v)) for -π/2 ≤ v ≤ π/2, acos(cos(v))
# for 0 ≤ v ≤ π, atan(tan(v)) for -π/2 < v < π/2, acosh(cosh(v)) for v ≥ 0
const _INVERSE_OF = Dict(:atanh => :tanh, :asinh => :sinh, :asin => :sin, :acos => :cos, :atan => :tan,
                         :acosh => :cosh)
const _INVERSE_RANGE = Dict(:asin => (-PI / 2, PI / 2, true), :acos => (ZERO, PI, true),
                            :atan => (-PI / 2, PI / 2, false), :acosh => (ZERO, nothing, true))
function _refine_inverse(f::Symbol, v::Sym, A::Assumptions)
    _isreal(v, A) === true || return nothing
    r = get(_INVERSE_RANGE, f, nothing)
    r === nothing && return v                      # atanh, asinh: all real v
    lo, hi, closed = r
    _holds(Relation(closed ? :ge : :gt, v, lo), A) === true || return nothing
    hi === nothing || _holds(Relation(closed ? :le : :lt, v, hi), A) === true || return nothing
    v
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
    e = _refine(Sym(ex), A)
    r = _simplify(e, A)
    # once more when the first pass changed it: a candidate of the result may reach
    # further (2cos(x)sin(x)sin(π cos x) is sin(2x)sin(π cos x); SymPy #18231)
    if !_eq(r, e) && _hassymbol(r) && _nodes(r) <= 40
        r2 = _simplify(r, A)
        _nodes(r2) < _nodes(r) && (r = r2)
    end
    r
end
# a conjunction of facts without those the others imply (x ≪ 1 & x ≪ 2 is x ≪ 1), a
# disjunction without those that imply another; one fact left is that fact (SymPy #23391)
function simplify(f::Union{FactAnd,FactOr}; assume = nothing)
    fs = Any[g isa Relation ? simplify(g; assume) : g for g in f.facts]
    any(g -> g === (f isa FactAnd ? false : true), fs) && return f isa FactAnd ? false : true
    fs = [g for g in fs if g !== (f isa FactAnd ? true : false)]
    keep = trues(length(fs))
    for i in eachindex(fs)
        others = [fs[j] for j in eachindex(fs) if j != i && keep[j]]
        isempty(others) && continue
        implied = if f isa FactAnd
            B = try Assumptions(others...) catch err; err isa ArgumentError || rethrow(); nothing end
            B !== nothing && _holds(fs[i], B) === true
        else
            B = try Assumptions(fs[i]) catch err; err isa ArgumentError || rethrow(); nothing end
            B !== nothing && any(h -> _holds(h, B) === true, others)
        end
        implied && (keep[i] = false)
    end
    fs = fs[keep]
    isempty(fs) && return f isa FactAnd
    length(fs) == 1 ? fs[1] : (f isa FactAnd ? FactAnd(fs) : FactOr(fs))
end

# each case's value (an integral by cases; SymPy #23909 extra)
simplify(cs::Cases; assume = nothing) =
    _map_cases(v -> v isa AbstractArray ? map(w -> simplify(w; assume), v) : simplify(v; assume), cs)
function simplify(r::Relation; assume = nothing)
    v = refine(r; assume)
    v isa Bool && return v
    l, rr = _drop_common_terms(simplify(v.lhs; assume), simplify(v.rhs; assume))
    A = _assumptions(assume)
    A === nothing && (A = Assumptions())
    d = l - rr
    # sides that simplify to the same constant: decided (cos² + sin² ~ 1 was 0 ~ 0;
    # SymPy #12522)
    if !_hassymbol(d)
        w = _holds(Relation(v.op, expand(d), ZERO), A)
        w === nothing || return w
    end
    if v.op === :eq
        # a + i b = 0 with a, b real: a = 0 and b = 0 (X = iY; SymPy #27759)
        parts = _real_imag_split(d, A)
        parts === nothing || return _fact_and([simplify(Relation(:eq, p, ZERO); assume) for p in parts])
        # terms of one sign summing to 0: each is 0 (X = -Y with X, Y ≥ 0; SymPy #27755)
        ts = _signed_terms(d, A)
        ts === nothing || return _fact_and([simplify(Relation(:eq, t, ZERO); assume) for t in ts])
    end
    # an equation's factors known not to vanish divided out: a(b - c) ~ 0 with a ≠ 0 is
    # b ~ c (SymPy #21108)
    if v.op in (:eq, :ne) && _hassymbol(d)
        f = try factor(simplify(d; assume)) catch err; err isa ArgumentError || rethrow(); d end
        if f.kind == MUL && length(f.args) >= 2
            keep = [!(_isnonzero(b, A) === true) for b in f.args]
            if !all(keep) && any(keep)
                rest = prod((b^p for (b, p, k) in zip(f.args, f.exps, keep) if k); init = ONE)
                return simplify(Relation(v.op, rest, ZERO); assume)
            end
        end
    end
    # one side, without numeric content, when shorter: -x ~ x is x ~ 0, y ≥ -y is y ≥ 0
    # (SymPy #27754, #27623)
    one = _one_sided(v.op, simplify(d; assume))
    one !== nothing && _nodes(one.lhs) + _nodes(one.rhs) < _nodes(l) + _nodes(rr) && return one
    Relation(v.op, l, rr)
end

_fact_and(fs) = any(f -> f === false, fs) ? false : (gs = [f for f in fs if f !== true];
    isempty(gs) ? true : length(gs) == 1 ? gs[1] : reduce(&, gs))

const _FLIP = Dict(:lt => :gt, :gt => :lt, :le => :ge, :ge => :le, :eq => :eq, :ne => :ne)

# d op 0 with d's numeric content divided out (a negative one flipping an inequality)
function _one_sided(op::Symbol, d::Sym)
    _hassymbol(d) || return nothing
    c = d.kind == MUL ? d.num : nothing
    if d.kind == ADD
        q = _content_free(d)
        r = cancel(d / q)
        r.kind == NUMBER || return Relation(op, d, ZERO)
        c = r.num
        d = q
    elseif c !== nothing && isexact(c)
        d = d / Sym(c)
    else
        return Relation(op, d, ZERO)
    end
    Relation(sign(c) < 0 ? _FLIP[op] : op, d, ZERO)
end

# (a, b) with d = a + i b, a and b real under A and b ≠ 0 structurally; else nothing
function _real_imag_split(d::Sym, A)
    _isreal(d, A) === true && return nothing
    a, b = try
        _refine(real(d), A), _refine(imag(d), A)
    catch err
        err isa InterruptException && rethrow()
        return nothing
    end
    (iszero(b) || _hasfn(a, :real) || _hasfn(a, :imag) || _hasfn(b, :real) || _hasfn(b, :imag)) && return nothing
    (_isreal(a, A) === true && _isreal(b, A) === true) || return nothing
    (a, b)
end

# the terms of a sum d when all are ≥ 0 under A, or all ≤ 0 (at least one with a symbol)
function _signed_terms(d::Sym, A)
    d.kind == ADD && iszero(d.num) || return nothing
    ts = [Sym(c) * t for (c, t) in zip(d.coefs, d.args)]
    ms = [_smask(t, A) for t in ts]
    any(isnothing, ms) && return nothing
    (all(m -> m & _NEG == 0, ms) || all(m -> m & _POS == 0, ms)) || return nothing
    ts
end

# the terms on both sides of a relation taken away: x² + 1 ~ y² + 1 is x² ~ y² (SymPy #28666)
function _drop_common_terms(l::Sym, r::Sym)
    terms(e) = e.kind == ADD ? [Sym(c) * t for (c, t) in zip(e.coefs, e.args)] : (e.kind == NUMBER ? Sym[] : Sym[e])
    konst(e) = e.kind == ADD ? Sym(e.num) : e.kind == NUMBER ? e : ZERO
    lt, rt = terms(l), terms(r)
    common = [t for t in lt if any(u -> _eq(u, t), rt)]
    c = konst(l)
    k = (!iszero(c) && _eq(c, konst(r))) ? c : ZERO
    (isempty(common) && iszero(k)) && return l, r
    d = _add_many(common) + k
    l - d, r - d
end

# float exponents with an integer value made exact: x^2.0 is x^2 (it then merges with
# x^2; SymPy #27706, #27796)
function _exact_int_exponents(e::Sym)
    _replace_walk(e) do t
        if t.kind == POW && t.args[2].kind == NUMBER && isfloat(t.args[2].num)
            v = Float64(t.args[2].num)
            isinteger(v) && abs(v) < 2^53 && return _exact_int_exponents(t.args[1])^Sym(Int(v))
        elseif t.kind == MUL && any(p -> p.kind == NUMBER && isfloat(p.num) && isinteger(Float64(p.num)), t.exps)
            return Sym(t.num) * prod((_exact_int_exponents(b)^(p.kind == NUMBER && isfloat(p.num) && isinteger(Float64(p.num)) ?
                                     Sym(Int(Float64(p.num))) : p) for (b, p) in zip(t.args, t.exps)); init = ONE)
        end
        nothing
    end
end

function _simplify(e::Sym, A::Assumptions)
    e = _exact_int_exponents(e)
    if !_hassymbol(e) && e.kind != NUMBER
        # function arguments first, as for expressions with symbols: acos(cos a cos b +
        # sin a sin b) has cos(a - b) inside (SymPy #15129)
        e = _map_fn_args(a -> _simplify(a, A), e)
        # a constant in radicals that is 0 though no form cancels (it numerically 0 at 512
        # bits; SymPy #20427)
        w = try _evalc(e) catch; nothing end
        w isa Number && abs(w) < 1e-10 && _const_zero(e) && return ZERO
        # until it stops shrinking (a product of square roots that a step leaves, √10√2√5,
        # merges in the next; found in the review of SymPy #24948)
        for _ in 1:3
            c = _simplify_const(e, A)
            (c.kind == NUMBER || _nodes(c) >= _nodes(e)) && return c
            e = c
        end
        return e
    end
    # function arguments first
    e = _map_fn_args(a -> _simplify(a, A), e)
    best, size = e, _nodes(e)
    powers(x) = simplify_powers(x; assume = A)
    logs(x) = simplify_log(x; assume = A)
    radicals(x) = simplify_radicals(x; assume = A)
    # an expression whose expansion would have millions of terms (powers of sums with
    # radicals, as antiderivatives through Euler substitutions have; SymPy #23967): not
    # expanded, nor put over a common denominator
    huge = _expansion_estimate(e) > _EXPANSION_MAX
    expansion(x) = _bigpower(x) || huge ? nothing : expand(x)   # (x + y)^30 expanded never wins
    simplify_rational(x) = huge ? nothing : SymJ.simplify_rational(x)
    function consider(g)
        c = try g(e) catch err; err isa InterruptException && rethrow(); nothing end
        c isa Sym || return
        c = _refine(c, A)
        n = _nodes(c)
        n < size && ((best, size) = (c, n))
    end
    logq(x) = _hasfn2log(x) ? simplify_log(cancel(_log_quotients(x)); assume = A) : nothing
    elog(x) = (y = _exp_log_powers(x); y === x ? nothing : y)
    # piecewise functions folded under the assumptions: shadowed or impossible pieces gone
    # (SymPy #22148)
    function pwfold(x)
        _hasfn(x, :piecewise) || return nothing
        vs = Set{Sym}()
        _pw_condsyms!(vs, x)
        local folded = nothing   # (not `best`, which this closure would overwrite)
        for v in _sort_by_name!(collect(vs))
            r = try piecewise_fold(x, v; assume = A) catch err; err isa InterruptException && rethrow(); nothing end
            r isa Sym && (folded === nothing || _nodes(r) < _nodes(folded)) && (folded = r)
        end
        folded
    end
    for g in (simplify_rational, powers, expansion, _factor_or_nothing, simplify_trig, simplify_hyperbolic, logs, radicals, simplify_gamma, logq,
              elog, simplify_rational ∘ (x -> _exp_log_powers(x)), pwfold,
              x -> huge ? nothing : _conjugate_radical_den(x),
              simplify_rational ∘ powers, expansion ∘ simplify_rational, expansion ∘ radicals,
              simplify_rational ∘ radicals,
              x -> (c = expansion(x); c === nothing ? nothing : simplify_rational(c)),
              x -> (c = expansion(x); c === nothing ? nothing : _pull_common_den(c)))
        consider(g)
    end
    # a sum's circular-function terms simplified apart from the rest (sin² + cos² + (x + y)³:
    # the whole is never smaller expanded; SymPy #27356)
    if e.kind == ADD && _has_fn(e, _CIRCFNS)
        ts = [Sym(c) * t for (c, t) in zip(e.coefs, e.args)]
        trig = [t for t in ts if _has_fn(t, _CIRCFNS)]
        if 1 < length(trig) < length(ts)
            rest = Sym(e.num) + _add_many(t for t in ts if !_has_fn(t, _CIRCFNS))
            consider(_ -> _simplify(_add_many(trig), A) + rest)
        end
    end
    # a product's circular factors of one argument simplified apart from the others
    # (2sin(a/2)cos(a/2) cos(b)cos(c) is sin(a) cos(b)cos(c); simplify_trig of the whole
    # turns it into sums; SymPy #26835)
    if e.kind == MUL && _has_fn(e, _CIRCFNS)
        fs = [(b, p) for (b, p) in zip(e.args, e.exps)]
        groups = Dict{Sym,Vector{Int}}()
        for (i, (b, _)) in enumerate(fs)
            b.kind == FN && b.name in _CIRCFNS && length(b.args) == 1 && push!(get!(groups, b.args[1], Int[]), i)
        end
        for idx in values(groups)
            (length(idx) >= 2 && length(idx) < length(fs)) || continue
            part = Sym(e.num) * prod((fs[i][1]^fs[i][2] for i in idx); init = ONE)
            rest = prod((fs[i][1]^fs[i][2] for i in eachindex(fs) if !(i in idx)); init = ONE)
            consider(_ -> simplify_trig(part) * rest)
        end
    end
    # several piecewise functions (as from a definite integral with a parameter): one
    _count_fn(e, :piecewise) >= 2 && consider(x -> piecewise_fold(x; assume = A))
    _has_fn(e, (:real,)) && _has_fn(e, (:imag,)) && consider(_reim_pairs)
    # hyperbolic functions with exponentials of a compound argument: all as exponentials
    # (e^(2w)/4 + e^(-2w)/4 - sinh(w)² - 1/2 is 0 for w = (x + b)/u; found in the review)
    _has_fn(e, _HYPMODE.fns) && _has_epow(e) &&
        consider(x -> expand(rewrite(x, (sinh, cosh, tanh, coth, sech, csch) => exp)))
    # max and min: a common positive factor, a common term, and a common sign out
    # (max(2a, 2b) = 2max(a, b), max(a + c, b + c) = c + max(a, b), max(-a, -b) =
    # -min(a, b); SymPy #24578)
    _has_fn(e, (:max, :min)) && consider(_extreme_common)
    # u(x) δ(x - a) = u(a) δ(x - a): x δ(x) = 0 (SymPy #13996)
    _has_fn(e, (:dirac,)) && consider(_dirac_sift)
    # erfc(u) = 1 - erf(u): erf(x) + erfc(x) = 1 (SymPy #28900)
    _has_fn(e, (:erfc,)) && consider(_erfc_to_erf)
    # conj(w) w = abs(w)^2
    _has_fn(e, (:conj,)) && consider(_conj_pairs_to_abs)
    # |u| sign(u) = u; |u|^m u^(2k) = |u|^(m + 2k) for real u (SymPy #20402, #19484, #20108)
    _has_fn(e, (:abs,)) && consider(x -> _abs_sign_pairs(x, A))
    # sin(iu) = i sinh(u), cos(iu) = cosh(u), tan(iu) = i tanh(u)
    _has_fn(e, (:sin, :cos, :tan)) && consider(x -> _replace_walk(x) do t
        (t.kind == FN && t.name in (:sin, :cos, :tan) && length(t.args) == 1) || return nothing
        u = _im_multiple(t.args[1])
        u === nothing ? nothing : t.name === :sin ? IM * sinh(u) : t.name === :cos ? cosh(u) : IM * tanh(u)
    end)
    # a power of a constant with a symbolic exponent beside exponentials: all as exp
    # (iⁿ - exp(iπn/2) = 0; SymPy #16797)
    _has_epow(e) && _has_node(e, t -> t.kind == POW && !_hassymbol(t.args[1]) && t.args[1] !== E && _hassymbol(t.args[2])) &&
        consider(x -> rewrite(x, exp))
    # positive factors out of radicals, sums factored first, so they can cancel:
    # √(ab)/√b = √a for b > 0, √(1 - x²)/√(1 - x) = √(1 + x) for -1 < x < 1 (SymPy #17634)
    !isempty(A.facts) && _has_node(e, t -> !isempty(_radical_bases(t))) && consider(x -> _split_positive_radicals(x, A))
    # the integer content of a product's sum factors in its coefficient: -√2(2μ - 2x)/4 is
    # √2(x - μ)/2 (SymPy #17233)
    e.kind == MUL && any(b -> b.kind == ADD, e.args) && consider(_sum_contents_out)
    # rational powers of integers in a product merged by their primes: 2^(4/5) 4^(3/5) is
    # 2^(4/5 + 6/5) = 4 (SymPy #15669)
    _has_node(e, t -> t.kind == MUL && count(i -> _int_rational_power(t.args[i], t.exps[i]), eachindex(t.args)) >= 2) &&
        consider(_merge_prime_powers)
    # |a||b| = |ab| (|t/x||x| = |t|; SymPy #18840)
    _count_fn(e, :abs) >= 2 && consider(_merge_abs)
    # conj(u) = re u - i im u, beside real and imag parts (SymPy #18271)
    _has_fn(e, (:conj,)) && _has_fn(e, (:real, :imag)) &&
        consider(x -> _replace_walk(x) do t
            t.kind == FN && t.name === :conj && length(t.args) == 1 ?
                apply_function(:real, t.args) - IM * apply_function(:imag, t.args) : nothing
        end)
    # a piece under s = c (s a symbol) has s = c in its value: piecewise(y = 0 => xy, …)
    # is piecewise(y = 0 => 0, …) (SymPy #16690)
    _has_fn(e, (:piecewise,)) && consider(_pw_eq_subs)
    # piecewise(u < 0 => -u, u ≥ 0 => u) is |u| (SymPy #19167)
    _has_fn(e, (:piecewise,)) && consider(_pw_as_abs)
    # Σ_a^b f - Σ_a^c f with b - c a small integer: the terms between (SymPy #19685)
    _has_fn(e, (:Sum,)) && consider(_combine_sums)
    # fractional powers of one base as powers of t = b^(1/q), cancelled as a rational
    # function of t, when smaller: (x - 1)/(√x - 1) = √x + 1 (cancel keeps such powers
    # independent; SymPy #19769)
    consider(_cancel_in_roots)
    # under a root, 1 ± cos(u) = 2cos²(u/2), 2sin²(u/2): √((1 + cos x)/2) = cos(x/2) for 0 < x < π
    # W(u) exp(W(u)) = u, on every branch
    _has_fn(e, (:lambertw,)) && consider(_lambertw_pairs)
    if _has_fn(e, (:cos,))
        # only when refine then takes a root (fewer fractional powers): √(cos²(x/2)) alone is no gain
        function halfangle(x)
            c = _refine(_half_angle_radicands(x), A)
            _root_of_trig_square(c) ? nothing : c
        end
        consider(halfangle)
        consider(x -> (c = halfangle(x); c === nothing ? nothing : simplify_trig(c)))   # 2cos(x/2)sin(x/2) = sin(x)
    end
    if _has_fn(e, _CIRCFNS)       # identities simplify_trig's normal form does not reach
        for g in (_simp_trig_rules, _expand_trig_small, simplify_rational ∘ _expand_trig_small)
            consider(g)
        end
    end
    best
end
# in each product, W(u)^k exp(k W(u) + r) as u^k exp(r), k an integer (W any branch of
# lambertw): W e^W = u
function _lambertw_pairs(e::Sym)
    _replace_walk(e) do t
        t.kind == MUL || return nothing
        bs, ps = collect(t.args), collect(t.exps)
        changed = false
        for i in eachindex(bs)
            w, k = bs[i], ps[i]
            (w.kind == FN && w.name === :lambertw && k.kind == NUMBER && isexact(k.num) && isone(_denominator(k.num))) || continue
            j = findfirst(j -> j != i && _eq(bs[j], E), eachindex(bs))
            j === nothing && continue
            ex = ps[j]
            rest = expand(ex - k * w)
            has(rest, w) && continue                        # the coefficient of W is not k
            bs[i], ps[i] = w.args[1], k
            ps[j] = rest
            changed = true
        end
        changed ? Sym(t.num) * prod(b^p for (b, p) in zip(bs, ps)) : nothing
    end
end

# whether e has a non-integer power of an expression with cos² or sin² in it (a root
# refine could not take)
function _root_of_trig_square(e::Sym)
    sq(b) = _has_node(b, t -> t.kind == POW && t.args[2] == Sym(2) && (_isfn1(t.args[1], :cos) || _isfn1(t.args[1], :sin))) ||
            (b.kind == MUL && any(i -> b.exps[i] == Sym(2) && (_isfn1(b.args[i], :cos) || _isfn1(b.args[i], :sin)), eachindex(b.args)))
    frac(q) = !(q.kind == NUMBER && isone(_denominator(q.num)))
    e.kind == POW && frac(e.args[2]) && sq(e.args[1]) && return true
    e.kind == MUL && any(i -> frac(e.exps[i]) && sq(e.args[i]), eachindex(e.args)) && return true
    any(_root_of_trig_square, e.args)
end
_has_node(e::Sym, p) = p(e) || any(a -> _has_node(a, p), e.args)

# n + n cos(u) as 2n cos²(u/2) and n - n cos(u) as 2n sin²(u/2), in the bases of non-integer
# powers only (where refine can then take the root)
function _half_angle_radicands(e::Sym)
    half(t::Sym) = _replace_walk(t) do s
        (s.kind == ADD && length(s.args) == 1 && _isfn1(s.args[1], :cos)) || return nothing
        n, c, u = s.num, s.coefs[1], s.args[1].args[1]
        iszero(n) && return nothing
        c == n && return 2Sym(n) * cos(u / 2)^2
        c == -n && return 2Sym(n) * sin(u / 2)^2
        nothing
    end
    _replace_walk(e) do t
        if t.kind == POW
            b, q = t.args
            (q.kind == NUMBER && isone(_denominator(q.num))) && return nothing
            nb = half(b)
            return _eq(nb, b) ? nothing : nb^q
        elseif t.kind == MUL
            changed = false
            fs = Sym[]
            for (b, q) in zip(t.args, t.exps)
                if !(q.kind == NUMBER && isone(_denominator(q.num)))
                    nb = half(b)
                    changed |= !_eq(nb, b)
                    push!(fs, nb^q)
                else
                    push!(fs, b^q)
                end
            end
            return changed ? Sym(t.num) * prod(fs) : nothing
        end
        nothing
    end
end

_hasfn2log(e::Sym) = (e.kind == FN && e.name === :log && length(e.args) == 2) || any(_hasfn2log, e.args) ||
    (e.kind == MUL && any(_hasfn2log, e.exps))

_has_epow(e::Sym) = (e.kind == POW && e.args[1] === E) || (e.kind == MUL && any(b -> b === E, e.args)) ||
    any(_has_epow, e.args) || (e.kind == MUL && any(_has_epow, e.exps))

# a sum with a denominator factor common to every term: (Σ tᵢ D)/D (the exponentials
# over k01 - k11 of SymPy #26329, rather than their cosh/sinh form)
function _pull_common_den(e::Sym)
    e.kind == ADD && iszero(e.num) || return nothing
    ts = [Sym(c) * a for (c, a) in zip(e.coefs, e.args)]
    dens(t) = t.kind == MUL ? Dict(b => -p for (b, p) in zip(t.args, t.exps) if p.kind == NUMBER && p.num < N0 && b.kind == ADD) :
              t.kind == POW && t.args[2].kind == NUMBER && t.args[2].num < N0 && t.args[1].kind == ADD ? Dict(t.args[1] => -t.args[2]) :
              Dict{Sym,Sym}()
    common = dens(ts[1])
    for t in ts[2:end]
        d = dens(t)
        common = Dict(b => (d[b].num < k.num ? d[b] : k) for (b, k) in common if haskey(d, b))
        isempty(common) && return nothing
    end
    D = prod((b^k for (b, k) in common); init = ONE)
    _add_many(t * D for t in ts) / D
end

# c real(u) + c i imag(u) = c u (SymPy #25564)
_reim_pairs(e::Sym) = _replace_walk(e) do t
    t.kind == ADD || return nothing
    cs, as = Sym[Sym(c) for c in t.coefs], collect(t.args)
    ts = [c * a for (c, a) in zip(cs, as)]
    used, out, changed = falses(length(as)), Sym(t.num), false
    for i in eachindex(as)
        used[i] && continue
        a = as[i]
        if a.kind == FN && a.name === :real && length(a.args) == 1
            target = cs[i] * IM * apply_function(:imag, Sym[a.args[1]])     # (i is (-1)^(1/2))
            j = findfirst(j -> !used[j] && j != i && _eq(ts[j], target), eachindex(as))
            if j !== nothing
                used[i] = used[j] = true
                out += cs[i] * a.args[1]
                changed = true
                continue
            end
        end
        used[i] = true
        out += ts[i]
    end
    changed ? out : nothing
end

# i^q, (-i)^q, (-1)^q for rational q as cos(θ) + i sin(θ) where exact: i^(3/2) + √(-i)
# is 0 (SymPy #26140)
function _unit_powers(e::Sym)
    function unit(b::Sym, q::Sym)
        (q.kind == NUMBER && isexact(q.num) && !_isint(q.num)) || return nothing
        f = _eq(b, IM) ? 1 // 2 : _eq(b, -IM) ? -1 // 2 : _eq(b, MINUS_ONE) ? 1 // 1 : nothing
        f === nothing && return nothing
        θ = Rational{BigInt}(f) * _rat(q.num)
        sn, cs = _sinpi(θ), _sinpi(θ + 1 // 2)
        (sn === nothing || cs === nothing) ? nothing : cs + IM * sn
    end
    _replace_walk(e) do t
        if t.kind == POW
            return unit(t.args[1], t.args[2])
        elseif t.kind == MUL && any(b -> b.kind == NUMBER, t.args)
            ch = false
            fs = map(zip(t.args, t.exps)) do (b, p)
                u = unit(b, p)
                u === nothing ? b^p : (ch = true; u)
            end
            return ch ? Sym(t.num) * prod(fs; init = ONE) : nothing
        end
        nothing
    end
end

# log(b, a) as log(a)/log(b), so that logarithms to one base combine
# (and the logarithm of an integer as its primes', so that log(10) = log(2) + log(5) relates
# the bases: log(10, 5)² + log(10, 5) log(10, 2) + log(10, 2) is 1; SymPy #25860)
function _log_quotients(x::Sym)
    lg(u) = (u.kind == NUMBER && isexact(u.num) && _isint(u.num) && u.num > N1) ? _log_primes(_toint(u.num), u) : log(u)
    _replace_walk(x) do t
        t.kind == FN && t.name === :log && length(t.args) == 2 ? lg(t.args[2]) / lg(t.args[1]) : nothing
    end
end

# log(n) as Σ k log(p) over n's prime factors (trial division below 10^6), else log(n)
function _log_primes(n::Integer, u::Sym)
    out, m, p = ZERO, n, 2
    while p * p <= m && p < 10^6
        k = 0
        while m % p == 0
            m ÷= p
            k += 1
        end
        k > 0 && (out += k * log(Sym(p)))
        p += p == 2 ? 1 : 2
    end
    if m > 1                                          # a prime left, or an unfactored part
        p * p > m || return log(u)
        out += log(Sym(m))
    end
    out
end

# a constant: powers, logarithms and radicals of numbers
_erfc_to_erf(x::Sym) = _replace_walk(x) do t
    t.kind == FN && t.name === :erfc && length(t.args) == 1 ? 1 - apply_function(:erf, Sym[t.args[1]]) : nothing
end

function _simplify_const(e::Sym, A::Assumptions)
    # large powers of sums of radicals expanded first: expand reduces √2² = 2 as it goes,
    # the rest (cancel over √2 as a generator) would not, at degree 2012 (SymPy #25109)
    if _bigpower(e)
        x = try expand(e) catch err; err isa InterruptException && rethrow(); e end
        _nodes(x) < _nodes(e) && (e = x)
    end
    # √(a + bi) with a² + b² a square and rational parts: √(3 + 4i) = 2 + i (SymPy #20024)
    c = _complex_sqrt_denest(e)
    _nodes(c) < _nodes(e) && (e = c)
    # an algebraic constant that is rational: its minimal polynomial of degree 1
    # (-48√(4 + 2√5)/(1 + √5)^(3/2) is -24; SymPy #18475)
    r = _rational_by_minpoly(e)
    r === nothing || return r
    best = e
    # (trigonometric and hyperbolic too: cos(2)² + sin(2)² is 1, SymPy #25705; and over
    # a common denominator, log(5)/log(10) + log(2)/log(10) is 1, #25860)
    for g in (simplify_powers, simplify_log, simplify_radicals, (x; assume) -> simplify_gamma(x),
              (x; assume) -> _erfc_to_erf(x), (x; assume) -> simplify_trig(x),
              (x; assume) -> simplify_hyperbolic(x), (x; assume) -> expand_trig(x),
              (x; assume) -> _unit_powers(x), (x; assume) -> cancel(simplify_radicals(x; assume)),
              (x; assume) -> simplify_log(cancel(_log_quotients(x)); assume),
              (x; assume) -> _trig_unity_zero(x), (x; assume) -> _merge_prime_powers(x),
              # inverse hyperbolic functions as logs, the logs combined (log(3/2) -
              # 2atanh(2) + 2atanh(3) is 0; SymPy #16730)
              (x; assume) -> _has_fn(x, (:atanh, :acoth, :asinh, :acosh)) ?
                  simplify_log(cancel(_log_quotients(expand(rewrite(x, log)))); assume) : x)
        c = try g(e; assume = A) catch err; err isa InterruptException && rethrow(); e end
        c = try expand(c) catch; c end
        _nodes(c) < _nodes(best) && (best = c)
    end
    _has_fn(best, _INVTRIG_FNS) || return best
    r = _invtrig_multiple_of_pi(best, A)
    r === nothing ? best : r
end

const _INVTRIG_FNS = (:acos, :asin, :atan, :acot, :asec, :acsc)

# 0 for a constant in sines and cosines of rational multiples of π that is 0: through
# exponentials, roots of unity reduced modulo the cyclotomic polynomial, a normal form
# (-32cos(π/14)cos(2π/7) - 32sin(3π/14)cos(3π/14) + … is 0; SymPy #23015); else e
function _trig_unity_zero(e::Sym)
    _has_fn(e, (:sin, :cos)) || return e
    w = try _unity_reduce(expand(rewrite(e, exp))) catch err; err isa InterruptException && rethrow(); return e end
    iszero(w) ? ZERO : e
end

# a constant with inverse trig functions whose value is q π, q rational with a small
# denominator: e - qπ is a multiple of 2π when cos(e - qπ) is 1 exactly, and so 0, its
# value being tiny (acos(2√5/5) - acos(-3√10/10) + 3π/4 is 0; SymPy #24948)
const _IN_INVTRIG_PI = ScopedValue(false)
function _invtrig_multiple_of_pi(e::Sym, A::Assumptions)
    _IN_INVTRIG_PI[] && return nothing
    v = try N(e) catch; return nothing end
    v isa Real || return nothing
    q = rationalize(v / π; tol = 1e-12)
    (denominator(q) <= 12 && abs(v - q * π) < 1e-10) || return nothing
    qπ = Sym(q) * PI
    c = with(_IN_INVTRIG_PI => true) do
        try simplify(simplify(cos(e - qπ))) catch; nothing end
    end
    isone(c) ? qπ : nothing
end

function _extreme_common(e::Sym)
    _replace_walk(e) do t
        (t.kind == FN && t.name in (:max, :min) && length(t.args) >= 2) || return nothing
        args = t.args
        other = t.name === :max ? :min : :max
        # a common term (same coefficient) in sums
        if all(a -> a.kind == ADD, args)
            common = [(c, u) for (c, u) in zip(args[1].coefs, args[1].args)
                      if all(a -> any(i -> _eq(a.args[i], u) && a.coefs[i] == c, eachindex(a.args)), args[2:end])]
            if !isempty(common)
                ct = _add_many(Sym(c) * u for (c, u) in common)
                return ct + apply_function(t.name, Sym[expand(a - ct) for a in args])
            end
        end
        # a common numeric factor: its absolute value out, and max ↔ min for a negative one
        cs = [a.kind == MUL ? a.num : a.kind == NUMBER ? a.num : N1 for a in args]
        all(isexact, cs) || return nothing
        allneg = all(c -> c < N0, cs)
        g = reduce(gcd, (Rational{BigInt}(_numerator(c), _denominator(c)) for c in cs))
        (isone(g) && !allneg) && return nothing
        k = Sym(SymNumber(allneg ? -g : g))
        apply_function(allneg ? other : t.name, Sym[expand(a / k) for a in args]) * k
    end
end

# about how many terms e has expanded, numerator and denominator together (saturating)
const _EXPANSION_MAX = 50_000.0
function _expansion_estimate(e::Sym)::Float64
    k = e.kind
    k == ADD && return min(1e15, sum(_expansion_estimate, e.args; init = iszero(e.num) ? 0.0 : 1.0))
    k == MUL && return min(1e15, prod((_power_estimate(b, p) for (b, p) in zip(e.args, e.exps)); init = 1.0))
    k == POW && return _power_estimate(e.args[1], e.args[2])
    1.0
end
function _power_estimate(b::Sym, p::Sym)
    (p.kind == NUMBER && _isint(p.num)) || return 1.0       # a radical or symbolic power: an atom
    n = abs(Float64(p.num))
    t = _expansion_estimate(b)
    t <= 1 && return 1.0
    min(1e15, t^n, exp(_lgamma_est(t + n) - _lgamma_est(t) - _lgamma_est(n + 1)))  # C(t+n-1, n)
end
_lgamma_est(x::Float64) = x < 1 ? 0.0 : (x - 0.5) * log(x) - x + 0.919   # Stirling, rough

_bigpower(e::Sym) = (e.kind == POW && e.args[1].kind == ADD && e.args[2].kind == NUMBER && abs(e.args[2].num) > SymNumber(10)) ||
    (e.kind == MUL && any(i -> e.args[i].kind == ADD && e.exps[i].kind == NUMBER && abs(e.exps[i].num) > SymNumber(10), eachindex(e.args))) ||
    any(_bigpower, e.args)

# (not tried when the expanded numerator would be large: factoring it can take minutes,
# as for the integral in SymPy #27957)
const _FACTOR_MAX_TERMS = 2000

function _factor_or_nothing(e::Sym)
    n, d = numer_denom(e)
    _expanded_terms(n) > _FACTOR_MAX_TERMS && return nothing
    _ispolylike(n) && (_ispolylike(d) || d.kind == NUMBER) || return nothing
    with(() -> factor(n) / factor(d), _FACTOR_EXPS => false)
end
# an upper bound on the number of terms of expand(e) (capped)
function _expanded_terms(e::Sym)
    cap = 10^9
    k = e.kind
    k == ADD && return min(cap, (iszero(e.num) ? 0 : 1) + sum(_expanded_terms, e.args; init = 0))
    if k == MUL || k == POW
        bs, ps = k == MUL ? (e.args, e.exps) : ((e.args[1],), (e.args[2],))
        t = 1
        for (b, p) in zip(bs, ps)
            m = _expanded_terms(b)
            m == 1 && continue
            q = p.kind == NUMBER && _isint(p.num) && p.num > N0 ? min(_toint_or(p.num, 64), 64) : 1
            t = min(cap, t * (q == 1 ? m : binomial(big(m + q - 1), q) > cap ? cap : Int(binomial(big(m + q - 1), q))))
        end
        return t
    end
    1
end
_toint_or(x::SymNumber, d::Int) = (x.tag == T_SMALL ? x.n : d)

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
# and a long exact number one more for every 4 digits past 6, so that 30 log(5) is
# smaller than log(931322574615478515625) (SymPy #23210)
_division_cost(p::Sym) = p.kind == NUMBER && p.num < N0 ? 1 : 0
function _digit_cost(q::SymNumber)
    q.tag == T_SMALL || q.tag == T_BIG || return 0
    d = q.tag == T_SMALL ? ndigits(q.n) + (q.d == 1 ? 0 : ndigits(q.d)) :
        ndigits(numerator(q.big)) + (isone(denominator(q.big)) ? 0 : ndigits(denominator(q.big)))
    max(0, d - 6) ÷ 4
end

function _nodes(e::Sym)
    k = e.kind
    k == NUMBER && return 1 + _digit_cost(e.num)
    isempty(e.args) && return 1
    n = 1
    if k == ADD
        n += iszero(e.num) ? 0 : 1 + _digit_cost(e.num)
        for (c, a) in zip(e.coefs, e.args)
            n += _nodes(a) + (abs(c) == N1 ? 0 : 1 + _digit_cost(c))
        end
    elseif k == MUL
        n += abs(e.num) == N1 ? 0 : 1 + _digit_cost(e.num)
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

_negated_exponent(p::Sym) = p.kind == MUL ? p.num < N0 :
    p.kind == ADD && p.num <= N0 && all(c -> c < N0, p.coefs)

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
                # b^(-p) joins the group of p as (1/b)^p, for b > 0 ((a + b)^(n+1) /
                # (a + 2b)^(n+1) = ((a + b)/(a + 2b))^(n+1); SymPy #22546)
                q = expand(-p2)
                if !haskey(groups, p2) && haskey(groups, q) && (b2.kind == NUMBER || _ispos(b2, A))
                    push!(groups[q], inv(b2))
                else
                    push!(get!(groups, p2, Sym[]), b2)
                end
            else
                other *= b2^p2
            end
        end
        for (p, bs) in groups
            b = prod(bs; init = ONE)
            # (a negated exponent as the power of the reciprocal, (u/v)^(n+1) rather than
            # (v/u)^(-n-1), when more than one base merged)
            length(bs) >= 2 && _negated_exponent(p) ? (other *= inv(b)^expand(-p)) : (other *= b^p)
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

# for simplify's candidates: e unchanged when a circular function's argument has an
# integer multiple beyond _CIRC_MAX_MULTIPLE (a polynomial of that degree; SymPy #27293)
function _expand_trig_small(e::Sym)
    big = false
    _replace_walk(e) do t
        if !big && t.kind == FN && length(t.args) == 1 && t.name in _CIRCFNS
            q = _trig_split(t.args[1])[1]
            abs(numerator(q)) > _CIRC_MAX_MULTIPLE && (big = true)
        end
        nothing
    end
    big ? e : _expand_trig(e)
end

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
    # e's own denominator 0 by the identities (before cancel takes the common factor out):
    # e is defined nowhere, and no form of it is simpler (sec² from (tan² - sec² + 1)/
    # (sin² + cos² - 1); SymPy #21340)
    d0 = numer_denom(e)[2]
    if _hassymbol(d0)
        pd = _circ_poly(d0, bases, vars, M)
        if pd !== nothing
            nd, dd = numer_denom(cancel(pd))
            any(o -> iszero(_circ_reduce(nd, vars, o, M)), Iterators.product(((true, false) for _ in bases)...)) && return e
        end
    end
    n, d = numer_denom(cancel(p))
    best, size = e, _nodes(e)
    consider(r) = (k = _nodes(r); k < size && ((best, size) = (r, k)))
    s, c = M.fns[1], M.fns[2]
    back(q) = _circ_tan(subs(q, (v => apply_function(f, Sym[g]) for ((S, C), g) in zip(vars, bases)
                                       for (v, f) in ((S, s), (C, c)))...), M)
    for order in Iterators.product(((true, false) for _ in bases)...)
        ns = _circ_reduce(n, vars, order, M)
        ds = _circ_reduce(d, vars, order, M)
        iszero(ds) && return e
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
# quotients in S, C; nothing if an argument is not a multiple of a base, or a multiple
# larger than _CIRC_MAX_MULTIPLE
const _CIRC_MAX_MULTIPLE = 32
function _circ_poly(e::Sym, bases::Vector{Sym}, vars, M::_Circ)
    ok = Ref(true)
    function find(u)
        for (i, g) in enumerate(bases)
            k = _circ_multiple(u, g)
            k === nothing && continue
            # a polynomial of degree k in S, C: sin(316801t) is beyond any use (SymPy #27293)
            abs(k) > _CIRC_MAX_MULTIPLE && return nothing
            return i, k
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
        # a sum p/q + r over its common denominator, log(n/d) = log(n) - log(d) for n > 0
        # and d ≥ 0 (d ≠ 0 where the logarithm is defined): log(1 + 1/x²) = log(x² + 1) -
        # log(x²) for real x (SymPy #23569)
        if u.kind == ADD && _hassymbol(u)
            n, d = numer_denom(u)
            if !(d.kind == NUMBER) && (pos(n) || n.kind != NUMBER && _ispos(n, A)) && (force || _nonneg(d, A))
                return _log_expand(log(n), A, force) - _log_expand(log(d), A, force)
            end
        end
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
            # (and the root's powers: 1/sqrt(2 + √3) too; SymPy #25915)
            if _denominator(q.num) == 2 && (d = _denest(b)) !== nothing
                return expand(d^Sym(2q.num))
            end
            # the principal cube root of a positive surd
            if _denominator(q.num) == 3 && (v = _realval(b)) !== nothing && v > 0 &&
               (d = _denest_cbrt(b)) !== nothing
                return expand(d^Sym(3q.num))
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
    if e.kind == FN && e.name === :cbrt && length(e.args) == 1      # the real cube root
        b = _radicals(e.args[1], A)
        d = _denest_cbrt(b)
        return d === nothing ? cbrt(b) : d
    end
    _rebuild_node(e, Sym[_radicals(a, A) for a in e.args])
end

_realval(b::Sym) = (v = try _evalc(b) catch; return nothing end; iszero(imag(v)) ? real(v) : nothing)

# cbrt(a + b sqrt(c)) = p + q sqrt(c) with p, q rational (SymPy #27082): then
# n = p^2 - q^2 c is the rational cube root of a^2 - b^2 c, and s = 2p a rational root of
# s^3 - 3n s - 2a (s the sum of the cube root and its conjugate, whose product is n)
function _denest_cbrt(u::Sym; recurse::Bool = true)
    t = _quadratic_surd(u)
    t === nothing && return nothing
    a, b, c = t
    Nm = a * a - b * b * c
    nn = cbrt(Sym(Nm))
    if !(nn.kind == NUMBER && isexact(nn.num))
        # the norm not a cube: cbrt(u) = cbrt(N u)/cbrt(N), N u of norm N³ (cbrt(22 + 9√6)
        # = (2 + √6)/2^(1/3); SymPy #26857)
        recurse || return nothing
        r = _denest_cbrt(Sym(Nm) * u; recurse = false)
        return r === nothing ? nothing : r / cbrt(Sym(Nm))
    end
    z = Sym(:_s)
    f = try factor(z^3 - 3nn * z - 2Sym(a)) catch; return nothing end
    for g in (f.kind == MUL ? f.args : [f])
        g.kind == ADD || continue
        k0, k1 = subs(g, z => 0), diff(g, z)
        (k1.kind == NUMBER && !iszero(k1)) || continue           # linear in z
        p = -k0 / k1 / 2
        q2 = (p^2 - nn) / Sym(c)
        q = sqrt(q2)
        q.kind == NUMBER || continue
        r = p + (b > N0 ? q : -q) * sqrt(Sym(c))
        iszero(expand(r^3 - u)) && return r
    end
    nothing
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

# n/d with d linear in a square root R = √Q of an expression with symbols, d = a + bR:
# times the conjugate a - bR over the norm a² - b²Q, which has no R, then cancelled
# (the derivative of log(b + 2√c√(a + bx + cx²) + 2cx) is √c/√(a + bx + cx²); SymPy
# #23641); nothing when there is no such root
function _conjugate_radical_den(e::Sym)
    # (expanded: a factor Q beside √Q would be taken as (√Q)² by the substitution)
    n, d = map(expand, numer_denom(cancel(e)))
    rs = Sym[]
    _replace_walk(d) do t
        t.kind == POW && t.args[2].kind == NUMBER && t.args[2].num == SymNumber(1 // 2) &&
            _hassymbol(t.args[1]) && !any(q -> _eq(q, t), rs) && push!(rs, t)
        nothing
    end
    isempty(rs) && return nothing
    best = nothing
    for R in rs
        u = Sym(Symbol("_conj_r"))
        du = expand(subs(d, R => u))
        P = try Poly{Sym}(du, [u]) catch; continue end
        degree(P, u) == 1 || continue
        conj = subs(expand(coeff(du, u, 0) - coeff(du, u, 1) * u), u => R)
        norm = expand(d * conj)
        has(norm, R) && continue
        c = try _cancel_halfpowers(expand(n * conj), norm) catch err; err isa InterruptException && rethrow(); continue end
        (best === nothing || _nodes(c) < _nodes(best)) && (best = c)
    end
    best
end

# p/q cancelled with √s of each symbol s in them a variable u (s = u², s^(3/2) = u³),
# so that c^(3/2) and √c share their factor
function _cancel_halfpowers(p::Sym, q::Sym)
    m = Dict{Sym,Sym}()
    function find(e)
        if e.kind == POW && _isnamed(e.args[1]) && e.args[2].kind == NUMBER && _denominator(e.args[2].num) == 2
            m[e.args[1]] = Sym(Symbol("_half_", e.args[1]))
        elseif e.kind == MUL
            for (b, x) in zip(e.args, e.exps)
                _isnamed(b) && x.kind == NUMBER && _denominator(x.num) == 2 && (m[b] = Sym(Symbol("_half_", b)))
            end
        end
        foreach(find, e.args)
    end
    find(p); find(q)
    isempty(m) && return cancel(p / q)
    function sub(e)
        _isnamed(e) && haskey(m, e) && return m[e]^2
        isempty(e.args) && return e
        if e.kind == POW && haskey(m, e.args[1]) && e.args[2].kind == NUMBER && isexact(e.args[2].num)
            return m[e.args[1]]^Sym(2 * e.args[2].num)
        elseif e.kind == MUL
            return Sym(e.num) * prod((haskey(m, b) && x.kind == NUMBER && isexact(x.num) ? m[b]^Sym(2 * x.num) : sub(b)^x
                                      for (b, x) in zip(e.args, e.exps)); init = ONE)
        end
        _rebuild_node(e, Sym[sub(a) for a in e.args])
    end
    r = cancel(sub(p) / sub(q))
    subs(r, (u => sqrt(b) for (b, u) in m)...)
end

# exp(w) with w linear in log(u) (a coefficient of any kind) as u^c exp(w - c log(u)):
# exp((a - 1) log x) is x^(a - 1) by the definition of the power (SymPy #21360)
function _exp_log_powers(e::Sym)
    _replace_walk(e) do t
        t.kind == POW && t.args[1] === E || return nothing
        w = expand(t.args[2])
        logs = Sym[]
        _replace_walk(w) do q
            q.kind == FN && q.name === :log && length(q.args) == 1 && !any(o -> _eq(o, q), logs) && push!(logs, q)
            nothing
        end
        isempty(logs) && return nothing
        out, rest = ONE, w
        for L in logs
            λ = Sym(Symbol("_elog"))
            r = expand(subs(rest, L => λ))
            c = expand(diff(r, λ))
            (has(c, λ) || iszero(c)) && continue
            out *= L.args[1]^c
            rest = expand(subs(r - c * λ, λ => L))
        end
        isone(out) ? nothing : out * exp(rest)
    end
end

# in each product, |u|^p sign(u)^q (p, q positive integers, q odd) as |u|^(p-1) u
# sign(u)^(q-1), and for real u, u^(2k) beside |u|^m as |u|^(m + 2k)
function _abs_sign_pairs(e::Sym, A::Assumptions)
    isempty(e.args) && return e
    e = _rebuild_node(e, Sym[_abs_sign_pairs(a, A) for a in e.args])
    e.kind == MUL || return e
    bs, ps = collect(e.args), collect(e.exps)
    posint(q) = q.kind == NUMBER && isexact(q.num) && _isint(q.num) && q.num > N0
    changed = false
    for i in eachindex(bs)
        b = bs[i]
        (b.kind == FN && b.name === :abs && length(b.args) == 1) || continue
        u = b.args[1]
        j = findfirst(j -> j != i && bs[j].kind == FN && bs[j].name === :sign && length(bs[j].args) == 1 &&
                           _eq(bs[j].args[1], u) && posint(ps[j]) && isodd(_toint(ps[j].num)) && posint(ps[i]), eachindex(bs))
        if j !== nothing
            ps[i] -= 1; ps[j] -= 1
            push!(bs, u); push!(ps, ONE)
            changed = true
            continue
        end
        _isreal(u, A) === true || continue
        k = findfirst(k -> k != i && _eq(bs[k], u) && ps[k].kind == NUMBER && isexact(ps[k].num) &&
                           _isint(ps[k].num) && iseven(_toint(ps[k].num)), eachindex(bs))
        if k !== nothing && ps[i].kind == NUMBER
            ps[i] += ps[k]; ps[k] = ZERO
            changed = true
        end
    end
    changed || return e
    Sym(e.num) * prod((b^p for (b, p) in zip(bs, ps)); init = ONE)
end

# each square root of a + bi (a, b rational, b ≠ 0) as x + iy when x = √((r + a)/2) and
# y = sign(b)√((r - a)/2), r = √(a² + b²), are rational: the principal root
function _complex_sqrt_denest(e::Sym)
    _replace_walk(e) do t
        if t.kind == POW && t.args[2].kind == NUMBER && t.args[2].num == SymNumber(1 // 2)
            r = _csqrt_parts(t.args[1])
            return r === nothing ? nothing : r
        elseif t.kind == MUL
            # √c √(a + bi), c a positive rational split off first (√(8 + 6i) is √2 √(4 + 3i)):
            # rejoined, √(c(a + bi)) (= 3 + i; SymPy #12420)
            half(q) = q.kind == NUMBER && q.num == SymNumber(1 // 2)
            ri = [i for (i, (b, q)) in enumerate(zip(t.args, t.exps)) if half(q) && b.kind == NUMBER && isexact(b.num) && b.num > N0]
            ci = [i for (i, (b, q)) in enumerate(zip(t.args, t.exps)) if half(q) && has(b, IM) && !_hassymbol(b) && b.kind != NUMBER]
            if length(ri) == 1 && length(ci) == 1
                r = _csqrt_parts(expand(t.args[ri[1]] * t.args[ci[1]]))
                if r !== nothing
                    return Sym(t.num) * r * prod((b^q for (i, (b, q)) in enumerate(zip(t.args, t.exps)) if i != ri[1] && i != ci[1]); init = ONE)
                end
            end
            changed = false
            fs = Sym[]
            for (b, q) in zip(t.args, t.exps)
                if q.kind == NUMBER && isexact(q.num) && _denominator(q.num) == 2 && (r = _csqrt_parts(b)) !== nothing
                    push!(fs, r^(Sym(_numerator(q.num))))
                    changed = true
                else
                    push!(fs, b^q)
                end
            end
            return changed ? Sym(t.num) * prod(fs) : nothing
        end
        nothing
    end
end
function _csqrt_parts(z::Sym)
    has(z, IM) && !_hassymbol(z) || return nothing
    a, b = real(z), imag(z)
    (a.kind == NUMBER && b.kind == NUMBER && isexact(a.num) && isexact(b.num) && !iszero(b.num)) || return nothing
    r = sqrt(a^2 + b^2)
    r.kind == NUMBER || return nothing
    x, y = sqrt((r + a) / 2), sqrt((r - a) / 2)
    (x.kind == NUMBER && y.kind == NUMBER) || return nothing
    x + IM * (b.num > N0 ? y : -y)
end

function _cancel_in_roots(e::Sym)
    exps = Dict{Sym,Vector{SymNumber}}()
    note(b, q) = (_hassymbol(b) && q.kind == NUMBER && isexact(q.num) && push!(get!(exps, b, SymNumber[]), q.num); nothing)
    _replace_walk(e) do t
        if t.kind == POW
            note(t.args[1], t.args[2])
        elseif t.kind == MUL
            foreach(note, t.args, t.exps)
        end
        nothing
    end
    bs = [b for (b, qs) in exps if any(q -> !_isint(q), qs)]
    length(bs) == 1 || return nothing
    b = only(bs)
    has(e, b) || return nothing
    q = reduce(lcm, (_denominator(r) for r in exps[b]); init = big(1))
    t = _fresh_dummy(e)
    intpow(v) = (w = v * SymNumber(q); t^Sym(SymNumber(_numerator(w))))      # (an integer exponent)
    # (the walk is bottom-up: a base b is already t^q when its power is reached)
    tq = t^Sym(q)
    isb(base) = _eq(base, b) || _eq(base, tq)
    sw(base, r) = (isb(base) && r.kind == NUMBER && isexact(r.num)) ? intpow(r.num) : nothing
    h = _replace_walk(e) do u
        _eq(u, b) && return tq
        if u.kind == POW
            return sw(u.args...)
        elseif u.kind == MUL && any(isb, u.args)
            return Sym(u.num) * prod((something(sw(a, r), a^r) for (a, r) in zip(u.args, u.exps)); init = ONE)
        end
        nothing
    end
    (has(h, b) || _has_node(h, u -> u.kind == POW && _eq(u.args[1], tq))) && return nothing                       # b still inside another way
    c = cancel(h)
    subs(c, t => b^(1 / Sym(q)))
end

# in a sum, c Sum(f, k, a, b) - c Sum(f, k, a', b') with the same f and k, and a = a',
# b - b' an integer d (|d| ≤ 20), as c times the terms k = b' + 1, …, b (or with b = b',
# the terms a, …, a' - 1); throughout e
function _combine_sums(e::Sym)
    isempty(e.args) && return e
    e = _rebuild_node(e, Sym[_combine_sums(a) for a in e.args])
    e.kind == ADD || return e
    ts = [(Sym(c), t) for (c, t) in zip(e.coefs, e.args)]
    issum(t) = t.kind == FN && t.name === :Sum && length(t.args) == 4
    smallint(d) = d.kind == NUMBER && isexact(d.num) && _isint(d.num) && abs(d.num) <= SymNumber(20)
    used = falses(length(ts))
    out = Sym[]
    for i in eachindex(ts), j in eachindex(ts)
        (i < j && !used[i] && !used[j]) || continue
        (ci, ti), (cj, tj) = ts[i], ts[j]
        (issum(ti) && issum(tj) && iszero(ci + cj)) || continue
        f, k, a, b = ti.args
        g, k2, a2, b2 = tj.args
        (_eq(f, subs(g, k2 => k))) || continue
        part = nothing
        if _eq(a, a2) && smallint(expand(b - b2))
            d = _toint(expand(b - b2).num)
            part = d >= 0 ? _add_many(subs(f, k => b2 + j2) for j2 in 1:d) :
                            -_add_many(subs(f, k => b + j2) for j2 in 1:(-d))
        elseif _eq(b, b2) && smallint(expand(a2 - a))
            d = _toint(expand(a2 - a).num)
            part = d >= 0 ? _add_many(subs(f, k => a + j2) for j2 in 0:(d - 1)) :
                            -_add_many(subs(f, k => a2 + j2) for j2 in 0:(-d - 1))
        end
        part === nothing && continue
        used[i] = used[j] = true
        push!(out, ci * part)
    end
    any(used) || return e
    Sym(e.num) + _add_many(c * t for (u, (c, t)) in zip(used, ts) if !u) + _add_many(out)
end

# in each product, the abs factors with one integer exponent merged: |a|^k |b|^k = |ab|^k
function _merge_abs(e::Sym)
    isempty(e.args) && return e
    e = _rebuild_node(e, Sym[_merge_abs(a) for a in e.args])
    e.kind == MUL || return e
    # |a|^p |b|^q = |a^p b^q| for integers p, q, cancelled inside: |x² - 1|/|x - 1| is
    # |x + 1| (SymPy #17292)
    ix = Int[]
    for (i, (b, q)) in enumerate(zip(e.args, e.exps))
        b.kind == FN && b.name === :abs && length(b.args) == 1 && q.kind == NUMBER && isexact(q.num) && _isint(q.num) &&
            push!(ix, i)
    end
    length(ix) >= 2 || return e
    inner = prod(e.args[i].args[1]^e.exps[i] for i in ix)
    any(i -> _toint(e.exps[i].num) < 0, ix) && (inner = cancel(inner))
    used = falses(length(e.args))
    used[ix] .= true
    # (an integer power inside as a power of |u|: |t|², not |t²|)
    ab = inner.kind == POW && inner.args[2].kind == NUMBER && _isint(inner.args[2].num) ?
         abs(inner.args[1])^inner.args[2] : abs(inner)
    Sym(e.num) * ab * prod((b^q for (u, b, q) in zip(used, e.args, e.exps) if !u); init = ONE)
end

function _pw_eq_subs(e::Sym)
    _replace_walk(e) do t
        (t.kind == FN && t.name === :piecewise) || return nothing
        args = copy(t.args)
        hit = false
        for i in 1:2:length(args)
            eqs = Pair{Sym,Sym}[]
            c = args[i]
            cs = c.kind == FN && c.name === :And ? c.args : Sym[c]
            for d in cs
                (d.kind == FN && get(_COND_RELS, d.name, nothing) === :eq) || continue
                l, r = d.args
                _isnamed(l) && !has(r, l) ? push!(eqs, l => r) : _isnamed(r) && !has(l, r) && push!(eqs, r => l)
            end
            isempty(eqs) && continue
            v = subs(args[i + 1], eqs...)
            _eq(v, args[i + 1]) || (args[i + 1] = v; hit = true)
        end
        hit ? apply_function(:piecewise, args) : nothing
    end
end

function _split_positive_radicals(e::Sym, A::Assumptions)
    function split(b::Sym, q::Sym)
        (q.kind == NUMBER && isexact(q.num) && !isone(_denominator(q.num))) || return nothing
        f = b.kind == ADD && _hassymbol(b) ? (try factor(b) catch err; err isa InterruptException && rethrow(); b end) : b
        f.kind == MUL || return nothing
        out, rest, hit = ONE, Sym(f.num), false
        for (c, p) in zip(f.args, f.exps)
            if _isreal(p, A) === true && _ispos(c, A)
                out *= c^expand(p * q); hit = true
            elseif p.kind == NUMBER && _isint(p.num) && _ispos(-c, A)
                out *= (-c)^expand(p * q); hit = true
                isodd(_toint(p.num)) && (rest = -rest)
            else
                rest *= c^p
            end
        end
        hit ? out * rest^q : nothing
    end
    _replace_walk(e) do t
        if t.kind == POW
            return split(t.args[1], t.args[2])
        elseif t.kind == MUL && any(q -> q.kind == NUMBER && isexact(q.num) && !isone(_denominator(q.num)), t.exps)
            r = Sym(t.num)
            for (b, q) in zip(t.args, t.exps)
                v = split(b, q)
                r *= v === nothing ? b^q : v
            end
            return r
        end
        nothing
    end
end

_int_rational_power(b::Sym, q::Sym) = b.kind == NUMBER && isexact(b.num) && _isint(b.num) && b.num > N1 &&
    q.kind == NUMBER && isexact(q.num) && !_isint(q.num)
# the prime factors of n by trial division to 1000, or nothing when a larger one remains
function _small_prime_factors(n::Integer)
    out = Dict{Int,Int}()
    m = BigInt(n)
    for p in 2:1000
        while iszero(mod(m, p))
            out[p] = get(out, p, 0) + 1
            m = div(m, p)
        end
        m == 1 && return out
    end
    nothing
end
_merge_prime_powers(e::Sym) = _replace_walk(e) do t
    t.kind == MUL || return nothing
    ix = [i for i in eachindex(t.args) if _int_rational_power(t.args[i], t.exps[i])]
    length(ix) >= 2 || return nothing
    acc = Dict{Int,SymNumber}()
    for i in ix
        f = _small_prime_factors(_toint(t.args[i].num))
        f === nothing && return nothing
        for (p, k) in f
            acc[p] = get(acc, p, N0) + SymNumber(k) * t.exps[i].num
        end
    end
    out = Sym(t.num)
    for (p, q) in acc
        out *= Sym(p)^Sym(q)
    end
    out * prod((t.args[i]^t.exps[i] for i in eachindex(t.args) if !(i in ix)); init = ONE)
end

function _sum_contents_out(e::Sym)
    e.kind == MUL || return e
    c = Sym(e.num)
    out = ONE
    for (b, p) in zip(e.args, e.exps)
        if b.kind == ADD && p.kind == NUMBER && _isint(p.num) && all(isexact, b.coefs) && isexact(b.num)
            ns = [_numerator(q) for q in b.coefs]
            iszero(b.num) || push!(ns, _numerator(b.num))
            g = gcd(reduce(gcd, ns), _denominator(c.num))       # (only what cancels)
            if g > 1
                c *= Sym(g)^p
                out *= expand(b / g)^p
                continue
            end
        end
        out *= b^p
    end
    c * out
end

# a piecewise function of two pieces, u < 0 (or u ≤ 0) => -u and the rest => u, as |u|
# (or the reverse, u > 0 => u); throughout e
function _pw_as_abs(e::Sym)
    _replace_walk(e) do t
        (t.kind == FN && t.name === :piecewise && length(t.args) == 4) || return nothing
        c1, v1, c2, v2 = t.args
        (c1.kind == FN && haskey(_COND_RELS, c1.name)) || return nothing
        r = _relation_of(c1)
        d = expand(r.lhs - r.rhs)
        u = r.op in (:lt, :le) ? d : r.op in (:gt, :ge) ? -d : return nothing
        # c1 says u < 0 (u ≤ 0): v1 = -w, v2 = w with w = ±u
        (iszero(expand(v1 + v2)) && (iszero(expand(v2 - u)) || iszero(expand(v2 + u)))) || return nothing
        w = iszero(expand(v2 - u)) ? u : -u
        # (c2 must be its complement or true)
        c2 == ONE || (c2.kind == FN && haskey(_COND_RELS, c2.name) &&
                      (q = _relation_of(c2); iszero(expand((q.lhs - q.rhs) - (r.lhs - r.rhs))) &&
                       q.op === (r.op === :lt ? :ge : r.op === :le ? :gt : r.op === :gt ? :le : :lt))) || return nothing
        _eq(w, u) ? abs(u) : -abs(u)
    end
end

_small_den(q) = (d = _denominator(q); d > 1000 ? 1001 : Int(d))
function _radical_degree_bound(e::Sym)
    seen = Set{Tuple{Sym,Int}}()
    _walk(u) = begin
        if u.kind == POW && u.args[2].kind == NUMBER && !_isint(u.args[2].num)
            push!(seen, (u.args[1], _small_den(u.args[2].num)))
        elseif u.kind == MUL
            for (b, x) in zip(u.args, u.exps)
                x.kind == NUMBER && !_isint(x.num) && push!(seen, (b, _small_den(x.num)))
            end
        end
        for a in u.args; _walk(a) end
    end
    _walk(e)
    d = 1
    for (_, q) in seen
        d = Base.checked_mul(d, q) ; d > 1000 && return d
    end
    d
end
function _rational_by_minpoly(e::Sym)
    (e.kind != NUMBER && _isexactex(e) && _nodes(e) <= 80) || return nothing
    _has_node(e, t -> t.kind == FN) && return nothing            # radicals only
    _has_node(e, t -> (t.kind == POW && t.args[2].kind == NUMBER && !_isint(t.args[2].num)) ||
                      (t.kind == MUL && any(q -> q.kind == NUMBER && !_isint(q.num), t.exps))) || return nothing
    # the minimal polynomial's degree is at most the product of the radical indices:
    # nested cube roots (Cardano) make it large and its coefficients huge
    _radical_degree_bound(e) <= 16 || return nothing
    v = try _evalc(e) catch; return nothing end
    abs(imag(v)) < 1e-12 * max(1, abs(v)) || return nothing
    t = Sym(:_mp_t)
    p = try minimal_polynomial(e, t) catch err; err isa InterruptException && rethrow(); return nothing end
    P = try Poly{SymNumber}(expand(p), [t]) catch; return nothing end
    degree(P, t) == 1 || return nothing
    c1, c0 = coeff(p, t, 1), coeff(p, t, 0)
    cancel(-c0 / c1)
end


# in each product with a factor δ(c(v - a)) (to the first power, v a symbol, c and a free
# of v), the other factors at v = a, where they are finite
function _dirac_sift(e::Sym)
    _replace_walk(e) do t
        t.kind == MUL || return nothing
        i = findfirst(j -> t.args[j].kind == FN && t.args[j].name === :dirac && length(t.args[j].args) == 1 &&
                          _isone(t.exps[j]), eachindex(t.args))
        i === nothing && return nothing
        g = t.args[i].args[1]
        vs = [v for v in _named(g) if _ispoly_in(g, v) && degree(Poly{Sym}(expand(g), [v]), v) == 1]
        isempty(vs) && return nothing
        v = vs[1]
        a = -coeff(expand(g), v, 0) / coeff(expand(g), v, 1)
        rest = t / t.args[i]
        has(rest, v) || return nothing
        r = try subs(rest, v => a) catch err; _isdivzero(err) || rethrow(); return nothing end
        _finite_value(r) || _hassymbol(r) || return nothing
        r * t.args[i]
    end
end
