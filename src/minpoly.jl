# Minimal polynomials of algebraic numbers (SymPy's minimal_polynomial, its "compose"
# method). An expression built from rationals by +, ×, rational powers (im = (-1)^(1/2)
# among them) and the trigonometric functions of rational multiples of π (through
# (-1)^r) is algebraic. Its minimal polynomial is found from the parts': for a root a of
# P and b of Q, a + b is a root of Res_y(P(y), Q(t - y)), a·b of Res_y(P(y), y^q Q(t/y))
# (q = deg Q) and a^(p/q) of Res_y(P(y), t^q - y^p). The resultant is factored over the
# rationals and the factor vanishing at the part's value, computed to high precision, is
# kept, which also picks the branch (cbrt(-8) is -2, (-8)^(1/3) is not).

const _MP_T = Sym(:_minpoly_t)
const _MP_Y = Sym(:_minpoly_y)

"""
    SymJ.minimal_polynomial(ex, x)

The minimal polynomial over the rationals of the algebraic number `ex`, as a polynomial
in `x` with integer coefficients, relatively prime, and a positive leading coefficient:
the irreducible polynomial of least degree with `ex` as a root. `ex` is built from
rationals by arithmetic and rational powers (`im`, `sqrt`, `cbrt` and other radicals),
and may contain `sin`, `cos`, `tan` of rational multiples of `π`. An `ArgumentError` for
an expression that is not of that form (with a free symbol, a float, `π`, `ℯ`, …). As
SymPy's `minimal_polynomial`. Not exported.

```julia
SymJ.minimal_polynomial(sqrt(Sym(2)) + sqrt(Sym(3)), x)    # x^4 - 10x^2 + 1
SymJ.minimal_polynomial((1 + sqrt(Sym(5))) / 2, x)         # x^2 - x - 1
SymJ.minimal_polynomial(cos(PI / 7), x)                    # 8x^3 - 4x^2 - 4x + 1
SymJ.minimal_polynomial(Sym(2)^(1 // 3) + im, x)          # degree 6
```
"""
function minimal_polynomial(ex, x::Sym)
    _isnamed(x) || throw(ArgumentError("minimal_polynomial: $x is not a symbol"))
    e = Sym(ex)
    _hassymbol(e) && throw(ArgumentError("minimal_polynomial: $e has free symbols; it is not a number"))
    _from_xcoeffs(_mp_normal(_mp(e)), x)
end

# the minimal polynomial of e, as coefficients (lowest degree first)
function _mp(e::Sym)::Vector{Sym}
    k = e.kind
    if k == NUMBER
        isexact(e.num) || throw(ArgumentError("minimal_polynomial: $e is a float; use nsimplify or rationalize first"))
        return Sym[Sym(-_numerator(e.num)), Sym(_denominator(e.num))]
    elseif k == ADD
        ts = Sym[Sym(c) * a for (c, a) in zip(e.coefs, e.args)]
        iszero(e.num) || push!(ts, Sym(e.num))
        return _mp_fold(ts, _mp_sum)
    elseif k == MUL
        fs = Sym[isone(p) ? b : b^p for (b, p) in zip(e.args, e.exps)]
        isone(e.num) || push!(fs, Sym(e.num))
        return _mp_fold(fs, _mp_product)
    elseif k == POW
        b, p = e.args
        (p.kind == NUMBER && isexact(p.num)) ||
            throw(ArgumentError("minimal_polynomial: $e is not algebraic (or not known to be): its exponent is not rational"))
        return _mp_select(_mp_power(_mp(b), Int(_numerator(p.num)), Int(_denominator(p.num))), e)
    elseif k == FN
        if e.name === :cbrt && length(e.args) == 1
            return _mp_select(_mp_power(_mp(e.args[1]), 1, 3), e)
        elseif e.name in (:sin, :cos, :tan) && length(e.args) == 1
            r = expand(e.args[1] / PI)
            if r.kind == NUMBER && isexact(r.num)
                z = MINUS_ONE^r                                       # exp(iπr)
                w = e.name === :cos ? (z + 1 / z) / 2 : e.name === :sin ? (z - 1 / z) / (2IM) :
                    (z - 1 / z) / (IM * (z + 1 / z))
                return _mp_select(_mp(w), e)
            end
        end
    end
    throw(ArgumentError("minimal_polynomial: $e is not algebraic (or not known to be)"))
end

# the parts combined two at a time, the minimal polynomial kept at each step
function _mp_fold(parts::Vector{Sym}, combine)
    acc, P = parts[1], _mp(parts[1])
    for q in parts[2:end]
        acc2 = combine === _mp_sum ? acc + q : acc * q
        P = _mp_select(combine(P, _mp(q)), acc2)
        acc = acc2
    end
    P
end

_mp_expr(c::Vector{Sym}, v::Sym) = _add_many(c[k + 1] * v^k for k in eachindex(c) .- 1 if !iszero(c[k + 1]))

# a polynomial (as an expression in _MP_T) with a + b as a root
_mp_sum(P, Q) = resultant(_mp_expr(P, _MP_Y), expand(_mp_expr(Q, _MP_T - _MP_Y)), _MP_Y)
# with a·b as a root: y^q Q(t/y)
function _mp_product(P, Q)
    q = length(Q) - 1
    R = _add_many(Q[k + 1] * _MP_T^k * _MP_Y^(q - k) for k in 0:q if !iszero(Q[k + 1]))
    resultant(_mp_expr(P, _MP_Y), R, _MP_Y)
end
# with a^(p/q) as a root (a reciprocal first for p < 0)
function _mp_power(P, p::Int, q::Int)
    p < 0 && (P = reverse(P); p = -p)
    iszero(p) && return Sym[MINUS_ONE, ONE]
    resultant(_mp_expr(P, _MP_Y), _MP_T^q - _MP_Y^p, _MP_Y)
end

# the irreducible factor of R (an expression in _MP_T) that vanishes at the value of e
function _mp_select(R::Sym, e::Sym)
    _, fs = factor(Poly{SymNumber}(expand(R), [_MP_T]))
    cands = [_xcoeffs(f) for (f, _) in fs if degree(f, _MP_T) >= 1]
    length(cands) == 1 && return cands[1]
    for prec in (256, 1024, 4096)
        sc = setprecision(BigFloat, prec) do
            v = evaluate(e, big(0.0))
            [_mp_residual(c, v) for c in cands]
        end
        i = argmin(sc)
        others = [sc[j] for j in eachindex(sc) if j != i]
        sc[i] < big(2.0)^(-prec ÷ 2) && minimum(others) > big(2.0)^(-prec ÷ 4) && return cands[i]
    end
    throw(ArgumentError("minimal_polynomial: could not tell which factor vanishes at $e"))
end
_mp_select(c::Vector{Sym}, e::Sym) = c                  # already minimal (a rational, say)

# |c(v)| relative to the size of its terms
function _mp_residual(c::Vector{Sym}, v)
    s, m = zero(v), zero(real(v))
    for k in eachindex(c)
        t = BigFloat(Rational{BigInt}(_rat(c[k].num))) * v^(k - 1)
        s += t
        m += abs(t)
    end
    abs(s) / m
end

# integer coefficients, relatively prime, the leading one positive
function _mp_normal(c::Vector{Sym})
    rs = [Rational{BigInt}(_rat(a.num)) for a in c]
    L = lcm(map(denominator, rs))
    ns = [numerator(r * L) for r in rs]
    g = gcd(ns)
    ns = ns .÷ g
    ns[end] < 0 && (ns = -ns)
    Sym[Sym(n) for n in ns]
end
