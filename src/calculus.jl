"""
    diff(ex::Sym, x::Sym, n=1)
    diff(ex::Sym, x::Sym, y::Sym, zs::Sym...)

Derivative of `ex` with respect to the symbol `x`, `n` times. With several symbols, the
mixed derivative: `ex` differentiated with respect to `x`, then `y`, then each of `zs`,
so `diff(ex, x, y)` is `∂²ex/∂y∂x` and `diff(ex, x, x, y)` is `∂³ex/∂y∂x²`.

Functions without a registered derivative rule give an unevaluated
`Derivative(f(...), x)`. An unknown function of an expression follows the chain rule,
its derivative at a point written as `H'(g)` writes it: `diff(H(x^2), x)` is
`2*H'(x^2)*x`.

```julia
diff(x^3 * y^2, x, 2)        # 6*x*y^2
diff(x^3 * y^2, x, y)        # 6*x^2*y
diff(sin(x * y), x, y)       # -x*y*sin(x*y) + cos(x*y)
```
"""
function Base.diff(ex::Sym, x::Sym, n::Integer = 1)
    _isnamed(x) || throw(ArgumentError("can only differentiate with respect to a symbol or parameter, got $x"))
    for _ in 1:n
        ex = _diff(ex, x)
    end
    ex
end
Base.diff(ex::Number, x::Sym, n::Integer = 1) = ZERO

function Base.diff(ex::Sym, x::Sym, y::Sym, zs::Sym...)
    for v in (x, y, zs...)
        ex = diff(ex, v)
    end
    ex
end
Base.diff(ex::Number, x::Sym, y::Sym, zs::Sym...) = ZERO

# Derivative(f, v) evaluates when f has a derivative rule (after a substitution:
# Derivative(y(x), x) with y(x) => sin(x) is cos(x)); Derivative(f, v, a), the derivative
# at a, when the derivative has no Derivative left
function _derivative_rule(args::Vector{Sym})
    length(args) in (2, 3) || return nothing
    f, v = args[1], args[2]
    _isnamed(v) || return nothing
    d = _diff(f, v)
    length(args) == 2 && return _eq(d, _mkfn_node(:Derivative, Sym[f, v])) ? nothing : d
    _hasfn(d, :Derivative) ? nothing : subs(d, v => args[3])
end
register_function!(:Derivative; eval = _derivative_rule)

# An integral SymJ could not do (no integrator, or the integrator failed), unevaluated;
# its derivative is the integrand. The RischIntegrate extension sets the integrator.
const _INTEGRATOR = Ref{Any}(nothing)
function _integrate(f::Sym, x::Sym)
    iszero(f) && return ZERO
    has(f, x) || return f * x
    g = _INTEGRATOR[]
    if g !== nothing
        r = try g(f, x) catch err; err isa InterruptException && rethrow(); nothing end
        r isa Sym && return r
    end
    r = _integrate_linear_fn(f, x)       # before expanding: (2x + 1)^-2 stays a power
    r === nothing || return r
    # polynomials in x, term by term
    p = expand(f)
    terms = p.kind == ADD ? [(Sym(c), t) for (c, t) in zip(p.coefs, p.args)] : [(ONE, p)]
    out = p.kind == ADD ? Sym(p.num) * x : ZERO
    for (c, t) in terms
        k = _monomial_power(t, x)
        if k === nothing
            r = _integrate_linear_fn(t, x)
            r === nothing && return _mkfn_node(:Integral, Sym[f, x])
            out += c * r
            continue
        end
        cf, e = k
        out += c * (e == -1 ? cf * log(x) : cf * x^(e + 1) / (e + 1))
    end
    out
end
# cf h(u) for h exp, sin, cos, sinh or cosh and u linear in x: its antiderivative, or nothing
function _integrate_linear_fn(t::Sym, x::Sym)
    cf, h = ONE, nothing
    fs = t.kind == MUL ? [(b, p) for (b, p) in zip(t.args, t.exps)] : [(t, ONE)]
    t.kind == MUL && (cf = Sym(t.num))
    for (b, p) in fs
        if !has(b, x) && !has(p, x)
            cf *= b^p
        elseif h === nothing
            h = b^p
        else
            return nothing
        end
    end
    h === nothing && return nothing
    lin(u) = (d = diff(u, x); has(u, x) && !has(d, x) ? d : nothing)
    if h.kind == POW && _eq(h.args[1], E)
        a = lin(h.args[2])
        return a === nothing ? nothing : cf * h / a
    end
    if h.kind == POW && h.args[2].kind == NUMBER     # (a x + b)^k
        u, k = h.args[1], h.args[2]
        a = lin(u)
        a === nothing && return nothing
        return k == -1 ? cf * log(u) / a : cf * u^(k + 1) / ((k + 1) * a)
    end
    (h.kind == FN && length(h.args) == 1) || return nothing
    a = lin(h.args[1])
    a === nothing && return nothing
    u = h.args[1]
    h.name === :sin && return -cf * cos(u) / a
    h.name === :cos && return cf * sin(u) / a
    h.name === :sinh && return cf * cosh(u) / a
    h.name === :cosh && return cf * sinh(u) / a
    nothing
end
# t = cf * x^e, cf free of x and e an integer, or nothing
function _monomial_power(t::Sym, x::Sym)
    _eq(t, x) && return (ONE, 1)
    has(t, x) || return (t, 0)
    if t.kind == POW && _eq(t.args[1], x) && t.args[2].kind == NUMBER && _isint(t.args[2].num)
        return (ONE, _toint(t.args[2].num))
    end
    t.kind == MUL || return nothing
    cf, e = Sym(t.num), 0
    for (b, q) in zip(t.args, t.exps)
        if _eq(b, x) && q.kind == NUMBER && _isint(q.num)
            e += _toint(q.num)
        elseif has(b, x) || has(q, x)
            return nothing
        else
            cf *= b^q
        end
    end
    cf, e
end

function _diff(ex::Sym, x::Sym)
    k = ex.kind
    _isnamed(ex) && return _eq(ex, x) ? ONE : ZERO
    k == NUMBER && return ex.num.tag == T_FLOAT && isnan(ex.num.f) ? ex : ZERO    # NaN stays NaN
    k == CONSTANT && return ZERO
    if k == ADD
        return _add_many(_scale(_diff(ex.args[i], x), ex.coefs[i]) for i in eachindex(ex.args))
    elseif k == MUL
        terms = Sym[]
        for i in eachindex(ex.args)
            df = _diff(_mkpow(ex.args[i], ex.exps[i]), x)
            iszero(df) && continue
            rest = _mkmul(ex.num, _deleteat(ex.args, i), _deleteat(ex.exps, i))
            push!(terms, mul(rest, df))
        end
        return _add_many(terms)
    elseif k == POW
        b, e = ex.args
        db, de = _diff(b, x), _diff(e, x)
        if iszero(de)
            iszero(db) && return ZERO
            return mul(mul(e, pow(b, add(e, MINUS_ONE))), db)
        end
        return mul(ex, add(mul(de, log(b)), mul(mul(e, db), pow(b, MINUS_ONE))))
    elseif ex.name === :Integral && length(ex.args) == 2 && _eq(ex.args[2], x)
        return ex.args[1]
    elseif ex.name === :Integral && length(ex.args) == 4 && !has(ex.args[1], x)
        # Integral(f, s, a, b) = ∫ₐᵇ f ds: f(b) b' - f(a) a'
        f, v, lo, hi = ex.args
        return add(mul(subs(f, v => hi), _diff(hi, x)), mul(MINUS_ONE, mul(subs(f, v => lo), _diff(lo, x))))
    elseif _isderivat(ex) && (_eq(ex.args[2], x) || !has(ex.args[1], x))
        # F′(a) for F = F(t): F″(a) a′ (a constant unless the point has x)
        F, t, a = ex.args
        da = _diff(a, x)
        iszero(da) && return ZERO
        return mul(apply_function(:Derivative, Sym[_mkfn_node(:Derivative, Sym[F, t]), t, a]), da)
    elseif ex.name === :piecewise                 # piece by piece
        return apply_function(:piecewise, Sym[isodd(i) ? a : _diff(a, x) for (i, a) in enumerate(ex.args)])
    else
        rule = get(FUNCTION_RULES, ex.name, nothing)
        if rule === nothing && length(ex.args) == 1 && ex.args[1].kind != SYMBOL
            # an unknown function of an expression, by the chain rule: H(g)′ = H′(g) g′,
            # with H′(g) the derivative at a point, as H'(g) writes it
            a = ex.args[1]
            da = _diff(a, x)
            iszero(da) && return ZERO
            return mul(SymDerivative(SymFunction(ex.name), 1)(a), da)
        end
        terms = Sym[]
        for i in eachindex(ex.args)
            da = _diff(ex.args[i], x)
            iszero(da) && continue
            (rule === nothing || rule.diff === nothing) && return _mkfn_node(:Derivative, Sym[ex, x])
            d = rule.diff(ex.args, i)
            d === nothing && return _mkfn_node(:Derivative, Sym[ex, x])   # no closed form
            push!(terms, mul(d::Sym, da))
        end
        return _add_many(terms)
    end
end
