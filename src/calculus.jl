"""
    diff(ex::Sym, x::Sym, n=1)
    diff(ex::Sym, x::Sym, y::Sym, zs::Sym...)

Derivative of `ex` with respect to the symbol `x`, `n` times. With several symbols, the
mixed derivative: `ex` differentiated with respect to `x`, then `y`, then each of `zs`,
so `diff(ex, x, y)` is `∂²ex/∂y∂x` and `diff(ex, x, x, y)` is `∂³ex/∂y∂x²`.

Functions without a registered derivative rule give an unevaluated
`Derivative(f(...), x)`.

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

function _diff(ex::Sym, x::Sym)
    k = ex.kind
    _isnamed(ex) && return _eq(ex, x) ? ONE : ZERO
    (k == NUMBER || k == CONSTANT) && return ZERO
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
    elseif ex.name === :piecewise                 # piece by piece
        return apply_function(:piecewise, Sym[isodd(i) ? a : _diff(a, x) for (i, a) in enumerate(ex.args)])
    else
        rule = get(FUNCTION_RULES, ex.name, nothing)
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
