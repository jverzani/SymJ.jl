"""
    diff(ex::Sym, x::Sym, n=1)

Derivative of `ex` with respect to the symbol `x`. Functions without a registered
derivative rule give an unevaluated `Derivative(f(...), x)`.
"""
function Base.diff(ex::Sym, x::Sym, n::Integer = 1)
    x.kind == SYMBOL || throw(ArgumentError("can only differentiate with respect to a symbol, got $x"))
    for _ in 1:n
        ex = _diff(ex, x)
    end
    ex
end
Base.diff(ex::Number, x::Sym, n::Integer = 1) = ZERO

function _diff(ex::Sym, x::Sym)
    k = ex.kind
    k == SYMBOL && return ex === x || ex.name === x.name ? ONE : ZERO
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
    else
        rule = get(FUNCTION_RULES, ex.name, nothing)
        terms = Sym[]
        for i in eachindex(ex.args)
            da = _diff(ex.args[i], x)
            iszero(da) && continue
            (rule === nothing || rule.diff === nothing) && return _mkfn_node(:Derivative, Sym[ex, x])
            push!(terms, mul(rule.diff(ex.args, i)::Sym, da))
        end
        return _add_many(terms)
    end
end
