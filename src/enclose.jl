# Interval enclosures: bounds lo ≤ e(x) ≤ hi for x in [a, b], by interval arithmetic in
# Float64, each step widened outward by two ulps (a library function's error is below
# one). Enough to show that 3 + cos(3x)^5/2 - (1 + sin(2x)^3/2) ≥ 1 or that
# 1 + t^2 sin(t)^2 > 0 has no real zeros (SymPy #27283, #27139). nothing when some part
# of e has no rule, a symbol other than x included.

const _Enc = Tuple{Float64,Float64}

_widen((l, h)::_Enc) = (prevfloat(prevfloat(l)), nextfloat(nextfloat(h)))
_encok((l, h)::_Enc) = !(isnan(l) || isnan(h)) && l <= h

function _enc_mul((a, b)::_Enc, (c, d)::_Enc)
    ps = (a * c, a * d, b * c, b * d)
    any(isnan, ps) && return (-Inf, Inf)               # 0 * Inf: no bound
    (minimum(ps), maximum(ps))
end

function _enc_ipow((l, h)::_Enc, k::Int)
    k == 0 && return (1.0, 1.0)
    if k < 0
        (l <= 0 <= h) && return nothing                # 1/0 inside
        l, h = 1 / h, 1 / l
        k = -k
    end
    lk, hk = l^k, h^k
    iseven(k) || return (lk, hk)
    l <= 0 <= h ? (0.0, max(lk, hk)) : (min(lk, hk), max(lk, hk))
end

# sin over [l, h]: its values at the ends, and ±1 where a peak or trough is inside
function _enc_sin((l, h)::_Enc)
    (isfinite(l) && isfinite(h) && h - l < 2π) || return (-1.0, 1.0)
    has(t) = ceil((l - t) / 2π) <= floor((h - t) / 2π)
    sl, sh = sin(l), sin(h)
    (has(-π / 2) ? -1.0 : min(sl, sh), has(π / 2) ? 1.0 : max(sl, sh))
end

_enc_mono(f, (l, h)::_Enc) = (f(l), f(h))

function _enclose(e::Sym, x::Sym, I::_Enc)
    r = _enc(e, x, I)
    r === nothing || _encok(r) || return nothing
    r
end

function _enc(e::Sym, x::Sym, I::_Enc)
    k = e.kind
    if k == NUMBER || k == CONSTANT
        v = try N(e) catch; return nothing end
        v isa Real && isfinite(v) || return nothing
        return _widen((Float64(v), Float64(v)))
    end
    _isnamed(e) && return _eq(e, x) ? I : nothing
    if k == ADD
        s = _widen((Float64(e.num), Float64(e.num)))
        for i in eachindex(e.args)
            t = _enc(e.args[i], x, I)
            t === nothing && return nothing
            c = _widen((Float64(e.coefs[i]), Float64(e.coefs[i])))
            t = _enc_mul(c, t)
            s = _widen((s[1] + t[1], s[2] + t[2]))
        end
        return s
    elseif k == MUL
        p = _widen((Float64(e.num), Float64(e.num)))
        for i in eachindex(e.args)
            t = _enc_pow(e.args[i], e.exps[i], x, I)
            t === nothing && return nothing
            p = _widen(_enc_mul(p, t))
        end
        return p
    elseif k == POW
        return _enc_pow(e.args[1], e.args[2], x, I)
    elseif k == FN && length(e.args) == 1
        t = _enc(e.args[1], x, I)
        t === nothing && return nothing
        l, h = t
        n = e.name
        n === :sin && return _widen(_enc_sin(t))
        n === :cos && return _widen(_enc_sin((l + π / 2, h + π / 2)))
        n in (:atan, :sinh, :tanh, :asinh, :cbrt) && return _widen(_enc_mono(_julia_fn(n), t))
        n === :abs && return l >= 0 ? t : h <= 0 ? (-h, -l) : (0.0, max(-l, h))
        n === :cosh && return _widen(l >= 0 ? (cosh(l), cosh(h)) : h <= 0 ? (cosh(h), cosh(l)) : (1.0, max(cosh(l), cosh(h))))
        n === :log && return l > 0 ? _widen((log(l), log(h))) : nothing
        n === :sqrt && return l >= 0 ? _widen((sqrt(l), sqrt(h))) : nothing
    end
    nothing
end

_julia_fn(n::Symbol) = n === :atan ? atan : n === :sinh ? sinh : n === :tanh ? tanh :
    n === :asinh ? asinh : cbrt

function _enc_pow(b::Sym, p::Sym, x::Sym, I::_Enc)
    b === E && return (t = _enc(p, x, I); t === nothing ? nothing : _widen((exp(t[1]), exp(t[2]))))
    p.kind == NUMBER && isexact(p.num) || return nothing
    t = _enc(b, x, I)
    t === nothing && return nothing
    if _isint(p.num)
        k = _toint(p.num)
        abs(k) <= 64 || return nothing
        r = _enc_ipow(t, k)
        return r === nothing ? nothing : _widen(r)
    end
    # a rational power of a positive base
    t[1] > 0 || return nothing
    q = Float64(_jlnum(p.num))
    _widen(q > 0 ? (t[1]^q, t[2]^q) : (t[2]^q, t[1]^q))
end

# the sign of u on [a, b], 1 or -1, when its enclosure excludes 0; else nothing
function _sign_on(u::Sym, x::Sym, a::Real, b::Real)
    r = _enclose(u, x, (Float64(a), Float64(b)))
    r === nothing && return nothing
    r[1] > 0 ? 1 : r[2] < 0 ? -1 : nothing
end

# g with every condition, heaviside and dirac in x whose argument has one sign on [a, b]
# replaced by its value there (a piecewise of such conditions collapses to a piece)
function _settle_conditions(g::Sym, x::Sym, a::Real, b::Real)
    has(g, x) || return g
    g.kind == SYMBOL && return g
    args = Sym[_settle_conditions(t, x, a, b) for t in g.args]
    if g.kind == MUL
        exps = Sym[_settle_conditions(t, x, a, b) for t in g.exps]
        g = Sym(g.num) * prod((b^p for (b, p) in zip(args, exps)); init = ONE)
        return g
    end
    g = _rebuild_node(g, args)
    g.kind == FN || return g
    if haskey(_COND_RELS, g.name) && length(g.args) == 2 && has(g, x)
        s = _sign_on(g.args[1] - g.args[2], x, a, b)
        s === nothing && return g
        r = _relation_of(g)
        return _REL_INFO[r.op][2](s, 0) ? ONE : ZERO
    elseif g.name in (:heaviside, :dirac) && has(g.args[1], x)
        s = _sign_on(g.args[1], x, a, b)
        s === nothing && return g
        return g.name === :dirac ? ZERO : s > 0 ? ONE : ZERO
    end
    g
end
