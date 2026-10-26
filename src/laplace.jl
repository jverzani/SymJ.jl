# Laplace transforms: SymJ.laplace_transform and SymJ.inverse_laplace_transform, and the
# Laplace method of SymJ.dsolve for linear equations with constant coefficients.
#
# Forward: a sum term by term. dirac(t - a) g(t) is e^(-as) g(a), heaviside(t - a) g(t)
# is e^(-as) L{g(t + a)} (a ≥ 0); otherwise sin, cos, sinh, cosh are written as
# exponentials and each term c tⁿ e^(αt) gives c n!/(s - α)ⁿ⁺¹, the sum cancelled back to
# a real rational function. A term tᵖ (p not an integer, p > -1) gives Γ(p + 1)/sᵖ⁺¹.
#
# Inverse: the terms grouped by a factor e^(-cs) (heaviside(t - c) g(t - c)); each group,
# a rational function, split by apart into N/qᵏ, q irreducible, and inverted by the
# residues of N e^(st)/qᵏ at the roots of q, then made real. A constant left over is a
# dirac(t).

"""
    SymJ.laplace_transform(f, t, s)

The Laplace transform `∫₀^∞ e^(-st) f(t) dt` of `f`: for sums of products of powers of
`t`, exponentials, sines, cosines, `sinh`, `cosh`, and their shifts by
[`SymJ.heaviside`](@ref)`(t - a)`, and [`SymJ.dirac`](@ref)`(t - a)` (`a ≥ 0`; a symbol `a`
is taken to be). Throws an `ArgumentError` for other terms. Not exported.

```julia
@syms t s
SymJ.laplace_transform(t * exp(-2t), t, s)            # 1/(s + 2)^2
SymJ.laplace_transform(sin(3t), t, s)                 # 3/(s^2 + 9)
SymJ.laplace_transform(SymJ.heaviside(t - 1), t, s)   # exp(-s)/s
```
"""
function laplace_transform(f, t::Sym, s::Sym)
    (_isnamed(t) && _isnamed(s)) || throw(ArgumentError("laplace_transform: t and s are symbols"))
    _lap(Sym(f), t, s)
end

function _lap(f::Sym, t::Sym, s::Sym)
    has(f, t) || return f / s
    f = expand(f)
    if f.kind == ADD
        out = iszero(f.num) ? ZERO : Sym(f.num) / s
        for (c, a) in zip(f.coefs, f.args)
            out += Sym(c) * _lap_term(a, t, s)
        end
        return _lap_tidy(out, s)
    end
    _lap_tidy(_lap_term(f, t, s), s)
end

# exponentials of s kept apart, the rest cancelled
_lap_tidy(F::Sym, s::Sym) = sum((_ode_tidy(cancel(R)) * exp(-c * s) for (c, R) in _lap_shifts(F, s)); init = ZERO)

function _lap_term(term::Sym, t::Sym, s::Sym)
    fs, cf = _lap_factors(term)
    for (i, (b, p)) in enumerate(fs)
        b.kind == FN && b.name in (:dirac, :heaviside) && length(b.args) == 1 && _eq(p, ONE) || continue
        u = b.args[1]
        d = diff(u, t)
        (has(d, t) || iszero(d)) && continue
        a = _ode_tidy(-subs(u, t => ZERO) / d)          # u = d (t - a)
        rest = cf * prod((bb^pp for (j, (bb, pp)) in enumerate(fs) if j != i); init = ONE)
        neg = a.kind == NUMBER && a.num < N0
        if b.name === :dirac                           # δ(d(t - a)) = δ(t - a)/|d|
            neg && return ZERO
            return subs(rest, t => a) * exp(-a * s) / abs(d)
        end
        if d.kind == NUMBER && d.num < N0              # H(-(t - a)) = 1 - H(t - a)
            return _lap(rest, t, s) - _lap(rest * heaviside(t - a), t, s)
        end
        neg && return _lap(rest, t, s)                 # H(t - a) = 1 on t ≥ 0
        return exp(-a * s) * _lap(subs(rest, t => t + a), t, s)
    end
    # tᵏ e^(αt) terms, sin and cos as exponentials
    e = expand(_ode_toexp(term))
    terms = e.kind == ADD ? [(Sym(c), x) for (c, x) in zip(e.coefs, e.args)] : [(ONE, e)]
    out = e.kind == ADD && !iszero(e.num) ? Sym(e.num) / s : ZERO
    for (c, x) in terms
        out += c * _lap_expoly(x, t, s)
    end
    _ode_tidy(expand(out))
end

# a term's factors as (base, exponent), and its coefficient
_lap_factors(e::Sym) = e.kind == MUL ? (collect(zip(e.args, e.exps)), Sym(e.num)) :
    e.kind == POW ? ([(e.args[1], e.args[2])], ONE) : ([(e, ONE)], ONE)

# L{cf tᵖ e^(αt)} = cf Γ(p + 1)/(s - α)ᵖ⁺¹
function _lap_expoly(term::Sym, t::Sym, s::Sym)
    p, α = ZERO, ZERO
    fs, cf = _lap_factors(term)
    for (b, q) in fs
        if !has(b, t) && !has(q, t)
            cf *= b^q
        elseif _eq(b, t) && q.kind == NUMBER
            p += q
        elseif b === E && (d = diff(q, t); !has(d, t))
            cf *= exp(subs(q, t => ZERO))
            α += d
        else
            throw(ArgumentError("laplace_transform: no transform found for $(term)"))
        end
    end
    p.kind == NUMBER && p.num > -N1 || throw(ArgumentError("laplace_transform: t^$p has no transform"))
    g = isinteger(p) ? Sym(factorial(big(_toint(p.num)))) : apply_function(:gamma, Sym[p + 1])
    cf * g / (s - α)^(p + 1)
end

# F as Σ R_c(s) e^(-cs): (c, R_c), the exponent of each e^(...) factor linear in s
function _lap_shifts(F::Sym, s::Sym)
    F = expand(F)
    groups = Dict{Sym,Sym}()
    terms = F.kind == ADD ? [Sym(c) * a for (c, a) in zip(F.coefs, F.args)] : [F]
    F.kind == ADD && !iszero(F.num) && push!(terms, Sym(F.num))
    for term in terms
        fs, cf = _lap_factors(term)
        c, rest = ZERO, cf
        for (b, q) in fs
            if b === E && has(q, s)
                d = diff(q, s)
                has(d, s) && throw(ArgumentError("inverse_laplace_transform: exp($q) is not a shift"))
                c -= d
                rest *= exp(subs(q, s => ZERO))
            else
                rest *= b^q
            end
        end
        c = _ode_tidy(c)
        groups[c] = get(groups, c, ZERO) + rest
    end
    sort!(collect(groups); by = p -> string(p[1]))
end

"""
    SymJ.inverse_laplace_transform(F, s, t)

The function of `t ≥ 0` whose Laplace transform is `F`: for `F` a sum of rational
functions of `s` times factors `exp(-c s)` (giving [`SymJ.heaviside`](@ref)`(t - c)`
shifts). A rational function is split into partial fractions and each inverted by
residues at the roots of its denominator, so the roots must be found; a constant term
gives [`SymJ.dirac`](@ref)`(t)`. Throws an `ArgumentError` otherwise. Not exported.

```julia
@syms t s
SymJ.inverse_laplace_transform(1 / (s^2 + 2s + 5), s, t)     # exp(-t)*sin(2t)/2
SymJ.inverse_laplace_transform(exp(-2s) / s^2, s, t)         # (t - 2)*heaviside(t - 2)
```
"""
function inverse_laplace_transform(F, s::Sym, t::Sym)
    (_isnamed(t) && _isnamed(s)) || throw(ArgumentError("inverse_laplace_transform: s and t are symbols"))
    out = ZERO
    for (c, R) in _lap_shifts(Sym(F), s)
        g = _ilap_rational(R, s, t)
        out += iszero(c) ? g : heaviside(t - c) * subs(g, t => t - c)
    end
    out
end

function _ilap_rational(R::Sym, s::Sym, t::Sym)
    R = cancel(R)
    n, d = numer_denom(R)
    has(d, s) || (has(n, s) ? throw(ArgumentError("inverse_laplace_transform: $R is a polynomial in $s, not a transform")) :
                  return R * dirac(t))
    A = apart(R, s)
    terms = A.kind == ADD ? [Sym(c) * a for (c, a) in zip(A.coefs, A.args)] : [A]
    A.kind == ADD && !iszero(A.num) && push!(terms, Sym(A.num))
    out = ZERO
    for term in terms
        tn, td = numer_denom(term)
        if !has(td, s)
            has(tn, s) && throw(ArgumentError("inverse_laplace_transform: $R is not proper in $s"))
            out += term * dirac(t)
            continue
        end
        out += _ilap_part(tn, td, s, t)
    end
    _lap_neat(expand(_ode_cis(expand(out))))
end

# a sum of tᵏ e^(αt) cos/sin terms: grouped by exponential, without simplify's
# rewriting of sin(2t)/2 as sin(t) cos(t)
function _lap_neat(e::Sym)
    e = _explog(e)
    for c in (cancel(e), simplify_hyperbolic(e))
        _nodes(c) < _nodes(e) && (e = c)
    end
    e
end

# N/D, D = lead qᵏ with q irreducible: Σ over the roots ρ of q of the residue of
# N e^(st)/D, (1/(k-1)!) d^(k-1)/ds^(k-1) [N e^(st) / (lead Π_(σ≠ρ) (s - σ))ᵏ] at ρ
function _ilap_part(N::Sym, D::Sym, s::Sym, t::Sym)
    q, k = nothing, 0
    for (b, p) in _ode_factors(factor(D))
        has(b, s) || continue
        q === nothing || throw(ArgumentError("inverse_laplace_transform: unexpected denominator $D"))
        q, k = b, _toint(p.num)
    end
    lead = cancel(D / q^k)                           # (with the numeric coefficient)
    P = Poly{Sym}(expand(q), [s])
    lq = Sym(coefficients(P)[1])                     # q = lq Π (s - ρ)
    ρs = solve(q ~ 0, s)
    length(ρs) == degree(P, s) || throw(ArgumentError("inverse_laplace_transform: the roots of $q were not found"))
    out = ZERO
    for ρ in ρs
        G = N * exp(s * t) / (lead * lq^k * prod(((s - σ)^k for σ in ρs if !_eq(σ, ρ)); init = ONE))
        out += subs(diff(G, s, k - 1), s => ρ) / factorial(big(k - 1))
    end
    out
end

## dsolve by Laplace transforms

# a linear equation with constant coefficients Σ aₖ y⁽ᵏ⁾ = g: Y(s) from the transformed
# equation, with y⁽ʲ⁾(a) the initial conditions (constants for those not given), and
# t shifted to the point a
function _ode_laplace(P::_ODE, Cs, conds)
    lc = _ode_linear_coefs(P)
    lc === nothing && return nothing
    A, g = lc
    x = P.x
    any(c -> has(c, x), A) && return nothing
    a = isempty(conds) ? ZERO : conds[1][2]
    all(c -> _eq(c[2], a), conds) || return nothing
    vals = Sym[Cs[j + 1] for j in 0:(P.n - 1)]
    for (o, _, v) in conds
        vals[o + 1] = v
    end
    s = Sym(:_ode_s)
    G = laplace_transform(subs(g, x => x + a), x, s)
    num = G + sum(A[k + 1] * sum((s^(k - 1 - j) * vals[j + 1] for j in 0:(k - 1)); init = ZERO) for k in 1:P.n)
    Y = num / sum(A[k + 1] * s^k for k in 0:P.n)
    y = inverse_laplace_transform(Y, s, x)
    _ode_tidy(subs(y, x => x - a))
end
