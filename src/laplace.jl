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
is taken to be); `sin(a√t)`, `sin(a√t)/√t`, `cos(a√t)/√t`; unknown functions `u(t)` and
their derivatives (held as `LaplaceTransform(u(t), t, s)`), and integrals `∫₀ᵗ`, of
which a held `∫ u(t) dt` is taken to be one. Throws an `ArgumentError` for other terms.
Not exported.

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
        held = ZERO                                    # unknown functions: kept apart
        for (c, a) in zip(f.coefs, f.args)
            r = _lap_unknown(a, t, s)
            r === nothing ? (out += Sym(c) * _lap_term(a, t, s)) : (held += Sym(c) * r)
        end
        return expand(held) + _lap_tidy(out, s)
    end
    r = _lap_unknown(f, t, s)
    r === nothing ? _lap_tidy(_lap_term(f, t, s), s) : expand(r)
end

# exponentials of s kept apart, the rest cancelled
# (an exponential not linear in s, exp(-1/(4s)) from sin(√t): the whole tidied as it is)
function _lap_tidy(F::Sym, s::Sym)
    sh = try _lap_shifts(F, s) catch err; err isa ArgumentError || rethrow(); return _ode_tidy(F) end
    sum((_ode_tidy(cancel(R)) * exp(-c * s) for (c, R) in sh); init = ZERO)
end

# an unknown function u(t) and its derivatives: L{u⁽ⁿ⁾(t)} = sⁿ U(s) - Σ sⁿ⁻¹⁻ᵏ u⁽ᵏ⁾(0),
# U(s) the held LaplaceTransform(u(t), t, s) (SymPy #28720)
register_function!(:LaplaceTransform)

function _lap_unknown(term::Sym, t::Sym, s::Sym)
    fs, cf = _lap_factors(term)
    i = findfirst(((b, p),) -> has(b, t) || has(p, t), fs)
    i === nothing && return nothing
    any(j -> j != i && (has(fs[j][1], t) || has(fs[j][2], t)), eachindex(fs)) && return nothing
    g, p = fs[i]
    _eq(p, ONE) || return nothing
    # ∫₀ᵗ h: H(s)/s; a held ∫ h(t) dt is taken from 0 (the antiderivative vanishing
    # there, as for initial values; SymPy #24545)
    if g.kind == FN && g.name === :Integral
        rest = cf * prod((b^q for (j, (b, q)) in enumerate(fs) if j != i); init = ONE)
        if length(g.args) == 2 && _eq(g.args[2], t)
            return rest * _lap(g.args[1], t, s) / s
        elseif length(g.args) == 4 && iszero(g.args[3]) && _eq(g.args[4], t) && !_eq(g.args[2], t)
            return rest * _lap(subs(g.args[1], g.args[2] => t), t, s) / s
        end
        return nothing
    end
    n = 0
    while g.kind == FN && g.name === :Derivative && length(g.args) == 2 && _eq(g.args[2], t)
        g, n = g.args[1], n + 1
    end
    (g.kind == FN && !haskey(FUNCTION_RULES, g.name) && length(g.args) == 1 && _eq(g.args[1], t)) || return nothing
    U = apply_function(:LaplaceTransform, Sym[g, t, s])
    out = s^n * U
    d = g
    for k in 0:(n - 1)
        out -= s^(n - 1 - k) * subs(d, t => ZERO)
        d = diff(d, t)
    end
    # the t-free factors too: k₁ u(t) is k₁ U(s) (#25815)
    cf * prod((b^q for (j, (b, q)) in enumerate(fs) if j != i); init = ONE) * out
end

function _lap_term(term::Sym, t::Sym, s::Sym)
    r = _lap_unknown(term, t, s)
    r === nothing || return r
    fs, cf = _lap_factors(term)
    # δ⁽ᵏ⁾(d(t - a)): L{f δ⁽ᵏ⁾(t - a)} = (-1)ᵏ dᵏ/dtᵏ [f e^(-st)] at t = a, over dᵏ|d| (SymPy #26775)
    for (i, (b, p)) in enumerate(fs)
        (b.kind == FN && b.name === :dirac && length(b.args) == 2 && _eq(p, ONE) &&
         b.args[2].kind == NUMBER && _isint(b.args[2].num)) || continue
        k = _toint(b.args[2].num)
        u = b.args[1]
        d = diff(u, t)
        (has(d, t) || iszero(d)) && continue
        a = _ode_tidy(-subs(u, t => ZERO) / d)
        a.kind == NUMBER && a.num < N0 && return ZERO
        rest = cf * prod((bb^pp for (j, (bb, pp)) in enumerate(fs) if j != i); init = ONE)
        v = subs(diff(rest * exp(-s * t), t, k), t => a)
        return expand((-1)^k * v / (d^k * abs(d)))
    end
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
            return _lap_value_at(rest, t, a) * exp(-a * s) / abs(d)
        end
        if d.kind == NUMBER && d.num < N0              # H(-(t - a)) = 1 - H(t - a)
            return _lap(rest, t, s) - _lap(rest * heaviside(t - a), t, s)
        end
        neg && return _lap(rest, t, s)                 # H(t - a) = 1 on t ≥ 0
        return exp(-a * s) * _lap(subs(rest, t => t + a), t, s)
    end
    # sin(a√t), sin(a√t)/√t and cos(a√t)/√t (SymPy #24309)
    r = _lap_sqrt_trig(fs, cf, t, s)
    r === nothing || return r
    # sin(a t) e^(αt)/t: atan(a/(s - α)) (L{f/t} = ∫ₛ^∞ F; SymPy #21372)
    r = _lap_sinc(fs, cf, t, s)
    r === nothing || return r
    # tᵏ e^(αt) terms, sin and cos as exponentials
    e = expand(_ode_toexp(term))
    terms = e.kind == ADD ? [(Sym(c), x) for (c, x) in zip(e.coefs, e.args)] : [(ONE, e)]
    out = e.kind == ADD && !iszero(e.num) ? Sym(e.num) / s : ZERO
    for (c, x) in terms
        out += c * _lap_expoly(x, t, s)
    end
    _ode_tidy(expand(out))
end

function _lap_sqrt_trig(fs, cf::Sym, t::Sym, s::Sym)
    i = findfirst(((b, p),) -> b.kind == FN && b.name in (:sin, :cos) && _eq(p, ONE) && has(b, t), fs)
    i === nothing && return nothing
    g = fs[i][1]
    a = cancel(g.args[1] / sqrt(t))
    has(a, t) && return nothing
    rest = cf * prod((b^p for (j, (b, p)) in enumerate(fs) if j != i); init = ONE)
    k = cancel(rest * sqrt(t))                          # rest = k/√t?
    q = exp(-a^2 / (4s))
    if g.name === :sin && !has(rest, t)
        return rest * a * sqrt(PI) / (2 * s^(Sym(3) // 2)) * q
    elseif !has(k, t)
        g.name === :cos && return k * sqrt(PI / s) * q
        return k * PI * apply_function(:erf, Sym[a / (2sqrt(s))])
    end
    nothing
end

function _lap_sinc(fs, cf::Sym, t::Sym, s::Sym)
    i = findfirst(((b, p),) -> _eq(b, t) && _eq(p, MINUS_ONE), fs)
    i === nothing && return nothing
    j = findfirst(((b, p),) -> b.kind == FN && b.name === :sin && _eq(p, ONE) && has(b, t), fs)
    j === nothing && return nothing
    a = diff(fs[j][1].args[1], t)
    (has(a, t) || !iszero(expand(subs(fs[j][1].args[1], t => ZERO)))) && return nothing
    c, α = cf, ZERO
    for (k, (b, p)) in enumerate(fs)
        (k == i || k == j) && continue
        if !has(b, t) && !has(p, t)
            c *= b^p
        elseif b === E && (d = diff(p, t); !has(d, t))
            c *= exp(subs(p, t => ZERO))
            α += d
        else
            return nothing
        end
    end
    c * atan(a / (s - α))
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
        elseif _eq(b, t) && !has(q, t)
            p += q
        elseif b === E && (d = diff(q, t); !has(d, t))
            cf *= exp(subs(q, t => ZERO))
            α += d
        elseif !has(b, t) && (d = diff(q, t); !has(d, t))   # bᵗ = e^(t log b) (SymPy #28136)
            cf *= b^subs(q, t => ZERO)
            α += d * log(b)
        else
            throw(ArgumentError("laplace_transform: no transform found for $(term)"))
        end
    end
    # (a symbolic power: Γ(p + 1)/(s - α)^(p + 1), for p > -1; SymPy #21372)
    p.kind == NUMBER && !(p.num > -N1) && throw(ArgumentError("laplace_transform: t^$p has no transform"))
    g = p.kind == NUMBER && isinteger(p) ? Sym(factorial(big(_toint(p.num)))) : apply_function(:gamma, Sym[p + 1])
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
gives [`SymJ.dirac`](@ref)`(t)`. A factor `exp(c s)` with `c > 0` (or `exp(a s)`, a
symbol taken to be positive) is a shift to `t < 0`, the transform of no function of
`t ≥ 0`. Throws an `ArgumentError` otherwise. Not exported.

```julia
@syms t s
SymJ.inverse_laplace_transform(1 / (s^2 + 2s + 5), s, t)     # exp(-t)*sin(2t)/2
SymJ.inverse_laplace_transform(exp(-2s) / s^2, s, t)         # (t - 2)*heaviside(t - 2)
```
"""
function inverse_laplace_transform(F, s::Sym, t::Sym)
    (_isnamed(t) && _isnamed(s)) || throw(ArgumentError("inverse_laplace_transform: s and t are symbols"))
    F = Sym(F)
    # floating point coefficients: the transform of the rationals they stand for, its
    # numbers floats again (apart needs exact coefficients; SymPy #26517)
    if !_isexactex(F) && !_float_exponent(F)
        # floats and no parameters: partial fractions numerically, simple poles only
        # (exactly, a cubic's roots in Cardano's form over 50-digit rationals took a
        # minute; SymPy #25352)
        r = _ilap_float(F, s, t)
        r === nothing || return r
        X = _exactfloats(F)
        _isexactex(X) && return _floatify(inverse_laplace_transform(X, s, t))
    end
    # held transforms of unknown functions: U(s) = LaplaceTransform(u(t), t, s) back to u(t)
    held = _held_transforms(F, s)
    if !isempty(held)
        out = ZERO
        rest = F
        for (U, g) in held
            c = coeff(expand(F), U)
            rest = expand(rest - c * U)
            # c(s) a polynomial in s: sᵏ U(s) = L{u⁽ᵏ⁾} + Σ sᵏ⁻¹⁻ʲ u⁽ʲ⁾(0)
            P = try Poly{Sym}(expand(c), [s]) catch err
                err isa ArgumentError || rethrow()
                throw(ArgumentError("inverse_laplace_transform: $U times $c is not supported"))
            end
            u = subs(g, g.args[1] => t)
            for (ev, ck) in zip(exponent_vectors(P), coefficients(P))
                k = ev[1]
                out += ck * (k == 0 ? u : diff(u, t, k))
                d = u
                for j in 0:(k - 1)
                    rest = expand(rest + ck * s^(k - 1 - j) * subs(d, t => ZERO))
                    d = diff(d, t)
                end
            end
        end
        has(rest, s) || iszero(rest) || throw(ArgumentError("inverse_laplace_transform: no inverse found for $rest"))
        return out + (iszero(rest) ? ZERO : inverse_laplace_transform(rest, s, t))
    end
    out = ZERO
    for (c, R) in _lap_shifts(F, s)
        # e^(-cs) with c < 0 is the transform of no function of t ≥ 0 (SymPy #24595); a
        # symbol c is taken to be ≥ 0, as in laplace_transform, so -a is refused too
        _ilap_negshift(c) && throw(ArgumentError("inverse_laplace_transform: exp($(-c)*s) is a shift to t < 0; no inverse for t ≥ 0"))
        g = _ilap_special(R, s, t)
        out += iszero(c) ? g : heaviside(t - c) * subs(g, t => t - c)
    end
    out
end

# terms c atan(a/(s - α)) (e^(αt) sin(at)/t) and c (s - α)^(-m), m symbolic
# (t^(m-1) e^(αt) c/Γ(m), the gamma functions simplified: factorial(n)/(s + a)^(n + 1)
# is tⁿ e^(-at); SymPy #21372); the rest a rational function
function _ilap_special(R::Sym, s::Sym, t::Sym)
    sympow(b, q) = has(b, s) && _hassymbol(q) && !has(q, s)
    special(x) = _hasfn(x, :atan) || _has_node(x, u -> (u.kind == POW && sympow(u.args...)) ||
                                                       (u.kind == MUL && any(i -> sympow(u.args[i], u.exps[i]), eachindex(u.args))))
    special(R) || return _ilap_rational(R, s, t)
    E1 = expand(R)
    terms = E1.kind == ADD ? [Sym(c) * x for (c, x) in zip(E1.coefs, E1.args)] : [E1]
    E1.kind == ADD && !iszero(E1.num) && push!(terms, Sym(E1.num))
    out, rest = ZERO, ZERO
    for x in terms
        g = special(x) ? _ilap_special_term(x, s, t) : nothing
        g === nothing ? (rest += x) : (out += g)
    end
    out + (iszero(rest) ? ZERO : _ilap_rational(rest, s, t))
end
function _ilap_special_term(x::Sym, s::Sym, t::Sym)
    fs, cf = _lap_factors(x)
    i = findfirst(((b, p),) -> has(b, s) || has(p, s), fs)
    i === nothing && return nothing
    any(k -> k != i && (has(fs[k][1], s) || has(fs[k][2], s)), eachindex(fs)) && return nothing
    c = cf * prod((b^p for (k, (b, p)) in enumerate(fs) if k != i); init = ONE)
    b, p = fs[i]
    lin(u) = (k1 = diff(u, s); has(k1, s) || iszero(k1) ? nothing : (k1, expand(subs(u, s => ZERO))))
    if b.kind == FN && b.name === :atan && length(b.args) == 1 && _eq(p, ONE)
        l = lin(cancel(1 / b.args[1]))                # 1/u = (s - α)/a
        l === nothing && return nothing
        k1, k0 = l
        a, α = 1 / k1, -k0 / k1
        return c * exp(α * t) * sin(a * t) / t
    end
    (has(p, s) || !_hassymbol(p)) && return nothing
    l = lin(b)
    l === nothing && return nothing
    k1, k0 = l                                        # b = k1 (s - α)
    α, m = -k0 / k1, -p
    c = c * k1^p
    q = c / apply_function(:gamma, Sym[m])
    q = try simplify_gamma(q) catch err; err isa InterruptException && rethrow(); q end
    q * t^(m - 1) * exp(α * t)
end

function _ilap_float(F::Sym, s::Sym, t::Sym)
    all(v -> _eq(v, s), free_symbols(F)) || return nothing
    _has_epow(F) && return nothing
    n, d = numer_denom(cancel(F))
    P = try (Poly{SymNumber}(_exactfloats(n), [s]), Poly{SymNumber}(_exactfloats(d), [s])) catch; return nothing end
    nc = [ComplexF64(_evalc(c)) for c in _xcoeffs(P[1])]
    dc = [ComplexF64(_evalc(c)) for c in _xcoeffs(P[2])]
    (length(dc) >= 2 && length(nc) < length(dc)) || return nothing
    ps = [ComplexF64(_evalc(r)) for r in _numeric_roots(Sym[Sym(real(c)) for c in dc])]
    scale = maximum(abs, ps; init = 1.0)
    for i in eachindex(ps), j in (i + 1):length(ps)
        abs(ps[i] - ps[j]) > 1e-6 * scale || return nothing            # a repeated pole
    end
    ev(c, z) = foldr((a, acc) -> acc * z + a, c; init = zero(ComplexF64))
    dd = [k * dc[k + 1] for k in 1:(length(dc) - 1)]
    out = ZERO
    for p in ps
        imag(p) < -1e-12 * scale && continue                             # its conjugate's term has it
        r = ev(nc, p) / ev(dd, p)
        if abs(imag(p)) <= 1e-12 * scale
            out += Sym(real(r)) * exp(Sym(real(p)) * t)
        else
            out += exp(Sym(real(p)) * t) * (Sym(2real(r)) * cos(Sym(imag(p)) * t) - Sym(2imag(r)) * sin(Sym(imag(p)) * t))
        end
    end
    out
end

# g(a), or its limit there for a removable singularity (sin(4π(t - 1))/(t - 1) at 1;
# SymPy #25293)
function _lap_value_at(g::Sym, t::Sym, a::Sym)
    v = try subs(g, t => a) catch err; _isdivzero(err) || rethrow(); nothing end
    v === nothing || !_finite_value(v) || return v
    l = limit(g, t => a)
    (l isa Sym && _finite_value(l)) || throw(ArgumentError("laplace_transform: $g is singular at the dirac's point $a"))
    l
end

function _ilap_negshift(c::Sym)
    iszero(c) && return false
    A = _ASSUMING[]
    if A !== nothing
        _nonneg(c, A) && return false
        _nonpos(c, A) && return true
    end
    _hassymbol(c) || return _isreal_const(c) === true && N(c) < 0
    _isneg(c)
end

function _held_transforms(F::Sym, s::Sym)
    out = Pair{Sym,Sym}[]
    _replace_walk(F) do e
        e.kind == FN && e.name === :LaplaceTransform && length(e.args) == 3 && _eq(e.args[3], s) &&
            !any(p -> _eq(p.first, e), out) && push!(out, e => e.args[1])
        nothing
    end
    out
end

function _ilap_rational(R::Sym, s::Sym, t::Sym)
    R = cancel(R)
    n, d = numer_denom(R)
    # a polynomial in s: Σ cₖ sᵏ is Σ cₖ δ⁽ᵏ⁾(t) (SymPy #26775)
    has(d, s) || return _ilap_poly(R, s, t)
    A = apart(R, s)
    terms = A.kind == ADD ? [Sym(c) * a for (c, a) in zip(A.coefs, A.args)] : [A]
    A.kind == ADD && !iszero(A.num) && push!(terms, Sym(A.num))
    out = ZERO
    for term in terms
        tn, td = numer_denom(term)
        if !has(td, s)
            out += _ilap_poly(term, s, t)
            continue
        end
        out += _ilap_part(tn, td, s, t)
    end
    _lap_neat(expand(_ode_cis(expand(out))))
end

function _ilap_poly(P::Sym, s::Sym, t::Sym)
    Q = Poly{Sym}(expand(P), [s])
    out = ZERO
    for (ev, c) in zip(exponent_vectors(Q), coefficients(Q))
        out += Sym(c) * (ev[1] == 0 ? dirac(t) : apply_function(:dirac, Sym[t, Sym(ev[1])]))
    end
    out
end

# a sum of tᵏ e^(αt) cos/sin terms: grouped by exponential, without simplify's
# rewriting of sin(2t)/2 as sin(t) cos(t)
function _lap_neat(e::Sym)
    e = _explog(e)
    # roots of a cubic or quartic (cube roots of radicals): cancel over all those
    # generators did not finish in minutes (SymPy #26560, #25796); left as they are
    _has_high_root(e) && return e
    for c in (cancel(e), simplify_hyperbolic(e))
        _nodes(c) < _nodes(e) && (e = c)
    end
    e
end

# a power with an exponent whose denominator is 3 or more (u^(1/3), u^(2/3))
_has_high_root(e::Sym) =
    (e.kind == POW && e.args[2].kind == NUMBER && isexact(e.args[2].num) && _denominator(e.args[2].num) >= 3) ||
    (e.kind == MUL && any(p -> p.kind == NUMBER && isexact(p.num) && _denominator(p.num) >= 3, e.exps)) ||
    any(_has_high_root, e.args) || (e.kind == MUL && any(_has_high_root, e.exps))

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
    ρs = degree(P, s) == 2 ? _ilap_quadratic_roots(q, s) : _solve_generic(q ~ 0, s)
    length(ρs) == degree(P, s) || throw(ArgumentError("inverse_laplace_transform: the roots of $q were not found"))
    out = ZERO
    for ρ in ρs
        G = N * exp(s * t) / (lead * lq^k * prod(((s - σ)^k for σ in ρs if !_eq(σ, ρ)); init = ONE))
        out += subs(diff(G, s, k - 1), s => ρ) / factorial(big(k - 1))
    end
    out
end

# the roots α ± iβ of a quadratic, β = sqrt(c/a - α²) with the parameters taken to be
# positive: s² + w² has roots ±iw, not ±sqrt(-w²). Either square root gives the same
# function of t (a conjugate pair contributes e^(αt)(B cos βt + (C + Bα)/β sin βt), even
# in β), and β = i sqrt(-β²) when that is simpler (real roots: s² - w², ±w).
function _ilap_quadratic_roots(q::Sym, s::Sym)
    c2, c1, c0 = (_ode_coeff(q, s, k) for k in 2:-1:0)
    α = _ode_tidy(-c1 / (2c2))
    β2 = _ode_tidy(c0 / c2 - α^2)
    ps = free_symbols(β2)
    A = isempty(ps) ? Assumptions() : Assumptions((p ≫ 0 for p in ps)...)
    r1 = refine(sqrt(β2); assume = A)
    r2 = IM * refine(sqrt(-β2); assume = A)
    β = _nodes(r2) < _nodes(r1) ? r2 : r1
    Sym[_ode_tidy(α + IM * β), _ode_tidy(α - IM * β)]
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
