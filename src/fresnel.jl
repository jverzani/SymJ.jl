# The Fresnel integrals S(z) = ∫₀ᶻ sin(πt²/2) dt and C(z) = ∫₀ᶻ cos(πt²/2) dt (SymPy's
# fresnels, fresnelc). Entire and odd; S(±∞) = C(±∞) = ±1/2. They give ∫ sin(ax² + bx + c)
# and ∫ cos(ax² + bx + c) (RischIntegrate's answers are erf of (-1)^(1/4) multiples).

import SpecialFunctions

"""
    SymJ.fresnels(z)
    SymJ.fresnelc(z)

The Fresnel integrals S(z) = ∫₀ᶻ sin(πt²/2) dt and C(z) = ∫₀ᶻ cos(πt²/2) dt: of a `Sym`, the
function (exact at 0 and ±∞, odd); of a number, its value (a complex number for a complex
`z`). Integrals of sin and cos of quadratics are given in them:
`integrate(sin(x^2), x)` is `√(π/2) S(√(2/π) x)`. Not exported.
"""
fresnels(z::Sym) = apply_function(:fresnels, Sym[z])
fresnelc(z::Sym) = apply_function(:fresnelc, Sym[z])
fresnels(z::Number) = _fresnel_num(true, z)
fresnelc(z::Number) = _fresnel_num(false, z)

# S (s = true) or C at z. Float64: by the power series for |z| ≤ 1 (the erf form loses
# S(z) ≈ πz³/6 to cancellation there), else by erf: C + iS = (1 + i)/2 erf((1 - i)√π z/2)
# for real z, and for complex z
#   S = (1 + i)/4 [erf(a z) - i erf(ā z)],  C = (1 - i)/4 [erf(a z) + i erf(ā z)],  a = (1 + i)√π/2.
# BigFloat: by the power series, with guard bits for the terms' growth, e^(π|z|²/2)
function _fresnel_num(s::Bool, z::Number)
    if z isa Union{BigFloat,Complex{BigFloat}}
        return _fresnel_series(s, z)
    end
    w = z isa Real ? Float64(z) : ComplexF64(z)
    isfinite(w) || (w isa Real && return copysign(0.5, w))
    abs(w) <= 1 && return _fresnel_series(s, w)
    a = (1 + im) * sqrt(π) / 2
    e1, e2 = SpecialFunctions.erf(a * w), SpecialFunctions.erf(conj(a) * w)
    v = s ? (1 + im) / 4 * (e1 - im * e2) : (1 - im) / 4 * (e1 + im * e2)
    w isa Real ? real(v) : v
end

# Σ (-1)ⁿ (π/2)^(2n+1) z^(4n+3) / ((2n+1)! (4n+3)) for S, Σ (-1)ⁿ (π/2)^(2n) z^(4n+1) / ((2n)! (4n+1)) for C
function _fresnel_series(s::Bool, z::T) where {T<:Number}
    big = real(T) === BigFloat
    big && abs(z) > 400 && throw(ArgumentError("fresnel: no BigFloat value for an argument this large ($(Float64(abs(z))))"))
    prec = big ? precision(BigFloat) + ceil(Int, π * Float64(abs2(z)) / 2 / log(2)) + 32 : 0
    f = function ()
        zz = big ? T(z) : z
        h = (big ? BigFloat(π) : π) / 2
        z4 = zz^4
        term = s ? h * zz^3 : zz                # (π/2)^(2n+1) z^(4n+3)/(2n+1)!, n = 0, or z for C
        k = s ? 3 : 1                            # 4n + 3, 4n + 1
        m = s ? 1 : 0                            # the factorial's index, 2n + 1 or 2n
        acc = term / k
        tol = big ? eps(BigFloat) : eps(Float64)
        for n in 1:10_000
            term *= -(h^2) * z4 / ((m + 1) * (m + 2))
            m += 2
            k += 4
            t = term / k
            acc += t
            abs(t) <= tol * abs(acc) && abs(term) <= abs(acc) && break
        end
        acc
    end
    big || return f()
    r = setprecision(f, BigFloat, prec)
    T(r)
end

_fresnel_value(s::Bool) = args -> begin
    a = args[1]
    iszero(a) && return ZERO
    a.kind == NUMBER && isfloat(a.num) && isinf(a.num.f) && return a.num.f > 0 ? HALF : -HALF
    _odd(s ? :fresnels : :fresnelc, a)
end

register_function!(:fresnels; julia = z -> _fresnel_num(true, z), eval = _fresnel_value(true),
    diff = (args, i) -> sin(PI * args[1]^2 / 2))
register_function!(:fresnelc; julia = z -> _fresnel_num(false, z), eval = _fresnel_value(false),
    diff = (args, i) -> cos(PI * args[1]^2 / 2))

# S(u), C(u) in erf, for limits (Gruntz takes erf of complex multiples): the formulas above
function _fresnel_as_erf(e::Sym)
    _replace_walk(e) do t
        (t.kind == FN && t.name in (:fresnels, :fresnelc) && length(t.args) == 1) || return nothing
        u = t.args[1]
        a = (1 + IM) * sqrt(PI) / 2
        e1 = apply_function(:erf, Sym[a * u])
        e2 = apply_function(:erf, Sym[(1 - IM) * sqrt(PI) / 2 * u])
        t.name === :fresnels ? (1 + IM) / 4 * (e1 - IM * e2) : (1 - IM) / 4 * (e1 + IM * e2)
    end
end

# ∫ k sin(q), ∫ k cos(q) dx for q = ax² + bx + c (k, a, b, c free of x, a's sign known): with
# t = √(2|a|/π)(x + b/(2a)) and θ = c - b²/(4a), q = θ ± πt²/2, so
#   ∫ sin q = √(π/(2|a|)) (sin θ C(t) ± cos θ S(t)),  ∫ cos q = √(π/(2|a|)) (cos θ C(t) ∓ sin θ S(t)).
# The terms of a sum that are such: (their antiderivative, the rest), or nothing for none
function _fresnel_split(ex::Sym, x::Sym, A)
    terms = ex.kind == ADD ? [(Sym(c), t) for (c, t) in zip(ex.coefs, ex.args)] : [(ONE, ex)]
    F, rest, any_ = ZERO, ex.kind == ADD ? Sym(ex.num) : ZERO, false
    for (c, t) in terms
        r = _fresnel_term(c * t, x, A)
        if r === nothing
            rest += c * t
        else
            F += r
            any_ = true
        end
    end
    any_ ? (F, rest) : nothing
end

function _fresnel_term(t::Sym, x::Sym, A)
    k, g = ONE, t
    if t.kind == MUL
        trig = [i for i in eachindex(t.args) if t.args[i].kind == FN && t.args[i].name in (:sin, :cos) && isone(t.exps[i])]
        length(trig) == 1 || return nothing
        i = only(trig)
        g = t.args[i]
        k = Sym(t.num) * prod((t.args[j]^t.exps[j] for j in eachindex(t.args) if j != i); init = ONE)
        has(k, x) && return nothing
    end
    (g.kind == FN && g.name in (:sin, :cos) && length(g.args) == 1) || return nothing
    q = expand(g.args[1])
    P = try Poly{Sym}(q, [x]) catch err; err isa ArgumentError || rethrow(); return nothing end
    degree(P, x) == 2 || return nothing
    a, b, c = coeff(q, x, 2), coeff(q, x, 1), coeff(q, x, 0)
    any(e -> has(e, x), (a, b, c)) && return nothing
    sg = if !_hassymbol(a)
        v = try _evalc(a) catch; return nothing end
        abs(imag(v)) > 0 ? 0 : real(v) > 0 ? 1 : real(v) < 0 ? -1 : 0
    elseif A !== nothing && _holds(Relation(:gt, a, ZERO), A) === true
        1
    elseif A !== nothing && _holds(Relation(:lt, a, ZERO), A) === true
        -1
    else
        0
    end
    sg == 0 && return nothing
    aa = sg * a
    tt = sqrt(2aa / PI) * (x + b / (2a))
    θ = c - b^2 / (4a)
    S, C = fresnels(tt), fresnelc(tt)
    m = sqrt(PI / (2aa))
    r = g.name === :sin ? sin(θ) * C + sg * cos(θ) * S : cos(θ) * C - sg * sin(θ) * S
    k * m * r
end

# S(u), C(u) for u → ±∞, with φ = πu²/2 and n terms after ±1/2 (by parts on ∫ᵤ^∞ e^(iπt²/2)):
#   S = sgn/2 - cos φ/(πu) - sin φ/(π²u³) + 3cos φ/(π³u⁵) + 15 sin φ/(π⁴u⁷) - …
#   C = sgn/2 + sin φ/(πu) - cos φ/(π²u³) - 3sin φ/(π³u⁵) + 15 cos φ/(π⁴u⁷) + …
function _fresnel_asymptotic(s::Bool, u::Sym, sg::Int, n::Int)
    φ = PI * u^2 / 2
    c, sn = cos(φ), sin(φ)
    termsS = (-c / (PI * u), -sn / (PI^2 * u^3), 3c / (PI^3 * u^5), 15sn / (PI^4 * u^7))
    termsC = (sn / (PI * u), -c / (PI^2 * u^3), -3sn / (PI^3 * u^5), 15c / (PI^4 * u^7))
    Sym(sg) / 2 + _add_many((s ? termsS : termsC)[1:n])
end

# the limit of f with its Fresnel integrals of x replaced by expansions: asymptotic where
# the argument tends to ±∞, Taylor polynomials about a finite limit L (S and C are entire),
# with fewer and with more terms, the answer when the two agree; nothing when an
# argument's limit is not found
function _fresnel_limit(f::Sym, x::Sym, x0, dir, assume)
    nodes = Sym[]
    _has_node(f, t -> (t.kind == FN && t.name in (:fresnels, :fresnelc) && has(t, x) &&
                       !any(n -> _eq(n, t), nodes) && push!(nodes, t); false))
    signs, Ls = Int[], Sym[]
    for t in nodes
        Lu = try limit(t.args[1], x => x0; dir, assume) catch err; err isa InterruptException && rethrow(); return nothing end
        Lu isa Sym || return nothing
        push!(signs, _isinfpt(Lu) ? _infsign(Lu) : 0)
        push!(Ls, Lu)
    end
    function with_terms(n)
        g = f
        for (i, (t, sg)) in enumerate(zip(nodes, signs))
            r = sg == 0 ? _fresnel_taylor(t.name === :fresnels, t.args[1], Ls[i], 4n + 2) :
                _fresnel_asymptotic(t.name === :fresnels, t.args[1], sg, n)
            g = _replace_walk(g) do w
                _eq(w, t) ? r : nothing
            end
        end
        limit(g, x => x0; dir, assume)
    end
    L2 = with_terms(2)
    L3 = with_terms(3)
    _sameanswer(L2, L3) || throw(ArgumentError("limit: the Fresnel integrals' expansions do not settle the limit"))
    L2
end

# S(L) + Σₖ₌₁ⁿ S⁽ᵏ⁾(L) (u - L)ᵏ/k! (or C's), S' = sin(πw²/2), C' = cos(πw²/2)
function _fresnel_taylor(s::Bool, u::Sym, L::Sym, n::Int)
    w = Sym(:_fresnel_w)
    d = s ? sin(PI * w^2 / 2) : cos(PI * w^2 / 2)
    v = s ? fresnels(L) : fresnelc(L)
    f = ONE
    for k in 1:n
        f *= k
        v += subs(d, w => L) * (u - L)^k / f
        d = diff(d, w)
    end
    v
end
