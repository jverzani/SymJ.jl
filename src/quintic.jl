# Solvable quintics: an irreducible quintic over ℚ whose Galois group is in F₂₀ (the
# Frobenius group of order 20, which holds D₅ and C₅) has its roots in radicals; the
# others (S₅, A₅) have none. Following Dummit ("Solving solvable quintics", 1991), with
# numerics only to pick an ordering and branches, the answer checked at the end:
#
# 1. f depressed (x = y - a₄/5) and scaled to a monic integer quintic in t = L y.
# 2. Its roots t₀, …, t₄ numerically (512 bits). F₂₀ is the group j ↦ aj + b of the
#    indices mod 5; θ = Σⱼ tⱼ² (tⱼ₊₁ tⱼ₊₄ + tⱼ₊₂ tⱼ₊₃) is invariant under it, and its six
#    values over the orderings are the roots of the resolvent sextic, an integer
#    polynomial. The group is in F₂₀ exactly when one of them is an integer (a simple
#    root); that ordering is the cyclic one.
# 3. The Lagrange resolvents Rₖ = Σⱼ ζ^(jk) tⱼ (ζ = e^(2πi/5), k = 1, …, 4; R₀ = 0),
#    kept as polynomials in a formal ζ reduced by ζ⁴ + ζ³ + ζ² + ζ + 1. Under the order
#    F₂₀ keeps, s₁ = R₁⁵ + R₄⁵ and P₁ = R₁R₄ are fixed by D₅ (their ζ-coefficients in
#    ℚ(√D)) and exchanged with s₂ = R₂⁵ + R₃⁵, P₂ = R₂R₃ by the rest; so s₁ + s₂ and
#    (s₁ - s₂)·V (V = ∏ᵢ<ⱼ (tᵢ - tⱼ) = ±√D) have integer coefficients, found by rounding.
# 4. R₁⁵ is a root of Y² - s₁Y + P₁⁵, R₁ its fifth root, the root and the branch the
#    ones the numeric R₁ has (picked by the exact expressions' values); R₄ = P₁/R₁. The
#    same for R₂, R₃.
# 5. tⱼ = (Σₖ ζ^(-jk) Rₖ)/5, and x = t/L - a₄/5; each checked against the numeric
#    roots at 256 bits and in Float64.
#
# nothing when f is not a solvable quintic, or anything fails to check.
function _solvable_quintic(c::Vector{Sym})
    length(c) == 6 || return nothing
    all(ci -> ci.kind == NUMBER && isexact(ci.num), c) || return nothing
    a = Rational{BigInt}[_rat(ci.num) for ci in c]
    a = a ./ a[6]
    h = -a[5] / 5
    # f(y + h), monic, its y⁴ coefficient 0 (synthetic division by y - h, repeatedly)
    e = copy(a)
    for i in 1:5, k in 5:-1:i
        e[k] += h * e[k + 1]
    end
    iszero(e[5]) || return nothing
    L = reduce(lcm, (denominator(e[k]) for k in 1:4); init = big(1))
    # the scaled quintic in t = L y has integer coefficients if L^(5-k) clears e[k]'s
    # denominator, which it does
    ti = BigInt[numerator(e[k] * L^(6 - k)) for k in 1:5]      # (constant first; t⁵ monic)
    all(k -> isone(denominator(e[k] * L^(6 - k))), 1:5) || return nothing
    setprecision(BigFloat, 512) do
        C = Complex{BigFloat}
        m = C[C.(ti); one(C)]
        z = _aberth(m, big(2.0)^-200)
        dm = [k * m[k + 1] for k in 1:5]
        for k in 1:5, _ in 1:20
            d = evalpoly(z[k], dm)
            iszero(d) && break
            δ = evalpoly(z[k], m) / d
            z[k] -= δ
            abs(δ) <= big(2.0)^-480 * max(1, abs(z[k])) && break
        end
        θ(t) = sum(t[j]^2 * (t[mod1(j + 1, 5)] * t[mod1(j + 4, 5)] + t[mod1(j + 2, 5)] * t[mod1(j + 3, 5)]) for j in 1:5)
        scale = 1 + maximum(abs, z)^4
        tol = big(10.0)^-100 * scale
        hits = Vector{Int}[]
        others = C[]
        for p in _perms4()
            o = [1; 1 .+ p]
            v = θ(z[o])
            if abs(imag(v)) <= tol && abs(real(v) - round(real(v))) <= tol
                push!(hits, o)
            else
                push!(others, v)
            end
        end
        # one class of four orderings (j ↦ aj), and its θ a simple root
        length(hits) == 4 || return nothing
        θ0 = θ(z[hits[1]])
        any(v -> abs(v - θ0) <= big(10.0)^-30 * scale, others) && return nothing
        t = z[hits[1]]
        ζn = cis(2 * big(π) / 5)
        # formal ζ: coefficient vectors mod ζ⁵ - 1, then mod Φ₅
        R(k) = (v = zeros(C, 5); for j in 0:4; v[mod(j * k, 5) + 1] += t[j + 1]; end; v)
        mul(u, w) = (r = zeros(C, 5); for i in 0:4, j in 0:4; r[mod(i + j, 5) + 1] += u[i + 1] * w[j + 1]; end; r)
        pw5(u) = mul(mul(mul(mul(u, u), u), u), u)
        red(u) = [u[i] - u[5] for i in 1:4]                     # ζ⁴ = -(1 + ζ + ζ² + ζ³)
        Rv = [R(k) for k in 1:4]
        s1 = red(pw5(Rv[1]) + pw5(Rv[4])); s2 = red(pw5(Rv[2]) + pw5(Rv[3]))
        P1 = red(mul(Rv[1], Rv[4]));       P2 = red(mul(Rv[2], Rv[3]))
        V = prod(t[i] - t[j] for i in 1:5 for j in (i + 1):5)
        D = round(BigInt, real(V^2))
        abs(V^2 - D) <= big(10.0)^-60 * max(1, abs(V)^2) || return nothing
        # exact V: ±√D, the sign the numeric V has
        Vx = sqrt(Sym(D))
        Vxn = C(_bigvalue(Vx, 512))
        Vx = abs(Vxn - V) < abs(Vxn + V) ? Vx : -Vx
        function exact(u1, u2)                    # u1's ζ-coefficients: (U + W V/D)/2
            out = Sym[]
            for i in 1:4
                U, W = u1[i] + u2[i], (u1[i] - u2[i]) * V
                (abs(imag(U)) <= tol && abs(imag(W)) <= tol * max(1, abs(V))) || return nothing
                Ui, Wi = round(BigInt, real(U)), round(BigInt, real(W))
                (abs(real(U) - Ui) <= big(10.0)^-60 * max(1, abs(U)) &&
                 abs(real(W) - Wi) <= big(10.0)^-60 * max(1, abs(W))) || return nothing
                push!(out, (Sym(Ui) + Sym(Wi) * Vx / Sym(D)) / 2)
            end
            out
        end
        S1, S2 = exact(s1, s2), exact(s2, s1)
        Q1, Q2 = exact(P1, P2), exact(P2, P1)
        any(isnothing, (S1, S2, Q1, Q2)) && return nothing
        # ζᵏ = cos(2πk/5) + i sin(2πk/5) in radicals; an element Σ cᵢ ζⁱ as a + ib, each
        # part expanded (rationals in √5, √D, √(10 ± 2√5))
        r5 = sqrt(Sym(5))
        cs5 = Sym[1, (r5 - 1) / 4, -(r5 + 1) / 4, -(r5 + 1) / 4, (r5 - 1) / 4]
        sn5 = Sym[0, sqrt(10 + 2r5) / 4, sqrt(10 - 2r5) / 4, -sqrt(10 - 2r5) / 4, -sqrt(10 + 2r5) / 4]
        ζp(k) = (k = mod(k, 5); k == 0 ? ONE : cs5[k + 1] + IM * sn5[k + 1])
        ev(cs) = expand(_add_many(Sym[cs[i] * cs5[i] for i in 1:4])) + IM * expand(_add_many(Sym[cs[i] * sn5[i] for i in 1:4]))
        smaller(e) = (f = expand(e); _nodes(f) < _nodes(e) ? f : e)
        Rn = [evalpoly(ζn, Rv[k]) for k in 1:4]
        # R₁ (or R₂) from s and P: the root of Y² - sY + P⁵ and the fifth root the numeric
        # R has; the partner P/R
        # the sign and the branch chosen by the exact expression's own value (Sym's
        # principal roots of its forms need not be those of the numbers they stand for)
        # (a radicand on the negative real axis, as for the real roots, sits on the
        # principal root's cut, where rounding picks the branch: N and 256 bits then
        # disagreed; i√(-w), -(-w)^(1/5) are roots too, away from it)
        oncut(v) = real(v) < 0 && abs(imag(v)) <= big(10.0)^-30 * abs(v)
        # a fifth root of y with the numeric value rn
        function root5(y, rn)
            vy = _bigvalue(y, 256)
            vy === nothing && return nothing
            r0 = oncut(C(vy)) ? -(-y)^(1 // 5) : y^(1 // 5)
            v0 = _bigvalue(r0, 256)
            v0 === nothing && return nothing
            k = argmin([abs(C(v0) * ζn^k - rn) for k in 0:4]) - 1
            abs(C(v0) * ζn^k - rn) <= big(10.0)^-40 * max(1, abs(rn)) || return nothing
            (k, r0)                                            # R = ζᵏ r0
        end
        function pair(S, Q, qv, rn, rpn)
            Se, Qe = ev(S), ev(Q)
            # (R = 0, then P = 0 and the partner's fifth power s: the roots ζʲu + ζ⁻ʲ/u of
            # x⁵ - 5x³ + 5x - 3 have R₂ = R₃ = 0)
            small = big(10.0)^-60 * max(1, maximum(abs, Rn))
            if abs(rn) <= small || abs(rpn) <= small
                abs(evalpoly(ζn, qv)) <= small * max(1, maximum(abs, Rn)) || return nothing
                z0 = (0, ZERO)
                abs(rn) <= small && abs(rpn) <= small && return (z0, z0)
                w = root5(Se, abs(rn) <= small ? rpn : rn)
                w === nothing && return nothing
                return abs(rn) <= small ? (z0, w) : (w, z0)
            end
            Δ = smaller(Se^2 - 4Qe^5)
            vΔ = _bigvalue(Δ, 256)
            vΔ === nothing && return nothing
            sq = oncut(C(vΔ)) ? IM * sqrt(-Δ) : sqrt(Δ)
            y5s = [(Se + sq) / 2, (Se - sq) / 2]
            vals = [_bigvalue(y, 256) for y in y5s]
            any(isnothing, vals) && return nothing
            i = argmin([abs(C(v) - rn^5) for v in vals])
            abs(C(vals[i]) - rn^5) <= big(10.0)^-40 * max(1, abs(rn)^5) || return nothing
            Rx = root5(y5s[i], rn)
            Rx === nothing && return nothing
            Rx, (-Rx[1], Qe / Rx[2])                            # P/R = ζ^(-k) P/r0
        end
        A = pair(S1, Q1, P1, Rn[1], Rn[4])
        B = pair(S2, Q2, P2, Rn[2], Rn[3])
        (A === nothing || B === nothing) && return nothing
        Rx = [A[1], B[1], B[2], A[2]]                    # R₁, R₂, R₃ = P₂/R₂, R₄ = P₁/R₁, each (m, ρ): ζᵐ ρ
        roots = Sym[]
        for j in 0:4
            tj = _add_many(Sym[ζp(Rx[k][1] - j * k) * Rx[k][2] for k in 1:4]) / 5
            push!(roots, tj / Sym(L) + Sym(SymNumber(h)))
        end
        # each root checked against the numeric ones, at 256 bits
        xs = [t[j] / L + h for j in 1:5]
        for (r, xn) in zip(roots, xs)
            v = _bigvalue(r, 256)
            (v !== nothing && abs(C(v) - xn) <= big(10.0)^-40 * max(1, abs(xn))) || return nothing
            # and in Float64, as N evaluates it
            w = try ComplexF64(_evalc(r)) catch; return nothing end
            abs(w - ComplexF64(xn)) <= 1e-8 * max(1, abs(xn)) || return nothing
        end
        roots
    end
end

# the 24 permutations of 1:4
function _perms4()
    out = Vector{Int}[]
    for a in 1:4, b in 1:4, c in 1:4
        length(unique((a, b, c))) == 3 || continue
        d = only(setdiff(1:4, (a, b, c)))
        push!(out, [a, b, c, d])
    end
    out
end
