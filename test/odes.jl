# Differential equations and Laplace transforms (included by runtests.jl)

@testset "differential equations" begin
    @syms y()
    D = SymJ.dsolve
    C1, C2, C3 = Sym(:C1), Sym(:C2), Sym(:C3)
    # initial values at an equilibrium: the constant solution (SymPy #29827)
    @test D(y'(x) ~ log(y(x)), y(x), y(1) => 1) == (y(x) ~ 1)
    @test D(y''(x) ~ y(x) * (y(x) - 1), y(x), y(0) => 1, y'(0) => 0) == (y(x) ~ 1)
    @test D(y'(x) ~ y(x), y(x), y(0) => 2) == (y(x) ~ 2exp(x))       # not an equilibrium
    # derivatives: bound variables, y' notation, evaluation after substitution
    @test y'(x) == diff(y(x), x) && y''(x) == diff(y(x), x, 2) && free_symbols(y'(0)) == Sym[]
    @test subs(diff(y(x), x), x => 0) == y'(0) && subs(y'(0), x => 5) == y'(0)
    @test subs(diff(y(x), x), y(x) => sin(x)) == cos(x) && subs(y''(0), y(x) => x^3 + x^2) == 2
    @test diff(y'(0), x) == 0 && diff(SymJ._integrate(exp(x^2), x), x) == exp(x^2)
    @test string(y'(0)) == "y'(0)" && string(y''(x + 1)) == "y''(x + 1)" && string(y'(0) + 2y(1)) == "y'(0) + 2*y(1)"
    # first order
    @test D(y'(x) ~ y(x), y(x)) == (y(x) ~ C1 * exp(x)) && D(y'(x) ~ x * y(x), y(x)) == (y(x) ~ C1 * exp(x^2 / 2))
    @test D(y'(x) + 2y(x) ~ exp(x), y(x)) == (y(x) ~ C1 * exp(-2x) + exp(x) / 3)
    @test D(y'(x) ~ x * y(x)^2, y(x)) == (y(x) ~ -2 / (x^2 + C1)) && D(y'(x) ~ y(x)^2, y(x)) == (y(x) ~ -1 / (x + C1))
    @test D(y'(x) ~ (x + y(x)) / x, y(x)) == (y(x) ~ (C1 + log(x)) * x)                    # homogeneous
    @test D(2x * y(x) + (x^2 + 3y(x)^2) * y'(x) ~ 0, y(x)) == (x^2 * y(x) + y(x)^3 ~ C1)   # exact, implicit
    @test D(y'(x) ~ sin(x) * y(x), y(x)) == (y(x) ~ C1 * exp(-cos(x)))
    # constant coefficients, any order
    @test D(y''(x) + y(x) ~ 0, y(x)) == (y(x) ~ C1 * cos(x) + C2 * sin(x))
    @test D(y''(x) - 3y'(x) + 2y(x) ~ 0, y(x)) == (y(x) ~ C1 * exp(x) + C2 * exp(2x))
    @test D(y''(x) - 2y'(x) + y(x) ~ 0, y(x)) == (y(x) ~ (C1 + C2 * x) * exp(x))         # a double root
    @test D(y''(x) + 2y'(x) + 5y(x) ~ 0, y(x)) == (y(x) ~ (C1 * cos(2x) + C2 * sin(2x)) * exp(-x))
    @test D(y''(x) + y(x) ~ x, y(x)) == (y(x) ~ C1 * cos(x) + C2 * sin(x) + x)
    @test D(y''(x) - y(x) ~ exp(x), y(x)) == (y(x) ~ C1 * exp(-x) + C2 * exp(x) + x * exp(x) / 2)   # resonance
    @test D(y''(x) + y(x) ~ sin(x), y(x)) == (y(x) ~ C1 * cos(x) + C2 * sin(x) - x * cos(x) / 2)
    @test D(y''(x) + 4y(x) ~ x * cos(x), y(x)) == (y(x) ~ C1 * cos(2x) + C2 * sin(2x) + x * cos(x) / 3 + 2sin(x) / 9)
    @test D(y''(x) ~ cos(2x) + 1, y(x)) == (y(x) ~ C1 * x + C2 + x^2 / 2 - cos(2x) / 4)
    r = D(diff(y(x), x, 3) - y(x) ~ 0, y(x))
    @test iszero(simplify(diff(r.rhs, x, 3) - r.rhs)) && all(C -> has(r.rhs, C), (C1, C2, C3))
    # Cauchy–Euler; y missing
    @test D(x^2 * y''(x) - 2y(x) ~ 0, y(x)) == (y(x) ~ C1 / x + C2 * x^2)
    @test D(x^2 * y''(x) + x * y'(x) + y(x) ~ 0, y(x)) == (y(x) ~ C1 * cos(log(x)) + C2 * sin(log(x)))
    r = D(x^2 * y''(x) - x * y'(x) + y(x) ~ x, y(x))                                       # a double root
    @test iszero(simplify(x^2 * diff(r.rhs, x, 2) - x * diff(r.rhs, x) + r.rhs - x)) && has(r.rhs, C1) && has(r.rhs, C2)
    @test D(x * y''(x) + y'(x) ~ 0, y(x)) == (y(x) ~ C1 * log(x) + C2)
    # initial conditions
    @test D(y''(x) + y(x) ~ 0, y(x), y(0) => 0, y'(0) => 1) == (y(x) ~ sin(x))
    @test D(y''(x) + 4y(x) ~ 0, y(x), y(0) => 1, y'(0) => 0) == (y(x) ~ cos(2x))
    @test D(y'(x) ~ y(x), y(x), y(0) => 3) == (y(x) ~ 3exp(x)) && D(y'(x) ~ x * y(x)^2, y(x), y(0) => 1) == (y(x) ~ -2 / (x^2 - 2))
    @test D(y'(x) ~ -x / y(x), y(x), y(0) => 2) == (y(x) ~ sqrt(4 - x^2))                 # the branch through (0, 2)
    @test D(y''(x) - y(x) ~ 1, y(x), y(0) => 0, y'(0) => 0) == (y(x) ~ cosh(x) - 1)
    # one constant, in exp: the equation's solve (in tan: test/extensions.jl)
    @test cancel(D(y'(x) ~ y(x) * (1 - y(x) / 100), y(x), y(0) => 20).rhs) == 100exp(x) / (exp(x) + 4)
    @syms a
    @test D(y'(x) ~ a * y(x) * (1 - y(x)), y(x)) == (y(x) ~ 1 / (C1 / exp(a * x) + 1))
    @test cancel(D(y'(x) ~ a * y(x) * (1 - y(x)), y(x), y(0) => Sym(1) / 2).rhs) == exp(a * x) / (exp(a * x) + 1)
    # an implicit solution quadratic in y: its two branches
    @test D(y'(x) ~ (1 - x) / y(x), y(x)) == [y(x) ~ -sqrt(-x^2 + C1 + 2x), y(x) ~ sqrt(-x^2 + C1 + 2x)]
    @test D(y'(x) ~ x / y(x)^2, y(x)) == (-x^2 / 2 + y(x)^3 / 3 ~ C1)
    # assumptions: complex characteristic roots give cosines and sines
    @syms k m
    @test D(y''(x) ~ -(k / m) * y(x), y(x)) == (y(x) ~ C1 / exp(sqrt(-k * m) * x / m) + C2 * exp(sqrt(-k * m) * x / m))
    @test D(y''(x) ~ -(k / m) * y(x), y(x); assume = (k ≫ 0, m ≫ 0)) == (y(x) ~ C1 * cos(sqrt(k * m) * x / m) + C2 * sin(sqrt(k * m) * x / m))
    @test assuming(() -> D(y''(x) ~ -(k / m) * y(x), y(x)), k ≫ 0, m ≫ 0) == (y(x) ~ C1 * cos(sqrt(k * m) * x / m) + C2 * sin(sqrt(k * m) * x / m))
    @test D(y''(x) + 2a * y'(x) + (a^2 + k^2) * y(x) ~ 0, y(x); assume = (a ∈ ℝ, k ≫ 0)) == (y(x) ~ (C1 * cos(k * x) + C2 * sin(k * x)) / exp(a * x))
    @test D(y''(x) ~ k * y(x), y(x); assume = (k ≫ 0,)) == (y(x) ~ C1 / exp(sqrt(k) * x) + C2 * exp(sqrt(k) * x))
    # autonomous, y'' = f(y, y'): p = y' as a function of y
    @test D(y(x) * y''(x) ~ y'(x)^2, y(x)).rhs == exp(C1 * (C2 + x))
    r = D(y(x) * y''(x) + y'(x)^2 ~ 0, y(x))                      # implicit, of second order
    @test r.lhs == y(x)^2 / (2C1) - x && r.rhs == C2
    @test D(y(x) * y''(x) + y'(x)^2 ~ 0, y(x), y(0) => 1, y'(0) => 1).rhs == sqrt(2x + 1)
    @test D(y''(x) ~ y'(x)^2 / y(x), y(x), y(0) => 1, y'(0) => 2).rhs == exp(2x)   # constants through the implicit form
    @test D(y''(x) + y'(x)^2 ~ 0, y(x)).rhs == C2 + log(C1 + x)
    # the built-in integrator: powers of a linear function
    @test SymJ._integrate(1 / (x + 2), x) == log(x + 2) && SymJ._integrate(3 / (2x + 1)^2, x) == -(3 // 2) / (2x + 1)
    # Riccati: linearized (u'' … = 0 by constant coefficients or Cauchy–Euler), or from a
    # particular solution
    @test D(y'(x) ~ y(x)^2 - 2 / x^2, y(x)).rhs == (-2x^3 + C1) / (x^4 + C1 * x)
    @test D(y'(x) ~ y(x)^2 + 4, y(x)).rhs == (2C1 * sin(2x) - 2cos(2x)) / (C1 * cos(2x) + sin(2x))   # not an Integral
    @test D(y'(x) ~ 1 + x^2 - 2x * y(x) + y(x)^2, y(x)).rhs == x + 1 / (C1 - x)                  # y₁ = x
    @test D(y'(x) ~ exp(2x) + (1 - 2exp(x)) * y(x) + y(x)^2, y(x)).rhs == exp(x) + 1 / (C1 / exp(x) - 1)   # y₁ = eˣ
    r = D(y'(x) ~ -y(x)^2 + 2 / x^2, y(x), y(1) => 1)
    @test iszero(simplify(diff(r.rhs, x) + r.rhs^2 - 2 / x^2)) && subs(r.rhs, x => 1) == 1
    # series solutions: Taylor at an ordinary point, Frobenius at a regular singular point
    S = (eq, ics...; kw...) -> D(eq, y(x), ics...; kw...).rhs
    @test S(y''(x) + y(x) ~ 0; series = 6) == C1 * (1 - x^2 / 2 + x^4 / 24) + C2 * (x - x^3 / 6 + x^5 / 120)
    @test S(y''(x) - x * y(x) ~ 0; series = 8) == C1 * (1 + x^3 / 6 + x^6 / 180) + C2 * (x + x^4 / 12 + x^7 / 504)   # Airy
    @test S((1 - x^2) * y''(x) - 2x * y'(x) + 6y(x) ~ 0; series = 6) == C1 * (1 - 3x^2) + C2 * (x - 2x^3 / 3 - x^5 / 5)
    @test S(y'(x) ~ x + y(x)^2, y(0) => 1; series = 5) == 1 + x + 3x^2 / 2 + 4x^3 / 3 + 17x^4 / 12   # nonlinear
    @test S(y''(x) + y(x) ~ x, y(1) => 0, y'(1) => 1; series = 4) == (x - 1) + (x - 1)^2 / 2
    J0 = 1 - x^2 / 4 + x^4 / 64
    @test S(x^2 * y''(x) + x * y'(x) + x^2 * y(x) ~ 0; series = 6) == C1 * J0 + C2 * (log(x) * J0 + x^2 / 4 - 3x^4 / 128)   # equal roots
    @test S(2x^2 * y''(x) + x * y'(x) - (x + 1) * y(x) ~ 0; series = 4) ==
          C1 * (x + x^2 / 5 + x^3 / 70 + x^4 / 1890) + C2 * (1 / sqrt(x) - sqrt(x) - x^(3 // 2) / 2 - x^(5 // 2) / 18)   # r = 1, -1/2
    @test has(S(x * y''(x) - y(x) ~ 0; series = 4), log(x))                                # roots 1, 0: a log term
    @test_throws ArgumentError S(x^3 * y''(x) + y(x) ~ 0; series = 4)                      # irregular singular point
    @test_throws ArgumentError S(x * y''(x) - y(x) ~ 0, y(0) => 1, y'(0) => 0; series = 4)
    @test_throws ArgumentError S(y'(x) ~ y(x); series = 0)
    # errors
    # two branches of y' = ±sqrt(…): one solution each, or the branch through the
    # initial conditions (C₁ fixed first, so the integral can be elementary)
    r = D(y''(x) ~ -sin(y(x)), y(x))
    @test r isa Vector && length(r) == 2 && all(q -> q.rhs == C2, r)
    @test D(y''(x) ~ 2y(x)^3, y(x), y(0) => 1, y'(0) => 1).rhs == -1 / (x - 1)
    @test D(y''(x) ~ 2y(x)^3, y(x), y(0) => 1, y'(0) => -1).rhs == 1 / (x + 1)
    @test D(y''(x) ~ exp(2y(x)), y(x), y(0) => 0, y'(0) => 1).rhs == log(-1 / (x - 1))
    r = D(y''(x) ~ -sin(y(x)), y(x), y(0) => 0, y'(0) => 2)                  # ∫₀^y … = x
    @test r.rhs == 0 && string(r.lhs) == "-x + Integral(sqrt(2)/(2*sqrt(cos(y) + 1)), (y, 0, y(x)))"
    # a definite integral prints as integrate takes it; its variable is bound
    @syms a b
    I4 = SymJ._mkfn_node(:Integral, Sym[cos(Sym(:u)) * a, Sym(:u), Sym(0), b])
    @test string(I4) == "Integral(a*cos(u), (u, 0, b))" && free_symbols(I4) == [a, b]
    @test subs(I4, Sym(:u) => 2, a => 3, b => x) == SymJ._mkfn_node(:Integral, Sym[3cos(Sym(:u)), Sym(:u), Sym(0), x])
    @test diff(I4, b) == a * cos(b)
    @test_throws ArgumentError D(y'(x) ~ y(x)^2 + x, y(x))                                  # Riccati to Airy: no method
    @test_throws ArgumentError D(y'(x) ~ y(x), y(x), y(0) => 1, y(1) => 2)
    @test_throws ArgumentError D(y(x) + x ~ 0, y(x))
    @test_throws ArgumentError D(y'(x) - y(x), y(x))
end

@testset "systems of differential equations" begin
    @syms t a u() v() w()
    D = SymJ.dsolve
    C1, C2, C3 = Sym(:C1), Sym(:C2), Sym(:C3)
    rhs_(r) = [s.rhs for s in r]
    # each solution satisfies its equations (the residuals simplify to 0)
    function ok(eqs, us, r)
        sb = [Pair(u, s) for (u, s) in zip(us, rhs_(r))]
        all(eq -> iszero(simplify(subs(eq.lhs - eq.rhs, sb...))), eqs)
    end
    @test rhs_(D([u'(t) ~ v(t), v'(t) ~ -u(t)], [u(t), v(t)])) == [C1 * sin(t) + C2 * cos(t), C1 * cos(t) - C2 * sin(t)]
    @test rhs_(D([u'(t) ~ 2u(t) + v(t), v'(t) ~ u(t) + 2v(t)], [u(t), v(t)])) == [-C1 * exp(t) + C2 * exp(3t), C1 * exp(t) + C2 * exp(3t)]
    @test rhs_(D([u'(t) ~ u(t) + v(t), v'(t) ~ v(t)], [u(t), v(t)])) == [(C2 * t + C1) * exp(t), C2 * exp(t)]   # defective
    @test rhs_(D([u'(t) ~ v(t), v'(t) ~ -u(t)], [u(t), v(t)], u(0) => 1, v(0) => 0)) == [cos(t), -sin(t)]
    @test rhs_(D([u'(t) ~ a * v(t), v'(t) ~ -a * u(t)], [u(t), v(t)])) == [C1 * sin(a * t) + C2 * cos(a * t), C1 * cos(a * t) - C2 * sin(a * t)]
    # complex eigenvalues, a 3×3 system, higher order, implicit form, forcing
    for (eqs, us) in (([u'(t) ~ u(t) - v(t), v'(t) ~ 5u(t) - 3v(t)], [u(t), v(t)]),
                      ([u'(t) ~ v(t), v'(t) ~ w(t), w'(t) ~ u(t)], [u(t), v(t), w(t)]),
                      ([u''(t) ~ -2u(t) + v(t), v''(t) ~ u(t) - 2v(t)], [u(t), v(t)]),
                      ([u''(t) ~ v(t), v'(t) ~ u(t)], [u(t), v(t)]),
                      ([u'(t) + v'(t) ~ u(t), u'(t) - v'(t) ~ v(t)], [u(t), v(t)]),
                      ([u'(t) ~ v(t) + t, v'(t) ~ -u(t)], [u(t), v(t)]),
                      ([u'(t) ~ v(t) + t * exp(t), v'(t) ~ u(t)], [u(t), v(t)]),
                      ([u'(t) ~ -v(t) + cos(t), v'(t) ~ u(t)], [u(t), v(t)]))
        r = D(eqs, us)
        @test ok(eqs, us, r) && all(s -> !has(s, IM), rhs_(r))
    end
    @test rhs_(D([u'(t) ~ u(t) + exp(2t), v'(t) ~ u(t) - v(t)], [u(t), v(t)]))[2] == C1 / exp(t) + C2 * exp(t) + exp(2t) / 3
    r = D([u''(t) + v(t) ~ 0, v'(t) ~ u(t)], [u(t), v(t)], u(0) => 1, u'(0) => 0, v(0) => 0)
    @test ok([u''(t) + v(t) ~ 0, v'(t) ~ u(t)], [u(t), v(t)], r) && subs(r[1].rhs, t => 0) == 1 && subs(r[2].rhs, t => 0) == 0
    # errors
    @test_throws ArgumentError D([u'(t) ~ v(t)^2, v'(t) ~ u(t)], [u(t), v(t)])           # nonlinear
    @test_throws ArgumentError D([u'(t) ~ t * v(t), v'(t) ~ u(t)], [u(t), v(t)])         # variable coefficients
    @test_throws ArgumentError D([u'(t) + v'(t) ~ u(t), u'(t) + v'(t) ~ v(t)], [u(t), v(t)])   # not in normal form
    @test_throws ArgumentError D([u'(t) ~ v(t)], [u(t), v(t)])
    @test_throws ArgumentError D([u'(t) ~ v(t), v'(t) ~ -u(t)], [u(t), v(t)], u(0) => 1)
    @test_throws ArgumentError D([u'(t) ~ v(t), v'(t) ~ -u(t)], [u(t), v(t)], u(0) => 1, w(0) => 1)
end

@testset "Laplace transforms" begin
    @syms t s w y() u() v()
    H, δ = SymJ.heaviside, SymJ.dirac
    L(f) = SymJ.laplace_transform(f, t, s)
    IL(F) = SymJ.inverse_laplace_transform(F, s, t)
    C1, C2 = Sym(:C1), Sym(:C2)
    @test L(Sym(1)) == 1 / s && L(t^3) == 6 / s^4 && L(exp(2t)) == 1 / (s - 2) && L(sin(3t)) == 3 / (s^2 + 9)
    @test L(cos(w * t)) == s / (s^2 + w^2) && L(sinh(t)) == 1 / (s^2 - 1) && L(t * sin(t)) == 2s / (s^2 + 1)^2
    # bᵗ = e^(t log b) (SymPy #28136); derivatives of an unknown u(t) (#28720)
    @test iszero(cancel(L(Sym(2)^t) - 1 / (s - log(Sym(2)))))     # (signs arranged differently on Julia 1.10)
    @test iszero(cancel(L(w^t) - 1 / (s - log(w))))
    let U = SymJ.apply_function(:LaplaceTransform, Sym[u(t), t, s])
        @test L(diff(u(t), t, 2)) == s^2 * U - s * u(0) - subs(diff(u(t), t), t => 0)
        @test IL(L(diff(u(t), t, 2))) == diff(u(t), t, 2)
        @test IL(L(diff(u(t), t) + 3u(t) + exp(t))) == diff(u(t), t) + 3u(t) + exp(t)
    end
    @test L(exp(-t) * cos(2t)) == (s + 1) / (s^2 + 2s + 5) && L(sqrt(t)) == SymJ.apply_function(:gamma, Sym[Sym(3 // 2)]) / s^(3 // 2)
    @test L(H(t - 1)) == exp(-s) / s && L(δ(t - 3)) == exp(-3s) && L(δ(t)) == 1 && L(H(t - 1) * exp(t)) == exp(1 - s) / (s - 1)
    @test IL(1 / (s^2 + 2s + 5)) == sin(2t) / (2exp(t)) && IL(1 / (s^2 * (s^2 + 1))) == t - sin(t)
    @test IL((s + 1) / (s^2 + 2s + 5)^2) == t * sin(2t) / (4exp(t)) && IL(1 / (s^2 - 2)) == sqrt(Sym(2)) * sinh(sqrt(Sym(2)) * t) / 2
    @test IL(exp(-2s) / s^2) == H(t - 2) * (t - 2) && IL(s / (s + 1)) == δ(t) - exp(-t)
    # parameters: the roots of a quadratic as α ± iβ, the parameters taken to be positive
    @syms a
    @test IL(1 / (s^2 + w^2)) == sin(w * t) / w && IL(s / (s^2 + w^2)) == cos(w * t) && IL(1 / (s^2 - w^2)) == sinh(w * t) / w
    @test IL(1 / ((s + a)^2 + w^2)) == sin(w * t) / (w * exp(a * t)) && IL(a / (s * (s + a))) == 1 - exp(-a * t)
    @test IL(1 / (s^2 + w^2)^2) == (sin(w * t) - w * t * cos(w * t)) / (2w^3)
    @test SymJ.dsolve(y''(t) + w^2 * y(t) ~ H(t - 1), y(t), y(0) => 0, y'(0) => 0).rhs == (1 - cos(w * (t - 1))) * H(t - 1) / w^2
    for f in (t^2 * exp(-t), cos(3t) - 2sin(t), exp(2t) * sin(t), t * cosh(t))
        @test iszero(simplify(IL(L(f)) - f))
    end
    @test_throws ArgumentError L(log(t))
    @test IL(s^2) == SymJ.dirac(t, 2)                                  # δ''(t) (SymPy #26775)
    # dsolve: chosen for a right side with heaviside or dirac, or with method = :laplace
    D = SymJ.dsolve
    @test D(y''(t) + y(t) ~ H(t - 1), y(t), y(0) => 0, y'(0) => 0).rhs == (1 - cos(t - 1)) * H(t - 1)
    @test D(y''(t) + 4y(t) ~ δ(t - 2), y(t), y(0) => 0, y'(0) => 0).rhs == H(t - 2) * sin(2t - 4) / 2
    @test D(y'(t) + y(t) ~ H(t - 1) - H(t - 2), y(t), y(0) => 0).rhs == (1 - exp(1 - t)) * H(t - 1) + (exp(2 - t) - 1) * H(t - 2)
    @test D(y''(t) + y(t) ~ H(t - 1), y(t)).rhs == C1 * cos(t) + C2 * sin(t) + (1 - cos(t - 1)) * H(t - 1)
    @test D(y''(t) + y(t) ~ H(t - 1), y(t), y(1) => 0, y'(1) => 0).rhs == 1 - cos(t - 1)          # shifted to the point
    @test D(y''(t) + 3y'(t) + 2y(t) ~ exp(-t), y(t), y(0) => 1, y'(0) => 0; method = :laplace).rhs == (t + 1) / exp(t)
    @test_throws ArgumentError D(y''(t) - t * y(t) ~ 0, y(t); method = :laplace)
    @test_throws ArgumentError D(y'(t) ~ y(t), y(t); method = :fourier)
    r = D([u'(t) ~ v(t), v'(t) ~ -u(t) + H(t - 1)], [u(t), v(t)], u(0) => 0, v(0) => 0)
    @test [q.rhs for q in r] == [(1 - cos(t - 1)) * H(t - 1), H(t - 1) * sin(t - 1)]
    @test [q.rhs for q in D([u'(t) ~ v(t), v'(t) ~ -u(t)], [u(t), v(t)], u(0) => 1, v(0) => 0; method = :laplace)] == [cos(t), -sin(t)]
end
