# Derivatives, series, limits, piecewise functions (included by runtests.jl)

@testset "diff" begin
    @test diff(x^3, x) == 3x^2
    @test isnan(N(diff(Sym(NaN), x)))                                 # NaN, not 0
    # an unknown function of an expression: the chain rule, H′ at a point
    let H = SymFunction(:H), z = Sym(:z)
        @test diff(H(x^2 + y), x) == 2x * H'(x^2 + y)
        @test diff(H(sin(x)), x, 2) == H''(sin(x)) * cos(x)^2 - H'(sin(x)) * sin(x)
        @test diff(H(2x), y) == 0
        @test string(diff(H(x), x)) == "Derivative(H(x), x)"            # of a symbol: as it was
        @test subs(diff(H(x^2), x), H(x) => sin(x)) == 2x * cos(x^2)
        # several arguments: Σ ∂ᵢH ∂argsᵢ/∂x, ∂ᵢH at the point in a fresh dummy (SymPy extra-11)
        t, D = Sym(:t), (e, v, a) -> SymJ.apply_function(:Derivative, Sym[e, v, a])
        @test diff(H(x + y, x * y), x) == D(H(t, x * y), t, x + y) + y * D(H(x + y, t), t, x * y)
        @test diff(H(x, x), x) == D(H(t, x), t, x) + D(H(x, t), t, x)
        @test string(diff(H(x, y), x)) == "Derivative(H(x, y), x)"      # symbols: as it was
        @test diff(H(x^2, t), x) == 2x * D(H(Sym(:s), t), Sym(:s), x^2)   # t is taken
        # shown as the call that makes it, not Derivative(H(t, x), t, x^2), SymPy's
        # order notation (#26958)
        @test string(D(H(t, x), t, x^2)) == "Derivative(H(t, x), t)(t => x^2)"
        @test diff(H(t, x), t)(t => x^2) == D(H(t, x), t, x^2)
        # the derivative of one at a point: in the point and in the other arguments
        H1, H2 = D(H(t, x), t, x^2), D(H(x^2, t), t, x)
        H11, H22 = D(diff(H(t, x), t), t, x^2), D(diff(H(x^2, t), t), t, x)
        H12 = D(diff(H(t, x), x), t, x^2)
        @test diff(H1, x) == 2x * H11 + H12
        @test diff(H2, x) == H22 + D(diff(H(x^2, t), x), t, x)
        @test diff(D(H(t, y), t, x^2), y) == D(diff(H(t, y), y), t, x^2)
        R = sqrt(x^2 + y^2 + z^2)                                         # curl of H(r) r̂ is 0
        M, N, P = H(R) * [x, y, z] / R
        @test diff(P, y) - diff(N, z) == 0
        @test diff(M, z) - diff(P, x) == 0
        @test diff(N, x) - diff(M, y) == 0
    end
    # mixed derivatives: in each variable in turn
    @test diff(x^3 * y^2, x, y) == 6x^2 * y == diff(diff(x^3 * y^2, x), y)
    @test diff(x^3 * y^2, y, x) == diff(x^3 * y^2, x, y)
    @test diff(x^3 * y^2, x, x, y) == diff(diff(x^3 * y^2, x, 2), y) == 6x * y * 2
    @test diff(sin(x * y), x, y) == cos(x * y) - x * y * sin(x * y)
    @test diff(x^2 * y * z, x, x, y, z) == 2
    @test diff(x^2 * y, x, z) == 0 && diff(3, x, y) == 0
    @test diff(exp(x * y), x, y) == exp(x * y) + x * y * exp(x * y)
    @test_throws ArgumentError diff(x * y, x, x + y)
    @test diff(x * y, x) == y
    @test diff(sin(x), x) == cos(x)
    @test diff(sin(x)^2, x) == 2sin(x) * cos(x)
    @test diff(exp(2x), x) == 2exp(2x)
    @test diff(log(x), x) == 1 / x
    @test diff(x^y, y) == log(x) * x^y
    @test diff(x^4, x, 2) == 12x^2
    @test diff(log(Sym(2), x), x) == 1 / (x * log(Sym(2)))
    @test kind(diff(f(x), x)) == FN
    @test diff(Sym(3), x) == 0
end

@testset "gradient, jacobian, hessian, divergence, curl, laplacian (SymPy #30186)" begin
    @syms r θ a
    grad, jac, hess = SymJ.gradient, SymJ.jacobian, SymJ.hessian
    @test grad(x^2 * y + sin(y), [x, y]) == [2x * y, x^2 + cos(y)]
    @test grad(x^2 * y, (x, y)) == [2x * y, x^2]                     # a tuple of variables
    @test grad(Sym(3), [x, y]) == [0, 0]
    @test jac([r * cos(θ), r * sin(θ)], [r, θ]) == Sym[cos(θ) -r*sin(θ); sin(θ) r*cos(θ)]
    @test simplify(SymJ.LinearAlgebra.det(jac([r * cos(θ), r * sin(θ)], [r, θ]))) == r
    @test size(jac([x, y, x * y], [x, y])) == (3, 2)                 # m × n
    @test jac(x^2 * y, [x, y]) == Sym[2x*y x^2]                      # a scalar: the row
    @test jac(Sym[x * y; x + y;;], [x, y]) == Sym[y x; 1 1]          # a column matrix
    @test jac((x * y, x + y), (x, y)) == Sym[y x; 1 1]               # tuples
    @test hess(x^3 + x * y^2, [x, y]) == Sym[6x 2y; 2y 2x]
    @test hess(exp(a * x * y), [x, y]) == jac(grad(exp(a * x * y), [x, y]), [x, y])
    @test hess(f(x, y), [x, y])[1, 2] == hess(f(x, y), [x, y])[2, 1]
    @test_throws ArgumentError grad(x^2, [x^2])                       # not a symbol
    @test_throws ArgumentError hess([x, y], [x, y])                  # not a scalar
    div, curl, lap = SymJ.divergence, SymJ.curl, SymJ.laplacian
    @test div([x^2, x * y, y * z], [x, y, z]) == 3x + y
    @test div((x, y), (x, y)) == 2
    @test_throws ArgumentError div([x, y], [x, y, z])               # lengths differ
    @test curl([-y, x, 0], [x, y, z]) == Sym[0, 0, 2]
    @test curl(grad(x^2 * y * z, [x, y, z]), [x, y, z]) == Sym[0, 0, 0]   # curl ∘ grad = 0
    @test div(curl([x * y, y * z^2, x^3], [x, y, z]), [x, y, z]) == 0      # div ∘ curl = 0
    @test_throws ArgumentError curl([x, y], [x, y])                  # three components
    @test lap(x^2 * y + y^3, [x, y]) == 8y
    @test simplify(lap(log(x^2 + y^2), [x, y])) == 0                # harmonic
    @test lap(x^2 * y * z, [x, y, z]) == div(grad(x^2 * y * z, [x, y, z]), [x, y, z])
    @test lap([x^2, y^3], [x, y]) == [2, 6y]                         # entrywise
end

@testset "series" begin
    @syms sa
    S(e, args...) = Sym(series(e, args...))
    @test S(sin(x), x) == x - x^3 / 6 + x^5 / 120
    @test S(sinc(x), x, 0, 4) == 1 - PI^2 * x^2 / 6                      # entire (SymPy #30545)
    @test S(exp(x), x, 0, 4) == 1 + x + x^2 / 2 + x^3 / 6
    @test S(log(1 + x), x, 0, 5) == x - x^2 / 2 + x^3 / 3 - x^4 / 4
    @test S(atan(x), x) == x - x^3 / 3 + x^5 / 5
    @test S(tan(x), x, 0, 8) == x + x^3 / 3 + 2x^5 / 15 + 17x^7 / 315
    @test S(1 / sin(x), x, 0, 3) == 1 / x + x / 6                      # Laurent
    @test S(cot(x), x, 0, 4) == 1 / x - x / 3 - x^3 / 45
    @test S((sin(x) - x) / x^3, x, 0, 4) == -1 // 6 + x^2 / 120          # cancellation
    @test S(sqrt(x + x^2), x, 0, 3) == sqrt(x) + x^(3 // 2) / 2 - x^(5 // 2) / 8   # Puiseux
    @test S(x^x, x, 0, 3) == 1 + log(x) * x + log(x)^2 * x^2 / 2       # logarithmic terms
    @test S(sqrt(x^2 + x), x, Inf, 3) == x + 1 // 2 - 1 / (8x) + 1 / (16x^2)   # at infinity
    @test S(sqrt(x^2 + 1), x, -Inf, 3) == -x - 1 / (2x)
    @test S(log(x + 1), x, Inf, 3) == log(x) + 1 / x - 1 / (2x^2)
    @test S(exp(x), x, 1, 3) == ℯ * (1 + (x - 1) + (x - 1)^2 / 2) |> expand ||
          iszero(expand(S(exp(x), x, 1, 3) - (Sym(ℯ) + Sym(ℯ) * (x - 1) + Sym(ℯ) * (x - 1)^2 / 2)))
    @test iszero(expand(S(tan(x), x, Sym(π) / 2, 3) - (-1 / (x - Sym(π) / 2) + (x - Sym(π) / 2) / 3)))
    @test S((1 + x)^sa, x, 0, 3) == 1 + sa * x + (sa^2 / 2 - sa / 2) * x^2   # symbolic exponent
    @test S(x^2 + 1, x) == x^2 + 1 && string(series(x^2 + 1, x)) == "1 + x^2"   # exact: no O
    # printing, and the error order
    @test string(series(sin(x), x)) == "x - x^3/6 + x^5/120 + O(x^6)"
    @test string(series(exp(x), x, 1, 2)) == "ℯ + ℯ*(x - 1) + O((x - 1)^2)"
    @test string(series(1 / (1 - x), x, Inf, 3)) == "-1/x - 1/x^2 + O(1/x^3)"
    # every term below the order is right: the error is O(h^n)
    for (f, x0, n) in ((exp(sin(x)), 0, 6), (1 / cos(x), 0, 8), (sqrt(1 + x) * log(1 + x), 0, 5),
                       (cot(x) - 1 / x, 0, 6), (exp(x) / (1 + x^2), 1, 5), (asin(x), 0, 7), (sinh(x) * cosh(x), 0, 7))
        T = Sym(series(f, x, x0, n))
        # in BigFloat (calling evaluates with Julia arithmetic), past Float64 rounding
        errs = [abs((f - T)(x0 + big(h))) for h in (1 // 100, 1 // 200)]
        @test errs[1] / errs[2] > 2.0^(n - 0.6)                    # halving h divides it by ~2^n
    end
    # arithmetic and derivatives of series
    s, c = series(sin(x), x), series(cos(x), x)
    @test Sym(s / c) == x + x^3 / 3 + 2x^5 / 15 && Sym(s * c) == x - 2x^3 / 3 + 2x^5 / 15
    @test Sym(diff(s, x)) == 1 - x^2 / 2 + x^4 / 24
    @test Sym(s^2 + c^2) == 1 && Sym(2 * s + 1) == 1 + 2x - x^3 / 3 + x^5 / 60
    @test s(0.1) ≈ sin(0.1) atol = 1e-6
    # no series
    @test_throws ArgumentError series(exp(1 / x), x)                  # essential singularity
    @test_throws ArgumentError series(SymFunction(:g2)(x, 2), x)       # undefined, of two arguments (f(x)'s series is in its derivatives; SymPy #3978)
    @test_throws ArgumentError series(sin(x), x + 1)
end

@testset "limit" begin
    # an integer variable at ∞: sin(πv) without v's integer polynomial part (SymPy #28313)
    let n = Sym(:n), Z = (n ∈ ℤ,)
        @test N(limit((1 + sin(PI * sqrt(4n^2 + 1)))^n, n => Inf; assume = Z)) ≈ exp(pi / 4)
        @test limit(sin(PI * n), n => Inf; assume = Z) == 0
        @test limit(cos(2PI * n) + 1 / n, n => Inf; assume = Z) == 1
        @test_throws SymJ.GruntzError limit(n * sin(PI * sqrt(n^2 + 1)), n => Inf; assume = Z)   # ±π/2
        @test_throws SymJ.GruntzError limit(sin(PI * n), n => Inf)              # n real: no limit
    end
    @syms a b
    @test limit(sin(x) / x, x => 0) == 1 && limit((1 + 1 / x)^x, x => Inf) == E
    # sinc and cosc through sin and cos (SymPy #30545)
    @test limit(sinc(x), x => 0) == 1
    @test limit(cosc(x), x => 0) == 0
    @test limit(sinc(x) / x^2 - 1 / x^2, x => 0) == -PI^2 / 6
    @test limit(1 / x, x => 0; dir = :+) == Sym(Inf) && limit(-x, x => Inf) == Sym(-Inf)
    @test limit(1 / x, x => 0; dir = :+) isa Sym && limit(x^2, x => Sym(Inf)) == Sym(Inf)
    @test limit(x * exp(1 / x), x => 0; dir = -) == 0
    @test_throws SymJ.GruntzError limit(sin(x), x => Inf)
    @test_throws ArgumentError limit(x, x => 0; dir = :up)
    # cases become Cases, with relations as conditions
    cs = limit(x^a, x => Inf)
    @test cs isa SymJ.Cases && length(cs) == 3 && cs[1].first == Any[a ≪ 0]
    @test cs(a => 2) == Sym(Inf) && cs(a => -1) == 0 && cs(a => 0) == 1
    # assumptions: signs for GruntzLimit, then cases filtered and conditions simplified
    @test limit(x^a, x => Inf; assume = @assume(a > 1)) == Sym(Inf)
    @test limit(x^a, x => Inf; assume = @assume(a == 3)) == Sym(Inf)          # substituted
    cs = limit(x^(a - 2), x => Inf; assume = @assume(a != 2))
    @test length(cs) == 2 && cs(a => 3) == Sym(Inf) && cs(a => 1) == 0
    @test limit((a^x + b^x)^(1 / x), x => Inf; assume = @assume(a > b > 0)) == a
    cs = limit(exp(a * x), x => Inf; assume = @assume(a < 0 || a > 1))   # a gap: filtered
    @test length(cs) == 2 && cs(a => 2) == Sym(Inf) && cs(a => -1) == 0
    @test (@assuming a > 0 begin limit(x^a, x => Inf) end) == Sym(Inf)
    @test_throws ArgumentError limit(x^a, x => Inf; assume = @assume(!(a in Reals)))
    @test limit((1 + a / x)^(b * x), x => Inf) == exp(a * b)
    # generic answers where there is one, as gruntz_limit gives them (no cases on
    # where the expression is defined, or on constants that may vanish)
    f = (sqrt(2 * a^3 * x - x^4) - a * cbrt(a^2 * x)) / (a - a^(1 // 4) * x^(3 // 4))
    @test limit(f, x => 1) == SymJ.GruntzLimit.gruntz_limit(f, x, 1) isa Sym
    @test limit(f, x => a) == 16a / 9                       # the classic L'Hôpital example
end

@testset "limit and subs of Cases" begin
    @syms a b
    # limit: the cases that hold near the point, each with its limit
    cs = SymJ.Cases{Sym}([Any[x ≪ 0] => -x, Any[x ~ 0] => Sym(0), Any[x ≫ 0] => x^2])
    @test limit(cs, x => 0) == 0
    @test limit(cs, x => 1) == 1
    @test limit(cs, x => -Inf) == Sym(Inf)
    jump = SymJ.Cases{Sym}([Any[x ≪ 1] => Sym(1), Any[x ≫ 1] => Sym(2)])
    @test limit(jump, x => 1; dir = :-) == 1
    @test limit(jump, x => 1; dir = :+) == 2
    @test_throws SymJ.GruntzError limit(jump, x => 1)
    @test limit(SymJ.Cases{Sym}([Any[x^2 ≫ 1] => Sym(5), Any[x^2 ≪ 1] => Sym(6)]), x => 1; dir = :+) == 5
    @test limit(jump, x => b; dir = :+, assume = @assume(b > 2)) == 2
    pc = SymJ.Cases{Sym}([Any[a ≪ 0] => x * a, Any[a ≫ 0] => x / a])       # cases on another symbol kept
    lc = limit(pc, x => 2)
    @test lc isa SymJ.Cases
    @test lc(a => 4) == Sym(1 // 2)
    # subs: into answers and conditions; decided conditions drop their case or themselves
    cb = SymJ.Cases{Sym}([Any[a ≦ 0] => Sym(Inf), Any[a ≫ 0] => 1 / a])
    @test subs(cb, a => 2) == Sym(1 // 2)
    @test subs(cb, a => -1) == Sym(Inf)
    @test subs(cb, a => b^2 + 1) == 1 / (b^2 + 1)                          # b real: b² + 1 > 0
    @test subs(jump, x => 2) == 2
    @test subs(jump, a => 3) isa SymJ.Cases
end

@testset "piecewise, heaviside, dirac" begin
    P, H, D = SymJ.piecewise, SymJ.heaviside, SymJ.dirac
    @syms a
    f = P(x ≪ 0 => -x, true => x^2)
    @test f(x => -2) == 2 && f(x => 3) == 9 && f(x => 0.5) == 0.25 && N(f(x => Sym(-1) // 3)) ≈ 1 / 3
    @test string(f) == "piecewise(x ≪ 0 => -x, true => x^2)"
    @test diff(f, x) == P(x ≪ 0 => -1, true => 2x)
    # built: false pieces dropped, a true one ends the list, equal neighbours merged
    @test P(x ≪ 0 => 1, true => 1) == 1 && P(Sym(-1) ≪ 0 => 3, true => 4) == 3
    @test P(Sym(1) ≪ 0 => 3, x ≪ 2 => x, true => 0) == P(x ≪ 2 => x, true => 0)
    @test isnan(N(P(x ≪ 0 => 1)(x => 1)))
    @test P(1 - x ≪ 0 => 1, true => 0) == P(x ≫ 1 => 1, true => 0)        # solved for x
    g = P((0 ≦ x, x ≪ 1) => x, true => 0)
    @test g(x => 1 // 2) == 1 // 2 && g(x => 2) == 0 && g(x => -1) == 0
    # lambdify: only the piece taken is evaluated
    l = SymJ.lambdify(P(x ≪ 0 => log(-x), true => sqrt(x)))
    @test l(-ℯ) ≈ 1 && l(4.0) ≈ 2
    # heaviside, dirac, sign
    @test H(Sym(-2)) == 0 && H(Sym(3)) == 1 && H(Sym(0)) == 1 // 2 && H(Sym(0), 1) == 1
    @test diff(H(x - 1), x) == D(x - 1) && diff(D(x), x) == D(x, 1) && diff(sign(x), x) == 2D(x)
    @test D(Sym(3)) == 0 && D(x, 0) == D(x)
    # rewriting, cases
    @test P(abs(x)) == P(x ≪ 0 => -x, true => x)
    @test all(v -> P(max(x, 1))(x => v) == max(v, 1) && P(min(x, 1))(x => v) == min(v, 1), (-2, 1, 3 // 2))
    @test P(H(x)) == P(x ≪ 0 => 0, (x ~ 0) => 1 // 2, true => 1)
    @test P(SymJ.summation(a^Sym(:k), (Sym(:k), 0, Inf))) == P((a ≫ -1, a ≪ 1) => 1 / (1 - a))
    # refine and limits
    @test refine(P(x ≪ 0 => -x, true => x); assume = @assume(x > 1)) == x
    @test refine(H(x - 1); assume = @assume(x > 2)) == 1
    @test limit(f, x => 0; dir = :+) == 0 && limit(P(x ≪ 1 => x, true => 2 - x), x => 1) == 1
    @test limit(H(x), x => 0; dir = :-) == 0 && limit(f, x => Inf) == Sym(Inf)
    @test_throws SymJ.GruntzError limit(H(x), x => 0)
    # conditions free of x split on the parameters (SymPy #30249)
    @test limit(P(Ne(a, 0) => 0, true => x), x => a; dir = :+) == P(Ne(a, 0) => 0, true => a)
    @test limit(P(a ≫ 0 => sin(x) / x, true => x^2), x => 0) == P(a ≫ 0 => 1, true => 0)
    # folding into one piecewise function; simplify does it when smaller
    F = SymJ.piecewise_fold
    @test F(abs(x) + abs(x - 1), x) == P(x ≦ 0 => 1 - 2x, x ≦ 1 => 1, true => 2x - 1)
    g = P(a ≫ 1 => a - 1 // 2, true => a^2 - a + 1 // 2) - P(a ≫ 0 => 0, true => a^2)
    @test F(g) == P(a ≦ 0 => 1 // 2 - a, a ≦ 1 => a^2 - a + 1 // 2, true => a - 1 // 2) == simplify(g)
    @test F(P(x ≪ 0 => -x, true => x) * P(x ≪ 1 => 0, true => x - 1), x) == P(x ≦ 1 => 0, true => x * (x - 1))
    @test F(H(x) - H(x - 1), x) == P(x ≪ 0 => 0, x ≦ 0 => 1 // 2, x ≪ 1 => 1, x ≦ 1 => 1 // 2, true => 0)
    @test F(max(x, 0) + min(x, 1), x) == P(x ≦ 0 => x, x ≦ 1 => 2x, true => x + 1)
    @syms b
    cs = F(H(x - a) - H(x - b), x)                                            # the order of a, b: by cases
    @test cs isa SymJ.Cases && length(cs) == 3 && cs[a => 1, b => 1].second == 0
    cs = F(abs(x - a) + abs(x - b), x)
    @test cs(a => 0, b => 1) == F(abs(x) + abs(x - 1), x) && cs(a => 2, b => 1) == P(x ≦ 1 => 3 - 2x, x ≦ 2 => 1, true => 2x - 3)
    @test F(abs(x - a) + abs(x), x; assume = @assume(a > 0)) == P(x ≦ 0 => a - 2x, x ≦ a => a, true => 2x - a)
    @test F(H(x - a) - H(x - b), x; assume = @assume(a < b)) == P(x ≪ a => 0, x ≦ a => 1 // 2, x ≪ b => 1, x ≦ b => 1 // 2, true => 0)
    @test_throws ArgumentError F(P(x ≪ 0 => 1, true => 2) + P(a ≪ 0 => 1, true => 2))   # which variable?
    @test F(x^2 + 1, x) == x^2 + 1
end

@testset "residue" begin
    @syms a
    R = SymJ.residue
    @test R((x^3 + 5) / ((x^4 - 1) * (x + 1)), x, -1) == -Sym(9) // 4      # SymPy #2519
    @test R(1 / sin(x), x, 0) == 1
    @test R(cot(x), x, PI) == 1
    @test R(exp(x) / x^3, x, 0) == Sym(1) // 2
    @test R(1 / (x^2 + 1), x, IM) == -IM / 2
    @test R(exp(x) / (x - a)^2, x, a) == exp(a)
    @test R(sin(x), x, 0) == 0                                               # analytic there
    # at infinity: -Res_{t=0} f(1/t)/t²
    @test R(1 / x, x, Inf) == -1
    @test R(x^2 / (x^3 + 1), x, Inf) == -1
    # no residue at a branch point
    @test_throws ArgumentError R(1 / sqrt(x), x, 0)
    @test_throws ArgumentError R(log(x) / x, x, 0)
end

@testset "unevaluated definite integrals: values by quadrature" begin
    @syms t
    Int4(f, v, a, b) = SymJ.Integral(f, (v, a, b))
    # SymJ.Integral builds the node, nothing integrated; several ranges nest, the first inside
    @test SymJ._isdefint(SymJ.Integral(x^2, (x, 0, 1)))
    @test SymJ.Integral(x * t, (x, 0, 1), (t, 0, 2)) == SymJ.Integral(SymJ.Integral(x * t, (x, 0, 1)), (t, 0, 2))
    @test length(SymJ.Integral(x^2, x).args) == 2
    @test_throws ArgumentError SymJ.Integral(x^2, (x^2, 0, 1))
    @test SymJ._quadgk(u -> exp(-u^2), -Inf, Inf)[1] ≈ sqrt(π) rtol = 1e-12
    @test SymJ._quadgk(u -> 1 / sqrt(u), 0, 1)[1] ≈ 2 rtol = 1e-8
    I1 = Int4(exp(x^3), x, 0, 1)
    @test N(I1) ≈ 1.3419044179774198 rtol = 1e-12
    @test N(Int4(exp(x^3), x, 1, 0)) ≈ -1.3419044179774198 rtol = 1e-12
    @test N(Int4(exp(-x^2), x, -Inf, Inf)) ≈ sqrt(π) rtol = 1e-12
    @test N(Int4(exp(x^3) / sqrt(x), x, 0, 1)) ≈ 2.3841376158 rtol = 1e-9   # Σ 1/(n!(3n + 1/2))
    # complex limits: along the segment (∫₀ⁱ e^(x²) = i ∫₀¹ e^(-t²))
    @test N(Int4(exp(x^2), x, 0, im)) ≈ 0.746824132812427im rtol = 1e-12
    # a parameter: by calling, by substituting; nested, the outer variable's value inside
    J = Int4(exp(x^3 * t), x, 0, 1)
    @test J(2.0) ≈ N(subs(J, t => 2)) rtol = 1e-14
    @test N(Int4(J, t, 0, 1)) ≈ 1.1537157127366073 rtol = 1e-10
    # a constant for float, calling, conversions; BigFloat to the precision asked
    G = Int4(exp(-x^2), x, -Inf, Inf)
    @test float(2G) ≈ 2sqrt(π) rtol = 1e-12
    @test G() ≈ sqrt(π) rtol = 1e-12
    @test Float64(G) ≈ sqrt(π) rtol = 1e-12
    g = SymJ.evalf(G, BigFloat, 128)
    @test precision(g) == 128
    @test abs(g - setprecision(() -> sqrt(BigFloat(π)), BigFloat, 128)) < BigFloat(2)^-120
    @test abs(SymJ.evalf(I1, BigFloat, 200) - setprecision(() -> sum(k -> 1 / (factorial(big(k)) * (3k + 1)), 0:80), BigFloat, 200)) < BigFloat(2)^-190
    @test_throws ArgumentError N(Int4(1 / x, x, -1, 1))
end

@testset "series with a power factor (exponent not rational)" begin
    @syms a
    s = series(x^a * exp(x), x, 0, 3)
    @test isequal(s.pre, x^a)
    @test Sym(s) == x^a * (1 + x + x^2 / 2)
    @test Sym(series(x^PI * sin(x), x, 0, 6)) == x^PI * (x - x^3 / 6 + x^5 / 120)
    @test Sym(series(sin(x)^PI, x, 0, 4)) == x^PI * (1 - PI * x^2 / 6)
    @test isequal(series(log(x)^a, x, 1, 3).pre, (x - 1)^a)          # SymPy #7168
    @test Sym(series(x^a / (1 + x), x, Inf, 3)) == x^a * (1 / x - 1 / x^2)
    @test isequal(Sym(series(x^sqrt(Sym(2)), x, 0, 3)), x^sqrt(Sym(2)))
    # diff, products, a rational exponent substituted into the ordinary series
    @test iszero(simplify(Sym(diff(s, x)) - expand(diff(x^a * (1 + x + x^2 / 2), x))))
    @test isequal((s * series(x^a, x, 0, 3)).pre, x^(2a))
    @test Sym(subs(s, a => Sym(1) // 2)) == sqrt(x) + x^(Sym(3) // 2) + x^(Sym(5) // 2) / 2
    @test_throws ArgumentError s + series(exp(x), x, 0, 3)
end

@testset "series with irrational exponents (SymPy #21358)" begin
    r2 = sqrt(Sym(2))
    @syms w a
    @test isequal(Sym(series(exp(w^(1 / PI)), w, 0, 1)), 1 + w^(1 / PI) + w^(2 / PI) / 2 + w^(3 / PI) / 6)
    @test string(series(x + x^r2, x)) == "x + x^sqrt(2) + O(x^6)"
    @test string(series(1 / (1 + x^PI), x, 0, 7)) == "1 - x^π + x^(2*π) + O(x^7)"
    @test string(series(x^2 + x^PI, x, 0, 3)) == "x^2 + O(x^3)"                   # π > 3
    @test isequal(Sym(series((1 + x^r2)^2, x, 0, 3)), 1 + 2x^r2 + x^(2r2))
    @test isequal(Sym(series(sin(x)^r2, x, 0, 4)), x^r2 * (1 - r2 * x^2 / 6))     # (a factor)
    # ordered by value: log 3/log 5 < 1
    @test isequal(series(x + x^(log(Sym(3)) / log(Sym(5))), x, 0, 2).ps.c[1], Sym(1))
    @test series(x + x^(log(Sym(3)) / log(Sym(5))), x, 0, 2).ps.e[2] == 1
    @test string(series(x^r2 + x, x, Inf, 3)) == "x^sqrt(2) + x + O(1/x^3)"
    @test string(series(x^r2 - x^(2 / r2) + x, x)) == "x"                          # 2/√2 = √2
    # a factor x^π folded into the exponents when added to a series without it
    @test isequal(Sym(series(x^PI * exp(x), x, 0, 2) + series(x, x, 0, 2)), x + x^PI + x^(PI + 1))
    @test isequal(Sym(diff(series(x + x^r2, x, 0, 3), x)), 1 + r2 * x^(r2 - 1))
    # numerically: 1/(1 + x^√2) at x = 1/100
    v = Sym(series(1 / (1 + x^r2), x, 0, 5))
    @test abs(N(subs(v, x => Sym(1 // 100))) - 1 / (1 + 0.01^sqrt(2))) < 1e-9
    # a symbolic exponent: which term leads depends on a
    @test_throws ArgumentError series(x^a + x^2, x)
end

