# Regressions from the SymPy issue reviews (included by runtests.jl)

@testset "crashes and hangs (SymPy review, round 3)" begin
    @syms t n x y()
    # NaN has no sign: compared as IEEE does, not an InexactError (SymPy #27440)
    @test refine(SymJ.Relation(:lt, Sym(NaN), Sym(0))) == false
    @test refine(SymJ.Relation(:gt, Sym(NaN), Sym(0))) == false
    # large multiples of an angle are not expanded (a polynomial of degree 316801; #27293)
    e = 4sin(316801t) - 4sin(316800t)
    @test simplify(e) == e
    @test simplify(sin(2x) / (2sin(x))) == cos(x)
    # constant coefficients: the characteristic roots, or the recurrence itself when they
    # are large radicals (Hyper took minutes; #27902)
    @test SymJ.rsolve(y(n + 2) - 2y(n + 1) + y(n) ~ 0, y(n), y(0) => 1, y(1) => 3) == 2n + 1
    r = SymJ.rsolve(y(n + 4) + 6y(n + 2) - y(n + 1) - y(n) ~ 0, y(n), y(0) => 1, y(1) => 0, y(2) => 0, y(3) => 0)
    @test [subs(r, n => k) for k in 6:9] == [-6, 1, 37, -12]
    # a coefficient's huge power stays a power; N past Float64's range (#27278)
    ex = (1 / (n * log((n + 3) / n)))^n * ((n + 1) * log((n + 4) / (n + 1)))^(n + 1) / 3
    v = subs(ex, n => Sym(10)^9)
    @test length(string(v)) < 200
    @test N(v) ≈ 1
    @test N(subs(ex, n => 1000)) ≈ 0.9999981335936851
end

@testset "missing cases (SymPy review, round 3)" begin
    @syms x y b c a X Y t θ() k
    # floats with symbols: solved with the rationals they stand for (#27528, #27712)
    @test solve(y ~ c + 0.5 * b, c) == [-0.5 * b + y]
    @test N(only(solve(175 ~ 1.5 * log(x), x))) ≈ exp(350 / 3)
    @test N(only(solve(1.5 ~ exp(x), x))) ≈ log(1.5)
    # sign(u) = 0 is u = 0 (#27322); products of exponentials by logarithms (#27330)
    @test solve(sign(1 - 1 / x^2) / x^3 ~ 0, x; assume = (x ∈ ℝ,)) == [-1, 1]
    @test sort(N.(solve(2^x * 3^(x^2) ~ 6, x))) ≈ [-log(6) / log(3), 1]
    @test solve(2^x * 5^x ~ 100, x) == [2]
    # a piecewise piece by piece, the roots its condition allows (#27037)
    @test solve(SymJ.piecewise(x ≪ 1 => x, true => 2x) ~ 4, x) == [2]
    @test solve(SymJ.piecewise(x ≪ 0 => -x, true => SymJ.piecewise(x ≪ 2 => x, true => x^2)) ~ 9, x) == [-9, 3]
    # sums from -∞ (#27824, #27893, #27256)
    @test SymJ.summation(1 / (k^2 + 1), k, -Inf, 0) == SymJ.summation(1 / (k^2 + 1), k, 0, Inf)
    @test N(SymJ.summation(1 / (k^2 + 1), k, -Inf, Inf)) ≈ π * coth(π)
    # relations: one side, real and imaginary parts, terms of one sign (#27754, #27623, #27759, #27755)
    @test simplify(y - x ~ y + x) == (x ~ 0)
    @test simplify(x + y ≧ x - y) == (y ≧ 0)
    @test string(simplify(X ~ im * Y; assume = (X ∈ ℝ, Y ∈ ℝ))) == "X ~ 0 & Y ~ 0"
    @test simplify(X ~ -X; assume = (X ≧ 0,)) == (X ~ 0)
    # signs by cases: floor, abs, max (#27564, #27364)
    @test refine(x - floor(x) ≧ 0; assume = (x ∈ ℝ,)) == true
    @test refine(1 - (x - floor(x)) ≫ 0; assume = (x ∈ ℝ,)) == true
    @test refine(a ≦ abs(a); assume = (a ∈ ℝ,)) == true
    @test refine(max(a, 1) ≧ a; assume = (a ∈ ℝ,)) == true
    # nested floors (#27400); x + ∞ for real x (#27605)
    @test floor(floor(k / 2) / 3) == floor(k / 6)
    @test ceil(ceil(x) / 3) == ceil(x / 3)
    @test refine((x + Sym(Inf))^2; assume = (x ∈ ℝ,)) == Inf
    # simplify: a sum's trig part alone (#27356); integer float exponents (#27706, #27796)
    @test simplify(sin(x)^2 + cos(x)^2 + (x + y)^3) == (x + y)^3 + 1
    @test simplify(x^2 + x^2.0) == 2x^2
    @test simplify((1.0x - y)^2.0 - (x - y)^2) == 0
    # lambdify with function calls as variables (#27006)
    let u = θ(t), du = diff(θ(t), t)
        @test lambdify(du^2 + u, u, du)(1.0, 2.0) == 5.0
    end
    # A' of a Matrix{Sym} is one of Sym
    @test collect(Sym[1 x; 0 2]') isa Matrix{Sym}
    # a piece at a point with the value there of the piece that takes over (#27627)
    @test SymJ.piecewise(Eq(x, 999) => 1, true => 1000 - x) == 1000 - x
    @test SymJ.piecewise(Eq(x, 1) => 0, x ≪ 0 => x, true => x - 1) == SymJ.piecewise(x ≪ 0 => x, true => x - 1)
    @test length(SymJ.piecewise(Eq(x, 2) => 5, true => 1000 - x).args) == 4
    @test length(SymJ.piecewise(Eq(x, 0) => 1, true => sin(x) / x).args) == 4     # needed at 0
    # a k-free denominator in rsolve's right side (#27901)
    let (n, k) = symbols("n k"), yf = SymFunction(:yf)
        rhs = binomial(n + 1, k) * 2^(-n - 1) - binomial(n, k) * 2^(-n)
        r = SymJ.rsolve(yf(k + 1) - yf(k) ~ rhs, yf(k))
        res = subs(r, k => k + 1) - r - rhs
        @test all(iszero, [subs(res, n => i, k => j) for i in 3:5 for j in 0:3])
    end
    # cube roots of quadratic surds denested (#27082)
    @test simplify(cbrt(85 + 62sqrt(Sym(7)))) == 1 + 2sqrt(Sym(7))
    @test simplify(cbrt(-85 - 62sqrt(Sym(7)))) == -1 - 2sqrt(Sym(7))
    @test SymJ.simplify_radicals(cbrt(2 + sqrt(Sym(5)))) == (1 + sqrt(Sym(5))) / 2
    @test simplify((26 + 15sqrt(Sym(3)))^(Sym(1) // 3)) == 2 + sqrt(Sym(3))
    @test simplify(cbrt(1 + sqrt(Sym(2)))) == cbrt(1 + sqrt(Sym(2)))                 # none
    # atan at ±∞ (#27109)
    @test Sym(series(atan(x), x, Inf, 4)) == PI / 2 - 1 / x + 1 / (3x^3)
    @test Sym(series(atan(x), x, -Inf, 4)) == -PI / 2 - 1 / x + 1 / (3x^3)
    @test Sym(series(atan(1 / x^2), x, 0, 4)) == PI / 2 - x^2
    @test_throws ArgumentError series(atan(1 / x), x, 0, 4)                       # by side
    # no zeros by bounds: 1 + x^2 sin(x)^2 ≥ 1 (#27139)
    @test isempty(SymJ._zeros_in(1 + x^2 * sin(x)^2, x, Sym(0), Sym(1)))
    # floor and mod: integer parts of their shifts out
    @test floor(x + Sym(7) // 4) == floor(x + Sym(3) // 4) + 1
    @test mod(x - PI / 2, 2PI) == mod(x + 3PI / 2, 2PI)
    # an unknown function of x: no limit; one free of x a constant (#27143)
    @test_throws ArgumentError limit(f(x), x => a)
    @test_throws ArgumentError limit(diff(f(x), x), x => a)
    @test limit(f(2) + x, x => a) == a + f(2)
end

@testset "wrong answers (SymPy review, round 4)" begin
    @syms y k1 k0 r G E1 E2 TT t s
    u = SymFunction(:u)
    # constant factors of an unknown function's transform kept (#25815)
    @test SymJ.laplace_transform(k1 * u(t), t, s) == k1 * SymJ.apply_function(:LaplaceTransform, Sym[u(t), t, s])
    # floor near an integer decided in BigFloat (#26368)
    @test floor(log(2, Sym(2)^360 - 1) + 1) == 360
    @test ceil(log(2, Sym(2)^360 - 1)) == 360
    # a float linear system: a rounding residue is not 0 = nonzero (#26223)
    @test N(only(solve([3.0y ~ 5, -y ~ -1.6666666666666667], [y]))[y]) ≈ 5 / 3
    @test N(only(solve([0.1x + 0.2y ~ 0.3, x - y ~ 0], [x, y]))[x]) ≈ 1
    # a tiny float coefficient: candidates checked at exact sample values (#25792)
    @test length(solve(G ~ 2.7^((E2 - E1) / (1.38e-23 * TT)), E1)) == 1
    # floats with parameters: cases of the exact equation (#26820)
    @test solve(r ~ (1.0k0 + 1.0k1) / (2.0(0.5k0 + 0.5k1)), r) == [1.0]
    # x outside sin and cos as well: not supported, not "no solution" (#26253)
    @test_throws ArgumentError solve(-2(15 - x) - 2(3 - sin(x)) * cos(x) ~ 0, x)
end

@testset "crashes and hangs (SymPy review, round 4)" begin
    @syms a s t df v1 v2 n1 n2
    # one real root with a negative Cardano radicand: the real cube root (#26069)
    S = solve(4x^3 - 4x^2 + x - 2 ≦ 0, x)
    @test S(1) == true
    @test S(Sym(6) // 5) == false
    @test solve(x^3 + 3x^2 + 3x - 1 ≧ 0, x) == solve(x ≧ cbrt(Sym(2)) - 1, x)
    @test N(solve(4x^3 - 4x^2 + x - 2 ~ 0, x)[1]) isa Real
    # Poly of zero (#25899)
    @test iszero(Poly(a - a))
    # inverse transforms through a cubic's roots: no cancel over its radicals (#26560, #25796)
    r = SymJ.inverse_laplace_transform((2s + 1) / (s * (s^3 + s^2 + 3s + 1)), s, t)
    @test real(N(subs(r, t => Sym(1) // 2))) ≈ 0.21788881804270455
    r = SymJ.inverse_laplace_transform(5 / ((s^3 + Sym(9) // 2 * s^2 + Sym(13) // 2 * s + 8) * s), s, t)
    @test real(N(subs(r, t => 1))) ≈ 0.27649212810980195
    # Welch's degrees of freedom for n1: seconds, not minutes (#25920)
    eq = df ~ (v1 / n1 + v2 / n2)^2 / (v1^2 / (n1^2 * (n1 - 1)) + v2^2 / (n2^2 * (n2 - 1)))
    @test @elapsed(solve(eq, n1)) < 30
end

@testset "missing cases (SymPy review, round 4)" begin
    # square factors above 1000²: 1009² and 51613² (#26903, #26317)
    @test sqrt(Sym(1009^2 * 1013)) == 1009sqrt(Sym(1013))
    @test sqrt(Sym(1001885125670250202210728125403128548905806811601)) == 51613sqrt(Sym(376096873138962280636080063132060838729))
    @test SymJ.minimal_polynomial(sqrt(Sym(1031316053)) - 1009sqrt(Sym(1013)), x) == x
end

@testset "smaller missing cases (SymPy review, round 4)" begin
    @syms a b c u s t k01 k11 Q L
    # constants: trigonometric identities, logarithms to one base, powers of i (#25705, #25860, #26140)
    @test simplify(cos(Sym(2))^2 + sin(Sym(2))^2) == 1
    @test simplify(log(Sym(10), Sym(5)) + log(Sym(10), Sym(2))) == 1
    @test simplify(log(Sym(10), Sym(5))^2 + log(Sym(10), Sym(5)) * log(Sym(10), Sym(2)) + log(Sym(10), Sym(2))) == 1
    @test simplify(sqrt(-IM * Sym(1)) + (IM * Sym(1))^(Sym(3) // 2)) == 0
    # trigonometric functions of inverse ones, and integer multiples of π (#26373, #25685)
    @test cos(asin(Sym(11) // 61)) == Sym(60) // 61
    @test simplify(cos(2acos(Sym(3) // 5) - 4PI)) == -Sym(7) // 25
    @test simplify(sin(2asin(Sym(11) // 61))) == Sym(1320) // 3721
    @test sin(x + 2PI) == sin(x)
    # a product's trigonometric factors of one argument (#26835)
    @test simplify(2sin(a / 2) * cos(a / 2) * cos(b) * cos(c)) == sin(a) * cos(b) * cos(c)
    # real and imaginary parts, hyperbolic functions as exponentials (#25564)
    @test simplify(real(a) + IM * imag(a)) == a
    w = (x + b) / u
    @test simplify(exp(2w) / 4 + exp(-2w) / 4 - sinh(w)^2 - Sym(1) // 2) == 0
    @test simplify((exp(k01 * t) - exp(k11 * t)) * exp(-t * (k01 + k11)) / (k01 - k11)) == (exp(-k11 * t) - exp(-k01 * t)) / (k01 - k11)
    # radicals: cube roots with a norm that is not a cube, roots in a denominator (#26857, #25915)
    @test simplify(cbrt(Sym(2)) * (cbrt(22 + 9sqrt(Sym(6))) - cbrt(9sqrt(Sym(6)) - 22))) == 4
    @test simplify((sqrt(Sym(6)) + sqrt(Sym(2))) / sqrt(sqrt(Sym(3)) + 2)) == 2
    # refine: powers, dirac, conditions given as facts (#26006, #25724, #26691, #25522)
    @test refine((b^x)^(log(a) / log(b)); assume = (a ≫ 0, b ≫ 0, x ∈ ℝ)) == a^x
    @test refine(3 * (-x^3)^(Sym(4) // 3) / (4x); assume = (x ≪ 0,)) == 3x^3 / 4
    @test refine((-x)^(Sym(4) // 3); assume = (x ≪ 0,)) == (-x)^(Sym(4) // 3)
    @test refine(SymJ.dirac(x - 1); assume = (x ≫ 2,)) == 0
    @test refine(SymJ.piecewise(Ne(x, 1) => x - 1, true => Sym(1)); assume = (x ≶ 1,)) == x - 1
    # |u^r| = |u|^r for real u, r: in solve and limit (#25779, #26513)
    @test solve(abs(sqrt(x)) ~ 2, x; assume = (x ∈ ℝ,)) == [-4, 4]
    @test limit(abs((-x / (x + 1))^x), x => Inf) == exp(-Sym(1))
    # an impossible piece is dropped (#25930)
    @test SymJ.piecewise((x ≦ 0, x ≫ PI) => x, (x ≧ 0, x ≦ PI) => 2x) == SymJ.piecewise((x ≧ 0, x ≦ PI) => 2x)
    # floats: with symbols in matrices, in an inverse transform (#26271, #26517)
    @test SymJ.LinearAlgebra.det(Sym[x-0.3 -0.06; -0.2 x-0.5]) == x^2 - 0.8x + 0.138
    r = SymJ.inverse_laplace_transform((20000.0s^2 + 1600.0s + 30.0) / (s * (20000.0s^3 + 5600.0s^2 + 266.0s + 3.0)), s, t)
    @test N(subs(r, t => 1000)) ≈ 9.999999770332835
    # a call form error with keywords: the message, not a MethodError (#25783)
    @test_throws ArgumentError solve(x ~ Q, [Q, L]; assume = (Q ≫ 0,))
end

@testset "bigger missing cases (SymPy review, round 4)" begin
    @syms y s t Ec
    u = SymFunction(:u)
    # Lambert W: k-th powers (y² = e^y as y = ±e^(y/2)), log(x) with x through x = e^y (#26461)
    @test N(solve(y^2 ~ exp(y), y)[1]) ≈ -0.7034674224983917
    @test any(r -> (v = N(r); v isa Real && v ≈ 0.4948664145165307), solve(1 / x - 1 / log(x)^2 ~ 0, x))
    @test N(only(solve(x * log(x) ~ 1, x))) ≈ 1.7632228343518968
    # numeric W of a complex argument (LambertW fails on -0.5 + 0im)
    w = N(lambertw(-Sym(1) // 2))
    @test w * exp(w) ≈ -0.5
    # Dirac derivatives in Laplace transforms (#26775)
    @test SymJ.inverse_laplace_transform(s, s, t) == SymJ.dirac(t, 1)
    @test SymJ.laplace_transform(SymJ.dirac(t - 1, 2), t, s) == s^2 * exp(-s)
    @test SymJ.laplace_transform(t^2 * SymJ.dirac(t - 2, 1), t, s) == (4s - 4) * exp(-2s)
    # real and imaginary parts of a quotient with real symbols (#26102)
    @test refine(imag(x^2 / (IM * y + Ec)); assume = (x ∈ ℝ, Ec ∈ ℝ, y ∈ ℝ)) == -x^2 * y / (y^2 + Ec^2)
    # a linear ODE that is an exact derivative, integrated once per order (#25479)
    sol = SymJ.dsolve(diff((1 + x)^3 * diff(u(x), x, 2), x, 2) ~ 0, u(x))
    @test iszero(simplify(diff((1 + x)^3 * diff(rhs(sol), x, 2), x, 2)))
end

@testset "wrong answers (SymPy review, round 5)" begin
    @syms a b s t
    # has sees the powers a product holds as factors (-im is a product; #25162)
    @test has(-IM, IM)
    @test has(x * sqrt(y), sqrt(y))
    sv = SymJ.LinearAlgebra.svdvals(Sym[0 IM; -IM IM])
    @test iszero(simplify(sv[1] - (1 + sqrt(Sym(5))) / 2))
    @test iszero(simplify(sv[2] - (sqrt(Sym(5)) - 1) / 2))
    # a float matrix: LAPACK's eigenvectors (inexact roots left empty null spaces; #24057)
    A = Sym[0.0 0 0; 0.92 0.86 0; 0 0.08 0.8]
    @test size(SymJ.LinearAlgebra.eigvecs(A), 2) == 3
    F = SymJ.LinearAlgebra.eigen(Sym.([0.0 0 0 0.45; 0.92 0.86 0 0; 0 0.08 0.8 0.83; 0 0.02 0.19 0]))
    @test N.(F.values) ≈ [-0.15445857152658957, -0.002629910195547378, 0.807969994528574, 1.009118487193563]
    # exp(c s) with c > 0 shifts to t < 0: no inverse (#24595)
    @test_throws ArgumentError SymJ.inverse_laplace_transform(exp(2s), s, t)
    @test_throws ArgumentError SymJ.inverse_laplace_transform(exp(a * s), s, t)
    @test SymJ.inverse_laplace_transform(exp(-2s) / s, s, t) == SymJ.heaviside(t - 2)
    # an infinite modulus with a constant: its float value, not sqrt(2) + NaN
    @test mod(sqrt(Sym(2)), Sym(Inf)) == Sym(sqrt(2))
    # a jump in a piecewise function is a dirac in its derivative, as for heaviside (#24127)
    @test diff(SymJ.piecewise(x ≪ 0 => 0, true => 1), x) == SymJ.dirac(x)
    @test diff(SymJ.piecewise(x ≪ 1 => x^2, true => 2x + 3), x) == 4SymJ.dirac(x - 1) + SymJ.piecewise(x ≪ 1 => 2x, true => 2)
    @test diff(SymJ.piecewise(x ≪ 0 => -x, true => x), x) == SymJ.piecewise(x ≪ 0 => -1, true => 1)
    # a vanishing coefficient with an irrational constant is a case (#24344)
    cs = solve((a - sqrt(Sym(3)) * b) * x^2 + x + 1 ~ 0, x)
    @test first(cs.cases).second == [-1]
    @test first(cs.cases).first == Any[a ~ sqrt(Sym(3)) * b]
end

@testset "crashes and hangs (SymPy review, round 5)" begin
    @syms a b k n s t w L1 L2 L3 C1 C2 C3 M23 RL U I1 I2 I3 I4
    LA = SymJ.LinearAlgebra
    # 0^(n - k) summed the other way (#23952)
    @test SymJ.summation(Sym(0)^(n - k), k, 0, n) == 1 - Sym(0)^(n + 1)
    # factorial at a float through Γ; an integer to a negative power (#24246)
    @test lambdify(factorial(x))(2.5) ≈ 3.3233509704478426
    @test lambdify(Sym(2)^x)(-3) == 0.125
    # a removable singularity at a dirac's point: the limit there (#25293)
    @test SymJ.laplace_transform(sin(4PI * (t - 1)) / (t - 1) * SymJ.dirac(t - 1), t, s) == 4PI * exp(-s)
    # a pattern's slot in a piecewise condition is shown, not decided (#25328)
    pw = SymJ.piecewise(x ≫ 7 => 15 - x, true => 7 - x)
    @test replace(pw, :(piecewise(~a, ~b, ~c, ~d)) => :(piecewise(~a, ~b, ~c, ~d))) == pw
    # sin of a number about 2^762451794: an error, not hours (#24279)
    @test_throws SymJ.HugeArgumentError N(sin(sin(exp(exp(exp(Sym(3))))) + 3))
    # large powers of sums of radicals expanded first (#25109)
    @test simplify(((-3 + 2sqrt(Sym(2)))^2012) * ((3 + 2sqrt(Sym(2)))^2011)) == 3 - 2sqrt(Sym(2))
    # a system of products split by factor (#24500)
    @syms c1 c3 g1
    sols = solve([x * (x - 1) * (c1 - y) ~ 0, y * (y - 1) * (c3 + g1 * x) ~ 0], [x, y])
    @test length(last(sols.cases).second) == 5
    # linear systems by Cramer's rule: a generic 4 × 4 (#24679), a circuit with im (#24135)
    M = Sym[Sym(Symbol(:m, i, j)) for i in 1:4, j in 1:4]
    bb = Sym[Sym(Symbol(:b, i)) for i in 1:4]
    X = M \ bb
    ps = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71)
    vals = [e => Sym(ps[i]) for (i, e) in enumerate([vec(M); bb])]
    ev(e) = N(subs(e, vals...))
    @test maximum(abs.(ev.(M) * ev.(X) - ev.(bb))) < 1e-10
    eqs = [I1 * (IM * L1 * w - IM / (C1 * w)) + I2 * IM / (C1 * w) ~ U,
           I2 * (IM * L2 * w - IM / (C2 * w) - IM / (C1 * w)) + I3 * (IM * M23 * w + IM / (C2 * w)) + I1 * IM / (C1 * w) ~ 0,
           I2 * (IM * M23 * w + IM / (C2 * w)) + I3 * (IM * L3 * w - IM / (C3 * w) - IM / (C2 * w)) + I4 * IM / (C3 * w) ~ 0,
           I4 * RL + I3 * IM / (C3 * w) - I4 * IM / (C3 * w) ~ 0]
    sol = only(last(solve(eqs, [I1, I2, I3, I4]).cases).second)
    pv = [L1 => 1, L2 => 2, L3 => 3, C1 => 5, C2 => 7, C3 => 11, M23 => 13, RL => 17, U => 19, w => 23]
    @test maximum(abs(N(subs(subs(e.lhs - e.rhs, sol...), pv...))) for e in eqs) < 1e-10
    # rational entries: the rows cleared of denominators (#25009)
    R = Sym[1/a 1; 1 1/b]
    @test simplify.(R * inv(R)) == Sym[1 0; 0 1]
    # a tower of radicals shown in linear time (#24565)
    r = sqrt(Sym(2))
    for _ in 1:12
        r = sqrt(r + sqrt(r^2 + 1) + 2)
    end
    @test length(string(r)) > 0
    # exact eigenvectors for quartic roots through the adjugate (#24057)
    A = Matrix(Sym.(rationalize.([0.0 0 0 0.45; 0.92 0.86 0 0; 0 0.08 0.8 0.83; 0 0.02 0.19 0]))')
    V, λ = LA.eigvecs(A), LA.eigvals(A)
    @test maximum(abs.(N.(A * V[:, 1] - λ[1] * V[:, 1]))) < 1e-12
end

@testset "missing cases (SymPy review, round 5)" begin
    @syms a b c k n s t z L
    LA = SymJ.LinearAlgebra
    # powers of positive constants merged; perfect-power bases reduced (#24090)
    @test (Sym(3)^(Sym(1) // 3))^(Sym(1) // 4) == Sym(3)^(Sym(1) // 12)
    @test refine(sqrt(Sym(256)^x); assume = (x ∈ ℝ,)) == Sym(16)^x
    # exp(iπr) reduced (#24082)
    @test exp(10 * IM * PI) == 1
    @test exp(2 * PI * IM / 5)^6 == exp(2 * PI * IM / 5)
    # max and min: common factors and terms out (#24578)
    @test simplify(max(2a, 2b)) == 2max(a, b)
    @test simplify(max(a + c, b + c)) == c + max(a, b)
    @test simplify(max(-a, -b)) == -min(a, b)
    # floats in apart (#24481)
    @test N(subs(apart((37s + 300) / ((s + 2) * (s^2 + 5.9999s + 5)), s) - (37s + 300) / ((s + 2) * (s^2 + 5.9999s + 5)), s => 3)) ≈ 0 atol = 1e-12
    # constant factors of factor's result combined (#24998)
    @test factor(x - x * exp(Sym(-10))) == (1 - exp(Sym(-10))) * x
    # exp kept in floating point inverse transforms
    @test !occursin("2.718", string(SymJ.inverse_laplace_transform((2.5e3 * s + 1.0e6) / (s * (s^2 + 1.2e3 * s + 4.0e6)), s, t)))
    # Laplace transforms of sin(√t), cos(√t)/√t and integrals (#24309, #24545)
    @test SymJ.laplace_transform(sin(sqrt(t)), t, s) == sqrt(PI) / (2 * s^(Sym(3) // 2)) * exp(-1 / (4s))
    @test SymJ.laplace_transform(cos(2sqrt(t)) / sqrt(t), t, s) == sqrt(PI / s) * exp(-1 / s)
    u = SymFunction(:u)
    @test SymJ.laplace_transform(SymJ.apply_function(:Integral, Sym[u(t), t]), t, s) == SymJ.apply_function(:LaplaceTransform, Sym[u(t), t, s]) / s
    # derivatives of real, imag, angle in a real variable (#24318)
    @test diff(real(x^2), x) == real(2x)
    @test simplify(refine(diff(angle(IM * x), x); assume = (x ∈ ℝ,))) == 0
    # a vanishing base with a positive exponent (#24067)
    @test limit(x^n, x => 0; assume = (n ∈ ℤ, n ≫ 0)) == 0
    # a piecewise limit at a point whose sign the assumptions give (#24042)
    @test limit(x * SymJ.heaviside(x), x => L; assume = (L ≫ 0,)) == L
    # x ≠ 0 in a tuple is the Bool true: an error that says so (#24183)
    @test_throws ArgumentError refine(x^2 ≫ 0; assume = (x ∈ ℝ, x ≠ 0))
    @test refine(x^2 ≫ 0; assume = (x ∈ ℝ, x ≶ 0)) == true
    # inverse trig constants equal to a rational multiple of π; simplify idempotent (#24948)
    @test simplify(acos(2sqrt(Sym(5)) / 5) - acos(-3sqrt(Sym(10)) / 10) + 3PI / 4) == 0
    @test simplify(cos(acos(2sqrt(Sym(5)) / 5) - acos(-3sqrt(Sym(10)) / 10) + 3PI / 4)) == 1
    # cyclotomic factors: roots as cos + i sin (#25192)
    rs = solve(x^7 + 1 ~ 0, x)
    @test !any(r -> occursin(".", string(r)), rs)
    @test maximum(abs(N(r^7 + 1)) for r in rs) < 1e-12
    # sums of squares of real expressions (#24207)
    @test solve((x - 5)^2 + (sin(x) - sin(Sym(5)))^2 ~ 0, x; assume = (x ∈ ℝ,)) == [5]
    # a complex constant in polar form (#24522)
    @test SymJ.rewrite((1 + IM * sqrt(Sym(3))) / 2, exp) == exp(IM * PI / 3)
    # series at branch points (#25175)
    @test Sym(series(acosh(1 + z), z, 0, 3)) == sqrt(Sym(2)) * sqrt(z) - sqrt(Sym(2)) * z^(Sym(3) // 2) / 12 + 3sqrt(Sym(2)) * z^(Sym(5) // 2) / 160
    # rank with cases (#24613)
    @test LA.rank(Sym[2 -1 7; 1 2 2; 0 5 -3; 5 0 5+k]) == 3
    cs = LA.rank(Sym[2 -1 7; 1 2 2; 0 5 -3; 5 0 5+k]; generic = false)
    @test first(cs.cases) == (Any[k ~ 11] => 2)
    # exponentials with a common ratio (#24138)
    @test solve(Sym(4)^x + Sym(6)^x ~ Sym(9)^x, x) == [log((1 + sqrt(Sym(5))) / 2) / log(Sym(3) / 2)]
    # floor: an interval (#25001)
    @test string(solve(floor(x / 2 - Sym(1) // 2) + 1 ~ 0, x)) == "[-1, 1)"
    # x through s² + x² only (#24418)
    @test length(solve(-5t^2 - 5x^2 + 5log(t^2 + x^2) + 7 ~ 0, t)) == 4
    # abs and square roots in systems (#24189, #24868)
    @test solve([abs(x) + abs(y) ~ 2, x - y ~ 1], [x, y]; assume = (x ∈ ℝ, y ∈ ℝ)) == [Dict(x => Sym(3) / 2, y => Sym(1) / 2), Dict(x => -Sym(1) / 2, y => -Sym(3) / 2)]
    f = sqrt((x + 1)^2 + (y - 1)^2 + (z + 1)^2)
    g = x^2 + y^2 + z^2 - 9
    @syms l
    @test length(solve([diff(f, x) ~ l * diff(g, x), diff(f, y) ~ l * diff(g, y), diff(f, z) ~ l * diff(g, z), g ~ 0], [x, y, z, l])) == 2
    # Jordan form with Cardano's trigonometric roots (#24942)
    @test LA.isdiag(SymJ.jordan_form(Sym[-3 1 2; 1 -1 0; 1 0 -2])[2])
    # exp of a scalar multiple: no sqrt(s^2)
    @test !occursin("sqrt(s^2)", string(exp(-Sym[-2 1; 1 -1] * s)))
end

@testset "remaining cases (SymPy review, round 5)" begin
    @syms Ω ζ ω ω1 r σ ϵ y s
    LA = SymJ.LinearAlgebra
    # a sum of squares in real parameters vanishes on its real components (#24518)
    A = SymJ.Assumptions(Ω ∈ ℝ, ζ ∈ ℝ, ω ∈ ℝ, ω1 ∈ ℝ)
    zs = SymJ._nonneg_zero_set(ζ^2 * ω^2 + Ω^2 - 2Ω * ω1 + ω1^2, [Ω, ζ, ω, ω1], A)
    @test all(d -> haskey(d, ω1) ? d[ω1] == Ω : d[Ω] == ω1, zs)
    @test any(d -> get(d, ζ, nothing) == 0, zs)
    # the DFT matrices' eigenvalues exactly (#24082)
    @test LA.eigvals(Sym[exp(2PI * IM / 5)^(i * j) for i in 0:4, j in 0:4]) == [-sqrt(Sym(5)), -IM * sqrt(Sym(5)), IM * sqrt(Sym(5)), sqrt(Sym(5)), sqrt(Sym(5))]
    @test LA.eigvals(Sym[exp(2PI * IM / 3)^(i * j) for i in 0:2, j in 0:2]) == [-sqrt(Sym(3)), IM * sqrt(Sym(3)), sqrt(Sym(3))]
    # nested fractions over the lcm of their denominators (#24537)
    cs = [Sym(Symbol(:c, i)) for i in 1:6]
    inner = sum(1 / (1 + c * s) for c in cs)
    e = inner^3 / (1 + inner)
    @test SymJ._expansion_estimate(SymJ.numer_denom(e)[1]) > 1e5
    ce = cancel(e)
    vals = [s => Sym(3) // 7; [c => Sym(i) // 5 for (i, c) in enumerate(cs)]]
    @test N(subs(ce, vals...)) ≈ N(subs(e, vals...))
    # a float inverse transform by numeric partial fractions: fast, the same values (#25352)
    @syms tt
    F = (1 / s) * (2.25214078181399s + 1.0) / (0.0935725677320455s^3 + 1.00201771101894s^2 + 2.27274666658392s + 1.0)
    tl = @elapsed f = SymJ.inverse_laplace_transform(F, s, tt)
    @test tl < 20
    @test N(subs(f, tt => 1)) ≈ 1.057626488214003
    # floats with parameters: the floats as constants, no unsupported cases (#24193)
    sol = solve([y ~ σ * 2^(1 / 6), ϵ * ((y / r)^12 - 2 * (y / r)^6) ~ s * 4 * ϵ * (-σ^6 / r^6 + σ^12 / r^12)], [y, s])
    @test !any(c -> c.second isa Exception, sol.cases)
    g = only(last(sol.cases).second)
    @test N(subs(g[s], r => 2, σ => 1)) ≈ 1
end

@testset "wrong answers (SymPy review, round 6)" begin
    @syms a b u v q i t
    LA = SymJ.LinearAlgebra
    # radicals in a linear system over the field they generate: √6√11 = √2√33 (#22852)
    x0, x1, x2 = Sym(:x0), Sym(:x1), Sym(:x2)
    sol = solve([sqrt(Sym(6)) * sqrt(Sym(11)) * x0 - sqrt(Sym(2)) * sqrt(Sym(33)) * x0 + x1 ~ 1, x1 - x2 ~ 0], [x0, x1, x2])
    @test !any(d -> any(w -> occursin("NaN", string(N(subs(w, x0 => 1, x1 => 1, x2 => 1)))), values(d)), sol)
    # exponentials of rational multiples of q seen together in rank and nullspace (#23531)
    m = Sym[-2cosh(q / 3) exp(-q) 1; exp(q) -2cosh(q / 3) 1; 1 1 -2cosh(q / 3)]
    @test LA.rank(m) == 2
    @test size(LA.nullspace(m), 2) == 1
    # a degenerate factor with complex roots is a case (#23459); derivative unknowns
    # keep their names in Cases
    cs = solve([a * u + b * v ~ 0, b^2 * u + a^2 * v ~ 0], [u, v])
    @test any(c -> any(r -> occursin("sqrt(3)", string(r)), c.first), cs.cases)
    f, g = SymFunction(:f), SymFunction(:g)
    cs = solve([a * diff(f(t), t) + b * diff(g(t), t) ~ 1, diff(f(t), t) - diff(g(t), t) ~ 2], [diff(f(t), t), diff(g(t), t)])
    @test haskey(only(last(cs.cases).second), diff(f(t), t))
    # telescoping infinite sums (#23119)
    @test SymJ.summation(1 / (i + 1) - 1 / (i + 2), i, 0, Inf) == 1
    @test SymJ.summation(1 / (5 * (5i + 2)) - 1 / (5 * (5i + 7)), i, 0, Inf) == Sym(1) / 10
    # log(0) = -∞, and 0*log(0) is NaN, not 0 (#22612)
    @test log(Sym(0)) == Sym(-Inf)
    @test isnan(N(subs(x * log(x), x => 0)))
    # an integer-valued float is an integer (#23731)
    @test solve(x ~ 1.0, x; assume = (x ∈ ℤ,)) == [Sym(1.0)]
    # binomial at negative integers, as binomial(n, n) = 1 (#23497)
    @test binomial(Sym(-1), Sym(-1)) == 1
    @test binomial(Sym(-3), Sym(-5)) == 6
    # floats with symbols and a quartic: refused, not formulas with no accuracy (#22685)
    @test_throws ArgumentError LA.eigvals(Sym[1.5 x 0 0; x 0.5 x 0; 0 x 2.5 x; 0 0 x 3.5])
    # an unknown function's argument is not floated (#22524)
    h = SymFunction(:h)
    @test solve(h(2) * 4.0 + x ~ 0, x) == [-4.0 * h(2)]
end

@testset "crashes and hangs (SymPy review, round 6)" begin
    @syms a b p q m kx ky
    LA = SymJ.LinearAlgebra
    # divrem of polynomials with coefficients of two types (#23828)
    @test divrem(Poly(a * x + b + p, a), Poly(a, a)) == (Poly(x, a), Poly(b + p, a))
    # an unknown function replaced by a function in subs (#23759)
    f, g = SymFunction(:f), SymFunction(:g)
    @test subs(f(g(x)), f => sin) == sin(g(x))
    @test subs(diff(f(x), x) + f(x), f => sin) == cos(x) + sin(x)
    # simplify of a conjunction of facts (#23391)
    @test simplify((x ≪ 1) & (x ≪ 2)) == (x ≪ 1)
    # subs into a piecewise: a branch whose condition is false is not evaluated (#23295)
    @test subs(SymJ.piecewise(abs(b) ≫ 0 => 1 / b, true => 0), b => 0) == 0
    # a zero pivot in nested radicals (#22604)
    J = SymJ.jordan_form(Sym[0 0 1 0; 0 0 0 1; -2 1 0 0; 1 -1 0 0])[2]
    @test LA.isdiag(J)
    # logs with a common coefficient: no degree-72 polynomial (#22462)
    @test length(solve(72log(x - Sym(1) // 288) / 5 - 72log(x) / 5 ~ q, x)) == 1
    # det with im and a symbolic denominator: fraction-free (#23408)
    I2, X, Y, Z = Sym[1 0; 0 1], Sym[0 1; 1 0], Sym[0 -IM; IM 0], Sym[1 0; 0 -1]
    H = ((kx^2 + ky^2) / (2m) - p) * kron(Z, I2) + a * kron(X, I2) + b * (kx * kron(I2, Y) - ky * kron(I2, X)) * kron(Z, I2) + q * kron(I2, X)
    vals = [kx => 1, ky => 2, m => 3, p => 5, a => 7, b => 11, q => 13]
    @test N(subs(LA.det(H), vals...)) ≈ real(LA.det(ComplexF64.(N.(subs.(H, Ref.(vals)...)))))
    # a quartic in 10 parameters: a clear error, not hours (#23604)
    ps = [Sym(Symbol(:c, i)) for i in 1:10]
    @syms u v
    sys = [u^2 * ps[1] + v^2 * ps[2] + u * v * ps[3] + ps[4] * u + ps[5] ~ 0, u^2 * ps[6] + v^2 * ps[7] + ps[8] * u * v + ps[9] * v + ps[10] ~ 0]
    @test_throws ArgumentError SymJ._parametric_degree_guard(Sym[e.lhs for e in sys], [u, v])
end

@testset "missing cases (SymPy review, round 6)" begin
    @syms x y a b c d n p v z R i
    PI = SymJ.PI
    LA = SymJ.LinearAlgebra
    # cse: terms that two sums share (#23297)
    reps, red = SymJ.cse([-a + d / 9 + c^2 / 2 + c / 3, -b + d / 9 + c^2 / 2 - c / 3])
    @test length(reps) == 1
    @test reps[1].second == c^2 / 2 + d / 9
    @test red == [-a + c / 3 + reps[1].first, -b - c / 3 + reps[1].first]
    # diff of max and min (#22726)
    @test diff(max(x, y), x) == SymJ.heaviside(x - y)
    @test diff(min(x, 2), x) == SymJ.heaviside(2 - x)
    # piecewise conditions with | and & (#22852)
    @test SymJ.piecewise(((x ≪ 0) | (x ≫ 1)) => 1, true => 0) == SymJ.piecewise(x ≪ 0 => 1, x ≫ 1 => 1, true => 0)
    @test SymJ.piecewise(((x ≪ 0) & (x ≫ -1)) => 1, true => 0) == SymJ.piecewise((x ≪ 0, x ≫ -1) => 1, true => 0)
    # powers with negated exponents combined (#22546)
    @test simplify((a + b)^(n + 1) / (a + 2b)^(n + 1); assume = (a ≫ 0, b ≫ 0)) == ((a + b) / (a + 2b))^(n + 1)
    # log of a sum over its denominator, x real (#23569)
    @test simplify(log(1 + 1 / x^2) / 2 - (log(x^2 + 1) - log(x^2)) / 2; assume = (x ∈ ℝ,)) == 0
    # a long integer counts in simplify's size (#23210)
    @test simplify(30log(Sym(5))) == 30log(Sym(5))
    # a conjugate radical in the denominator (#23641)
    @test simplify(diff(log(b + 2sqrt(c) * sqrt(a + b * x + c * x^2) + 2c * x), x)) == sqrt(c) / sqrt(a + b * x + c * x^2)
    # a trigonometric constant that is 0 (#23015)
    @test simplify(-32cos(PI / 14) * cos(2PI / 7) - 32sin(3PI / 14) * cos(3PI / 14) +
                   32cos(3PI / 14) * cos(2PI / 7) + 32sin(3PI / 14) * cos(PI / 14)) == 0
    # float exponents (#23284, #22584)
    r = solve(1 - 0.8 ~ exp(-6.7784905546796565e-6 * x^1.1), x)
    @test length(r) == 1
    @test N(r[1]) ≈ 77064.73080118217 rtol = 1e-10
    # every x except the excluded point (#23254)
    @test_throws "except x = 1" solve(x / (x - 1) - 1 / (x - 1) - 1 ~ 0, x)
    # the log first, then Lambert W (#23253)
    r = solve(4 / log(3 / x + 1)^2 - x^2 ~ 0, x)
    @test !isempty(r)
    @test all(s -> abs(N(subs(4 / log(3 / x + 1)^2 - x^2, x => s))) < 1e-8, r)
    # cos² + sin² in the denominators of a linear system (#23140)
    q2, q3 = Sym(:q2), Sym(:q3)
    @syms u1 u2 u3 d1 d2 d3
    s = solve([-(u1 * sin(q3) + u2 * cos(q3)) / sin(q2) + d1 ~ 0, -u1 * cos(q3) + u2 * sin(q3) + d2 ~ 0,
               (u1 * sin(q3) + u2 * cos(q3)) * cos(q2) / sin(q2) - u3 + d3 ~ 0], [u1, u2, u3])
    @test s[1][u1] == cos(q3) * d2 + d1 * sin(q2) * sin(q3)
    # log(A) = log(B) as A = B; polynomial equations first; floats with a log (#23662, #22723, #23318)
    @test solve([-y^2 + log(y^2 / x) + 2 ~ 0, -2x * y + 2x / y ~ 0], [x, y]) == [Dict(x => E, y => 1), Dict(x => E, y => -1)]
    s = solve([53.5780461486929 ~ x * log(y / (5.0 - y) + 1) / y, x ~ 0.0015 * z, 0.0015 ~ 7845.32 * y / z], [x, y, z])
    @test length(s) == 1
    @test N(s[1][y]) ≈ 0.03403016806816325 rtol = 1e-8
    s = solve([log(x * y) ~ log(Sym(2)), x + y ~ 3], [x, y])
    @test length(s) == 2
    # roots through x^(1/6) simplified and once (#23028)
    @test solve(2x^(1 // 3) - 1 - sqrt(Sym(3)) * sqrt(x) ~ 0, x) == [-1]
    # a small polynomial system with large radicals, quickly (#23165)
    @syms vc vs xn
    s = solve([vc / vs - xn / (170 - xn) ~ 0, 1 // 2 * 30000 * vc * 120 * xn - 102500vs - 240(200 - xn) ~ 0,
               19760400 - 102500vs * (170 - xn / 3) - 240(200 - xn) * (100 + xn / 6) ~ 0], [vc, vs, xn])
    @test sort([real(N(t[xn])) for t in s]) ≈ [-200.0201888275084, -200.0201888275084, -3.9343295061690924, 3.8038738278525512]
    # the case a² + b² = 0 of a circle and a line (#23889)
    s = solve([x^2 + y^2 ~ R^2, a * x + b * y ~ c], [x, y])
    @test !any(cv -> cv.second isa Exception, s)
    # special solutions within a family dropped (#23110)
    @test length(solve([x * y ~ 0, x * z ~ 0, x * v ~ 0], [x, y, z, v])) == 2
    # unknowns with an equation of their own solved first (#23322)
    s = solve([d + p - x ~ 0, d^2 + 1 - 9 ~ 0, -1 - p^2 + 4 ~ 0], [d, p, x])
    @test Set(t[x] for t in s) == Set([2sqrt(Sym(2)) + sqrt(Sym(3)), 2sqrt(Sym(2)) - sqrt(Sym(3)),
                                        -2sqrt(Sym(2)) + sqrt(Sym(3)), -2sqrt(Sym(2)) - sqrt(Sym(3))])
    # cosines summed over a full period, through exponentials (#23576)
    @test SymJ.summation(cos(2PI * i / v), i, 1, v; assume = (v ∈ ℤ, v ≫ 1)) == 0
    # minimal polynomial of Cardano's roots in the casus irreducibilis (#23677)
    rs = solve(4000000x^3 - 239960000x^2 + 4782399900x - 31663998001 ~ 0, x)
    m = SymJ.minimal_polynomial(rs[1] + rs[2], y)
    @test abs(N(subs(m, y => N(rs[1] + rs[2])))) < 1e-3
    @test SymJ.degree(Poly(m, y)) == 3
    # a series solution with float coefficients (#22525)
    ρ = SymFunction(:ρ)
    @syms N0 s0
    @test SymJ.dsolve(4.0 * N0 * sqrt(N0 - 1) * ρ(s0) + (4s0^2 * (N0 - 1) + (N0 - 2s0 * (N0 - 1))^2) * diff(ρ(s0), s0, 2) ~ 0,
                      ρ(s0); series = 5) isa SymJ.Relation
end

@testset "wrong answers (SymPy review, round 7)" begin
    @syms x y z k a n
    f = SymFunction(:f)
    # sums with an infinite part beside others: decided as a whole (#21434, #21315, #21151)
    @test SymJ.summation(1 - cos(1 / k), k, 1, Inf) == SymJ._summation_unevaluated(1 - cos(1 / k), k, Sym(1), Sym(Inf))
    @test SymJ.summation(1 - 1 / k, k, 1, Inf) == Inf
    @test SymJ.summation(1 / k - log(1 + 1 / k), k, 1, Inf) == -SymJ.SpecialFunctionsSupport.SF.digamma(Sym(1))
    @test N(SymJ.summation(3 / (4k) - log(1 + 3 / (4k)), k, 1, Inf)) ≈ 3 * Base.MathConstants.eulergamma / 4 + log(SymJ.SpecialFunctionsSupport.SF.gamma(1.75)) rtol = 1e-12
    # a quotient whose denominator is 0 by the identities is defined nowhere (#21340)
    e = (tan(x)^2 - sec(x)^2 + 1) / (sin(x)^2 + cos(x)^2 - 1)
    @test simplify(e) == e
    # a piece of a piecewise that holds identically: its interval (#21079)
    S = solve(SymJ.piecewise(x ≪ -1 => 0, x ≦ 1 => x^2, true => log(x)) ~ 0, x)
    @test S isa SymJ.SymSet
    @test S(-3)
    @test S(0)
    @test !S(Sym(1) // 2)
    # x^n = 1 for an integer n: refused, not the one root 1 (#21014)
    @test_throws ArgumentError solve(x^n ~ 1, x; assume = (n ∈ ℤ, n ≫ 0))
    # |u| branches that are identities: the half-lines (#22259)
    @test solve(abs(x) ~ -x, x; assume = (x ∈ ℝ,)) == solve(x ≦ 0, x)
    @test solve(abs(x) / x ~ 1, x; assume = (x ∈ ℝ,)) == solve(x ≫ 0, x)
    # f(x) => y keeps f'(x) (#21658)
    @test subs(diff(f(x), x)^2 + f(x), f(x) => y) == y + diff(f(x), x)^2
    # a BigFloat beyond Float64's range: an error, not Inf (#21860)
    @test_throws ArgumentError Sym(big"1.5e400")
    # |x| y = 0: the case x = 0 (#21943)
    @test solve(abs(x) * y ~ 0, y) isa SymJ.Cases
    # numbers as subs keys leave exponents alone (#21373)
    @test subs(x * y, 1 => 2) == x * y
    @test subs(x / z, 1 => 1.0) == x / z
end

@testset "crashes and hangs (SymPy review, round 7)" begin
    @syms x t ka k21 ke k12
    # rational powers with huge denominators: no 2^k, no loop over k (#21762, #21860, #21647)
    @test Sym(3)^(Sym(16666666667) // 50000000000) isa Sym
    @test subs((x^2 + 6)^(Sym(33333333333333333) // 50000000000000000), x => Sym(1) // 2) isa Sym
    @test subs(log((Sym(567) // 500)^(811 * (Sym(567) // 500)^x / 100)), x => 6) isa Sym
    # a 3×3 linear system in 11 parameters: its cases quickly (#21694)
    @syms A B C n m nu l R h K D E F
    PI = SymJ.PI
    eq1 = -PI^2 * A * m^2 / l^2 + A * n^2 * nu / (2R^2) - A * n^2 / (2R^2) + PI * B * m * n * nu / (2R * l) + PI * B * m * n / (2R * l) + PI * C * m * nu / (R * l) + D / K
    eq2 = PI * A * m * n * nu / (2R * l) + PI * A * m * n / (2R * l) + PI^2 * B * m^2 * nu / (2l^2) - PI^2 * B * m^2 / (2l^2) - B * n^2 / R^2 - PI^2 * C * h^2 * m^2 * n / (12R^2 * l^2) - C * n / R^2 - C * h^2 * n^3 / (12R^4) + E / K
    eq3 = -PI * A * m * nu / (R * l) + PI^2 * B * h^2 * m^2 * n / (12R^2 * l^2) + B * n / R^2 + B * h^2 * n^3 / (12R^4) + PI^4 * C * h^2 * m^4 / (12l^4) + PI^2 * C * h^2 * m^2 * n^2 / (6R^2 * l^2) + C / R^2 + C * h^2 * n^4 / (12R^4) - F / K
    @test solve([eq1 ~ 0, eq2 ~ 0, eq3 ~ 0], [A, B, C]) isa SymJ.Cases
    # ODE systems: eigenvalues of a cubic in parameters refused; Cardano's complex pair
    # recognized, not more columns than unknowns (#21867)
    Ac, Ap, Ad = SymFunction(:Ac), SymFunction(:Ap), SymFunction(:Ad)
    @test_throws ArgumentError SymJ.dsolve([diff(Ac(t), t) ~ Ad(t) * ka + Ap(t) * k21 - Ac(t) * (ke + k12),
        diff(Ap(t), t) ~ Ac(t) * k12 - Ap(t) * k21, diff(Ad(t), t) ~ -Ac(t) * ka], [Ac(t), Ap(t), Ad(t)])
    r = SymJ.dsolve([diff(Ac(t), t) ~ 2Ad(t) + 3Ap(t) - 12Ac(t), diff(Ap(t), t) ~ 7Ac(t) - 3Ap(t), diff(Ad(t), t) ~ -2Ac(t)],
                    [Ac(t), Ap(t), Ad(t)])
    @test length(r) == 3
    @test !any(e -> has(e.rhs, IM), r)
end

@testset "missing cases (SymPy review, round 7)" begin
    @syms h p A
    as = (h ≫ 0, p ≫ 0, A ∈ ℝ, sin(A) ≫ 0)
    # h cos A + h + p sin A = h(cos A + 1) + p sin A > 0 (#21811)
    @test refine(abs(h * cos(A) + h + p * sin(A)); assume = as) == h * cos(A) + h + p * sin(A)
    q = p + p^2 * sin(A) / (2abs(h * cos(A) + h + p * sin(A))) - p^2 * sin(A) / (2(h * cos(A) + h + p * sin(A)))
    @test simplify(q; assume = as) == p
    @test SymJ._hasfn(refine(abs(h * cos(A) - h + p * sin(A)); assume = as), :abs)   # undecided: kept
    # complex floats in a linear system, solved in complex floating point (#21445)
    @syms IR y x
    r = solve([IR + y ~ 1.0im, IR - y ~ 2.0], [IR, y])
    @test isapprox(ComplexF64(N(only(r)[IR])), 1.0 + 0.5im)
    @test isapprox(ComplexF64(N(only(solve([IR + 1.5 + 2.0im ~ 0], [IR]))[IR])), -1.5 - 2.0im)
    # a polynomial in exponentials e^(c x), c rational, factored in e^(x/q) (#21511)
    e = y^2 - y * exp(x / 7) - y * exp(x / 2) + exp(9x / 14)
    @test factor(e).kind == SymJ.MUL
    @test iszero(expand(factor(e) - e))
    @test factor(exp(2x) - 1) == (exp(x) - 1) * (exp(x) + 1)
    # the minimal polynomial of cos(2π/11): its roots cos(2πk/11) (#21631)
    PI = SymJ.PI
    @test Set(solve(32x^5 + 16x^4 - 32x^3 - 12x^2 + 6x + 1 ~ 0, x)) == Set(cos(2k * PI / 11) for k in 1:5)
    C11 = Sym[(abs(i - j) == 1 || abs(i - j) == 10) ? Sym(1) // 2 : Sym(0) for i in 1:11, j in 1:11]
    @test count(==(cos(2PI / 11)), SymJ.LinearAlgebra.eigvals(C11)) == 2
    # Laplace transforms: f/t, symbolic powers of t, and their inverses (#21372)
    @syms t s a n
    L, IL = SymJ.laplace_transform, SymJ.inverse_laplace_transform
    @test L(sin(a * t) / t, t, s) == atan(a / s)
    @test L(exp(-2t) * sin(3t) / t, t, s) == atan(3 / (s + 2))
    @test IL(atan(a / s), s, t) == sin(a * t) / t
    @test IL(factorial(n) / (s + a)^(n + 1), s, t) == t^n * exp(-a * t)
    @test IL(L(t^n * exp(-a * t), t, s), s, t) == t^n * exp(-a * t)
end

@testset "crashes and hangs (SymPy review, round 8)" begin
    @syms x z L Ca0 t a
    LA = SymJ.LinearAlgebra
    # boundary-value problems: two value conditions, a derivative at a symbolic point (#20547)
    C = SymFunction(:C)
    @test SymJ.dsolve(diff(C(z), z, z) ~ 0, C(z), C(0) => 1, C(1) => 0) == (C(z) ~ 1 - z)
    r = SymJ.dsolve(diff(C(z), z, z) - C(z) ~ 0, C(z), C(0) => 1, C(1) => 0)
    @test iszero(simplify(subs(r.rhs, z => 1)))
    r = SymJ.dsolve(diff(C(z), z, z) - C(z) ~ 0, C(z), C(0) => Ca0, subs(diff(C(z), z), z => L) => 0)
    @test iszero(simplify(subs(diff(r.rhs, z), z => L)))
    # isposdef and Hermitian for symbolic matrices (#19547)
    @test LA.isposdef(Sym[2 0; 0 1])
    @test !LA.isposdef(Sym[1 2; 2 1])
    @test LA.isposdef(LA.Hermitian(Sym[2 1; 1 2]))
    # Poly: higher derivatives; conj, transpose (#19590, #19531)
    p = SymJ.Poly(x^2 + 2x + 1, x)
    @test diff(p, x, 2) == SymJ.Poly(Sym(2), x)
    @test transpose([p, p])[1] == p
    @test conj(SymJ.Poly(x^2 + SymJ.IM * x, x)) == SymJ.Poly(x^2 - SymJ.IM * x, x)
    # a combination of relations (#20080)
    @test solve((x ≧ 0) & (x ≦ 1), x) == solve(x ≧ 0, x) ∩ solve(x ≦ 1, x)
    @test isempty(solve((x ≫ 0) & (x ≪ 0), x))
    # a matrix function with eigenvalues in a quartic's radicals: refused, not a hang (#19565)
    I_, s2 = SymJ.IM, sqrt(Sym(2))
    M = Sym[1 -I_ I_*s2 0; -I_ Sym(1)//2 0 I_*s2; I_*s2 0 Sym(1)//2 -I_; 0 I_*s2 -I_ Sym(1)//2]
    @test_throws ArgumentError exp(-M * t)
    @test exp(t * Sym[0 1; -1 0]) == Sym[cos(t) sin(t); -sin(t) cos(t)]
end

@testset "wrong answers (SymPy review, round 8)" begin
    @syms phi x y k a
    IM, PI = SymJ.IM, SymJ.PI
    # exp(iφ) with bounds on φ is not a positive real (#20453)
    as = (phi ≫ 0, phi ≪ 1)
    @test refine(imag(exp(IM * phi)); assume = as) == sin(phi)
    @test refine(abs(exp(IM * phi)); assume = as) == 1
    @test !(refine(exp(IM * phi) ≫ 0; assume = as) isa Bool)
    @test SymJ._isreal(IM * phi, SymJ._assumptions(as)) !== true
    # inverse trig inequalities; acot jumps at 0 (#20755, #19977)
    @test solve(acot(x) ≪ 0, x) == solve(x ≪ 0, x)
    @test solve(atan(log(x)) ≫ 0, x) == solve(x ≫ 1, x)
    @test SymJ._signs(acot(x), SymJ._assumptions((x ∈ ℝ,))) === nothing
    # trig roots in every period of an assumed interval (#20798)
    @test solve(cos(2x) ~ 1 // 2, x; assume = @assume(0 ≤ x ≤ 2PI)) == [PI / 6, 5PI / 6, 7PI / 6, 11PI / 6]
    @test length(solve(sin(x) ~ 1 // 2, x; assume = @assume(0 ≤ x ≤ 4PI))) == 4
    @test solve(cos(2x) ~ 1 // 2, x) == [PI / 6, 5PI / 6]
    # a disguised zero: sign 0, floor 0 (log evaluated at the working precision; #20897, #20865)
    L(n) = log(Sym(n))
    z = 5L(5) / (-6L(6) - 5L(5) + 6L(2) + 6L(3)) + 1
    @test sign(z) == 0
    @test floor(z) == 0
    s3 = sqrt(Sym(3)) + 1 - 1 / (sqrt(Sym(3)) / 2 - Sym(1) // 2)
    @test ceil(s3) == 0
    @test sign(log(Sym(2)) - 1) == -1
    # Poly coefficients that are 0 are dropped (#20389, #20427)
    @test degree(Poly(y * (x * (x + 1) - x^2 - x) + 1, y)) == 0
    q = Sym(1) // 3
    r = sqrt(Sym(9165)); u = (24201 + 253r)^q; D = 217603955769048 * u + 2273005839412 * r * u
    c = -117968192370600 * Sym(18)^q / D - 15720318185 * Sym(2)^(2q) * Sym(3)^q * u^2 / D +
        15720318185 * Sym(12)^q * u^2 / D + 117968192370600 * Sym(2)^q * Sym(3)^(2q) / D
    @test degree(Poly(c * x + 1, x)) == 0
    @test simplify(c) == 0
    @test degree(Poly((sqrt(Sym(2)) - 1) * x + 1, x)) == 1
    # a negative base: real only at integers (#19816); exact integer logs
    @test solve((-2)^x ~ -8, x; assume = (x ∈ ℤ,)) == [3]
    @test solve((-2)^x ~ 8, x; assume = (x ∈ ℤ,)) == Sym[]
    @test solve(2^x ~ 8, x; assume = (x ∈ ℤ,)) == [3]
    @test solve(Sym(3)^x ~ 5, x) == [log(Sym(5)) / log(Sym(3))]
    # an equation free of the unknowns: conditions on the parameters (#19564)
    r = solve([(k + 1)^2 - 4(k + 4) ~ 0, x^2 * (k + 4) + x * (k + 1) + 1 ~ 0], [x])
    @test r isa SymJ.Cases
    @test any(cv -> cv.first == Any[k ~ 5] && cv.second == [Dict(x => -Sym(1) // 3)], r)
    @test any(cv -> cv.first == Any[k ~ -3] && cv.second == [Dict(x => Sym(1))], r)
    # k = a, a + 1, … up to a non-integer bound (#19745)
    @test SymJ.summation(x, x, 1, Sym(3) // 2) == 1
    @test SymJ.summation(x, x, Sym(1) // 2, Sym(7) // 2) == 8
    # a root at a parametric pole: a case (a = 0, where 0 is the pole)
    r = solve(x^2 / (a + x) - 2 ~ 0, x)
    @test r isa SymJ.Cases
    @test any(cv -> cv.first == Any[a ~ 0] && cv.second == [Sym(2)], r)
end

@testset "missing cases (SymPy review, round 8)" begin
    @syms x y z k n h m
    IM, PI = SymJ.IM, SymJ.PI
    u = SymFunction(:u)
    # |u| sign(u) = u; |x|^m x^(2k) for real x (#20402, #19484, #20108)
    @test simplify(abs(x) * sign(x)) == x
    @test simplify(x^2 / abs(x); assume = (x ∈ ℝ,)) == abs(x)
    @test simplify(x^2 / abs(x)) == x^2 / abs(x)
    # √(a + bi) with rational parts (#20024)
    @test simplify(sqrt(3 + 4IM)) == 2 + IM
    @test simplify(sqrt(-3 + 4IM)) == 1 + 2IM
    # factor in powers of an integer base; simplify over a root (#19776, #19769)
    @test factor(4^x - 2^x - 2) == (2^x - 2) * (2^x + 1)
    @test simplify((x - 1) / (sqrt(x) - 1)) == sqrt(x) + 1
    # solve: conj, b^e = 0, offsets in exponents, reciprocal bases, sinh + cosh (#20305,
    # #20205, #20199, #20709, #19489)
    @test solve(conj(x) - 3 + IM ~ 0, x) == [3 + IM]
    @test solve((x - 1)^x ~ 0, x) == [Sym(1)]
    @test length(solve(exp(2x) + 2exp(x + 1) ~ 3, x)) == 2
    @test Set(solve(sqrt(2 + sqrt(Sym(3)))^x + sqrt(2 - sqrt(Sym(3)))^x ~ 4, x)) == Set(Sym[-2, 2])
    @test solve(sinh(8x) + cosh(8x) ~ 2, x) == [log(Sym(2)) / 8]
    # an exponential system through logs (#19587)
    @test solve([32 * 2^m * 2^n - 4^n ~ 0, 27 * 3^m - 3^(-n) ~ 0], [m, n]) == [Dict(m => Sym(-4), n => Sym(1))]
    # inequalities in piecewise functions and abs (#19718)
    @test solve(SymJ.piecewise((x ≫ 10) => 1, (x ≦ 10) => 0) ≫ 0, x) == solve(x ≫ 10, x)
    @test string(solve(abs(x) - 1 ≪ 0, x)) == "(-1, 1)"
    # sums: a shifted start, the summand refined on the range, shifted limits (#20464,
    # #19687, #19685)
    @test SymJ.summation(2^k / factorial(k), k, 1, Inf) == exp(Sym(2)) - 1
    @test SymJ.summation(SymJ.piecewise((k ≦ n) => k * SymJ.binomial(n, k) / 2^n, true => 0), k, 0, n) == n / 2
    @test SymJ.summation(1 / abs(x + 9), x, 7, Inf) == Inf
    @test simplify(SymJ.summation(u(k), k, 0, n) - SymJ.summation(u(k), k, 0, n - 1)) == u(n)
    # subs of a key with a numeric coefficient (#19422)
    @test subs(2 * x * y, 2x => z) == y * z
    @test subs(2x + 4, 2x => z) == z + 4
    @test subs(3 * x * y, 2x => z) == 3 * x * y
    # minimal polynomial of e^(iπ/3) (#19744)
    @test SymJ.minimal_polynomial(exp(IM * PI / 3), x) == x^2 - x + 1
    # eigenvectors of an exact matrix with numeric eigenvalues (#20351)
    M = Sym[0 1 1 1 1 0 1; 1 0 0 0 0 0 0; 1 1 0 0 0 0 0; 0 1 1 0 1 0 0; 1 0 1 1 0 1 0; 1 0 0 0 1 0 0; 0 0 0 0 1 0 0]
    T = transpose(M) * SymJ.LinearAlgebra.diagm(Sym[1 // 5, 1, 1 // 2, 1 // 3, 1 // 4, 1 // 2, 1])
    @test size(SymJ.LinearAlgebra.eigvecs(T)) == (7, 7)
    @test any(==(Sym(1)), SymJ.LinearAlgebra.eigen(T).values)
end

@testset "crashes and hangs (SymPy review, round 9)" begin
    @syms x y z d n
    # subs and the simplify_* functions on a relation (#18412)
    @test subs(x ~ y^3 + 1, y^3 => z) == (x ~ z + 1)
    @test SymJ.simplify_powers(x^2 * x ~ 1) == (x^3 ~ 1)
    # subs into a series: x => x^2 about 0, its order doubled; a parameter (#19120)
    s = subs(series(cos(x), x, 0, 8), x => x^2)
    @test Sym(s) == 1 - x^4 / 2 + x^8 / 24 - x^12 / 720
    @test Sym(subs(series(exp(d * x), x, 0, 3), d => 2)) == 1 + 2x + 2x^2
    @test_throws ArgumentError subs(series(cos(x), x, 0, 8), x => x + 1)
    # refine of an Or of relations, and of a membership (#18389, #19145)
    @test refine((n ≫ -1) | (n ≪ -5); assume = (n ≫ 0,)) == true
    @test refine((n ~ 0) | (n ≧ 1)) isa SymJ.FactOr
    @test refine(1 / (1 + d) ∈ ℤ; assume = (d ≫ 0,)) == false
    @test refine(d ∈ ℤ; assume = (d ∈ ℤ,)) == true
    # a root index past Int: the power as it is (#18526)
    @test Sym(2)^(1 // big(10)^50) isa Sym
    @test isapprox(Float64(N(Sym(2)^(1 // big(10)^50))), 1.0)
    # an array of polynomials times a polynomial
    p = SymJ.Poly(x + 1, x)
    @test [p p] * SymJ.Poly(x - 1, x) == [SymJ.Poly(x^2 - 1, x) SymJ.Poly(x^2 - 1, x)]
end

@testset "wrong answers (SymPy review, round 9)" begin
    @syms n x y
    IM = SymJ.IM
    # an alternating series whose terms do not tend to 0 diverges (#17942)
    @test_throws ArgumentError SymJ.summation((-1)^n * (1 / n + 1), n, 1, Inf)
    @test_throws ArgumentError SymJ.summation((-1)^n * n / (n + 1), n, 1, Inf)
    @test SymJ.summation((-1)^n / n, n, 1, Inf) == log(Sym(1) // 2)
    # infinite sides compared as they are, not through ∞ - ∞ (#19245)
    @test simplify(Sym(Inf) ~ Sym(Inf)) == true
    @test simplify(Sym(Inf) ≦ Sym(Inf)) == true
    @test simplify(Sym(Inf) ≪ Sym(Inf)) == false
    @test simplify(-Sym(Inf) ≪ Sym(Inf)) == true
    # real, imag, conj of the unknowns: through real and imaginary parts (#19349)
    @test solve([real(x) ~ 1, imag(x) ~ 2], [x]) == [Dict(x => 1 + 2IM)]
    @test solve([x + conj(x) ~ 4, x - conj(x) ~ 2IM], [x]) == [Dict(x => 2 + IM)]
    r = solve([x * y ~ 220, real(x) ~ 2], [x, y])
    @test length(r) == 2
    yi = Sym(:y_im)
    @test all(s -> iszero(simplify(subs(s[x] * s[y] - 220, yi => 7))), r)
    @test all(s -> iszero(simplify(subs(real(s[x]), yi => 7) - 2)), r)
    # abs in a system: over ℂ (a family), or with the unknowns assumed real (±2)
    @syms I2 R2
    @test solve([I2 * R2 ~ 220, abs(I2) ~ 2], [I2, R2]; assume = (I2 ∈ ℝ, R2 ∈ ℝ)) == [Dict(I2 => Sym(2), R2 => Sym(110)), Dict(I2 => Sym(-2), R2 => Sym(-110))]
    r = solve([I2 * R2 ~ 220, abs(I2) ~ 2], [I2, R2])
    @test length(r) == 2
    @test any(s -> has(s[I2], IM), r)
    # assumptions on the unknowns filter the solutions
    @test solve([x^2 ~ 4, y ~ x], [x, y]; assume = (x ≫ 0,)) == [Dict(x => Sym(2), y => Sym(2))]
    @test length(solve([x^2 ~ -4, y ~ x], [x, y]; assume = (x ∈ ℝ,))) == 0
end

@testset "missing cases (SymPy review, round 9)" begin
    @syms x y t z n i b
    IM, PI = SymJ.IM, SymJ.PI
    # limits through factorial, binomial, lambertw by name (#19379, #18992, #18969)
    @test SymJ.summation(factorial(n) / factorial(n + 2), n, 1, Inf) == Sym(1) // 2
    @test limit(x / factorial(x)^(1 / x), x => Inf) == SymJ.E
    @test limit(lambertw(x) / log(x), x => Inf) == 1
    # cot inverted; heaviside's set refused clearly (#19115, #18939)
    @test solve(cot(x) ~ 0, x) == [PI / 2]
    @test_throws ArgumentError solve(SymJ.heaviside(x) ~ 1, x)
    @test solve(SymJ.heaviside(x) ~ Sym(1) // 2, x) == [Sym(0)]
    # roots of 3cos b = sin b tidied; multiple angles (#19102, #18315)
    @test solve(3cos(b) ~ sin(b), b) == [atan(Sym(3)), atan(Sym(3)) + PI]
    @test all(r -> abs(N(sin(r) + cos(r / 2))) < 1e-12, solve(sin(x) + cos(x / 2) ~ 0, x))
    @test length(solve(sin(x) + cos(x / 2) ~ 0, x)) == 4
    # exponential systems (#17941)
    @test solve([3exp(x) - 2exp(y) ~ 17, 5exp(x) + 3exp(y) ~ 41], [x, y]) == [Dict(x => log(Sym(7)), y => log(Sym(2)))]
    # real/imag of log of a complex constant; conj recombined (#17948, #18271)
    @test real(log(2 + 3IM)) == log(sqrt(Sym(13)))
    @test imag(log(1 + IM)) == PI / 4
    @test simplify(conj(z) - (real(z) - IM * imag(z))) == 0
    # abs merged; a piecewise as abs; floor of floor + 1 (#18840, #19167, #18689)
    @test simplify(abs(t * abs(x) / x)^2) == abs(t)^2
    @test simplify(SymJ.piecewise(x ≪ 0 => -x, x ≧ 0 => x)) == abs(x)
    @test floor(floor(x) + 1) == floor(x) + 1
    # factorial's poles; sin of an imaginary multiple under simplify (extras)
    @test_throws SymJ.DivisionByZero factorial(Sym(-3))
    @test simplify(sin(IM * t) / sinh(t)) == IM
    @test simplify(cos(IM * t) - cosh(t)) == 0
    # a sum with a point condition (#18388)
    @test SymJ.summation(SymJ.piecewise(i ⩵ 0 => 1, true => 0), i, 0, n) == 1
    @test SymJ.summation(SymJ.piecewise(i ⩵ 0 => 2, true => 1), i, 0, n) == n + 2
    # a rational radical constant; a second simplify pass (#18475, #18231)
    @test simplify(-48 * sqrt(4 + 2sqrt(Sym(5))) / (1 + sqrt(Sym(5)))^(Sym(3) // 2)) == -24
    @test simplify((cos(x - PI * cos(x)) - cos(x + PI * cos(x))) * cos(x)) == sin(2x) * sin(PI * cos(x))
end

@testset "wrong answers (SymPy review, round 10)" begin
    @syms x y z n a c0 c1 c2 c3 c4
    IM, PI = SymJ.IM, SymJ.PI
    # tiny constants are not 0; hidden zeros still are (#17450, #22604)
    @test simplify(exp(-Sym(300))) == 1 / exp(Sym(300))
    @test simplify((sqrt(Sym(2)) - 1)^300) == (sqrt(Sym(2)) - 1)^300
    @test simplify(sqrt((3 + sqrt(Sym(5))) / 2) - (1 + sqrt(Sym(5))) / 2) == 0
    @test inv(Sym[sqrt(Sym(2)) exp(-Sym(400)); 0 1])[1, 2] == -sqrt(Sym(2)) / (2exp(Sym(400)))
    @test N(imag((1 - exp(-Sym(40)))^IM)) ≈ sin(log1p(-exp(-40.0))) rtol = 1e-10   # ≈ -4.2e-18, not 0
    # the generic quartic at a degenerate point: q = 0 is a case (#17563)
    r = solve(c4 * x^4 + c3 * x^3 + c2 * x^2 + c1 * x + c0 ~ 0, x)
    s = r(c4 => 1 // 2, c3 => 0, c2 => 1, c1 => 0, c0 => -1)
    @test maximum(abs(N(v^4 / 2 + v^2 - 1)) for v in s) < 1e-12
    @test length(unique(N.(s))) == 4
    # series at a symbolic infinity (#16777)
    @test series(1 / (n + 1), n, Sym(Inf), 3) == series(1 / (n + 1), n, Inf, 3)
    @test series(exp(1 / n), n, limit(x, x => Inf), 3) == series(exp(1 / n), n, Inf, 3)
    # float matrices: a pivot tolerance (#17218)
    @test SymJ.LinearAlgebra.rank(Sym[0.8 -0.8 -0.4; -0.3 0.9 -0.4; -0.5 -0.1 0.8]) == 2
    @test size(SymJ.LinearAlgebra.nullspace(Sym[0.8 -0.8 -0.4; -0.3 0.9 -0.4; -0.5 -0.1 0.8]), 2) == 1
    # |x|/x is undefined at 0: no root there (#17438)
    @test solve(x^2 - 1 + abs(x) / x ~ 0, x; assume = (x ∈ ℝ,)) == [-sqrt(Sym(2))]
    @test solve(x + abs(x) / x ~ 1, x; assume = (x ∈ ℝ,)) == Sym[]
    # trig of a nonlinear argument on an assumed interval: every root, or an error (#17579)
    @test_throws ArgumentError solve(sin(log(x)) ~ 0, x; assume = (x ≫ 0, x ≦ 1))
    @test length(solve(cos(x^2) ~ 0, x; assume = (x ≧ 0, x ≦ 3))) == 3
    @test length(solve(cos(x^2) ~ 0, x; assume = (x ≧ -3, x ≦ 3))) == 6
    # a residual vanishing for some parameters: a case, not a dropped root
    r = solve(sqrt(z - x) + sqrt(y - x) ~ 0, z)
    @test r(y => x) == [x]
    @test r(y => 2, x => 1) == Sym[]
    # float poles raise as exact ones; N of an infinity is real when it can be (#16847)
    @test_throws SymJ.DivisionByZero SymJ.SpecialFunctionsSupport.SF.gamma(Sym(-1.0))
    @test N(SymJ.apply_function(:expinti, Sym[Sym(0)])) == -Inf
end

@testset "crashes and hangs (SymPy review, round 10)" begin
    @syms x y z a b Δ ω G
    LA = SymJ.LinearAlgebra
    # inequalities against ±∞ (#16582)
    @test isempty(solve(x ≪ -Inf, x).intervals)
    @test solve(x ≪ Inf, x) == solve(x ≫ -Inf, x)
    @test length(solve(1 / x ≪ Inf, x).intervals) == 2
    # subs on And, Or and memberships (#17530, #16627)
    @test string(subs((x ≪ 0) | (y ≪ 0), y => z)) == "x ≪ 0 | z ≪ 0"
    @test subs((x ∈ ℤ) & (y ≫ 0), x => 3) == (y ≫ 0)
    @test subs((x ∈ ℤ) & (y ≫ 0), x => Sym(1) // 2) === false
    # symbolic eigenvectors through nested radicals: the adjugate, quickly (#17732)
    M = Sym[-Δ G 0 G; G ω G 0; 0 -G Δ -G; -G 0 -G -ω]
    V = LA.eigvecs(M)
    @test size(V) == (4, 4)
    vals = (Δ => Sym(3) // 10, ω => Sym(7) // 10, G => Sym(1) // 5)
    ev = LA.eigvals(M)
    @test maximum(abs(N(subs(e, vals...))) for e in (M * V[:, 1] - ev[1] * V[:, 1])) < 1e-12
    # exponentials of a + b with exp(a), exp(b): one generator each (#17130)
    @test LA.rank(Sym[exp(a) exp(a + b); 1 exp(b)]) == 1
    @test SymJ.jordan_form(Sym[exp(a)+exp(b) exp(a)-exp(b); exp(a)-exp(b) exp(a)+exp(b)]) isa Tuple
    # a Poly in a relation; irrational constants; numeric values without a method
    @test solve(SymJ.Poly(x^2 - 1, x) ~ 0, x) == [-1, 1]
    @test Sym(Base.MathConstants.golden) == (1 + sqrt(Sym(5))) / 2
    @test N(Sym(Base.MathConstants.catalan)) ≈ Float64(Base.MathConstants.catalan) rtol = 1e-14
    @test N(Sym(Base.MathConstants.eulergamma)) ≈ Float64(Base.MathConstants.eulergamma) rtol = 1e-14
    @test_throws ArgumentError N(SymJ.apply_function(:polygamma, Sym[IM, PI]))
    # a power of a polynomial = 0: its base, not its expansion (#17867)
    @test solve((x + 2)^(10^6) ~ 0, x) == [-2]
    @test solve((x - 1)^50 * (x^2 + 1)^30 ~ 0, x) == [1, -IM, IM]
end

@testset "missing cases (SymPy review, round 10)" begin
    @syms x y z t a b n m L i mu σ
    IM, PI = SymJ.IM, SymJ.PI
    LA = SymJ.LinearAlgebra
    # poles of inverse functions (#17143)
    @test_throws SymJ.DivisionByZero acsc(Sym(0))
    @test_throws SymJ.DivisionByZero acoth(Sym(1))
    @test_throws SymJ.DivisionByZero atan(IM)
    # csc and sec inverted (#17536)
    @test solve(csc(2PI / x) ~ 0, x) == Sym[]
    @test solve(sec(x) ~ 2, x) == [PI / 3, 5PI / 3]
    # a real x is finite (#17113)
    @test refine(x ≪ Inf; assume = (x ∈ ℝ,)) === true
    @test refine(x ≫ Inf; assume = (x ∈ ℝ,)) === false
    # float matrices: singular values numerically (#17332)
    @test N(LA.svdvals(Sym[1.0 2.0; 3.0 4.0])[1]) ≈ 5.464985704219043
    # a scaled sum factor matches (#17225)
    @test subs(7 * (x + y) / z, (x + y) / z => t) == 7t
    # |a|^p |b|^q merged and cancelled (#17292)
    @test simplify(abs(x^2 - 1) / abs(x - 1)) == abs(x + 1)
    # positive factors out of radicals (#17634)
    @test simplify(sqrt(a * b) / sqrt(b); assume = (b ≫ 0,)) == sqrt(a)
    @test simplify(sqrt(1 - x^2) / sqrt(1 - x); assume = (x ≪ 1, x ≫ -1)) == sqrt(x + 1)
    # exponentials on two sides of a difference (#17604)
    @test solve(-2^(3x / 11) * exp(x / 11) + PI^(x / 11) ~ 0, x) == [0]
    # inverse of a periodic function at a constant (#17300, #17208)
    @test atan(tan(Sym(5))) == 5 - 2PI
    @test asin(sin(Sym(2))) == PI - 2
    @test log(exp(IM * PI / 3)) == IM * PI / 3
    # max and min decided (#16926)
    @test max(x, x + 1) == x + 1
    @test refine(max(n, log(n) / log(Sym(2))); assume = (n ∈ ℤ, n ≫ 0)) == n
    # a symbolic point condition in a sum (#16532)
    f = SymFunction(:f)
    @test SymJ.summation(SymJ.piecewise(Ne(m, n) => 0, true => L) * f(n), n, 0, Inf; assume = (m ∈ ℤ, m ≧ 0)) == L * f(m)
    # f(x) f(conj x) is real (#16500)
    @test simplify(imag(cot(x) * cot(conj(x)))) == 0
    # sin/cos as tan (#16499)
    @test SymJ.rewrite(sin(x) / cos(x), tan) == tan(x)
    # iⁿ as an exponential (#16797)
    @test simplify(IM^n - exp(IM * PI * n / 2)) == 0
    # inverse hyperbolic constants through logs (#16730)
    @test simplify(log(Sym(3) / 2) - 2atanh(Sym(2)) + 2atanh(Sym(3))) == 0
    # a common factor free of the unknown function (#17397)
    u = SymFunction(:u)
    @test SymJ.dsolve(x * (diff(u(x), x, 2) + u(x)) ~ 0, u(x)) == SymJ.dsolve(diff(u(x), x, 2) + u(x) ~ 0, u(x))
    # an initial value problem: the candidate valid near the point (#16880)
    @test SymJ.dsolve(diff(f(x), x) ~ sqrt(f(x)), f(x), f(0) => 1) == (f(x) ~ (x / 2 + 1)^2)
    # resultants of float polynomials in several variables (#16646)
    @syms cc
    @test N(subs(SymJ.resultant(0.5y + cc^2 - 1.0, m - cc + 3, cc), y => 1, m => 2)) ≈ 24.5
    # the content of sum factors (#17233)
    @test simplify(diff(exp(-(x - mu)^2 / (2σ^2)) / (sqrt(2PI) * σ^2), mu)) ==
          -sqrt(Sym(2)) * (mu - x) / (2sqrt(PI) * σ^4 * exp((x - mu)^2 / (2σ^2)))
    # an equality condition in its piece's value (#16690)
    @test simplify(SymJ.piecewise(y ⩵ 0 => x * y, true => x)) == SymJ.piecewise(y ⩵ 0 => 0, true => x)
end

@testset "wrong answers (SymPy review, round 11)" begin
    @syms x y c
    IM = SymJ.IM
    LA = SymJ.LinearAlgebra
    # a fact on real(c) or imag(c) does not make c real (#14978)
    @test SymJ._hasfn(refine(imag(c); assume = (real(c) ≫ 0,)), :imag)
    @test refine(real(c) ≫ 0; assume = (real(c) ≫ 0,)) === true
    @test refine(imag(c) ≫ 0; assume = (imag(c) ≫ 0,)) === true
    @test solve(c^2 ~ -1, c; assume = (real(c) ≫ -1,)) == [-IM, IM]
    @test refine(imag(x); assume = (x ≫ 0,)) == 0
    # huge coefficients: roots in BigFloat, not NaN (#15355)
    r = solve(x^5 - big(10)^400 * x + 1 ~ 0, x)
    @test all(v -> !isnan(abs(N(v))), r)
    @test length(solve(x^5 - big(10)^400 * x + 1 ~ 0, x; assume = (x ∈ ℝ,))) == 3
    # (u)^n = c through its base, exactly; numeric roots polished (#15296)
    @test solve((1 + x)^5 ~ 3 // 2, x; assume = (x ∈ ℝ,)) == [(Sym(3) / 2)^(1 // 5) - 1]
    @test maximum(abs(N((1 + v)^30 - 3 // 2)) for v in solve((1 + x)^30 ~ 3 // 2, x)) < 1e-12
    q = solve(x^5 - x + 1 ~ 0, x)
    @test count(v -> imag(N(v)) != 0 && any(w -> N(w) == conj(N(v)), q), q) == 4    # exact conjugate pairs
    # relations among cosines in a rank (#15585)
    @syms z
    M = Sym[1 cos(x/4+y/4) cos(x/4+z/4) cos(y/4+z/4); cos(x/4+y/4) 1 cos(y/4-z/4) cos(x/4-z/4);
            cos(x/4+z/4) cos(y/4-z/4) 1 cos(x/4-y/4); cos(y/4+z/4) cos(x/4-z/4) cos(x/4-y/4) 1]
    @test LA.rank(M) == 2
    @test size(LA.nullspace(M), 2) == 2
    # has sees the numbers a node stores (#16054)
    @test has(x - 1, Sym(-1))
    @test has(3x, Sym(3))
end

@testset "crashes and hangs (SymPy review, round 11)" begin
    @syms x n a
    LA = SymJ.LinearAlgebra
    # Σ xⁿ/n: no DivisionByZero from the retry at 0; -log(1 - x) for |x| < 1 (#14909)
    @test any(cv -> cv.second == -log(1 - x), SymJ.summation(x^n / n, n, 1, Inf))
    @test SymJ.summation(1 / (n * 2^n), n, 1, Inf) == log(Sym(2))
    # ℤ ⊆ ℝ (#16268)
    @test ℤ ⊆ ℝ
    @test !(ℂ ⊆ ℝ)
    # a triangular matrix: its diagonal (#15395)
    @syms A B C D
    @syms E F
    @test LA.eigvals(Sym[-(A+B) 0 0 0 0; A -C 0 0 0; B C -D 0 0; 0 0 D -E 0; 0 0 0 E -F]) == [-A - B, -C, -D, -E, -F]
    # a linear system with Cardano's roots: quickly, and right (#16245)
    M = Sym[15 -4 1; -10 12 -6; 20 4 -2]
    @syms x1 x2 x3
    λ = LA.eigvals(M)[1]
    s = only(solve(collect(M * [x1, x2, x3] .~ λ * [x1, x2, x3]), [x1, x2, x3]))
    v = [N(subs(s[q], x3 => 1)) for q in (x1, x2, x3)]
    @test maximum(abs.(N.(M) * v - N(λ) * v)) < 1e-10
    # eigenvalues too large for symbolic eigenvectors: an error, not a hang (#15076)
    @syms y z
    Mc = Sym[1 cos(x/4+y/4) cos(x/4+z/4) cos(y/4+z/4); cos(x/4+y/4) 1 cos(y/4-z/4) cos(x/4-z/4);
             cos(x/4+z/4) cos(y/4-z/4) 1 cos(x/4-y/4); cos(y/4+z/4) cos(x/4-z/4) cos(x/4-y/4) 1]
    @test_throws ArgumentError LA.eigen(Mc)
end

@testset "missing cases (SymPy review, round 11)" begin
    @syms x y z n a F t c
    IM, PI = SymJ.IM, SymJ.PI
    SF = SymJ.SpecialFunctionsSupport.SF
    # log of a power of |u|; mod of integer coefficients (#15854, #15493)
    @test simplify(log(abs(x)^2) / log(abs(x))) == 2
    @test refine(mod(3n, 2); assume = (n ∈ ℤ,)) == mod(n, 2)
    # Hurwitz ζ and ψ at half-integers (round-11 review)
    @test SF.zeta(Sym(2), Sym(1) // 2) == PI^2 / 2
    @test SF.zeta(Sym(2), Sym(3) // 2) == PI^2 / 2 - 4
    @test N(SF.digamma(Sym(1) // 2)) ≈ -0.5772156649015329 - 2log(2)
    # √ of a negative constant sum: i √(-u)
    @test sqrt(-sqrt(Sym(3)) - 2) == IM * sqrt(sqrt(Sym(3)) + 2)
    # dsolve's constants avoid function names (#15056)
    C1 = SymFunction(:C1)
    @test SymJ.dsolve(diff(C1(t), t) ~ c * C1(t), C1(t)) == (C1(t) ~ Sym(:C2) * exp(c * t))
    # |u| = 0 over ℂ; a power equal to 1; a base with parameters (#15731)
    @test solve(abs(x - 2) ~ 0, x) == [2]
    @test solve((x^2 - 7x + 11)^(x^2 - 13x + 42) ~ 1, x) == [2, 3, 4, 5, 6, 7]
    @test solve(x^x ~ 1, x) == [1]
    @test solve(y^x ~ z, x) isa SymJ.Cases
    # angle with a positive real part (#15945)
    @test refine(angle(a + IM * F); assume = (a ≫ 0, F ∈ ℝ)) == atan(F / a)
    # constants inside function arguments simplified (#15129)
    @test simplify(acos(cos(Sym(1) // 50) * cos(Sym(1) // 25) + sin(Sym(1) // 50) * sin(Sym(1) // 25))) == Sym(1) // 50
    # rational powers of integers merged by primes (#15669)
    @test simplify(Sym(2)^(Sym(4) // 5) * Sym(4)^(Sym(3) // 5)) == 4
    # parts of constant complex powers (#16348)
    @test real(IM^IM) == exp(-PI / 2)
    @test imag(sqrt(IM)) == sqrt(Sym(2)) / 2
    # conjugate roots of unity paired (#15594)
    u = SymFunction(:u)
    @test SymJ.dsolve(diff(u(x), x, 5) - u(x) ~ 0, u(x)) isa Relation
    # gather inside denominators (#15394)
    @test SymJ.gather(1 + 1 / (x * y + x), x) == 1 + 1 / (x * (y + 1))
end


@testset "wrong answers (SymPy review, round 12)" begin
    @syms x y w n s t
    f = SymFunction(:f)
    # power series of radius 1 by cases: past |z| = 1 they diverge (#14213)
    @test SymJ.summation(Sym(3)^n / n, n, 1, Inf) == Inf
    @test_throws ArgumentError SymJ.summation((-1)^(n - 1) * Sym(3)^n / n, n, 1, Inf)
    @test_throws ArgumentError SymJ.summation(x^n / n, n, 1, Inf; assume = (x ≫ 2,))
    @test SymJ.summation(x^n / n, n, 1, Inf) isa SymJ.Cases
    @test SymJ.summation(1 / (n * Sym(2)^n), n, 1, Inf) == log(Sym(2))
    @test SymJ.summation((-1)^(n - 1) / n, n, 1, Inf) == log(Sym(2))
    @test SymJ.summation(binomial(2n, n), n, 0, Inf) == Inf
    @test SymJ.summation(binomial(2n, n) / Sym(8)^n, n, 0, Inf) == sqrt(Sym(2))
    @test SymJ.summation((-1)^n * binomial(2n, n) / Sym(4)^n, n, 0, Inf) == sqrt(Sym(2)) / 2
    @test SymJ.summation(x^n * binomial(2n, n), n, 0, Inf) isa SymJ.Cases
    # a term undefined at the start of the range: no sum (#14410)
    @test_throws SymJ.DivisionByZero SymJ.summation(3 - 1 / sqrt(x), x, 0, Inf)
    @test SymJ.summation(3 - 1 / sqrt(x), x, 1, Inf) == Inf
    # cos(∞) is no pole: left as it is (#13499)
    @test SymJ._hasfn(subs(cos(x), x => Inf), :cos)
    # a negative order of differentiation (#13873)
    @test_throws ArgumentError diff(sin(x), x, -1)
    # 1^∞ is indeterminate, not the float 1.0 (#13914)
    @test kind(Sym(1)^Sym(Inf)) == POW
    # substitution into a derivative captures no variable (#13795)
    @test subs(diff(f(x, y), y), y => x) == SymJ.apply_function(:Derivative, Sym[f(x, y), y, x])
    yp = Sym(Symbol("y′"))
    @test subs(diff(f(x, y), y), x => y) == SymJ.apply_function(:Derivative, Sym[f(y, yp), yp, y])
    @test subs(diff(f(x), x), x => y) == diff(f(y), y)
    # a repeated pole with float coefficients (#14326)
    @test SymJ.inverse_laplace_transform(6 / (s - 0.5)^4, s, t) == t^3 * exp(t / 2)
    # integer and rational unknowns of a system (#14170)
    @test isempty(solve([x + y ~ 1, x - y ~ 0], [x, y]; assume = (x ∈ ℤ, y ∈ ℤ)))
    @test solve([x + y ~ 2, x - y ~ 0], [x, y]; assume = (x ∈ ℤ, y ∈ ℤ)) == [Dict(x => 1, y => 1)]
    @test isempty(solve([x^2 ~ 2, y ~ 1], [x, y]; assume = (x ∈ ℚ,)))
    # back-substitution at 256 bits: no points with residuals to 1e34 (#14799)
    S = solve([x^2 + y^2 + w^2 ~ 1, x^5 + y - 3 // 10 ~ 0, x * y * w ~ 1 // 100], [x, y, w])
    @test length(S) == 22
    res(d) = (z = [ComplexF64(N(d[v])) for v in (x, y, w)];
              maximum(abs, (z[1]^2 + z[2]^2 + z[3]^2 - 1, z[1]^5 + z[2] - 0.3, z[1] * z[2] * z[3] - 0.01)))
    @test maximum(res, S) < 1e-12
end

@testset "crashes and hangs (SymPy review, round 12)" begin
    @syms x y z a b k n th
    LA = SymJ.LinearAlgebra
    # sums and derivatives of Cases (#14235, #14211)
    r = SymJ.summation(SymJ.summation(a^n * b^k, k, 0, Inf), n, 0, Inf)
    @test any(cv -> cv.second isa Sym && iszero(simplify(cv.second - 1 / ((1 - a) * (1 - b)))), r)
    d = diff(SymJ.summation(x^k, k, 0, Inf), x)
    @test any(cv -> cv.second isa Sym && iszero(simplify(cv.second - 1 / (1 - x)^2)), d)
    # cse of a matrix and of a tuple of matrices, in their shapes (#14118)
    reps, red = SymJ.cse(Sym[x+y (x+y)^2; z x+y])
    @test length(reps) == 1
    @test size(red) == (2, 2)
    reps2, red2 = SymJ.cse((Sym[x+y (x+y)^2; z x+y], Sym[x + y, 2]))
    @test red2[2] == Sym[first(reps2[1]), 2]
    # the 1- and ∞-norms of a symbolic matrix (#13985)
    @test LA.opnorm(Sym[1 -2; 3 4], 1) == 6
    @test LA.opnorm(Sym[1 -2; 3 4], Inf) == 7
    # a series simplified keeps its order term; numer_denom of one is an error (#14622)
    @test simplify(series(sin(x) / x, x, 0, 5)) isa SymJ.SymSeries
    @test_throws ArgumentError numer_denom(series(sin(x) / x, x, 0, 5))
    # e^(-18000000x) is not factored as a polynomial of that degree (#14300)
    s = solve(1 - exp(-18000000x) ~ y, x)
    @test s isa SymJ.Cases
    @test factor(exp(2x) - 1) == (exp(x) + 1) * (exp(x) - 1)
    # constant residuals of radical roots decided numerically (#13625): two real roots
    # and two complex ones
    rs = solve(sin(th) + cos(th) + sin(th) * cos(th) ~ 0, th)
    @test length(rs) == 4
    @test maximum(r -> (v = ComplexF64(N(r)); abs(sin(v) + cos(v) + sin(v) * cos(v))), rs) < 1e-12
end

@testset "missing cases (SymPy review, round 12)" begin
    @syms x y z n m k t
    LA = SymJ.LinearAlgebra
    # the principal root is never negative (#14815)
    @test refine(sqrt(x) ≪ 0; assume = (x ∈ ℝ,)) === false
    @test refine(sqrt(x) ≪ 0) === false
    # divergence of series with no closed form (#14579, #14112, #14111)
    @test SymJ.summation(sqrt(x), x, 1, Inf) == Inf
    @test SymJ.summation(1 / sqrt(n), n, 1, Inf) == Inf
    @test SymJ.summation(1 / log(log(n)), n, 22, Inf) == Inf
    @test SymJ._hasfn(SymJ.summation((-1)^n / sqrt(n), n, 1, Inf), :Sum)
    # \ for a non-square matrix: least squares, least norm (#14548)
    @test Sym[2 1; 1 0; 1 0] \ Sym[3, 1, 1] == Sym[1, 1]
    @test Sym[1 1 0; 0 1 1] \ Sym[1, 1] == Sym[1 // 3, 2 // 3, 1 // 3]
    # im stays im when numbers are floated (#14371)
    @test !occursin("sqrt(-1.0)", string(apart(z^2 / (z + 2 + 1.0im)^2, z)))
    # exp of a real matrix with complex eigenvalues: real entries (#14330)
    @test exp(Sym[1 -1; 1 1]) == Sym[cos(Sym(1))*SymJ.E -sin(Sym(1))*SymJ.E; sin(Sym(1))*SymJ.E cos(Sym(1))*SymJ.E]
    @test exp(t * Sym[1 -1; 1 1])[1, 1] == cos(t) * exp(t)
    # |b^q| for constant b and q, i^q (#14274, #14260, #14262)
    @test abs((4 + 5IM)^(6 + 7IM)) == 68921 * exp(-7atan(Sym(5) // 4))
    @test refine(abs(IM^n); assume = (n ∈ ℝ,)) == 1
    @test IM^IM == exp(-PI / 2)
    @test simplify(log(IM^IM)) == -PI / 2
    # harmonic numbers in limits, as digamma (#14204, #14125); an unevaluated sum refused (#14196)
    @test limit(SymJ.summation(1 / k, k, n, 2n), n => Inf) == log(Sym(2))
    @test limit(n^2 * (SymJ.harmonic(n) - log(n + 1 // 2) - Sym(MathConstants.eulergamma)), n => Inf) == Sym(1) // 24
    @test_throws ArgumentError limit(SymJ.summation(m^k, m, 1, n) / n^(k + 1), n => Inf; assume = (k ∈ ℤ, k ≫ 0))
    # telescoping sums; harmonic numbers of shifted arguments combined (#14195)
    @test SymJ.summation(log(m + 1) - log(m), m, 1, n) == log(n + 1)
    @test SymJ.summation(sqrt(m) - sqrt(m - 1), m, 1, n) == sqrt(n)
    @test SymJ.summation(1 / m - 1 / (m + 1), m, 1, n) == 1 - 1 / (n + 1)
    # an empty harmonic sum (#14185)
    @test SymJ.harmonic(0, n) == 0
    # the binomial theorem in two variables (#14142)
    @test SymJ.summation(binomial(n, m) * x^m * y^(n - m), m, 0, n) == (x + y)^n
    # u(x) δ(x - a) = u(a) δ(x - a) (#13996)
    @test simplify(x * SymJ.dirac(x)) == 0
    @test simplify(x^2 * SymJ.dirac(x - 2)) == 4SymJ.dirac(x - 2)
    # signs of binomial and factorial (#13617, #13616)
    @test refine(binomial(n, k) ≧ 0; assume = (n ∈ ℤ, n ≧ 0, k ∈ ℤ, k ≧ 0)) === true
    @test refine(factorial(n) ≫ 0; assume = (n ∈ ℤ, n ≧ 0)) === true
    # a mod by the same modulus inside a mod (#13543)
    @test mod(mod(x + 1, 2) + 1, 2) == mod(x, 2)
end

@testset "wrong answers (SymPy review, round 13)" begin
    @syms x y z a a1 b1 c1
    # zeros at the poles of tan, cot, sec, csc in a denominator (#12527)
    @test solve(1 / tan(x) ~ 0, x) == [PI / 2, 3PI / 2]
    @test solve(1 / sec(x) ~ 0, x) == [PI / 2, 3PI / 2]
    @test solve(sin(x) / tan(x) ~ 0, x) == [PI / 2, 3PI / 2]
    # an unknown call with other functions of its variable (#12129)
    @test solve(cos(2a) ~ 0, cos(a)) == [-sqrt(Sym(2)) / 2, sqrt(Sym(2)) / 2]
    @test solve(cos(3a) ~ 0, cos(a)) == [-sqrt(Sym(3)) / 2, 0, sqrt(Sym(3)) / 2]
    @test_throws ArgumentError solve(sin(2a) ~ 0, sin(a))
    # a family narrowed by the facts about its free unknowns (#12932)
    @test isempty(solve([a1 + b1 ~ 0], [a1, b1]; assume = (a1 ≫ 0, b1 ≫ 0)))
    @test solve([a1 + b1 ~ 0, c1 ~ 1], [a1, b1, c1]; assume = (a1 ≧ 0, b1 ≧ 0)) == [Dict(a1 => 0, b1 => 0, c1 => 1)]
    # integer unknowns: a family with no integer member (#12546)
    @test isempty(solve([2x + 4y ~ -3], [x, y]; assume = (x ∈ ℤ, y ∈ ℤ)))
    @test isempty(solve([2x + 4y + 6z ~ 3], [x, y, z]; assume = (x ∈ ℤ, y ∈ ℤ, z ∈ ℤ)))
    @test solve([2x + 4y ~ 6], [x, y]; assume = (x ∈ ℤ, y ∈ ℤ)) == [Dict(x => 3 - 2y, y => y)]
    # |u| in a system whose solutions are families or have parameters (#13306)
    @test_throws ArgumentError solve([abs(x) + y ~ 1], [x, y]; assume = (x ∈ ℝ, y ∈ ℝ))
    @test solve([abs(x) + y ~ 1], [x]; assume = (x ∈ ℝ, y ∈ ℝ)) isa SymJ.Cases
    # complex coefficients in an inequality: real values only where the imaginary part
    # vanishes (#13105)
    @test solve((x + im) * (x + 2im) ≪ 0, x) == SymJ.SymSet([SymJ.SymInterval(Sym(0), Sym(0), true, true)])
    @test isempty(solve(x + im ≪ 0, x))
    # a symbolic exponent that may vanish: cases (#12774)
    @test solve(x^a ~ 1, x) isa SymJ.Cases
    @test solve(x^a ~ 1, x)[2] == ([a ≶ 0] => [1])
    @test solve(x^a ~ 1, x; assume = (a ≶ 0,)) == [1]
    # the poles of a condition are breakpoints (#12095)
    @test solve(abs((x - 1) / (x - 5)) ≦ 1 // 3, x) == SymJ.SymSet([SymJ.SymInterval(Sym(-1), Sym(2), true, true)])
    @test 9 in solve(abs((x - 1) / (x - 5)) ≦ 2, x)
    @test solve(abs(x - 1) / abs(x - 5) ≦ 1 // 3, x) == SymJ.SymSet([SymJ.SymInterval(Sym(-1), Sym(2), true, true)])
end

@testset "crashes and hangs (SymPy review, round 13)" begin
    @syms x y a b c d e f g w root
    # a zero polynomial to a power (#12114)
    @test iszero(SymJ.Poly{SymJ.SymNumber}(Sym(0), [x])^2)
    # a symbol is no rule for rewrite (#12639)
    @test_throws ArgumentError SymJ.rewrite(max(a, b), :piecewise)
    # NaN in a piecewise branch of a system (#12414)
    @test length(solve([SymJ.piecewise((x ≫ 0) => y^2 - x + 1, true => NaN) ~ 0, x ~ 2], [x, y])) == 2
    # a symbolic Cholesky factorization (#12500)
    @test SymJ.LinearAlgebra.cholesky(Sym[4 2; 2 3]).U == Sym[2 1; 0 sqrt(Sym(2))]
    @test SymJ.LinearAlgebra.cholesky(Sym[4 12 -16; 12 37 -43; -16 -43 98]).U == Sym[2 6 -8; 0 1 5; 0 0 3]
    @test_throws SymJ.LinearAlgebra.PosDefException SymJ.LinearAlgebra.cholesky(Sym[1 2; 2 1])
    # a gcd that is a constant because one side lacks a variable (#13126)
    @test gcd(x * (a + b) + 1, (a + b) * (y + 2)) == 1
    @test iszero(expand(resultant(x + c, root + x * (a^2 + b * c) + a * b * c, x) - (root - c * (a^2 + b * c) + a * b * c)))
    # back-substitution through a linear basis element, unexpanded (#12345)
    @test length(solve([x^2 - y ~ sqrt(Sym(2)), x^2 + x * y - y^2 ~ 0], [x, y])) == 4
    # a positive-dimensional system: no empty reduce (#12114)
    @test solve([1 + a * b + d * e ~ 0, 1 + a * c + d * f ~ 0, 1 + b * c + e * f ~ 0, g - a^2 - d^2 ~ 0, g - b^2 - e^2 ~ 0, g - c^2 - f^2 ~ 0], [a, b, c, d, e, f, g]) isa Vector
end

@testset "missing cases (SymPy review, round 13)" begin
    @syms x y z a b n k v
    # a relation whose sides simplify to one constant is decided (#12522)
    @test simplify(cos(x)^2 + sin(x)^2 ~ 1) === true
    # δ(c u) = δ(u)/|c| (#12257)
    @test simplify(SymJ.dirac(-x) - SymJ.dirac(x)) == 0
    @test SymJ.dirac(2x - 2) == SymJ.dirac(x - 1) / 2
    # complex floats multiplied and raised to integer powers (#12154)
    @test Sym(2.0 + 1.0im)^2 == Sym(3.0 + 4.0im)
    @test Sym(2.0 + 1.0im) * Sym(1.0 + 3.0im) == Sym(-1.0 + 7.0im)
    # min/max absorption (#12638)
    @test min(a, b, max(a, y)) == min(a, b)
    @test min(1, y, max(4, x)) == min(1, y)
    @test max(a, min(a, y)) == a
    # integer multiples of the modulus dropped (#12592)
    @test mod(3n + z, n) == mod(z, n)
    @test refine(mod(n * y + z, n); assume = (y ∈ ℤ, n ∈ ℤ)) == mod(z, n)
    # angle(x + iy) = atan(y, x) for real x, y (#11988)
    @test simplify(angle(x + y * im) - atan(y, x); assume = (x ∈ ℝ, y ∈ ℝ)) == 0
    # ψ and ζ(s, a) at negative half-integers (#12753)
    @test SymJ.SpecialFunctionsSupport.SF.digamma(Sym(-1) // 2) == SymJ.SpecialFunctionsSupport.SF.digamma(Sym(1) // 2) + 2
    @test simplify(SymJ.summation(1 / ((2v + 1)^2 - 4)^2, v, 0, Inf)) == PI^2 / 64
    # Σ 1/k^p for symbolic p: ζ(p) (#12994)
    @test SymJ.summation(1 / k^a, k, 1, Inf; assume = (a ≫ 1,)) == SymJ.SpecialFunctionsSupport.SF.zeta(a)
    @test SymJ.summation(1 / k^a, k, 1, Inf; assume = (a ≦ 1,)) == Inf
    # √(c(a + bi)) denested with the rational content split off (#12420)
    @test simplify(sqrt(Sym(8) + 6IM)) == 3 + IM
    # solve inverts lambertw (#12783)
    @test solve(SymJ.lambertw(n) ~ 2, n) == [2exp(Sym(2))]
    @test isempty(solve(SymJ.lambertw(n) ~ -2, n))
    # cancel with powers of one exponential as one generator (#12787)
    @test cancel((exp(2x) - 1) / (exp(x) - 1)) == exp(x) + 1
    @test cancel((x^(2a) - 1) / (x^a - 1)) == x^a + 1
    # finite sums of sin and cos of arithmetic arguments (#12107)
    @test SymJ.summation(sin(k * x), k, 1, n) isa SymJ.Cases
    @test N(subs(SymJ.summation(cos(2k + 1), k, 0, n), n => 6)) ≈ sum(cos(2j + 1) for j in 0:6)
end

@testset "wrong answers (SymPy review, round 14)" begin
    @syms x y
    # float matrices factored with pivoting by magnitude (#11549, #10780)
    A = Sym[1e-17 1 1 1; 1 1e-17 -1 1; -1 1 1 1; -1 -1 1 1]
    @test N(SymJ.LinearAlgebra.det(A)) ≈ 4.0
    B = Sym[6e-17 1 4; -1 0 8; 0 0 1]
    @test maximum(abs.(N.(B * inv(B) - Sym[1 0 0; 0 1 0; 0 0 1]))) < 1e-12
    @test N.(Sym[1e-20 1; 1 1] \ Sym[1, 2]) ≈ [1.0, 1.0]
    # lu: no pivot that is 0 by sin² + cos² = 1 (#11905)
    F = SymJ.LinearAlgebra.lu(Sym[sin(x)^2 + cos(x)^2 - 1 1; 1 1])
    @test F.p == [2, 1]
    # float results outside Float64's range: an error, not 0 or Inf (#11004)
    @test_throws ArgumentError Sym(1e300) * Sym(1e300)
    @test_throws ArgumentError Sym(2.0)^Sym(-2000)
    # ζ(s, a) for a < 0 numerically as Hurwitz's (#10592)
    Z = SymJ.SpecialFunctionsSupport.SF.zeta
    @test N(Z(Sym(0), Sym(-3))) ≈ 3.5
    @test N(Z(Sym(-1), Sym(-3))) ≈ -73 / 12
    @test Z(Sym(0), Sym(-3)) == Sym(7) // 2
    # a point piece's derivative is the limit of its neighbour's (#11402)
    p = SymJ.piecewise(x ≶ 0 => sin(x) / x, true => 1)
    @test diff(p, x, 2).args[4] == -Sym(1) // 3
    # roots with |y| checked by y's sign (#11534)
    @test solve(-y + x / sqrt(1 - x^2) ~ 0, x; assume = (x ∈ ℝ, y ≫ 0)) == [sqrt((4y^2 + 4) * y^2) / (2y^2 + 2)]
    @test length(solve(-y + x / sqrt(1 - x^2) ~ 0, x; assume = (x ∈ ℝ, y ∈ ℝ))) == 3
    # 0^x = 1 at x = 0 (#11536)
    @test solve(Sym(0)^x ~ 1, x) == [0]
end

@testset "crashes and hangs (SymPy review, round 14)" begin
    @syms x y n
    # factor of 0 (#11151)
    @test factor(Sym(0)) == 0
    @test SymJ.factor_squarefree(Sym(0)) == 0
    # lambdify of a relation (#11248)
    @test lambdify(x ~ y, [x, y])(1, 1) === true
    @test lambdify(x ≪ y, x, y)(2, 1) === false
    # piecewise_fold of a relation (#11120)
    @test SymJ.piecewise_fold(Gt(SymJ.piecewise(x ≫ 0 => x, true => 0) - 1, 2), x) isa Relation
    # Bools among a system's equations (#11659)
    @test solve(Any[true, x ~ y], [x, y]) == [Dict(x => y, y => y)]
    @test isempty(solve(Any[false, x ~ y], [x, y]))
    # cos x + i sin x = 1: tan(x/2) = i is no point (#10809)
    @test solve(cos(x) + im * sin(x) ~ 1, x) == [0]
    # svd of an empty matrix (#10768)
    @test isempty(SymJ.LinearAlgebra.svdvals(Matrix{Sym}(undef, 0, 0)))
    # huge fractional powers kept as powers; their value through logs (#11004)
    v = subs((Sym(2) // 3)^(n + Sym(1) // 2), n => Sym(big(602214179000000000000000)))
    @test v.kind == SymJ.POW
    st(m) = sqrt(2PI * m) * (m / SymJ.E)^m
    q = simplify(st(n) / st(n / 2)^2 / 2^n / (st(n) / st(n / 3)^3 / 3^n))
    @test N(subs(q, n => Sym(big(602214179000000000000000)))) ≈ 7.487098316702582e11 rtol = 1e-8
    # simplify of a constant with sin(12) and erf of complex arguments: no polar rewriting
    # of everything (the endpoint value of ∫₁³ cos(2x)²√x dx; #11354)
    erf_ = SymJ.SpecialFunctionsSupport.SF.erf
    c = Sym(3)^(Sym(3) // 2) / 3 + sin(Sym(12)) * sqrt(Sym(3)) / 8 +
        IM * erf_(2sqrt(-IM) * sqrt(Sym(3))) * sqrt(PI) / (64sqrt(-IM)) - IM * erf_(2sqrt(IM) * sqrt(Sym(3))) * sqrt(PI) / (64sqrt(IM))
    @test simplify(c) isa Sym
end
