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
