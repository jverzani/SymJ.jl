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
