# Assumptions, refine, simplify and rewriting (included by runtests.jl)

@testset "assumptions" begin
    @syms a b n
    J = SymJ
    sg(e, A) = (m = J._signs(e, A); m === nothing ? nothing :
                Set(s for (bit, s) in ((0x01, -1), (0x02, 0), (0x04, 1)) if m & bit != 0))
    A = Assumptions(x ≫ 1, n ∈ ℤ, n ≫ 0, a ≦ b, y ∈ ℝ)
    @test sg(x, A) == sg(x - 1, A) == sg(x^2 - 1, A) == sg(1 / x, A) == sg(log(x), A) == Set([1])
    @test sg(x^3 - x, A) == Set([1])                       # through its factors
    @test sg(sqrt(x - 1), A) == sg(exp(y), A) == sg(y^2 + 1, A) == Set([1])
    @test sg(y^2, A) == sg(b - a, A) == Set([0, 1])
    @test sg(a - b, A) == Set([-1, 0])
    @test sg(a, A) == Set([-1, 0, 1]) && sg(z, A) === nothing   # z complex
    @test sg(n - 1, A) == Set([0, 1])                     # n ≥ 1
    @test J._isinteger(2n + 1, A) === true && J._isinteger(x, A) === nothing
    @test J._isreal(y, A) === true && J._isreal(z, A) === nothing
    @test J._isreal(z, Assumptions(z ∉ ℝ)) === false
    @test J._isinteger(x, Assumptions(x ∉ ℤ)) === false
    # linear relations among several symbols
    B = Assumptions(x ≫ 0, y ≫ 0, x + y ≪ 1)
    @test sg(x + y - 1, B) == sg(x - 1, B) == sg(x + 2y - 2, B) == Set([-1])
    @test sg(x * y, B) == Set([1])
    # equalities are substituted, the lone symbol on the left
    C = Assumptions(x ~ 2y + 1, y ≫ 0)
    @test C.subst == Dict(x => 2y + 1) && sg(x - 1, C) == Set([1])
    # sets: Not and Or about one symbol; integers narrowed
    @test string(Assumptions(!(Sym(0) ≦ x) | (x ≫ 1)).sets[x]) == "(-∞, 0) ∪ (1, ∞)"
    @test string(Assumptions(n ∈ ℤ, n ≫ 0, n ≪ 3).sets[n]) == "[1, 2]"
    @test_throws ArgumentError Assumptions((x ≫ 0) | (y ≫ 0))       # Or across symbols
    # contradictions
    for bad in ((x ≫ 1, x ≪ 0), (n ∈ ℤ, n ≫ 0, n ≪ 1), (x ≫ 0, y ≫ 0, x + y ≪ 0),
                (x ~ 2, x ≶ 2), (n ∈ ℤ, n ∉ ℚ), (exp(x) ≫ 2, x ≪ 0), (z ∉ ℝ, z ≫ 0))
        @test_throws ArgumentError Assumptions(bad...)
    end
    # ≠ on a complex symbol
    F = Assumptions(a ≶ 0)
    @test J._isnonzero(a, F) === true && J._isnonzero(a^2, F) === true && sg(abs(a), F) == Set([1])
    G = Assumptions(x ∈ ℝ)
    @test sg(abs(x), G) == Set([0, 1]) && sg(exp(x) - 1, G) == Set([-1, 0, 1])
    # membership on numbers is a Bool
    @test 2 ∈ ℤ && !(Sym(1) / 2 ∈ ℤ) && Sym(π) ∈ ℝ && 1.5 ∈ ℚ && !(π ∈ ℚ) && im ∉ ℝ
    # signs from the range on an interval: critical points and end values (SymPy #29342)
    let y = Sym(:y)
        @test refine(2^y - y ≫ 0; assume = (y ≫ 0,)) == true
        @test refine(log(x) - x ≪ 0; assume = (x ≫ 0,)) == true
        @test refine(exp(x) - x - 1 ≧ 0; assume = (x ∈ ℝ,)) == true
        @test refine(exp(x) - x - 1 ≫ 0; assume = (x ≫ 0,)) == true
        @test refine(exp(x) - x - 1 ≫ 0; assume = (x ∈ ℝ,)) isa Relation   # 0 at x = 0
        @test refine(x - log(x) - 2 ≫ 0; assume = (x ≫ 0,)) isa Relation   # negative at 1
        @test limit(x^y / x^(2^y), x => Inf; assume = (y ≫ 0,)) == 0
    end
    # a polynomial in π with algebraic coefficients is irrational (SymPy #27993)
    @test refine(PI + PI^2 ∈ ℚ) == 0
    @test refine(sqrt(Sym(2)) * PI + PI^2 + 1 ∈ ℚ) == 0
    @test_throws ArgumentError refine(PI + E ∈ ℚ)
    @test solve(x - PI - PI^2 ~ 0, x; assume = (x ∈ ℚ,)) == Sym[]
    @test (x ∈ ℤ) isa SymJ.Membership && string(x ∉ ℤ) == "x ∉ ℤ"
    # blocks nest and add to the keyword
    assuming(x ≫ 0) do
        @test J._assumptions(nothing).facts == Any[x ≫ 0]
        assuming(y ≫ 0) do
            @test length(J._assumptions((n ∈ ℤ,)).facts) == 3
        end
    end
    @test J._assumptions(nothing) === nothing
    @test Assumptions(A, z ≫ 0).facts[end] == (z ≫ 0)
    # @assume: Julia's operators; ASCII names of the sets
    B = @assume(x > 0 && 1 < y <= 2, n in Integers, a != 0)
    @test B.facts == Any[x ≫ 0, 1 ≪ y, y ≦ 2, n ∈ ℤ, a ≶ 0]
    @test string(B) == "Assumptions(x > 0, 1 < y, y <= 2, n ∈ ℤ, a != 0)"
    @test string(@assume(!(0 <= x < 1)).sets[x]) == "(-∞, 0) ∪ [1, ∞)"
    @test @assume(x ≥ 1, x ≠ 2, x == y + 1).subst == Dict(x => y + 1)
    @test_throws ArgumentError @assume(x > 1 && x < 0)
    @test Integers === ℤ && Rationals === ℚ && Reals === ℝ && Complexes === ℂ
    @test Sym(2) in Integers && (n in Reals) isa SymJ.Membership
end

@testset "refine and simplify" begin
    @syms a n
    R(e, as...) = refine(e; assume = Assumptions(as...))
    @test R(sqrt(x^2), x ≫ 0) == x
    # inverse compositions where the inner value is in the range (SymPy #28376)
    @test R(atanh(tanh(x)), x ∈ ℝ) == x
    @test R(asinh(sinh(x)), x ∈ ℝ) == x
    @test R(atanh(tanh(x))) == atanh(tanh(x))
    @test R(atan(tan(x)), x ≫ -1, x ≪ 1) == x
    @test R(atan(tan(x)), x ≫ 0) == atan(tan(x))
    @test R(acos(cos(x)), x ≧ 0, x ≦ 1) == x
    # acos, asin outside [-1, 1] (SymPy #28541)
    @test simplify(real(acos(Sym(2)))) == 0
    @test N(simplify(imag(acos(Sym(2))))) ≈ imag(acos(complex(2.0)))
    @test R(real(acos(x)), x ≫ 1) == 0
    @test R(real(asin(x)), x ≪ -1) == -PI / 2
    # an equation: common terms taken off both sides (SymPy #28666)
    @test simplify(x^2 + 1 ~ a^2 + 1) == (x^2 ~ a^2)
    @test simplify(x + 1 ~ a + 1) == (x ~ a)
    # contradictory conditions of a piece (SymPy #28323)
    @test simplify(SymJ.piecewise((x ≧ 0, x^3 ≪ -8) => 1, true => 0)) == 0
    @test simplify(SymJ.piecewise((x ≧ 0, x ≪ 2) => 1, true => 0)) == SymJ.piecewise((x ≧ 0, x ≪ 2) => 1, true => 0) && R(sqrt(x^2), x ∈ ℝ) == abs(x) && R(sqrt(x^2)) == sqrt(x^2)
    @test R(sqrt(x^4), x ∈ ℝ) == x^2 && R((x^3)^(1 // 3), x ≫ 0) == x
    @test R((x^3)^(1 // 3), x ∈ ℝ) == (x^3)^(1 // 3)          # a real cube root is not x^(1/3)'s branch
    @test R(abs(x), x ≪ 0) == -x && R(abs(2x * y), x ≫ 0) == 2x * abs(y) && R(sign(x - 1), x ≫ 2) == 1
    @test R(log(exp(x)), x ∈ ℝ) == x && R(log(exp(x))) == log(exp(x)) && R(log(x^2), x ≫ 0) == 2log(x)
    # exp(iπ w): integer multiples of 2πi drop out (SymPy #29196)
    @test R(exp(2PI * IM * n), n ∈ ℤ) == 1
    @test R(exp(PI * IM * n), n ∈ ℤ) == (-1)^n
    @test R(exp(2PI * IM * (n + 1 // 2)), n ∈ ℤ) == -1
    @test R(exp(3PI * IM * n + IM * x), n ∈ ℤ) == (-1)^n * exp(IM * x)
    @test R(exp(PI * IM * n / 2), n ∈ ℤ) == exp(PI * IM * n / 2)
    @test R(exp(2PI * IM * x)) == exp(2PI * IM * x)
    # half angles under a root, and the signs of sin, cos on an interval (SymPy #29371)
    @test simplify(sqrt((1 + cos(x)) / 2); assume = (x ≫ 0, x ≪ PI)) == cos(x / 2)
    @test simplify(sqrt(1 + cos(x)) * sqrt(1 - cos(x)); assume = (x ≫ 0, x ≪ PI)) == sin(x)
    @test R(sin(x) ≫ 0, x ≫ 0, x ≪ PI) === true
    @test R(cos(x) ≪ 0, x ≫ PI / 2, x ≪ 3PI / 2) === true
    @test R(cos(x) ≫ 0, x ≫ 0, x ≪ PI) isa Relation                      # changes sign
    @test R(sin(x) ≫ 0, x ≧ 0, x ≦ PI) isa Relation                       # 0 at the ends
    @test R((-1)^(2n), n ∈ ℤ) == 1 && R((-1)^(2n + 1), n ∈ ℤ) == -1
    @test R(sin(n * PI), n ∈ ℤ) == 0 && R(cos(n * PI), n ∈ ℤ) == (-1)^n
    @test R(sin(x + 2n * PI), n ∈ ℤ) == sin(x) && R(cos(x + n * PI), n ∈ ℤ) == (-1)^n * cos(x)
    @test R(floor(n), n ∈ ℤ) == n && R(floor(x + n), n ∈ ℤ) == floor(x) + n
    @test R(sqrt(exp(x)), x ∈ ℝ) == exp(x / 2) && R(x^2 + a, a ~ 3) == x^2 + 3
    @test R(x ≫ 0, x ≫ 1) === true && R(x ≫ 2, x ≫ 1) == (x ≫ 2) && R(x^2 ≧ 0, x ∈ ℝ) === true
    @test (@assuming x > 0 begin refine(sqrt(x^2)) end) == x
    # a radical's real factors to even powers come out: (uv)^q = u^q v^q for u ≥ 0
    @test R(sqrt(x^2 * (a^2 + 1)), x ≫ 0) == x * sqrt(a^2 + 1)
    @test R(sqrt(x^2 * (a^2 + 1)), x ∈ ℝ) == abs(x) * sqrt(a^2 + 1)
    @test R(sqrt(4x^2 * a), x ∈ ℝ) == 2abs(x) * sqrt(a)
    @test R(sqrt(-x^2 * a), x ∈ ℝ) == abs(x) * sqrt(-a)
    @test R(sqrt(x * a), x ≫ 0, a ≫ 0) == sqrt(x * a)                  # not split
    # refine = true: the block's value refined with the facts
    @test (@assuming x > 0 begin sqrt(x^2) end) == sqrt(x^2)
    @test (@assuming refine = true x > 0 begin sqrt(x^2) end) == x
    @test (@assuming refine = false x > 0 begin sqrt(x^2) end) == sqrt(x^2)
    @test (@assuming refine = true x > 0 begin Sym[sqrt(x^2), x^2] end) == [x, x^2]
    @test (@assuming refine = true x > 0 begin (sqrt(x^2), 3) end) == (x, 3)
    @test (@assuming refine = true x > 0 begin sqrt(x^2) ≫ 0 end) === true
    @test assuming(x ≫ 0; refine = true) do
        Dict(:a => sqrt(x^2))
    end == Dict(:a => x)
    @test (@assuming refine = true x > 0 begin "text" end) == "text"
    @test_throws ArgumentError macroexpand(@__MODULE__, :(@assuming foo = 1 x > 0 begin x end))
    S(e, as...) = simplify(e; assume = Assumptions(as...))
    @test S((x^2 - 1) / (x - 1)) == x + 1 && S((x + 1)^2 - x^2) == 2x + 1
    @test S(sqrt(Sym(2)) * sqrt(Sym(3))) == sqrt(Sym(6)) && S(Sym(2)^(1 // 3) * Sym(4)^(1 // 3)) == 2
    @test S(sqrt(x^2) / x, x ≫ 0) == 1 && S(sin((x^2 - 1) / (x - 1))) == sin(x + 1)
    @test S(x^2 + 2x + 1) == (x + 1)^2 && S(x / (x + 1) + 1 / (x + 1)) == 1
    @test S(sqrt(x) * sqrt(y), x ≫ 0, y ≫ 0) == sqrt(x * y) && S(sqrt(x) * sqrt(y)) == sqrt(x) * sqrt(y)
    @test S((x^3 - 1) / (x^2 + x + 1)) == x - 1 && S(x^2 - 1) == x^2 - 1    # not factored: no smaller
    @test simplify(x ≫ 0; assume = @assume(x > 1)) === true
    @test SymJ.simplify_powers(sqrt(Sym(2)) * sqrt(Sym(8))) == 4 && SymJ.simplify_rational(1 / x + 1 / y) == (x + y) / (x * y)
    @test SymJ._nodes(sqrt(x * y)) < SymJ._nodes(sqrt(x) * sqrt(y))
    @test simplify((x + y)^30 - (x + y)^30) == 0 && simplify((x + y)^30) == (x + y)^30   # no expansion tried
    # trigonometric: a normal form modulo sin^2 + cos^2 = 1, and linearization
    T = SymJ.simplify_trig
    @test T(sin(x)^2 + cos(x)^2) == 1 && T(1 - cos(x)^2) == sin(x)^2 && T(sec(x)^2 - tan(x)^2) == 1
    @test T(2sin(x) * cos(x)) == sin(2x) && T(cos(x)^2 - sin(x)^2) == cos(2x) && T(sin(2x) / sin(x)) == 2cos(x)
    @test T(tan(x) * cos(x)) == sin(x) && T((1 - cos(x)^2) / sin(x)) == sin(x) && T(sin(x) / cos(x)) == tan(x)
    @test T(1 / (1 + tan(x)^2)) == cos(x)^2 && T(sin(x)^4 - cos(x)^4) == -cos(2x)
    @test T(sin(x / 2) * cos(x / 2)) == sin(x) / 2 && T(cos(2x) + 2sin(x)^2) == 1
    @test T(sin(3x) - 3sin(x) + 4sin(x)^3) == 0 && T(x * sin(x)^2 + x * cos(x)^2) == x
    @test T(sin(x) * cos(y) + cos(x) * sin(y)) == sin(x + y) && T(cos(x) * cos(y) - sin(x) * sin(y)) == cos(x + y)
    @test T(cos(x) * cos(y) + sin(x) * sin(y)) == cos(x - y) && T(sin(x) * cos(y) - cos(x) * sin(y)) == sin(x - y)
    @test T(sin(x + y) - sin(x) * cos(y) - cos(x) * sin(y)) == 0
    @test T(sin(x) + 1) == sin(x) + 1 && T(x^2) == x^2                  # nothing to do
    @test simplify(sin(x)^2 + cos(x)^2 + x) == x + 1 && simplify((1 - cos(x)^2) / sin(x)^2) == 1
    # logarithms: combined or expanded where the arguments are positive
    @test S(log(Sym(2)) + log(Sym(3))) == log(Sym(6)) && S(2log(Sym(2)) - log(Sym(4))) == 0
    @test S(log(x) + log(y)) == log(x) + log(y) && S(log(x) + log(y), x ≫ 0, y ≫ 0) == log(x * y)
    @test S(log(x * y) - log(x), x ≫ 0) == log(y) && S(log(x^3) - 3log(x), x ≫ 0) == 0
    @test S(log(2x) - log(Sym(2))) == log(x)                                 # 2 > 0: always
    @test S(exp(log(x) + log(y)), x ≫ 0, y ≫ 0) == x * y
    @test SymJ.simplify_log(log(x * y); assume = Assumptions(x ≫ 0, y ≫ 0)) == log(x * y)  # the smaller
    # expansion: trigonometric and hyperbolic sums and multiples, logarithms
    X = SymJ.expand_trig
    @test X(sin(2x)) == 2sin(x) * cos(x) && X(cos(2x)) == 2cos(x)^2 - 1 && X(cos(3x)) == 4cos(x)^3 - 3cos(x)
    @test X(sin(x + y)) == sin(x) * cos(y) + cos(x) * sin(y) && X(cos(x - y)) == cos(x) * cos(y) + sin(x) * sin(y)
    @test X(sinh(x + y)) == sinh(x) * cosh(y) + cosh(x) * sinh(y) && X(cosh(2x)) == 2cosh(x)^2 - 1
    @test cancel(X(tan(2x)) - 2tan(x) / (1 - tan(x)^2)) == 0 && cancel(X(tanh(x + y)) - (tanh(x) + tanh(y)) / (1 + tanh(x) * tanh(y))) == 0
    @test X(sin(x + Sym(π) / 2)) == cos(x) && X(sin(x / 2)) == sin(x / 2) && X(x^2 + 1) == x^2 + 1
    @test X(sin(2sin(2x))) == 2sin(2sin(x) * cos(x)) * cos(2sin(x) * cos(x))        # outside in
    @test T(X(sin(3x))) == sin(3x) && T(X(cos(x + y + z))) == cos(x + y + z)
    for e in (sin(3x), cos(4x), sin(x + y + z), tan(3x), tan(x - y), cot(2x), cot(x + y), sec(2x), csc(x + y),
              sinh(3x), cosh(x - 2y), tanh(3x), coth(x + y), sech(2x), csch(x + y), sin(2x + 1) * exp(cos(2y)))
        @test isapprox(N(subs(e - X(e), x => 0.3, y => 0.7, z => -1.1)), 0; atol = 1e-12)
    end
    LX = SymJ.expand_log
    @test LX(log(x^2 * y); assume = Assumptions(x ≫ 0, y ≫ 0)) == 2log(x) + log(y)
    @test LX(log(x * y)) == log(x * y) && LX(log(x * y); assume = Assumptions(x ≫ 0)) == log(x) + log(y)
    @test LX(log(2x)) == log(Sym(2)) + log(x) && LX(log(x^2); assume = @assume(x in Reals)) == log(x^2)
    @test LX(log(x / y); force = true) == log(x) - log(y) && LX(log((x * y)^2); force = true) == 2log(x) + 2log(y)
    @test LX(log(sqrt(x) * y^3); force = true) == log(x) / 2 + 3log(y) && LX(log(x^y); force = true) == y * log(x)
    @test LX(exp(x) * log(x * y) + 1; force = true) == exp(x) * (log(x) + log(y)) + 1
    # factorials and binomials: integer shifts combined, ratios as binomials
    GS = SymJ.simplify_gamma
    @syms nn kk
    @test GS(factorial(nn + 2) / factorial(nn)) == (nn + 1) * (nn + 2) && GS(factorial(nn) / factorial(nn - 1)) == nn
    @test GS(factorial(nn) / factorial(nn + 3)) == 1 / ((nn + 1) * (nn + 2) * (nn + 3))
    @test GS(binomial(nn + 1, kk) / binomial(nn, kk)) == (nn + 1) / (nn + 1 - kk)
    @test GS(binomial(nn, kk + 1) / binomial(nn, kk)) == (nn - kk) / (kk + 1)
    @test GS(factorial(nn) / (factorial(kk) * factorial(nn - kk))) == binomial(nn, kk)
    @test GS(factorial(2nn) / factorial(nn)^2) == binomial(2nn, nn) && GS(factorial(nn) * (nn + 1)) == factorial(nn + 1)
    @test GS(factorial(nn + 1) - nn * factorial(nn)) == factorial(nn) && GS(factorial(nn + 1) + factorial(nn)) == (nn + 2) * factorial(nn)
    @test GS((nn + 1) * binomial(nn, kk) / (nn - kk + 1)) == binomial(nn + 1, kk)
    @test GS(sin(factorial(nn + 1) / factorial(nn))) == sin(nn + 1) && GS(x^2 + 1) == x^2 + 1
    @test simplify(factorial(nn + 2) / factorial(nn)) == (nn + 1) * (nn + 2) && simplify(factorial(nn) * (nn + 1)) == factorial(nn + 1)
    for e in (binomial(2nn + 2, nn + 1) / binomial(2nn, nn), factorial(nn + 1) / factorial(nn - 2) + binomial(nn, kk) / binomial(nn - 1, kk))
        @test N(subs(e - GS(e), nn => 9, kk => 4)) == 0
    end
    # radicals: denested, denominators rationalized, radicands factored
    @test S(sqrt(3 + 2sqrt(Sym(2)))) == 1 + sqrt(Sym(2)) && S(sqrt(5 - 2sqrt(Sym(6)))) == sqrt(Sym(3)) - sqrt(Sym(2))
    @test S(1 / (1 + sqrt(Sym(2)))) == sqrt(Sym(2)) - 1 && S(1 / (sqrt(Sym(3)) - 1)) == (sqrt(Sym(3)) + 1) / 2
    @test S(2 / (1 + sqrt(Sym(5)))^2) == (3 - sqrt(Sym(5))) / 4
    @test S(sqrt(x^2 + 2x + 1), x ∈ ℝ) == abs(x + 1) && S(sqrt(x^2 + 2x + 1), x ≫ -1) == x + 1
    @test S(sqrt(4x^2 + 4x + 1) / (2x + 1), x ≫ 0) == 1
    @test SymJ.simplify_radicals(sqrt(2 + sqrt(Sym(3)))) == sqrt(Sym(6)) / 2 + sqrt(Sym(2)) / 2
    @test S(sqrt(2 + sqrt(Sym(3)))) == sqrt(2 + sqrt(Sym(3)))               # denested is larger
    # hyperbolic, and exponentials with trigonometric and hyperbolic functions
    H = SymJ.simplify_hyperbolic
    @test H(cosh(x)^2 - sinh(x)^2) == 1 && H(2sinh(x) * cosh(x)) == sinh(2x) && H(cosh(x)^2 + sinh(x)^2) == cosh(2x)
    @test H(cosh(x) + sinh(x)) == exp(x) && H(cosh(x) - sinh(x)) == exp(-x) && H(exp(x) - cosh(x)) == sinh(x)
    @test H(exp(x) + exp(-x)) == 2cosh(x) && H((exp(x) - exp(-x)) / 2) == sinh(x) && H(tanh(x) * cosh(x)) == sinh(x)
    @test H(1 - tanh(x)^2) == sech(x)^2 && H(sinh(2x) / sinh(x)) == 2cosh(x) && H(coth(x)^2 - 1) == csch(x)^2
    # reciprocal powers as sec, csc, sech, csch; tan and cot from what pairs
    @test T(1 + tan(x)^2) == sec(x)^2 && S(1 / sinh(x)^2) == csch(x)^2 && S(cosh(x) / sinh(x)) == coth(x)
    @test S(sin(x) / cos(x)^2) == tan(x) * sec(x) && S(1 / (sin(x) * cos(x))) == csc(x) * sec(x)
    @test H(sinh(x) * cosh(y) + cosh(x) * sinh(y)) == sinh(x + y)
    @test H((exp(2x) - 1) / (exp(x) - exp(-x))) == exp(x) && H(exp(x)) == exp(x) && H(sinh(x) + 1) == sinh(x) + 1
    @test T((exp(IM * x) + exp(-IM * x)) / 2) == cos(x) && T((exp(IM * x) - exp(-IM * x)) / (2IM)) == sin(x)
    @test T(cos(x) + IM * sin(x)) == exp(IM * x) && T(cos(x) - IM * sin(x)) == exp(-IM * x)
    @test T(exp(2IM * x) - 2IM * sin(x) * cos(x)) == cos(2x)
    @test simplify(cosh(x)^2 - sinh(x)^2 + x) == x + 1 && simplify(exp(x) + exp(-x)) == 2cosh(x)
    # every simplification agrees with its input numerically
    for e in (sin(x)^2 + cos(x)^2, 1 - cos(x)^2, 2sin(x) * cos(x), sin(2x) / sin(x), sin(x)^4 - cos(x)^4,
              sin(x)^2 * cos(x)^2, cos(x) * cos(y) + sin(x) * sin(y), tan(x) + tan(y), 1 / (1 + tan(x)^2),
              cosh(x) - sinh(x), cosh(x)^2 + sinh(x)^2, exp(x) - cosh(x), 1 - tanh(x)^2, exp(2x) - 1,
              cos(x) - IM * sin(x), exp(IM * x) - cos(x), log(2x) - log(Sym(2)), (x^3 - 1) / (x - 1),
              x / (x + 1) + 1 / (x + 1), (x + 1)^2 - x^2, sqrt(Sym(2)) * sqrt(Sym(6)), sin(x / 2) * cos(x / 2))
        s = simplify(e)
        for (vx, vy) in ((0.37, 1.13), (1.9, -0.6))
            @test isapprox(N(subs(s, x => vx, y => vy)), N(subs(e, x => vx, y => vy)); atol = 1e-12)
        end
    end
end

@testset "simplify with rules" begin
    @syms t a b
    @test simplify((1 - cos(2t)) / 2) == sin(t)^2 && simplify((1 + cos(2t)) / 2) == cos(t)^2
    @test simplify(3 - 3cos(2t) + a) == a + 6sin(t)^2 && simplify((cosh(2t) - 1) / 2) == sinh(t)^2
    @test simplify(sin(t) / (1 + cos(t))) == tan(t / 2) && simplify((1 - cos(t)) / sin(t)) == tan(t / 2)
    @test simplify(a * sin(t) / (1 + cos(t))) == a * tan(t / 2)
    @test simplify(1 - cos(t)) == 1 - cos(t)                 # not 2sin(t/2)^2, which is larger
    @test simplify(sin(t + a) - cos(a) * sin(t)) == cos(t) * sin(a)
    # exp(u) and radicands are simplified as function arguments are
    @test simplify(exp(sin(t)^2 + cos(t)^2)) == Sym(ℯ) && simplify(sqrt(1 - sin(t)^2)) == sqrt(cos(t)^2)
    @test simplify(t * exp(t^2 + 2t + 1 - (t + 1)^2)) == t
    @test SymJ._halvable(2t) && SymJ._halvable(2t + 4a) && !SymJ._halvable(t) && !SymJ._halvable(2t + 1)
end

@testset "pattern matching" begin
    P = SymJ.pattern
    M(p, ex) = match(P(p), ex)
    Ms(p, ex) = eachmatch(P(p), ex)
    D(ps...) = Dict{Symbol,Any}(ps...)
    @syms t a b c f()
    # literals, function arguments, powers, repeated slots
    @test M(:(~x), sin(x)) == D(:x => sin(x)) && M(:(sin(~x)), sin(x + 1)) == D(:x => x + 1)
    @test M(:(sin(~x)), cos(x)) === nothing && M(:((~x)^2), x^3) === nothing
    @test M(:((~x)^(~x)), (x + a)^(x + a)) == D(:x => x + a) && M(:((~x)^(~x)), (x + a)^(x + b)) === nothing
    @test M(:(sin(~x) * cos(~x)), sin(t) * cos(t)) == D(:x => t)
    @test M(:(sin(~x)^2 + cos(~x)^2), sin(t)^2 + cos(t)^2) == D(:x => t)
    @test M(:(f(~x, ~~r)), f(a, b, c)) == D(:x => a, :r => Sym[b, c]) && M(:(f(~~~r)), f()) === nothing
    @test M(:((~u)^(~m)), 1 / sqrt(x)) == D(:u => x, :m => -1 // 2) && M(:(exp(~u)), exp(2x)) == D(:u => 2x)
    @test M(:((~u)^(~m) * (~v)^(~n)), a^2 / b^3) in Ms(:((~u)^(~m) * (~v)^(~n)), a^2 / b^3)
    # sums and products: slots take groups of terms (factors)
    @test length(Ms(:(~x + ~y), a + b + c)) == 6 && length(Ms(:(~x * ~y), a * b * c)) == 6
    @test D(:x => a + b, :y => c) in Ms(:(~x + ~y), a + b + c)
    @test M(:(~x + sin(~y)), a + b + sin(c)) == D(:x => a + b, :y => c)
    @test M(:(~x * (~y)^2), a * b * c^2) == D(:x => a * b, :y => c)
    @test M(:(sin(~x)^2 + cos(~x)^2 + ~~r), sin(t)^2 + cos(t)^2 + a) == D(:x => t, :r => a)
    @test M(:(sin(~x)^2 + cos(~x)^2 + ~~r), sin(t)^2 + cos(t)^2) == D(:x => t, :r => 0)
    @test M(:(sin(~x)^2 + cos(~x)^2), sin(t)^2 + cos(t)^2 + a) === nothing
    @test M(:(~x + ~~~y), a) === nothing && length(Ms(:(~x + ~~y), a + b + c)) == 7
    @test M(:(~x * ~y), a + b) === nothing
    # defaults, and which match comes first: numbers to ~! slots, else none
    @test M(:(~!a * ~x + ~!b), 3y + 5) == D(:a => 3, :x => y, :b => 5)
    @test M(:(~!a * ~x + ~!b), sin(t)) == D(:a => 1, :x => sin(t), :b => 0)
    @test M(:(~!a * (~x)^(~!n)), 2y^3) == D(:a => 2, :x => y, :n => 3)
    @test M(:(~!a * (~x)^(~!n)), y) == D(:a => 1, :x => y, :n => 1)
    @test M(:(~x + ~!r), a + b) == D(:x => a + b, :r => 0)
    @test M(:(~a * ~x + ~b), 2x + 3) == D(:a => 2, :x => x, :b => 3)
    @test M(:((~a + ~b * ~x)^(~m)), (2 + 3x)^5) == D(:a => 2, :b => 3, :x => x, :m => 5)
    @test M(:((~a + ~b * ~x)^(~m) * (~c + ~d * ~x)^(~n)), (1 + 2x)^3 * (4 + 5x)^6) in
          (D(:a => 1, :b => 2, :m => 3, :c => 4, :d => 5, :n => 6, :x => x), D(:a => 4, :b => 5, :m => 6, :c => 1, :d => 2, :n => 3, :x => x))
    @test M(:(~x * (~!a + ~!b * ~x)^(~!m)), x * (2 + 3x)^4) == D(:x => x, :a => 2, :b => 3, :m => 4)
    # numbers divide into slots
    @test M(:(~x - ~y), a - 2b) == D(:x => a, :y => 2b) && M(:(2 * ~x), 6a) == D(:x => 3a)
    @test D(:x => a, :y => -b) in Ms(:(~x - ~y), a + b)
    @test M(:(-sin(~x)), -sin(a)) == D(:x => a) && M(:(-sin(~x)), sin(a)) === nothing
    @test M(:(~x + 1), a + 3) === nothing && M(:(~x + 1), a + b + 1) == D(:x => a + b)
    # the pattern is its canonical form
    @test M(:(~x + ~x), 2a) == D(:x => a) && M(:(sin(-~x)), sin(-a)) == D(:x => a)
    @test M(:(exp(log(~x))), sin(t)) == D(:x => sin(t))
    @test_throws ArgumentError P(:(~x - ~x))
    @test_throws ArgumentError P(:(~x^2))
    @test_throws ArgumentError P(:(~x + ~~x))
    # ~n and ~!n may name one slot; a slot bound to 1 lets b^(~n) match b
    @test M(:((~!b * t^(~!n)) * (~c + ~!d * t^(~n))), t * (1 + 3t)) == D(:b => 1, :n => 1, :d => 3, :c => 1)
    @test M(:((~!b * t^(~!n)) * (~c + ~!d * t^(~n))), t^2 * (1 + 3t^2)) == D(:b => 1, :n => 2, :d => 3, :c => 1)
    @test M(:((~!b * t^(~!n)) * (~c + ~!d * t^(~n))), t^2 * (1 + 3t^3)) === nothing
    @test M(:(f(~!n, t^(~n))), f(1, t)) == D(:n => 1) && M(:(f(~n, t^(~n))), f(2, t)) === nothing
    # head slots and predicates
    @test M(:((~f)(~x)), sin(t)) == D(:f => sin, :x => t) && M(:((~f)(~x)), y) === nothing
    @test M(:((~f)(~x)^2), cos(x)^2) == D(:f => cos, :x => x)
    @test M(:((~f)(~x) + (~f)(~y)), cos(b) + sin(a)) === nothing
    @test M(:((~g)(~x)), f(t)) == D(:g => f, :x => t)
    @test M(:((~x)^(~n::iseven)), x^2) == D(:x => x, :n => 2) && M(:((~x)^(~n::iseven)), x^3) === nothing
    @test M(:(~a::isinteger * ~x), 3x) == D(:a => 3, :x => x) && M(:((~x)^(~n::isodd) + ~y), x^3 + a) == D(:x => x, :n => 3, :y => a)
    @test isinteger(Sym(3)) && isinteger(Sym(2.0)) && !isinteger(Sym(3 // 2)) && !isinteger(x) && !isinteger(Sym(2.5))
    @test iseven(Sym(4)) && !iseven(Sym(3)) && isodd(Sym(-3)) && !isodd(x) && iseven(Sym(big(2)^70)) && isodd(Sym(5.0))
    isnum(v) = v isa Sym && v.kind == SymJ.NUMBER
    @test match(P(:(~n::$isnum * ~x)), 3a) == D(:n => 3, :x => a)
    @test match(P(:(~n::$isnum * sin(~x))), a * sin(t)) === nothing
    # occursin
    @test occursin(P(:(sin(~x)^2)), 1 + sin(t)^2) && !occursin(P(:(cos(~x))), 1 + sin(t)^2)
    @test occursin(P(:(exp(~u))), x * exp(2x)) && string(P(:(~!a * sin(~x) + ~~r))) == "pattern(sin(~x)*~!a + ~~r)"
end

@testset "rewrite rules" begin
    @syms t a b c f() h()
    R = SymJ.rewrite
    r1 = SymJ.@rule sin(~x)^2 + cos(~x)^2 + ~~r => 1 + ~r
    @test r1(sin(t)^2 + cos(t)^2 + a) == 1 + a && r1(sin(t)^2 + cos(t)^2) == 1 && r1(sin(t)) === nothing
    @test string(r1) == "rule(cos(~x)^2 + sin(~x)^2 + ~~r => 1 + ~r)"
    # predicates and values from the calling scope; a right side of nothing rejects the match
    issym(v) = v isa Sym && v.kind == SymJ.SYMBOL
    @test (SymJ.@rule ~x::issym * sin(~y) => ~x)(a * sin(t)) == a
    @test (SymJ.@rule ~x * ~y => (isinteger(~x) ? nothing : ~x))(2a) == a
    c0 = Sym(3)
    @test match(SymJ.@pattern(~x + $c0), a + 3) == Dict{Symbol,Any}(:x => a) && match(SymJ.@pattern(~x + $c0), a + 2) === nothing
    @test (SymJ.@rule (~x)^(~n::isinteger) => ~n * (~x)^(~n - 1))(x^3) == 3x^2
    @test (SymJ.@rule (~x)^(~n::isinteger) => ~n * (~x)^(~n - 1))(x^y) === nothing
    @test (SymJ.@rule (~f)(-~x) => (~f)(~x))(h(-t)) == h(t)
    @test_throws ArgumentError SymJ.@rule(sin(~x) => cos(~y))
    # rules at run time: a template or a function of the bindings
    @test SymJ.rule(:(sin(~x)^2) => :(1 - cos(~x)^2))(sin(t)^2) == 1 - cos(t)^2
    @test SymJ.rule(:((~x)^(~n)) => d -> d[:n] * d[:x]^(d[:n] - 1))(x^4) == 4x^3
    @test SymJ.rule(:(f(~x, ~~r)) => :(f(~~r, ~x)))(f(a, b, c)) == f(b, c, a)
    @test SymJ.rule(:((~f)(~x)^2) => :((~f)(2 * ~x)))(cos(t)^2) == cos(2t)
    @test_throws ArgumentError SymJ.rule(:(sin(~x)) => :(cos(~y)))
    # rewriting every part, to a fixpoint or once
    r6 = SymJ.@rule sin(~x)^2 => 1 - cos(~x)^2
    @test R(2sin(t)^2 + cos(t)^2, r6) == 2 - cos(t)^2 && R(sin(sin(t)^2), r6) == sin(1 - cos(t)^2)
    r7 = SymJ.@rule log(~x * ~y) => log(~x) + log(~y)
    @test R(log(a * b * c), r7) == log(a) + log(b) + log(c) && R(log(a * b * c), r7; order = :pre) == log(a) + log(b) + log(c)
    @test length(replace(log(a * b * c), r7).args) == 2 && replace(sin(t), r7) == sin(t)
    @test R(x * sin(t)^2, [r6, r7]) == x * (1 - cos(t)^2)
    @test (@test_logs (:warn,) R(a, SymJ.@rule(~x::issym => ~x + 1); maxsteps = 3)) == a + 3
    @test_throws ArgumentError R(a, r6; order = :up)
    # indexed rule sets keep the rules' order
    rh = SymJ.@rule (~f)(~x) => (~f === sin ? cos(~x) : nothing)
    rs = SymJ.@rule sin(~x) => tan(~x)
    @test replace(sin(t), rh, rs) == cos(t) && replace(sin(t), rs, rh) == tan(t)
    @test replace(sin(t) + h(t), rs, rh) == tan(t) + h(t) && R(sin(t)^3, SymJ.@rule((~x)^3 => ~x)) == sin(t)
    # a bound slot whose value is 1 takes no factors in a product
    @test match(SymJ.pattern(:(h(~c, ~c * a))), h(1, a)) == Dict(:c => 1)
    @test match(SymJ.pattern(:(~!c - ~!c * a)), 1 - a) == Dict(:c => 1)
end

@testset "rewrite to a target" begin
    @syms a b n k
    R = SymJ.rewrite
    same(e, r, v = 0.37) = abs(N(subs(e - r, x => v, y => 0.81))) < 1e-12
    for f in (sin, cos, tan, cot, sec, csc, sinh, cosh, tanh, coth, sech, csch)
        @test same(f(x), R(f(x), exp))
        @test same(f(x), R(f(x), f in (sinh, cosh, tanh, coth, sech, csch) ? tanh : tan))
        @test !SymJ._hasfn(R(f(x), exp), nameof(f))
    end
    for (f, v) in ((asin, 0.37), (acos, 0.37), (atan, 0.37), (acot, 0.37), (asec, 1.7), (acsc, 1.7),
                   (asinh, 0.37), (acosh, 1.7), (atanh, 0.37), (acoth, 1.7), (asech, 0.6), (acsch, 0.37))
        @test same(f(x), R(f(x), log), v)
    end
    @test R(sin(x), exp) == (exp(im * x) - exp(-im * x)) / (2im)
    @test R(exp(im * x), cos) == cos(x) + im * sin(x)
    @test R(exp(2 + 3im * x + y), sin) == exp(y + 2) * (cos(3x) + im * sin(3x))
    @test R(exp(x), cos) == exp(x)                          # no imaginary part
    @test R(exp(x), cosh) == cosh(x) + sinh(x)
    @test R(tan(x)^2 + sec(x), sin) == sin(x)^2 / cos(x)^2 + 1 / cos(x)
    @test R(tan(x)^2 + sec(x), cos) == R(tan(x)^2 + sec(x), sin)
    @test R(tanh(x), cosh) == sinh(x) / cosh(x)
    @test R(sin(x), tan) == 2tan(x / 2) / (1 + tan(x / 2)^2)
    @test R(asinh(x), log) == log(x + sqrt(x^2 + 1))
    @test R(sin(x) + cos(x), sin => exp) == cos(x) + (exp(im * x) - exp(-im * x)) / (2im)
    @test R(sin(x) + cosh(x) + cos(y), (sin, cosh) => exp) == R(sin(x) + cosh(x), exp) + cos(y)
    @test R(binomial(n, k), factorial) == factorial(n) / (factorial(k) * factorial(n - k))
    @test SymJ.simplify(R(sin(x)^2 + cos(x)^2, exp)) == 1
    @test R(x^2 + y, exp) == x^2 + y                         # nothing to rewrite
    @test_throws ArgumentError R(sin(x), y)
    @test_throws ArgumentError R(sin(x), sqrt)
    @test_throws ArgumentError R(sin(x), x => exp)
end

@testset "replace" begin
    @syms t a u()
    isnum = SymJ.isnumber
    @test replace(sin(x) + sin(y)^2, :(sin(~u)) => :(cos(~u))) == cos(x) + cos(y)^2
    @test replace(exp(x) + exp(2y), :(exp(~u)) => :(1 + ~u)) == x + 2y + 2
    @test replace(exp(x) * sin(x), :(exp(~u)) => :((~u)^2)) == x^2 * sin(x)    # a factor of a product
    @test replace(sin(sin(x)), :(sin(~u)) => :(2 * ~u)) == 4x                    # from the leaves up
    @test replace(u(x) + u(2x), :(u(~s)) => :((~s)^2)) == 5x^2
    @test replace(sin(x) + cos(y), :(sin(~u)) => :(cos(~u)), :(cos(~u)) => :(sin(~u))) == cos(x) + sin(y)
    @test replace(sin(x), SymJ.pattern(:(sin(~u))) => d -> d[:u]^2) == x^2
    # predicates: numeric coefficients and constants are parts too
    @test replace(x^2 + 3x + 5, (e -> isnum(e) && e > 2) => (e -> e - 1)) == x^2 + 2x + 4
    @test replace(x - y, isnum => (e -> 5e)) == x - 5y
    @test replace(e -> isnum(e) ? 2e : e, x^2 + 3x + 1) == x^4 + 6x + 2
    @test replace(x + y, x => a) == a + y && replace(x + 2, 2 => a) == a + x
    r = SymJ.@rule log(~x * ~y) => log(~x) + log(~y)
    @test replace(log(a * x) + sin(x), r, :(sin(~u)) => :(cos(~u))) == log(a) + log(x) + cos(x)
    let e0 = x^2 * sin(y) + 3
        @test replace(e0, a => x) === e0                                         # unchanged: the same node
    end
    @test_throws ArgumentError replace(sin(x), sin => cos)                       # use a pattern
end
