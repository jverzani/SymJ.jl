# Elementary and other functions (included by runtests.jl)

@testset "functions" begin
    # poles are errors, as 1/0 is (tan(π/2) left as a value made a wrong integral; SymPy #27185)
    @test_throws SymJ.DivisionByZero tan(PI / 2)
    @test_throws SymJ.DivisionByZero sec(PI / 2)
    @test_throws SymJ.DivisionByZero cot(Sym(0))
    @test_throws SymJ.DivisionByZero csc(PI)
    @test tan(PI / 4) == 1
    @test sin(Sym(0)) == 0
    @test sin(PI) == 0
    @test cos(PI) == -1
    @test sin(PI / 6) == Sym(1 // 2)
    @test sin(PI / 3) == sqrt(Sym(3)) / 2
    @test cos(PI / 4) == sqrt(Sym(2)) / 2
    @test sin(7PI / 6) == Sym(-1 // 2)
    @test tan(PI / 4) == 1
    @test sin(-x) == -sin(x)
    @test cos(-x) == cos(x)
    @test sin(y - x) == -sin(x - y)
    # the minus comes out when the term printed first is negative
    @test string(sin(x - y)) == "sin(x - y)" && string(sin(y - x)) == "-sin(x - y)" && string(cos(y - x)) == "cos(x - y)"
    @test string(tan(1 - x)) == "-tan(x - 1)" && string(sinh(2 - x^2)) == "-sinh(x^2 - 2)" && string(sin(x - y + 1)) == "sin(x - y + 1)"
    @test string(simplify(cos(y) * sin(x) - cos(x) * sin(y))) == "sin(x - y)"
    @test log(Sym(1)) == 0
    @test log(E) == 1
    @test log(exp(Sym(3))) == 3
    # the principal value at a negative number or ±im (SymPy's log(-1) = iπ)
    @test log(Sym(-1)) == IM * PI
    @test log(Sym(-2)) == log(Sym(2)) + IM * PI
    @test log(IM) == IM * PI / 2
    @test log(-IM) == -IM * PI / 2
    @test kind(log(Sym(2), x)) == FN && length(log(Sym(2), x).args) == 2
    @test log(x, x) == 1
    @test log(x, 1) == 0
    @test kind(sin(Sym(2))) == FN          # exact argument stays symbolic
    @test sin(Sym(2.0)) == Sym(sin(2.0))   # float argument evaluates
    @test exp(Sym(0.5)) == Sym(exp(0.5))   # a constant to a float power too
    @test PI^1.5 == Sym(Float64(π)^1.5)
    @test kind(exp(x)^1.0) == POW            # (with a symbol, not)
    @test kind((-PI)^0.5) == POW
    # floats evaluate before symbolic rules: Julia's principal value on a branch cut
    @test N(asin(Sym(-17.0))) ≈ asin(complex(-17.0))
    @test N(acosh(Sym(0.5))) ≈ acosh(complex(0.5))
    @test f(x) == f(x)
    @test f(x) != g(x)
    @test f(x) != f(y)
    @test f(x) + f(x) == 2f(x)
    @test string(f(x, 2)) == "f(x, 2)"
end

@testset "more elementary functions" begin
    # reciprocal trig: exact values through cos/sin/tan, poles unevaluated, parity
    @test sec(Sym(0)) == 1 && sec(PI / 3) == 2 && csc(PI / 6) == 2
    @test cot(PI / 4) == 1 && cot(PI / 2) == 0
    @test_throws SymJ.DivisionByZero sec(PI / 2)              # a pole: an error, as 1/0
    @test_throws SymJ.DivisionByZero csc(Sym(0))
    @test sec(-x) == sec(x) && csc(-x) == -csc(x) && cot(-x) == -cot(x)
    @test sec(acos(x)) == 1 / x
    @test kind(sec(x)) == FN && operation(sec(x)) === sec
    @test asec(Sym(1)) == 0 && acsc(Sym(1)) == PI / 2
    @test acot(Sym(0)) == PI / 2 && acot(Sym(1)) == PI / 4
    # odd only away from 0 (acot(0) = π/2): the sign comes out of a nonzero number only
    @test acot(-Sym(3)) == -acot(Sym(3))
    @test string(acot(-x)) == "acot(-x)"
    @test subs(acot(-x), x => 0) == PI / 2
    @test sech(Sym(0)) == 1 && coth(-x) == -coth(x) && csch(-x) == -csch(x)
    @test asech(Sym(1)) == 0
    @test acoth(-Sym(3)) == -acoth(Sym(3))
    @test string(acoth(-x)) == "acoth(-x)"
    # real cube root
    @test cbrt(Sym(-8)) == -2 && cbrt(Sym(27)) == 3
    @test cbrt(Sym(2)) == Sym(2)^(1 // 3) && cbrt(Sym(-16)) == -2 * Sym(2)^(1 // 3)
    @test cbrt(x)(-8.0) == -2.0 && cbrt(-x) == -cbrt(x)
    # piecewise constant
    @test sign(Sym(-3)) == -1 && sign(sqrt(Sym(2)) - 1) == 1 && sign(-x) == -sign(x)
    @test floor(Sym(7 // 2)) == 3 && floor(Sym(-7 // 2)) == -4 && ceil(Sym(π)) == 4
    @test trunc(Sym(-7 // 2)) == -3 && round(Sym(5 // 2)) == 2 && round(Sym(7 // 2)) == 4
    @test round(sqrt(Sym(2))) == 1 && kind(floor(x)) == FN && floor(floor(x)) == floor(x)
    # floor(-v) = -ceil(v), ceil(-v) = -floor(v); trunc and round are odd (SymPy #30300)
    @test floor(x) + ceil(-x) == 0
    @test floor(-x) == -ceil(x)
    @test trunc(-x) == -trunc(x)
    @test round(-x) == -round(x)
    @test round(x)(2.5) == 2.0
    # rewrites
    @test log2(x) == log(Sym(2), x) && log10(x) == log(Sym(10), x)
    @test log1p(x) == log(1 + x) && expm1(x) == exp(x) - 1
    @test exp2(x) == Sym(2)^x && exp10(x) == Sym(10)^x
    @test sinpi(Sym(1 // 6)) == Sym(1 // 2) && cospi(Sym(1)) == -1 && tanpi(Sym(1 // 4)) == 1
    @test sincos(x) == (sin(x), cos(x))
    @test sind(Sym(30)) == Sym(1 // 2) && cosd(Sym(60)) == Sym(1 // 2) && asind(Sym(1)) == 90
    @test sincosd(Sym(30)) == (Sym(1 // 2), sqrt(Sym(3)) / 2)
    @test deg2rad(Sym(180)) == PI && rad2deg(PI) == 180
    @test hypot(Sym(3), 4) == 5 && hypot(x, y) == sqrt(x^2 + y^2)
    @test abs2(x) == abs(x)^2 && fourthroot(Sym(16)) == 2
    @test (2 \ x) == x / 2
    @test numerator(x / y + 1) == x + y && denominator(x / y + 1) == y
    @test muladd(x, 2, 1) == 2x + 1 && fma(2, x, y) == 2x + y
    @test expand(evalpoly(x, (1, 2, 3))) == 3x^2 + 2x + 1
    # sinc and cosc
    @test sinc(Sym(0)) == 1 && sinc(Sym(1 // 2)) == 2 / PI && sinc(Sym(3)) == 0
    @test cosc(Sym(0)) == 0 && sinc(-x) == sinc(x) && cosc(-x) == -cosc(x)
    # mod, rem, div, fld, cld
    @test mod(Sym(-7), 3) == 2 && rem(Sym(-7), 3) == -1 && mod(Sym(7 // 2), 1) == Sym(1 // 2)
    @test mod2pi(Sym(7)) == 7 - 2PI && kind(mod(x, 3)) == FN && mod(x, 3)(7.5) == 1.5
    @test div(Sym(7), 2) == 3 && fld(Sym(-7), 2) == -4 && cld(Sym(7), 2) == 4
    @test divrem(Sym(7), Sym(2)) == (3, 1)
    # max and min
    @test max(x, 1, 3, x) == max(x, 3) && max(y, x) == max(x, y)
    @test max(max(x, y), 2) == max(x, y, 2) && min(Sym(2), sqrt(Sym(3))) == sqrt(Sym(3))
    @test max(Sym(1), 2) == 2 && minmax(Sym(3), 2) == (2, 3) && max(x, 2)(5) == 5
    @test clamp(x, 0, 1) == min(max(x, 0), 1) && clamp(x, 0, 1)(2.5) == 1
    # factorial and binomial
    @test factorial(Sym(20)) == factorial(20) && factorial(Sym(30)) == Sym(factorial(big(30)))
    @test kind(factorial(x)) == FN && factorial(x)(5) == 120
    @test binomial(x, 3) == x * (x - 1) * (x - 2) / 6 && binomial(Sym(10), 3) == 120
    @test kind(binomial(x, y)) == FN && binomial(x, -1) == 0

    # numeric values agree with Base, through evaluation and N
    for (f, vals) in ((sec, (0.4, -1.2)), (csc, (0.4, -1.2)), (cot, (0.4, 2.0)), (asec, (1.5, -3.0)),
                      (acsc, (1.5, -3.0)), (acot, (0.4, -2.0)), (sech, (0.4, -1.2)), (csch, (0.4, -1.2)),
                      (coth, (0.4, -1.2)), (asech, (0.3, 0.9)), (acsch, (0.4, -2.0)), (acoth, (1.5, -3.0)),
                      (cbrt, (5.0, -5.0)), (sinc, (0.3, -1.7)), (cosc, (0.3, -1.7)), (sind, (30.5, -100.0)),
                      (cosd, (30.5, -100.0)), (tand, (30.5, -100.0)), (secd, (30.5,)), (cscd, (30.5,)),
                      (cotd, (30.5,)), (asind, (0.3,)), (acosd, (0.3,)), (atand, (2.0,)), (asecd, (2.0,)),
                      (acscd, (2.0,)), (acotd, (2.0,)), (tanpi, (0.1,)), (log2, (3.0,)), (log10, (3.0,)),
                      (log1p, (0.3,)), (exp2, (1.3,)), (exp10, (0.3,)), (expm1, (0.3,)), (deg2rad, (45.0,)),
                      (rad2deg, (1.0,)), (fourthroot, (5.0,)), (abs2, (-3.0,)), (sign, (-2.5,)),
                      (floor, (2.5, -2.5)), (ceil, (2.5,)), (trunc, (-2.5,)), (round, (3.5, 2.5)))
        for v in vals
            @test f(x)(v) ≈ f(v)
            @test N(f(Sym(v))) ≈ f(v)
        end
    end

    # derivative rules agree with central differences
    for (f, v) in ((sec, 0.4), (csc, 0.4), (cot, 0.4), (asec, 1.5), (acsc, 1.5), (acot, 0.4), (sech, 0.4),
                   (csch, 0.4), (coth, 0.4), (asech, 0.3), (acsch, 0.4), (acoth, 1.5), (cbrt, 5.0),
                   (sinc, 0.3), (cosc, 0.3), (sind, 30.5), (asind, 0.3), (log2, 3.0), (log1p, 0.3), (expm1, 0.3))
        h = big(1) / 10^20
        fd = (f(x)(big(v) + h) - f(x)(big(v) - h)) / (2h)
        @test Float64(diff(f(x), x)(big(v))) ≈ Float64(fd) rtol = 1e-12
    end
end

@testset "atan(y, x), predicates, rounding and conversions" begin
    @syms x y
    # atan(y, x): exact for constants of known sign, in (-π, π]
    @test atan(Sym(1), Sym(-1)) == 3PI / 4
    @test atan(Sym(-1), -1) == -3PI / 4
    @test atan(Sym(0), -2) == PI
    @test atan(Sym(0), Sym(0)) == 0
    @test atan(-3, Sym(0)) == -PI / 2
    @test atan(sqrt(Sym(3)), 1) == PI / 3
    @test atand(1, Sym(-1)) == 135
    @test string(atan(y, x)) == "atan(y, x)"
    @test kind(atan(y, x)) == FN
    @test N(atan(1.0, Sym(-1))) ≈ atan(1.0, -1.0)
    for (a, b) in ((1.0, -1.0), (-0.5, -2.0), (2.0, 0.3), (-1.0, 0.0))
        @test atan(y, x)(y => a, x => b) == atan(Sym(a), Sym(b))
        @test lambdify(atan(y, x), [y, x])(a, b) ≈ atan(a, b)
    end
    @test diff(atan(y, x), x) == -y / (x^2 + y^2)
    @test diff(atan(y, x), y) == x / (x^2 + y^2)
    # placed in its quadrant under assumptions
    @test refine(atan(y, x); assume = @assume(x > 0, y in Reals)) == atan(y / x)
    @test refine(atan(y, x); assume = @assume(x < 0, y >= 0)) == atan(y / x) + PI
    @test refine(atan(y, x); assume = @assume(x < 0, y < 0)) == atan(y / x) - PI
    @test refine(atan(y, x); assume = @assume(x ~ 0, y > 0)) == PI / 2
    @test refine(atan(y, x); assume = @assume(x < 0)) == atan(y, x)        # y may be complex

    # predicates: true when known
    @test !isreal(x)
    @test isreal(PI + 1)
    @test !isreal(IM)
    @test (@assuming x in Reals begin isreal(x^2 + 1) end)
    @test !isfinite(x)
    @test isfinite(sqrt(Sym(2)))
    @test !isfinite(Sym(Inf))
    @test (@assuming x > 0 begin isfinite(x) end)
    @test !isfinite(x + Sym(Inf))
    @test isinf(Sym(-Inf))
    @test !isinf(x)
    @test !isinf(PI)
    @test !isnan(x)
    @test !isnan(Sym(1))
    @test modf(Sym(7) / 2) == (Sym(1) / 2, 3)

    # rounding: modes, digits, and to an integer type
    h = Sym(7) / 2
    @test floor(Int, h) === 3
    @test ceil(Int, h) === 4
    @test trunc(Int, -h) === -3
    @test round(Int, h) === 4
    @test round(Int, Sym(5) / 2) === 2                        # ties to even, as Base
    @test round(Sym(5) / 2, RoundNearestTiesAway) == 3
    @test round(Sym(-5) / 2, RoundNearestTiesAway) == -3
    @test round(Sym(5) / 2, RoundNearestTiesUp) == 3
    @test round(h, RoundDown) == 3
    @test round(PI; digits = 2) == Sym(157) / 50
    @test round(PI; sigdigits = 3) == Sym(157) / 50
    @test round(Sym(1234); sigdigits = 2) == 1200
    @test floor(PI; digits = 3) == Sym(3141) / 1000
    @test floor(x; digits = 1) == floor(10x) / 10
    @test kind(round(x)) == FN
    @test_throws ArgumentError round(x; sigdigits = 2)

    # conversions
    @test Int(Sym(3)) === 3
    @test Integer(Sym(3)) === 3
    @test Integer(Sym(2)^70) == big(2)^70
    @test Int(Sym(3.0)) === 3
    @test_throws InexactError Int(h)
    @test_throws InexactError Int(x)
    @test_throws InexactError Int(sqrt(Sym(2)))
    @test Rational(h) === 7 // 2
    @test Rational{BigInt}(h) == big(7) // 2
    @test Float32(h) === 3.5f0
    @test BigFloat(PI) == big(π)
    @test_throws InexactError BigFloat(IM)

    # hypot of one and of three or more arguments
    @test hypot(x) == abs(x)
    @test hypot(Sym(2), 3, 6) == 7
    @test hypot(1, 2, x) == sqrt(x^2 + 5)
    @test hypot(x, y, 1, 1) == sqrt(x^2 + y^2 + 2)
end

@testset "inverse trig values" begin
    @test asin(Sym(1 // 2)) == PI / 6 && asin(-sqrt(Sym(2)) / 2) == -PI / 4 && asin(sqrt(Sym(3)) / 2) == PI / 3
    @test acos(Sym(-1 // 2)) == 2PI / 3 && acos(sqrt(Sym(2)) / 2) == PI / 4 && acos(Sym(-1)) == PI
    @test atan(sqrt(Sym(3))) == PI / 3 && atan(-sqrt(Sym(3)) / 3) == -PI / 6 && atan(Sym(-1)) == -PI / 4
    @test kind(asin(Sym(1 // 3))) == FN && asin(-x) == -asin(x)
    # positive constants come out of radicals
    @test sqrt(4x) == 2sqrt(x) && sqrt(2x + 2y) == sqrt(Sym(2)) * sqrt(x + y)
    @test sqrt(-8x) == 2sqrt(Sym(2)) * sqrt(-x)
end

@testset "lambertw" begin
    @test lambertw(Sym(0)) == 0 && lambertw(SymJ.E) == 1 && lambertw(-1 / SymJ.E) == -1
    @test lambertw(-1 / SymJ.E, -1) == -1 && kind(lambertw(x)) == FN
    @test diff(lambertw(x), x) == lambertw(x) / (x * (1 + lambertw(x)))
    @test N(lambertw(Sym(1))) ≈ 0.5671432904097838 && lambertw(x)(1.0) ≈ 0.5671432904097838
    @test lambdify(lambertw(x))(1.0) ≈ 0.5671432904097838
    @test N(lambertw(Sym(-1 // 5), -1)) ≈ -2.542641357773527
    # W(u) exp(W(u)) = u, on every branch (SymPy #29737)
    @test simplify(exp(lambertw(x)) * lambertw(x) - x) == 0
    @test simplify(y * lambertw(x) * exp(lambertw(x) + y)) == x * y * exp(y)
    @test simplify(lambertw(x, -1) * exp(lambertw(x, -1))) == x
end

@testset "real, imag, conj, angle" begin
    I = IM
    # constants, and identities for every complex value
    @test real(2 + 3I) == 2 && imag(2 + 3I) == 3 && conj(2 + 3I) == 2 - 3I && reim(1 + 2I) == (1, 2)
    @test angle(Sym(-2)) == PI && angle(1 + I) == PI / 4 && angle(-1 - I) == -3PI / 4 && abs(3 + 4I) == 5
    @test real(PI) == PI && conj(sqrt(Sym(2))) == sqrt(Sym(2)) && real(exp(I * PI / 3)) == 1 // 2 && imag((1 + I)^2) == 2
    @test real(x + 2) == real(x) + 2 && real(3x) == 3real(x) && real(I * x) == -imag(x) && imag(I * x) == real(x)
    @test conj(x + I) == conj(x) - I && conj(x * y) == conj(x) * conj(y) && conj(x^2) == conj(x)^2
    @test conj(sin(x)) == sin(conj(x)) && conj(conj(x)) == x && conj(exp(I * x)) == exp(-I * conj(x))
    @test conj(log(x)) isa Sym && SymJ._isfn1(conj(log(x)), :conj)         # a branch cut: not log(conj(x))
    @test SymJ._isfn1(conj(sqrt(x)), :conj) && SymJ._isfn1(real(x), :real)
    @test real(abs(x)) == abs(x) && imag(abs(x)) == 0 && real(conj(x)) == real(x) && imag(conj(x)) == -imag(x)
    # real for every value: conj(w)^k w^k, real-valued functions (SymPy #29801)
    @test imag(x * conj(x)) == 0
    @test real(x * conj(x)) == conj(x) * x
    @test imag(abs(x)^2) == 0
    @test imag(conj(x)^2 * x^2 * conj(y) * y) == 0
    @test imag(conj(sqrt(x)) * sqrt(x)) == 0
    @test SymJ._isfn1(imag(conj(x) * x^2), :imag)                         # not real
    @test isreal(x * conj(x))
    @test iszero(simplify(abs(x)^2 - conj(x) * x))
    # numeric evaluation, lambdify, adjoints
    @test N(real(x)(x => 2 + 3im)) == 2 && lambdify(real(x) + imag(x))(2 + 3im) == 5
    @test [1+I x; y 2]' == [1-I conj(y); conj(x) 2] && transpose([x y]) == [x; y;;]
    # under assumptions
    R(e, as...) = refine(e; assume = Assumptions(as...))
    @test R(real(x), x ∈ ℝ) == x && R(imag(x), x ∈ ℝ) == 0 && R(conj(x), x ≫ 0) == x
    @test R(angle(x), x ≫ 0) == 0 && R(angle(x), x ≪ 0) == PI && R(real(z), y ≫ 0) == real(z)
    @test R(real(exp(I * x)), x ∈ ℝ) == cos(x) && R(imag(exp(I * x)), x ∈ ℝ) == sin(x)
    @test R(real(x + I * y), x ∈ ℝ, y ∈ ℝ) == x && R(conj(x + I * y), x ∈ ℝ, y ∈ ℝ) == x - I * y
    @test R(real((x + I * y)^2), x ∈ ℝ, y ∈ ℝ) == x^2 - y^2 && R(imag(sin(x + I * y)), x ∈ ℝ, y ∈ ℝ) == cos(x) * sinh(y)
    @test R(angle(x + I * y), x ≫ 0, y ≫ 0) == atan(y / x) && R(real(x * conj(x)), x ∈ ℝ) == x^2
    @test R(abs(exp(I * x)), x ∈ ℝ) == 1 && R(abs(x + I * y), x ∈ ℝ, y ∈ ℝ) == sqrt(x^2 + y^2)
    @test R(x ≫ real(x) - 1, x ∈ ℝ) === true
end
