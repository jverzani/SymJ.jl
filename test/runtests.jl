using SymJ
using SymJ: PI, E, IM, ADD, MUL, POW, FN, NUMBER, kind
using TermInterface
using Test
import ForwardDiff

@syms x y z f() g()

@testset "SymJ.jl" begin
    @testset "SymNumber" begin
        a = SymNumber(typemax(Int))
        @test (a + a) == SymNumber(2 * big(typemax(Int)))
        @test (a + a) - a == a
        @test (a + a - a).tag == SymJ.T_SMALL        # demotes back to Int
        @test SymNumber(1 // 3) + SymNumber(2 // 3) == SymNumber(1)
        @test SymNumber(2 // 3)^-2 == SymNumber(9 // 4)
        @test SymNumber(3)^100 == SymNumber(big(3)^100)
        @test !isequal(SymNumber(1), SymNumber(1.0))
        @test SymNumber(1) == SymNumber(1.0)
        @test SymNumber(1 // 2) < SymNumber(2 // 3)
    end

    @testset "canonical sums" begin
        @test x + x == 2x
        @test x - x == 0
        @test x + y == y + x
        @test (x + y) + z == x + (y + z)
        @test 2x + 3y - 2x == 3y
        @test kind(x + y) == ADD
        @test 2 * (x + y) == 2x + 2y          # numbers distribute over a sum
        @test (x + y) - (x + y) == 0
        @test x + 1 - 1 == x
        @test sum([x, y, x, 1, -1]) == 2x + y
    end

    @testset "canonical products and powers" begin
        @test x * x == x^2
        @test x * y == y * x
        @test x^2 * x^3 == x^5
        @test x * y / y == x
        @test x^2 / x^2 == 1
        @test (x * y)^2 == x^2 * y^2
        @test (x^2)^3 == x^6
        @test kind((x^2)^(1 // 2)) == POW      # not simplified to x
        @test (x^(1 // 2))^2 == x
        @test kind((x + y)^2) == POW          # no automatic expansion
        @test x^0 == 1
        @test x^1 == x
        @test 1^x == 1
        @test exp(x) * exp(y) == exp(x + y)
        @test exp(log(x)) == x
    end

    @testset "exact numbers" begin
        @test sqrt(Sym(4)) == 2
        @test sqrt(Sym(8)) == 2 * sqrt(Sym(2))
        @test sqrt(Sym(2)) * sqrt(Sym(2)) == 2
        @test sqrt(Sym(2))^3 == 2 * sqrt(Sym(2))
        @test Sym(4)^(1 // 4) == sqrt(Sym(2))
        @test Sym(27)^(2 // 3) == 9
        @test Sym(2)^(-1 // 2) == sqrt(Sym(2)) / 2
        @test sqrt(Sym(1 // 4)) == Sym(1 // 2)
        @test sqrt(Sym(-1)) == IM
        @test IM * IM == -1
        @test IM^4 == 1
        @test sqrt(Sym(-4)) == 2IM
        @test Sym(2)^100 == Sym(big(2)^100)
        @test Sym(1) / 3 == Sym(1 // 3)
        @test N(sqrt(Sym(2))) ≈ sqrt(2)
        @test N(IM) == im
        @test Sym(π) == PI
        @test N(Sym(π) / 4) ≈ π / 4
        @test Sym(0.5) + Sym(0.25) == Sym(0.75)
    end

    @testset "functions" begin
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
        @test log(Sym(1)) == 0
        @test log(E) == 1
        @test log(exp(Sym(3))) == 3
        @test kind(log(Sym(2), x)) == FN && length(log(Sym(2), x).args) == 2
        @test log(x, x) == 1
        @test log(x, 1) == 0
        @test kind(sin(Sym(2))) == FN          # exact argument stays symbolic
        @test sin(Sym(2.0)) == Sym(sin(2.0))   # float argument evaluates
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
        @test kind(sec(PI / 2)) == FN && kind(csc(Sym(0))) == FN
        @test sec(-x) == sec(x) && csc(-x) == -csc(x) && cot(-x) == -cot(x)
        @test sec(acos(x)) == 1 / x
        @test kind(sec(x)) == FN && operation(sec(x)) === sec
        @test asec(Sym(1)) == 0 && acsc(Sym(1)) == PI / 2
        @test acot(Sym(0)) == PI / 2 && acot(Sym(1)) == PI / 4 && acot(-x) == -acot(x)
        @test sech(Sym(0)) == 1 && coth(-x) == -coth(x) && csch(-x) == -csch(x)
        @test asech(Sym(1)) == 0 && acoth(-x) == -acoth(x)
        # real cube root
        @test cbrt(Sym(-8)) == -2 && cbrt(Sym(27)) == 3
        @test cbrt(Sym(2)) == Sym(2)^(1 // 3) && cbrt(Sym(-16)) == -2 * Sym(2)^(1 // 3)
        @test cbrt(x)(-8.0) == -2.0 && cbrt(-x) == -cbrt(x)
        # piecewise constant
        @test sign(Sym(-3)) == -1 && sign(sqrt(Sym(2)) - 1) == 1 && sign(-x) == -sign(x)
        @test floor(Sym(7 // 2)) == 3 && floor(Sym(-7 // 2)) == -4 && ceil(Sym(π)) == 4
        @test trunc(Sym(-7 // 2)) == -3 && round(Sym(5 // 2)) == 2 && round(Sym(7 // 2)) == 4
        @test round(sqrt(Sym(2))) == 1 && kind(floor(x)) == FN && floor(floor(x)) == floor(x)
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
    end

    @testset "subs inside products" begin
        @test subs(x * exp(x), exp(x) => y) == x * y
        @test subs(3x^2 * y, x^2 => z) == 3y * z
    end

    @testset "diff" begin
        @test diff(x^3, x) == 3x^2
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

    @testset "expand" begin
        @test expand((x + 1)^2) == x^2 + 2x + 1
        @test expand((x + y) * (x - y)) == x^2 - y^2
        @test expand(x * (y + z)) == x * y + x * z
        @test expand((sqrt(Sym(2)) + 1)^2) == 3 + 2sqrt(Sym(2))
        @test length(expand((x + y + z + 1)^10).args) == 285
        @test expand(sin((x + 1)^2)) == sin(x^2 + 2x + 1)
    end

    @testset "Poly" begin
        @syms a b
        p = Poly((x + y + 1)^3, x, y)
        @test p isa Poly{SymNumber}
        @test nterms(p) == 10
        @test Sym(p) == SymJ._expand_generic((x + y + 1)^3)
        @test gens(p) == [x, y]
        @test degree(p) == 3 && degree(p, y) == 3
        @test exponent_vectors(Poly(x^2 * y + 3, x, y)) == [[2, 1], [0, 0]]
        @test coefficients(Poly(2x^2 - 1 // 2, x)) == [SymNumber(2), SymNumber(-1 // 2)]

        # arithmetic agrees with expression arithmetic
        q = Poly(x - 2y + 3, x, y)
        for r in (p + q, p - q, p * q, q^4, -q, 2 * p, p + 1, p * x)
            @test r isa Poly
        end
        @test Sym(p * q) == expand(Sym(p) * Sym(q))
        @test Sym(q^4) == SymJ._expand_generic((x - 2y + 3)^4)
        @test p - p == zero(p) && iszero(p - p)
        @test isone(q^0)
        @test p * q == q * p
        @test (p + q) * q == p * q + q * q

        # coefficient paths: Int128 fast path, BigInt fallback, rationals, floats
        for ex in ((x + 2^40 * y + 1)^5, (x + typemax(Int) * y)^3, (x / 3 + y / 2 - 1)^6, (0.5x + y)^3)
            @test Sym(Poly(ex, x, y)) == SymJ._expand_generic(ex)
        end
        @test coefficients(Poly((x + 2^40)^4, x))[end] == SymNumber(big(2)^160)

        # powers: multinomial theorem (sparse bases) and repeated multiplication (dense)
        rep(p, k) = foldl(*, fill(p, k))
        for (b, k) in ((Poly(x + y, x, y), 100), (Poly(x - 3y, x, y), 7), (Poly(x / 2 + y - 1, x, y), 9),
                       (Poly(big(2)^70 * x + 1, x), 5), (Poly(sum(x^i for i in 0:9), x), 20),
                       (Poly(x * y + z + 1, x, y, z), 12))
            @test b^k == rep(b, k)
        end
        @test coefficients(Poly(x + y, x, y)^100)[51] == SymNumber(binomial(big(100), 50))
        @test Poly(Sym(:a) * x + 1, x)^3 == rep(Poly(Sym(:a) * x + 1, x), 3)   # Sym coefficients

        # differing generators are unified
        @test Poly(x + 1, x) * Poly(y - 1, y) == Poly(x * y - x + y - 1, x, y)

        # symbolic coefficients
        s = Poly(a * x^2 + b * x + 3, x)
        @test s isa Poly{Sym}
        @test coefficients(s) == [a, b, Sym(3)]
        @test expand(Sym(s * s)) == expand((a * x^2 + b * x + 3)^2)
        @test_throws ArgumentError Poly{SymNumber}(a * x, x)
        @test_throws ArgumentError Poly(sin(x), x)
        @test_throws ArgumentError Poly(sqrt(x), x)

        @test diff(p, x) == Poly(diff(Sym(p), x), x, y)
        @test diff(p, z) == zero(p)
        @test p(x => 1, y => 1) == 27
        @test string(Poly(x^2 - 2x * y + 1, x, y)) == "Poly(x^2 - 2*x*y + 1, x, y)"

        # exponent overflow is detected (16 generators -> 4 bits each, max exponent 14)
        g = [Sym(Symbol(:v, i)) for i in 1:16]
        h = Poly(g[1] + 1, g...)
        @test_throws OverflowError h^15
        @test nterms(h^14) == 15
        @test expand((g[1] + 1)^20) == SymJ._expand_generic((g[1] + 1)^20)  # unaffected (1 generator)
    end

    @testset "subs and calling" begin
        ex = x^2 + 2x + 1
        @test subs(ex, x => 1) == 4
        @test ex(x => y) == y^2 + 2y + 1
        @test (x * y)(x => 2, y => 3) == 6
        @test sin(x)(x => PI) == 0
        @test f(x)(x => 2) == f(Sym(2))
        @test f(2) == f(Sym(2))
        @test free_symbols(x * sin(y) + z) == [x, y, z]
        @test has(sin(x) + y, x)
        @test !has(sin(x) + y, z)
        @test N((x^2)(x => 2.5)) == 6.25
    end

    @testset "calling evaluates numerically" begin
        u = x^2
        @test u(3) === 9
        @test u(3.0) === 9.0
        @test u(1 // 2) === 1 // 4
        @test u(big(3)) isa BigInt && u(big(3)) == 9
        @test u(2im) == -4
        @test u(y) == y^2                       # symbolic argument substitutes
        @test u(Sym(3)) == Sym(9)
        @test u(x => 3) == Sym(9)               # pairs substitute
        # no arguments: the value of a constant expression, as a Julia number; an
        # expression with free symbols is returned as is (nothing to substitute)
        @test u() === u && x() === x
        @test Sym(2)() === 2 && Sym(-1)() === -1
        @test Sym(3 // 4)() === 3 // 4
        @test Sym(π)() === π && SymJ.E() === ℯ
        @test sqrt(Sym(2))() ≈ sqrt(2)
        @test (Sym(π) / 2)() ≈ π / 2
        @test IM() == im
        @test Poly(Sym(5), x)() === 5
        @test (1 / x)(2) === 0.5
        @test (x^-2)(2 // 1) === 1 // 4
        @test (x / 3)(3) == 1 && (x / 3)(3.0) === 1.0
        @test sqrt(x)(4) === 2.0
        @test sqrt(x)(-4) == 2im
        @test (sin(x) + x / 3)(1.0) ≈ sin(1.0) + 1 / 3
        @test exp(x)(0.0) === 1.0
        @test exp(x)(1) ≈ ℯ
        @test log(Sym(2), x)(8.0) ≈ 3.0
        @test (x^y)(x => 2)(3) === 8
        @test (IM * x)(2) == 2im
        @test (x^(1 // 3))(8.0) ≈ 2.0
        @test abs(x)(-3) === 3

        # constants follow the precision of a floating point argument
        @test (Sym(π) * x)(big(1.0)) isa BigFloat
        @test (Sym(π) * x)(big(1.0)) == big(π)
        @test (sqrt(Sym(2)) * x)(big(1.0)) == sqrt(big(2.0))
        @test (x + Sym(1 // 3))(1.0f0) isa Float32
        @test sqrt(Sym(2))(1.0f0) isa Float32   # no free symbols: a constant function
        @test Sym(2)(5) === 2

        # broadcasting and use as an ordinary function
        @test u.([1, 2, 3]) == [1, 4, 9]
        @test map(u, [1.0, 2.0]) == [1.0, 4.0]
        g0 = x^2 - 2
        dg = diff(g0, x)
        z0 = 1.0
        for _ in 1:6
            z0 -= g0(z0) / dg(z0)              # Newton's method with expressions as functions
        end
        @test z0 ≈ sqrt(2)

        # polynomials call the same way
        @test Poly(x^2 + 1, x)(2) === 5
        @test Poly(x^2 + 1, x)(x => y) == y^2 + 1

        # more than one free symbol, or more than one value, needs pairs
        @test_throws ArgumentError (x * y)(2)
        @test_throws ArgumentError (x * y)(Sym(2))
        @test_throws ArgumentError u(1, 2)
        @test (x * y)(x => 2) == 2y
        @test_throws ArgumentError f(x)(2)      # no numeric method for f
        @test f(x)(Sym(2)) == f(Sym(2))
    end

    @testset "lambdify" begin
        w = sin(x)^2 + 3x^5 - 2x + exp(-x / 2) + sqrt(Sym(2)) * x / 3
        fw = lambdify(w)
        for v in (1.5, -0.25, 2)
            @test fw(v) ≈ w(v)
        end
        @test fw(1.5) == w(1.5)
        @test @inferred(fw(1.5)) isa Float64
        @test lambdify(x^2 + 1)(3) === 10                  # same semantics as calling
        @test lambdify(x / 3)(3) == 1 // 1 && lambdify(x / 3)(3.0) === 1.0
        @test lambdify(1 / x)(2) === 0.5 && lambdify(x^-2)(2 // 1) === 1 // 4
        @test lambdify(sqrt(x))(-4.0) == 2im && lambdify(x^(1 // 3))(8.0) ≈ 2
        @test lambdify(IM * x)(2) == 2im
        @test lambdify(sqrt(Sym(2)) * x)(big(1.0)) == sqrt(big(2.0))
        @test lambdify(Sym(π) * x)(1.0f0) isa Float32
        @test lambdify(x^y, x, y)(2, 10) === 1024
        @test lambdify(log(Sym(2), x))(8.0) ≈ 3
        @test lambdify(Sym(7))() === 7                     # no variables
        @test lambdify(Poly(x^2 + 1, x))(2) === 5

        # variable order is explicit; several outputs share subexpressions
        h = lambdify(x - 2y, y, x)
        @test h(1, 10) == 8 && lambdify(x - 2y, (x, y))(1, 10) == -19
        J = lambdify([diff(w, x) x * y; y 1], x, y)
        @test J(1.0, 2.0) ≈ [diff(w, x)(1.0) 2.0; 2.0 1.0] && size(J(1.0, 2.0)) == (2, 2)
        @test lambdify((x + 1, x * y), (x, y))(1, 2) == (2, 2)
        code = string(lambdify_expr(sin(x) * y + sin(x)^2 + cos(sin(x)), (x, y)))
        @test count("Base.sin(arg1)", code) == 1           # computed once

        # errors
        @test_throws ArgumentError lambdify(x * y)          # order needed
        @test_throws ArgumentError lambdify(x * y, x)       # y not a variable
        @test_throws ArgumentError lambdify(f(x), x)        # no numeric method for f

        # generic arithmetic: automatic differentiation through both paths
        @test ForwardDiff.derivative(fw, 1.5) ≈ diff(w, x)(1.5)
        @test ForwardDiff.derivative(w, 1.5) ≈ diff(w, x)(1.5)
        @test ForwardDiff.gradient(v -> lambdify(x^2 * y, x, y)(v...), [2.0, 3.0]) ≈ [12.0, 4.0]
    end

    @testset "parameters" begin
        @parameters a b
        e = a * x^2 + b * x + 1
        @test isparameter(a) && !isparameter(x) && !isparameter(e) && kind(a) == SymJ.PARAMETER
        @test a != Sym(:a) && parameter(:a) == a && hash(a) != hash(Sym(:a))
        @test free_symbols(e) == [x] && free_parameters(e) == [a, b] && free_parameters(x^2) == Sym[]
        # calling substitutes the variable and keeps the parameters
        @test e(3) == 9a + 3b + 1 && e(x => 2) == 4a + 2b + 1 && e(y) == a * y^2 + b * y + 1
        @test (2a)(3) == 2a && a() === a && e() === e
        @test (x^2)(3) === 9                                  # no parameters: numeric, as before
        @test_throws ArgumentError evaluate(e, 2.0)           # parameters have no value
        @test e(x => 2, a => 1, b => 0)() === 5               # all substituted: a number
        # calculus and algebra treat parameters as symbols
        @test diff(e, x) == 2a * x + b && diff(e, a) == x^2 && diff(sin(a * x), a) == x * cos(a * x)
        @test expand((a + x)^2) == x^2 + 2a * x + a^2
        @test factor(a^2 * x^2 - b^2) == (a * x + b) * (a * x - b)
        @test cancel((x^2 - a^2) / (x - a)) == x + a
        @test gcd(x^2 - a^2, x^2 + 2a * x + a^2) == x + a
        @test Poly(e, x) isa Poly{Sym} && coefficients(Poly(e, x)) == [a, b, Sym(1)]
        @test string(a * x^2 + b * x + 1) == "a*x^2 + b*x + 1" && string(x + a) == "x + a"
        # lambdify: the variable, then the parameters
        @test lambdify(e)(2.0, 3.0, 4.0) == 21.0 && lambdify(e, x, a, b)(2, 3, 4) == 21
        @test lambdify(e, (b, a, x))(4, 3, 2) == 21
        # constant tests see parameters as unknown
        @test kind(floor(a)) == FN && kind(sign(a)) == FN && kind(max(a, 1)) == FN
        @test kind(mod(a, 3)) == FN
        @test sin(a)(x => 1) == sin(a)                        # nothing to substitute
        @test N(subs(sin(a) + x, a => 1, x => 2)) ≈ sin(1) + 2
    end

    @testset "relations" begin
        @parameters a
        r = x^2 + 1 ~ 3
        @test r isa Relation && lhs(r) == x^2 + 1 && rhs(r) == 3 && (2 ~ x) == Eq(2, x)
        @test (x ⩵ y) == Eq(x, y) && (x ≪ 1) == Lt(x, 1) && (x ≦ 1) == Le(x, 1)
        @test (x ≫ 1) == Gt(x, 1) && (x ≧ 1) == Ge(x, 1) && (x ≶ 1) == Ne(x, 1) && (x ≷ 1) == Ne(x, 1)
        @test (x ~ 1) != (x ~ 2) && (x ~ 1) != (x ≪ 1) && hash(x ~ 1) == hash(x ~ 1)
        @test string(r) == "x^2 + 1 ~ 3" && string(x ≪ 1) == "x ≪ 1" && string(x ≧ a) == "x ≧ a"
        @test (x + 1 ≪ y) == Lt(x + 1, y)                   # comparison precedence
        @test (x == y) === false && (x ~ x) isa Relation    # == stays structural
        # calling: substitute; decided when no symbols remain
        @test (x^2 ~ 4)(2) === true && (x^2 ~ 4)(3) === false && (x^2 ~ 4)(x => y) == (y^2 ~ 4)
        @test (x ≪ 1)(0.5) === true && (x ≧ 1)(1) === true && (x ≶ 1)(1) === false
        @test (sqrt(Sym(2)) ≪ Sym(3) // 2)() === true && (Sym(π) ~ 22 // 7)() === false
        @test (x ~ a)(x => a) === true && (x ≪ a)(x => a) === false   # identical sides
        @test (x ≪ 1)(x => a) == Lt(a, 1)
        # TermInterface
        @test iscall(r) && operation(r) === (~) && arguments(r) == [x^2 + 1, Sym(3)]
        @test operation(x ≦ 1) === Le && maketerm(Relation, Le, [x, 1], nothing) == (x ≦ 1)
    end

    @testset "solve" begin
        @parameters a b c
        satisfies(eq, rs) = all(r -> abs(N(subs(lhs(eq) - rhs(eq), x => r))) < 1e-8 * max(1, abs(N(r))), rs)
        cases = [(2x + 3 ~ 0, 1), (x^2 - 2 ~ 0, 2), (x^2 ~ -1, 2), (x^2 + x + 1 ~ 0, 2), (x^4 ~ 1, 4),
                 (x^5 ~ 2, 5), (x^3 - 2x - 5 ~ 0, 3), (x^3 - 3x + 1 ~ 0, 3), (x^4 + x + 1 ~ 0, 4),
                 (x^4 - 4x^3 + x + 7 ~ 0, 4), (x^4 + 1 ~ 0, 4), (x^4 - 10x^2 + 1 ~ 0, 4), (x^5 - x - 1 ~ 0, 5),
                 ((x^2 - 1) / (x - 1) ~ 0, 1), (1 / x ~ 0, 0), (x^6 - 1 ~ 0, 6), (0.5x^2 - 1 ~ 0, 2),
                 ((x - 1)^3 * (x + 2) ~ 0, 2), (x^3 ~ x, 3), (2x^2 ~ x + 1, 2)]
        for (eq, nroots) in cases
            r = solve(eq, x)
            @test length(r) == nroots
            @test satisfies(eq, r)
        end
        @test solve(x^2 - 2 ~ 0, x) == [-sqrt(Sym(2)), sqrt(Sym(2))]      # exact, sorted
        @test solve(x^2 ~ -1, x) == [-IM, IM]
        @test solve(2x^2 ~ x + 1, x) == [Sym(-1 // 2), Sym(1)]
        @test solve(x^3 - 3x + 1 ~ 0, x) == [2cos(10PI / 9), 2cos(4PI / 9), 2cos(2PI / 9)]  # real form
        @test all(r -> imag(N(r)) == 0, solve(x^3 - 3x + 1 ~ 0, x))
        @test all(r -> !SymJ._hassymbol(r) && N(r) isa Number, solve(x^5 - x - 1 ~ 0, x))  # numeric (a + b*im)
        # symbolic coefficients
        @test solve(a * x + b ~ 0, x) == [-b / a]
        @test solve(x^2 ~ a, x) == [-sqrt(a), sqrt(a)]
        q = solve(a * x^2 + b * x + c ~ 0, x)
        @test length(q) == 2 && all(r -> cancel(subs(a * x^2 + b * x + c, x => r)) == 0 ||
                                     abs(N(subs(a * x^2 + b * x + c, x => r, a => 2, b => 3, c => -5))) < 1e-12, q)
        @test solve(x^2 - a^2 ~ 0, x) == [-a, a] || Set(solve(x^2 - a^2 ~ 0, x)) == Set([-a, a])
        # errors
        @test_throws ArgumentError solve(x^2 + 1, x)          # an equation is required
        @test string(solve(x ≪ 1, x)) == "(-∞, 1)"            # inequalities: sets
        @test_throws ArgumentError solve(x ~ x, x)
        # through one function of x: substitute, solve, invert; extraneous roots dropped
        @test solve(sin(x) ~ 1 // 2, x) == [PI / 6, 5PI / 6]
        @test solve(cos(x) ~ 0, x) == [-PI / 2, PI / 2] && solve(tan(x) ~ 1, x) == [PI / 4]
        @test solve(sin(2x) ~ 1, x) == [PI / 4]
        @test solve(sin(x)^2 ~ 1 // 4, x) == [-PI / 6, PI / 6, 5PI / 6, 7PI / 6]
        @test solve(exp(x) ~ 2, x) == [log(Sym(2))] && solve(log(x) ~ 2, x) == [exp(Sym(2))]
        @test solve(exp(x) ~ a, x) == [log(a)]
        @test solve(sqrt(x) ~ 3, x) == [Sym(9)] && solve(sqrt(x) ~ -1, x) == Sym[]
        @test solve(2^x ~ 8, x) == [log(Sym(8)) / log(Sym(2))]
        # powers of one exponential, also with integer bases that are powers of one integer
        @test solve(exp(2x) - 3exp(x) + 2 ~ 0, x) == [Sym(0), log(Sym(2))]
        @test solve(exp(x) + exp(-x) ~ 2, x) == [Sym(0)]
        @test solve(4^x - 3 * 2^x + 2 ~ 0, x) == [Sym(0), Sym(1)] && solve(8^x ~ 2 * 4^x, x) == [Sym(1)]
        # a square root mixed with x: squared, extraneous roots dropped
        @test solve(sqrt(x + 1) ~ x - 1, x) == [Sym(3)]                # 0 is extraneous
        @test solve(sqrt(x) + x ~ 6, x) == [Sym(4)] && solve(x - sqrt(x) ~ 0, x) == [Sym(0), Sym(1)]
        @test_throws ArgumentError solve(sin(x) + cos(x) ~ 1, x)     # two functions: not yet
        # Lambert W: (αx + β) e^(γx + δ) = K, both real branches when there are two
        @test solve(x * exp(x) ~ 1, x) == [lambertw(Sym(1))] && solve(x ~ exp(-x), x) == [lambertw(Sym(1))]
        @test Set(solve(x * exp(x) ~ -1 // 5, x)) == Set([lambertw(Sym(-1 // 5)), lambertw(Sym(-1 // 5), -1)])
        @test solve(x * exp(x) ~ 0, x) == [Sym(0)]
        for eq in (2^x ~ 3x, exp(x) ~ x + 2, (2x + 1) * exp(3x) ~ 4)
            rs = solve(eq, x)
            @test !isempty(rs) && satisfies(eq, rs)
        end
        @test_throws ArgumentError solve(x * exp(x^2) ~ 1, x)       # not a Lambert form

        # linear systems: a Dict; free variables omitted; nothing when inconsistent
        holds(eqs, d) = all(eq -> cancel(subs(lhs(eq) - rhs(eq), d)) == 0 ||
                                  abs(N(subs(lhs(eq) - rhs(eq), d))) < 1e-12, eqs)
        @test solve([x + y ~ 3, x - y ~ 1], [x, y]) == Dict(x => 2, y => 1)
        sys3 = [2x + 3y - z ~ 1, x - y + 2z ~ 4, 3x + y + z ~ 6]
        @test solve(sys3, [x, y, z]) == Dict(x => Sym(11 // 5), y => Sym(-1), z => Sym(2 // 5))
        @test solve([x + y + z ~ 1, x - z ~ 0], [x, y, z]) == Dict(x => z, y => 1 - 2z)
        @test solve([x + y ~ 1, x + y ~ 2], [x, y]) === nothing
        @test solve([x + y ~ 1, 2x + 2y ~ 2], [x, y]) == Dict(x => 1 - y)
        @test solve((x + y ~ 3, x - y ~ 1), (x, y)) == Dict(x => 2, y => 1)       # tuples
        for (eqs, vs) in (([a * x + y ~ 1, x - y ~ b], [x, y]), ([a * x + b * y ~ 1, b * x - a * y ~ 0], [x, y]),
                          ([x / 2 + y / 3 ~ 1, x ~ y], [x, y]), ([x + y ~ sin(a), x - y ~ cos(a)], [x, y]),
                          ([x + 2y + z ~ a, y - z ~ b, x + z ~ 0], [x, y, z]))
            d = solve(eqs, vs)
            @test d !== nothing && length(d) == length(vs) && holds(eqs, d)
        end
        @test solve([a * x + b * y ~ 1, b * x - a * y ~ 0], [x, y]) == Dict(x => a / (a^2 + b^2), y => b / (a^2 + b^2))
        d = solve([0.5x + y ~ 1, x - y ~ 0.25], [x, y])
        @test N(d[x]) ≈ 5 / 6 && N(d[y]) ≈ 7 / 12
        @test_throws ArgumentError solve([x * y ~ 1, x + y ~ 2], [x, y])     # nonlinear: not yet
        @test_throws ArgumentError solve([1 / x ~ 1], [x])
        @test string(solve([x ≪ 1], [x])) == "(-∞, 1)"                # a one-unknown system

        # < and friends compare constants as numbers; with symbols, a pointer to ≪
        @test sqrt(Sym(2)) < 2 && Sym(π) > 3 && !(2 ≥ sqrt(Sym(5))) && Sym(1 // 3) <= 0.5
        @test sort([Sym(3), sqrt(Sym(2)), Sym(π)]) == [sqrt(Sym(2)), Sym(3), Sym(π)]
        @test_throws ArgumentError x + y < 2
        @test_throws ArgumentError x >= 1
        @test occursin("≪", (try x < 1 catch e; sprint(showerror, e) end))
        @test occursin("≧", (try x >= 1 catch e; sprint(showerror, e) end))

        # systems of inequalities: one unknown (sets), several (cylindrical, linear)
        @test string(solve([x ≫ 1, x ≪ 3], [x])) == "(1, 3)"
        @test string(solve([x^2 ≪ 4, x ≧ 0], [x])) == "[0, 2)"
        @test string(solve([x^2 ~ 4, x ≫ 0], [x])) == "{2}"
        @test string(solve([x ≶ 2, x^2 ≦ 9], [x])) == "[-3, 2) ∪ (2, 3]"
        C2 = solve([x + y ≪ 2, y ≪ 3], [x, y])
        @test C2 isa SymJ.CylindricalSet && string(C2) == "y ∈ (-∞, 3), x ∈ (-∞, -y + 2)"
        T2 = solve([x ≧ 0, y ≧ 0, x + y ≦ 1], [x, y])
        @test string(T2) == "y ∈ [0, 1], x ∈ [0, -y + 1]"
        @test string(solve([x ≧ 0, y ≧ 0, x + y ≦ 1], [y, x])) == "x ∈ [0, 1], y ∈ [0, -x + 1]"
        @test Dict(x => 0.2, y => 0.3) in T2 && !(Dict(x => 0.8, y => 0.3) in T2) && Dict(x => 0, y => 1) in T2
        @test isempty(solve([x + y ≪ 1, x + y ≫ 2], [x, y]))
        S3 = solve([x + y + z ~ 1, x ≧ 0, y ≧ 0, z ≧ 0], [x, y, z])
        @test string(S3) == "z ∈ [0, 1], y ∈ [0, -z + 1], x = -y - z + 1"
        @test Dict(x => 0.5, y => 0.25, z => 0.25) in S3 && !(Dict(x => 0.5, y => 0.5, z => 0.25) in S3)
        @test string(solve([x ≧ 0, y ≧ 0, z ≧ 0, x + 2y + 3z ≦ 6], [x, y, z])) ==
              "z ∈ [0, 2], y ∈ [0, -3*z/2 + 3], x ∈ [0, -2*y - 3*z + 6]"
        @test_throws ArgumentError solve([x * y ≪ 1, x ≫ 0], [x, y])        # nonlinear: not yet
        @test_throws ArgumentError solve([x ≶ y, x ≫ 0], [x, y])           # not convex

        # inequalities: unions of intervals with exact endpoints
        S(r) = string(solve(r, x))
        @test S(x^2 ≫ 1) == "(-∞, -1) ∪ (1, ∞)" && S(x^2 ≦ 1) == "[-1, 1]"
        @test S(x^2 ≦ 0) == "{0}" && S(x^2 ≪ 0) == "∅" && S(x^2 + 1 ≫ 0) == "(-∞, ∞)"
        @test S((x - 1) / (x + 2) ≧ 0) == "(-∞, -2) ∪ [1, ∞)" && S(1 / x ≪ 0) == "(-∞, 0)"
        @test S(x^3 - x ≧ 0) == "[-1, 0] ∪ [1, ∞)" && S(x^2 - 2 ≪ 0) == "(-sqrt(2), sqrt(2))"
        @test S(x ≶ 3) == "(-∞, 3) ∪ (3, ∞)" && S((x - 1)^2 ≫ 0) == "(-∞, 1) ∪ (1, ∞)"
        @test S(x^2 ≧ x) == "(-∞, 0] ∪ [1, ∞)"
        P2s, T = solve(x^2 ≫ 1, x), solve(x ≦ 3, x)
        @test P2s isa SymJ.SymSet && !(0 in P2s) && 2 in P2s && -1.5 in P2s && !(-1 in P2s)
        @test string(P2s ∪ T) == "(-∞, ∞)" && string(P2s ∩ T) == "(-∞, -1) ∪ (1, 3]"
        @test string(setdiff(T, P2s)) == "[-1, 1]" && isempty(setdiff(T, T))
        @test sqrt(Sym(2)) in solve(x^2 ≦ 2, x) && !(sqrt(Sym(2)) in solve(x^2 ≪ 2, x))
        @test_throws ArgumentError solve(a * x ≫ 1, x)               # parameters: not yet
        @test_throws ArgumentError solve(sin(x) ≫ 0, x)
        @test_throws ArgumentError solve(x^5 + a * x + 1 ~ 0, x)  # degree 5, symbolic
    end

    @testset "TermInterface" begin
        ex = 2x + sin(y)
        @test iscall(ex) && isexpr(ex)
        @test !iscall(x) && !isexpr(x)
        @test operation(ex) === (+)
        @test operation(sin(y)) === sin
        @test operation(f(y)) isa SymFunction
        @test operation(x^y) === (^)
        @test operation(2x * y) === (*)
        for e in (ex, 3x * y^2 + 1, sin(x)^2, f(x, y), log(Sym(2), x), x^y)
            @test maketerm(Sym, operation(e), arguments(e), nothing) == e
            @test maketerm(Sym, head(e), children(e), nothing) == e
            @test arguments(e) isa Vector{Sym}
        end
        @test Set(arguments(x + y + 1)) == Set([x, y, Sym(1)])
        @test sorted_arguments(x^2 + x + 1) == [x^2, x, Sym(1)]
        @test arity(f(x, y)) == 2
    end

    @testset "printing" begin
        @test string(x^2 + 2x + 1) == "x^2 + 2*x + 1"
        @test string(x - y) == "x - y"
        @test string(-x / (2y)) == "-x/(2*y)"
        @test string(x^-2) == "1/x^2"
        @test string(sqrt(Sym(8))) == "2*sqrt(2)"
        @test string(Sym(2)^(1 // 3)) == "2^(1/3)"
        @test string(exp(x)) == "exp(x)"
        @test string(2IM) == "2*im"
        @test string((x + y)^2) == "(x + y)^2"
        @test string(PI) == "π"
        @test sprint(show, MIME"text/plain"(), x + 1) == "x + 1"
    end

    @testset "callable / Function" begin
        @test x isa Function
        @test f isa Function
        @test (x .+ [1, 2]) == [x + 1, x + 2]
    end

    @testset "Poly division and gcd" begin
        P(e) = Poly(e, x, y, z)
        # small deterministic pseudo-random polynomials (LCG; no extra test deps)
        state = Ref(UInt64(12345))
        rnd(k) = (state[] = state[] * 6364136223846793005 + 1442695040888963407; Int((state[] >> 33) % k))
        randpoly(nt, d) = P(sum((rnd(7) - 3) * x^rnd(d) * y^rnd(d) * z^rnd(2) for _ in 1:nt) + rnd(5) + 1)

        # long division
        q, r = divrem(Poly(x^3 - 1, x), Poly(x - 1, x))
        @test q == Poly(x^2 + x + 1, x) && iszero(r)
        q, r = divrem(Poly(x^2 + 1, x), Poly(2x, x))
        @test q == Poly(x / 2, x) && r == Poly(Sym(1), x)
        @test_throws DivideError divrem(P(x), zero(P(x)))
        @test divexact(P((x + y)^3), P(x + y)) == P((x + y)^2)
        @test_throws ArgumentError divexact(P(x + 1), P(y))
        @test divides(P(x * y + x), P(y + 1)) == (true, P(x))

        for _ in 1:30
            a, b = randpoly(6, 4), randpoly(3, 3)
            q, r = divrem(a, b)
            @test q * b + r == a
            @test !any(m -> SymJ._mono_divides(b.monos[1], m, b.bits, 3), r.monos)
        end

        # gcd
        @test gcd(P((x + 1) * (x - 2)^2 * (x + y)), P((x + 1) * (x + y) * (y - 3))) == P((x + 1) * (x + y))
        @test gcd(P(2x + 2), P(4x + 4)) == P(2x + 2)              # integer content kept
        @test gcd(P(x / 2 + 1 // 2), P(x^2 - 1)) == P(x + 1)      # rational input: monic
        @test gcd(P(-x - 1), P(x + 1)) == P(x + 1)                # positive leading coefficient
        @test isone(gcd(P(x^2 + y^2), P(x + y)))
        @test gcd(zero(P(x)), P(3x + 6)) == P(3x + 6)
        @test gcd(P(Sym(6)), P(Sym(4))) == P(Sym(2))
        @test gcd(Poly(x^2 - 1, x), Poly(y * x + y, y, x)) == Poly(x + 1, x, y)  # gens unified
        @test_throws ArgumentError gcd(P(0.5x), P(x))
        @test lcm(P(x^2 - 1), P(x^2 + 2x + 1)) == P((x - 1) * (x + 1)^2)

        for _ in 1:20
            a, b, g = randpoly(3, 3), randpoly(3, 3), randpoly(2, 2)
            A, B = a * g, b * g
            G = gcd(A, B)
            @test G == gcd(B, A)
            @test divides(A, G)[1] && divides(B, G)[1]
            @test divides(G, g)[1]
            @test isone(gcd(divexact(A, G), divexact(B, G)))
        end

        # modular gcd agrees with the PRS reference implementation
        for _ in 1:40
            a, b, g = randpoly(4, 4), randpoly(4, 3), randpoly(3, 3)
            @test gcd(a * g, b * g) == SymJ._gcd(a * g, b * g)
            @test gcd(a, b) == SymJ._gcd(a, b)
        end
        # large coefficients: several primes and Chinese remaindering
        g = P(big(3)^40 * x^2 * y - big(7)^25 * z + 2^40)
        A, B = g * P(x^3 + 5y - big(11)^20), g * P(y^2 * z - 3x + 1)
        @test gcd(A, B) == g
        @test gcd(P(big(2)^100 * (x + 1)), P(big(6)^80 * (x + 1) * (y + 2))) == P(big(2)^80 * (x + 1))
        # the first prime used, 2^31 - 1, is unlucky here: the cofactors agree mod p
        @test gcd(P((x + y) * x), P((x + y) * (x + (2^31 - 1)))) == P(x + y)
        @test gcd(P((x + y) * (z + 1) * x), P((x + y) * (z + 1) * (x + (2^31 - 1) * y))) == P((x + y) * (z + 1))
        # univariate, high degree
        u, v = Poly((x + 1)^30 * (x - 2)^20, x), Poly((x + 1)^25 * (x + 3)^20, x)
        @test gcd(u, v) == Poly((x + 1)^25, x)

        # squarefree factorization
        prodf(u, fs) = foldl((r, fm) -> r * fm.first^fm.second, fs; init = P(Sym(u)))
        u, fs = factor_squarefree(P((x - 1)^3 * (x + 2)^2 * (x + 3)))
        @test u == 1 && fs == [P(x + 3) => 1, P(x + 2) => 2, P(x - 1) => 3]
        @test factor_squarefree(P(6 * (x + 1)^2)) == (SymNumber(6), [P(x + 1) => 2])
        @test factor_squarefree(P(-(x - 1)^2)) == (SymNumber(-1), [P(x - 1) => 2])
        @test factor_squarefree(P((x / 2 + 1 // 2)^2)) == (SymNumber(1 // 4), [P(x + 1) => 2])
        @test factor_squarefree(P(x * (y + 1)^2)) == (SymNumber(1), [P(x) => 1, P(y + 1) => 2])
        @test factor_squarefree(P(Sym(5))) == (SymNumber(5), Pair{Poly{SymNumber},Int}[])
        @test_throws ArgumentError factor_squarefree(zero(P(x)))
        @test squarefree_part(P((x + y)^2 * (x - y) * (y + 1)^3)) == P((x + y) * (x - y) * (y + 1))
        for _ in 1:15
            fs0 = [randpoly(2, 2) for _ in 1:3]
            p = P(Sym(rnd(5) + 1)) * fs0[1] * fs0[2]^2 * fs0[3]^3
            u, fs = factor_squarefree(p)
            @test prodf(u, fs) == p
            @test issorted(last.(fs)) && allunique(last.(fs))
            for (f, _) in fs                     # squarefree: gcd(f, ∂f/∂x, ∂f/∂y, ∂f/∂z) = 1
                @test isone(foldl(gcd, (diff(f, v) for v in (x, y, z)); init = f))
                @test sign(f.coefs[1]) > 0
            end
            for i in eachindex(fs), j in (i + 1):length(fs)   # pairwise coprime
                @test isone(gcd(fs[i].first, fs[j].first))
            end
        end
        # expressions
        @test factor_squarefree(x^3 - x^2 - x + 1) == (x + 1) * (x - 1)^2
        @test factor_squarefree(2x^2 + 4x + 2) == 2 * (x + 1)^2
        @test factor_squarefree((x^2 - 1) / (x^2 + 2x + 1)^2) == (x - 1) / (x + 1)^3
        @test factor_squarefree(sin(x)^2 + 2sin(x) + 1) == (sin(x) + 1)^2
        @test factor_squarefree(x / (x + 0.5)) == x / (x + 0.5)   # inexact: unchanged
        @test squarefree_part(x^3 - x^2 - x + 1) == (x + 1) * (x - 1)

        # on expressions
        @test gcd(x^2 - 1, x^2 + 2x + 1) == x + 1
        @test gcd(x^2 * y - y, x * y + y) == x * y + y
        @test lcm(x^2 - 1, x + 1) == x^2 - 1
        @test_throws ArgumentError gcd(sin(x), x)
    end

    @testset "numer_denom and cancel" begin
        @test numer_denom(x / y + 1) == (x + y, y)
        @test numer_denom(3x / 4) == (3x, Sym(4))
        @test numer_denom(1 / sqrt(x)) == (Sym(1), sqrt(x))
        @test numer_denom(x) == (x, Sym(1))
        @test numer_denom(Sym(2) / 3) == (Sym(2), Sym(3))

        @test cancel((x^2 - 1) / (x + 1)) == x - 1
        @test cancel((x^2 - 1) / (x^2 + 2x + 1)) == (x - 1) / (x + 1)
        @test cancel(1 / x + 1 / x^2) == (x + 1) / x^2
        @test cancel((2x + 2) / (4x)) == (x + 1) / (2x)
        @test cancel((x * y + x) / (y^2 - 1)) == x / (y - 1)
        @test cancel((x - 1) / (1 - x)) == -1
        @test cancel(1 / (x + 1) - 1 / (x - 1)) == -2 / (x^2 - 1)
        @test cancel((x / 2 + 1) / (x + 2)) == Sym(1 // 2)
        @test cancel(x * (x + 1)) == x^2 + x                   # polynomials come back expanded
        @test cancel((sin(x)^2 - 1) / (sin(x) + 1)) == sin(x) - 1
        @test cancel((x - 1) / (sqrt(x) - 1)) == (x - 1) / (sqrt(x) - 1)  # radicals are independent
        @test cancel(x / (x + 0.5)) == x / (x + 0.5)           # floats: unchanged
        @test cancel(Sym(3) / 6) == Sym(1 // 2)
        @test cancel(2) == Sym(2)

        # cancel never changes the value
        for e in ((x^3 - y^3) / (x - y) + 1 / (x + y), (x^2 * y - y) / (x * y + y) - x / (y + 1),
                  (1 / x - 1 / y) / (x - y), ((x + 1)^3 - 1) / x)
            c = cancel(e)
            for (a, b) in ((2, 3), (-1 // 2, 5), (7, -4))
                @test e(x => a, y => b) == c(x => a, y => b)
            end
        end
        @test cancel((x^3 - y^3) / (x - y)) == x^2 + x * y + y^2
        @test cancel((1 / x - 1 / y) / (x - y)) == -1 / (x * y)
    end

    @testset "factor and apart" begin
        P1(e) = Poly(e, x)
        P3(e) = Poly(e, x, y, z)
        prodf(p, (u, fs)) = foldl((r, fm) -> r * fm.first^fm.second, fs; init = Poly(Sym(u), gens(p)...))
        @test factor(P1(x^4 - 1)) == (SymNumber(1), [P1(x + 1) => 1, P1(x - 1) => 1, P1(x^2 + 1) => 1])
        @test factor(P1(6x^2 + 5x + 1)) == (SymNumber(1), [P1(2x + 1) => 1, P1(3x + 1) => 1])
        @test factor(P1(2x^3 - 6x)) == (SymNumber(2), [P1(x) => 1, P1(x^2 - 3) => 1])
        @test factor(P1(-(x^2 + 1)^2 * (x - 3) / 4)) == (SymNumber(-1 // 4), [P1(x - 3) => 1, P1(x^2 + 1) => 2])
        @test length(factor(P1(x^5 - x - 1))[2]) == 1                 # irreducible
        @test length(factor(P1(x^4 + 1))[2]) == 1                     # irreducible over Q, splits mod every p
        @test Set(factor(P1(x^4 + 4))[2]) == Set([P1(x^2 - 2x + 2) => 1, P1(x^2 + 2x + 2) => 1])
        @test [degree(g) for (g, _) in factor(P1(x^60 - 1))[2]] == [1, 1, 2, 2, 2, 4, 4, 4, 8, 8, 8, 16]
        @test length(factor(P1(prod(x - i for i in -15:15)))[2]) == 31
        big2 = P1((big(10)^20 * x^3 + 7x + big(3)^40) * (x^4 - big(2)^60 * x + 1))
        @test prodf(big2, factor(big2)) == big2 && length(factor(big2)[2]) == 2
        # Swinnerton-Dyer: irreducible of degree 16, many factors mod every prime
        s2, s3, s5, s7 = sqrt(Sym(2)), sqrt(Sym(3)), sqrt(Sym(5)), sqrt(Sym(7))
        sd = expand(prod(x + a + b + c + d for a in (s2, -s2), b in (s3, -s3), c in (s5, -s5), d in (s7, -s7)))
        @test length(factor(P1(sd))[2]) == 1
        # multivariate
        @test Set(factor(P3(x^2 * y^2 - 1))[2]) == Set([P3(x * y - 1) => 1, P3(x * y + 1) => 1])
        @test factor(P3(x^3 + y^3))[2] == [P3(x + y) => 1, P3(x^2 - x * y + y^2) => 1]
        @test factor(P3(x * y + x + y + 1))[2] == [P3(x + 1) => 1, P3(y + 1) => 1]
        @test length(factor(P3(x^6 - y^6))[2]) == 4
        @test factor(P3((x + y + z + 1)^2 * (x - y) * (x * y * z - 1)))[2] ==
              [P3(x - y) => 1, P3(x * y * z - 1) => 1, P3(x + y + z + 1) => 2]
        @test_throws ArgumentError factor(zero(P1(x)))

        # EEZ: against Kronecker substitution (the reference) on random products
        w = Sym(:w)
        P4(e) = Poly(e, x, y, z, w)
        state4 = Ref(UInt64(77))
        rnd4(k) = (state4[] = state4[] * 6364136223846793005 + 1442695040888963407; Int((state4[] >> 33) % k))
        rp4() = P4(sum((rnd4(5) - 2) * x^rnd4(3) * y^rnd4(2) * z^rnd4(2) for _ in 1:3) + rnd4(3) + 1)
        for _ in 1:12
            p = rp4() * rp4() * rp4()
            _, sq = factor_squarefree(p)
            for (q, _) in sq
                act = [i for i in 1:4 if SymJ._deg_in_all(q, i) > 0]
                length(act) > 1 || continue
                @test Set(SymJ._factor_eez(q, act)) == Set(SymJ._factor_kronecker(q, act))
            end
        end
        # a point that splits an irreducible factor: recombination recovers the factors
        P2(e) = Poly(e, x, y)
        for (e, us, want) in ((x^2 - y, [P2(x - 2), P2(x + 2)], [x^2 - y]),
                              ((x^2 - y) * (x + y), [P2(x - 2), P2(x + 2), P2(x + 4)], [x + y, x^2 - y]))
            qs = SymJ._shift_var(P2(expand(e)), 2, 4)
            gs = [SymJ._shift_var(g, 2, -4) for g in SymJ._eez_split(qs, us, 1, [2])]
            @test Set(Sym.(gs)) == Set(want)
        end
        # van Hoeij recombination: the same factors as subset search, and the
        # Swinnerton-Dyer polynomials (irreducible, many factors mod every prime)
        densez(S) = (v = zeros(BigInt, degree(S) + 1);
                     for (e, c) in zip(exponent_vectors(S), coefficients(S)); v[e[1]+1] = big(SymJ._numerator(c)); end; v)
        function vh_vs_subsets(e)
            f = densez(P1(expand(e)))
            p, mf = SymJ._choose_prime(f)
            vh = SymJ._vanhoeij(f, mf, p)
            B = big(2)^(length(f) - 1) * (isqrt(sum(c -> c^2, f)) + 1)
            k, pk = 1, big(p)
            while pk <= 2 * abs(f[end]) * B; k += 1; pk *= p; end
            zs = SymJ._recombine(f, SymJ._hensel_lift(SymJ._zmod(f, pk), mf, p, k), pk)
            length(mf), vh, zs
        end
        for e in (prod(x - i for i in 1:12), prod(x^2 + i for i in 1:6), x^30 - 1,
                  (x^4 - 10x^2 + 1) * (x^2 - 2) * (x - 1) * (x^3 - 2))
            nmod, vh, zs = vh_vs_subsets(e)
            @test vh !== nothing && Set(vh) == Set(zs)
        end
        function swinnerton_dyer(k)
            S = x
            for q in (2, 3, 5, 7, 11)[1:k]
                r = sqrt(Sym(q))
                S = expand(expand(S(x => x + r)) * expand(S(x => x - r)))
            end
            P1(S)
        end
        for k in (4, 5)
            f = densez(swinnerton_dyer(k))
            p, mf = SymJ._choose_prime(f)
            k == 5 && @test length(mf) >= SymJ._VANHOEIJ_MIN_FACTORS  # the lattice path
            @test SymJ._vanhoeij(f, mf, p) == [f]                     # irreducible
        end
        @test length(factor(swinnerton_dyer(5))[2]) == 1

        # Wang's leading-coefficient distribution, at the point y = 0 where its test holds
        @test SymJ._wang_nondivisors(BigInt[2, 5], big(1), big(1)) == BigInt[2, 5]
        @test SymJ._wang_nondivisors(BigInt[3, 6], big(2), big(1)) === nothing   # 6 shares 2 and 3
        qw = P2(expand(((y + 2) * x + 1) * ((y + 5) * x^2 + 3)))
        gw = SymJ._eez_wang(qw, [P2(2x + 1), P2(5x^2 + 3)], 1, [2])
        @test gw !== nothing && Set(Sym.(gw)) == Set(expand.([(y + 2) * x + 1, (y + 5) * x^2 + 3]))
        # beyond Kronecker substitution's reach (image degree ~ 10*9*6*6)
        big4 = P4(expand((x^5 * y^3 + z^4 * w^2 - 7) * (x^4 - y^5 * z + w^3 + 2) * (x * y * z * w + 1)))
        @test_throws ArgumentError SymJ._factor_kronecker(big4, [1, 2, 3, 4])
        u4, f4 = factor(big4)
        @test prodf(big4, (u4, f4)) == big4 && length(f4) == 3
        # content in the main variable is factored separately
        @test Set(factor(P4((y^2 - 1) * (x * y + z)))[2]) == Set([P4(y - 1) => 1, P4(y + 1) => 1, P4(x * y + z) => 1])

        # random products: exact product, irreducible factors, right multiplicities
        state = Ref(UInt64(2024))
        rnd(k) = (state[] = state[] * 6364136223846793005 + 1442695040888963407; Int((state[] >> 33) % k))
        rp() = P3(sum((rnd(5) - 2) * x^rnd(3) * y^rnd(2) for _ in 1:3) + rnd(4) + 1)
        for _ in 1:12
            p = rp() * rp()^2 * rp()
            iszero(p) && continue
            u, fs = factor(p)
            @test prodf(p, (u, fs)) == p
            @test all(fm -> length(factor(fm.first)[2]) == 1 && sign(fm.first.coefs[1]) > 0, fs)
        end
        for _ in 1:12
            p = P1(sum((rnd(9) - 4) * x^i for i in 0:rnd(6)) + 1) * P1(sum((rnd(9) - 4) * x^i for i in 0:rnd(6)) + 1)
            iszero(p) && continue
            u, fs = factor(p)
            @test prodf(p, (u, fs)) == p
            @test all(fm -> length(factor(fm.first)[2]) == 1, fs)
        end

        # expressions
        @test factor(x^4 - 1) == (x + 1) * (x - 1) * (x^2 + 1)
        @test factor(x^2 * y^2 - 1) == (x * y + 1) * (x * y - 1)
        @test factor((x^2 - 1) / (x^2 + 2x + 1)) == (x - 1) / (x + 1)
        @test factor(sin(x)^2 - 1) == (sin(x) - 1) * (sin(x) + 1)
        @test factor(x / (x + 0.5)) == x / (x + 0.5)
        @test factor(12) == 12

        # partial fractions
        @test apart(1 / (x^2 - 1)) == (1 // 2) / (x - 1) - (1 // 2) / (x + 1)
        @test apart((x^3 + 1) / (x^2 - 1)) == x + 1 / (x - 1)
        @test apart(1 / (x^3 + x)) == 1 / x - x / (x^2 + 1)
        @test apart((3x + 5) / (x + 1)^2) == 3 / (x + 1) + 2 / (x + 1)^2
        @test apart((2x^2 + 3) / (x^3 - 2)) == (2x^2 + 3) / (x^3 - 2)   # irreducible denominator
        @test apart(x^2 + 1) == x^2 + 1
        for e in (1 / ((x - 1)^2 * (x + 1)), (x^4 + 1) / (x^2 * (x^2 + 1)^2), 1 / (x^4 - 1),
                  (x / 2 + 1) / (3x^2 - 3), (x^5 + 2) / ((x^2 + x + 1)^2 * (x - 2)^3 * x))
            @test cancel(apart(e) - e) == 0
        end
        # symbolic coefficients: over ℚ(params), with variables or parameters
        @parameters a b c
        @test apart(1 / (x^2 - a^2), x) == 1 / (2a * (x - a)) - 1 / (2a * (x + a))
        @test apart(1 / (x^2 + a), x) == 1 / (x^2 + a)                 # irreducible over ℚ(a)
        for e in (1 / ((x - a) * (x - b)), (x + 1) / ((x - a)^2 * (x + b)), x^3 / ((x - a) * (x - 1)),
                  1 / (a * x * (x + b)), (x^2 + c) / ((x - a)^2 * (x^2 + b)), 1 / ((x - y) * (x + y)),
                  (a * x + 1) / ((x - a) * (x^2 - b)), 1 / ((x - a) * (x - 2a) * (x - 3a)))
            r = apart(e, x)
            @test cancel(r - e) == 0
        end
        @test length(SymJ.arguments(apart(1 / ((x - a) * (x - b) * (x - c)), x))) == 3
        @test_throws ArgumentError apart(x / y)                         # two symbols: give x
    end

    @testset "RatFun" begin
        R = SymJ.RatFun
        P(e) = Poly(e, x, y)
        r1 = SymJ._ratfun(P(x^2 - y^2), P(2x - 2y))            # reduces to (x + y)/2
        @test SymJ._ispoly(r1) && r1.num == P(x / 2 + y / 2)
        r2 = SymJ._ratfun(P(x), P(-2y))                         # integer num/den, positive den
        @test r2.den == P(2y) && r2.num == P(-x)
        @test r1 + r2 == R(P(x / 2 + y / 2)) + SymJ._ratfun(P(-x), P(2y))
        @test iszero(r1 - r1) && isone(r1 / r1) && r1 * inv(r1) == one(R)
        @test (r2^2) == r2 * r2 && r2^-1 == inv(r2)
        @test R(SymNumber(3)) + R(P(x)) == R(P(x + 3))           # constants mix with polynomials
        @test Sym(r2) == -x / (2y) && zero(R) == R(SymNumber(0))
        @test_throws DivideError inv(zero(R))
    end

    @testset "iteration (like Number)" begin
        for ex in (x, x + y, sin(x), Sym(2), Poly(x + 1, x))
            @test collect(ex) == fill(ex)
            @test length(ex) == 1 && size(ex) == () && ndims(ex) == 0 && !isempty(ex)
            @test eltype(ex) == typeof(ex)
            @test [e for e in ex] == fill(ex)
            @test first(ex) == ex && last(ex) == ex && only(ex) == ex
            @test ex[] == ex && ex[1] == ex && ex[1, 1] == ex
            @test_throws BoundsError ex[2]
            @test ex in ex
            (a,) = ex
            @test a == ex
        end
        @test sum(x + y) == x + y        # sum over a one-element iterable
        @test [x; y] == [x, y]
        @test x .+ [1, 2] == [x + 1, x + 2]
    end

    @testset "ambiguities" begin
        @test isempty(Test.detect_ambiguities(SymJ; recursive = true))
    end

    include("extensions.jl")
end
