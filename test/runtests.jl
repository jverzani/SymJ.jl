using SymJ
using SymJ: PI, E, IM, ADD, MUL, POW, FN, NUMBER, kind
using SymJ: @parameters, parameter, isparameter, free_parameters   # internal for now
using TermInterface
using Test
import ForwardDiff
using Latexify: latexify

@syms x y z f() g()
_add_many_test(f, k, n) = sum(subs(f, k => i) for i in 1:n)

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
        # the minus comes out when the term printed first is negative
        @test string(sin(x - y)) == "sin(x - y)" && string(sin(y - x)) == "-sin(x - y)" && string(cos(y - x)) == "cos(x - y)"
        @test string(tan(1 - x)) == "-tan(x - 1)" && string(sinh(2 - x^2)) == "-sinh(x^2 - 2)" && string(sin(x - y + 1)) == "sin(x - y + 1)"
        @test string(simplify(cos(y) * sin(x) - cos(x) * sin(y))) == "sin(x - y)"
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

    @testset "Gröbner bases" begin
        @syms gx gy gz gr ga
        same(A, B) = length(A) == length(B) && all(iszero(expand(a - b)) for (a, b) in zip(A, B))
        @test same(groebner([gx^2 + gy^2 - 1, gx - gy], gx, gy), [gx - gy, gy^2 - 1 // 2])
        @test same(groebner([gx^2 + gy^2 - 1, gx - gy], gx, gy; order = :grevlex), [gy^2 - 1 // 2, gx - gy])
        @test same(groebner([gx * gy - 1, gx^2 + gy^2 - 4], gx, gy), [gx + gy^3 - 4gy, gy^4 - 4gy^2 + 1])
        @test same(groebner([gx * gy - 1, gx^2 + gy^2 - 4], gx, gy; order = :grevlex),
                   [gy^3 + gx - 4gy, gx^2 + gy^2 - 4, gx * gy - 1])
        @test same(groebner([gx + gy + gz, gx * gy + gy * gz + gz * gx, gx * gy * gz - 1], gx, gy, gz),
                   [gx + gy + gz, gy^2 + gy * gz + gz^2, gz^3 - 1])                        # cyclic-3
        @test same(groebner([gx + 2gy + 2gz - 1, gx^2 + 2gy^2 + 2gz^2 - gx, 2gx * gy + 2gy * gz - gy], gx, gy, gz),
                   [gx - 60gz^3 + 158gz^2 / 7 + 8gz / 7 - 1, gy + 30gz^3 - 79gz^2 / 7 + 3gz / 7,
                    gz^4 - 10gz^3 / 21 + gz^2 / 84 + gz / 84])                              # katsura-2
        @test groebner([gx - 1, gx - 2], gx) == [1]                                          # no common zero
        @test same(groebner([gy - gx^2, gz - gx^3], gz, gy, gx), [gz - gx^3, gy - gx^2])     # infinitely many
        @test same(groebner([gy - gx^2, gz - gx^3], gx, gy, gz; order = :grevlex), [gx^2 - gy, gx * gy - gz, gy^2 - gx * gz])
        @test same(groebner([gx^2 + gy^2 - gr^2, gx - gy], gx, gy), [gx - gy, gy^2 - gr^2 / 2])   # a parameter
        @test same(groebner([ga * gx - 1, gx * gy - ga], gx, gy), [gx - 1 / ga, gy - ga^2])
        @test same(groebner([gx^2 - 2, gy - gx]), [gx - gy, gy^2 - 2])                      # all the symbols
        @test same(groebner([gx^2 + gy^2 ~ 1, gx ~ gy], gx, gy), [gx - gy, gy^2 - 1 // 2])   # equations
        @test same(groebner([gx^300 - 1, gx^200 - 1], gx), [gx^100 - 1])
        vs = [Sym(Symbol("gv$i")) for i in 1:7]                                            # wider words
        @test same(groebner([vs[i]^200 - vs[i+1] for i in 1:6], vs...; order = :grevlex)[end:end], [vs[6]^200 - vs[7]])
        # FGLM agrees with lex computed directly, and the inputs are in the ideal
        for (F, xs) in (([gx^2 + 2gy * gz - 3, gy^2 - gx * gz + 1, gz^2 + gx * gy - 2], [gx, gy, gz]),
                        ([gx + gy + gz, gx * gy + gy * gz + gz * gx, gx * gy * gz - 1], [gx, gy, gz]))
            @test same(groebner(F, xs...), SymJ._groebner(F, xs, Sym[], :lex; fglm = false))
            for ord in (:lex, :grlex, :grevlex)
                G = groebner(F, xs...; order = ord)
                @test all(f -> iszero(reduced(f, G, xs...; order = ord)[2]), F)
            end
        end
        # division
        G = groebner([gx^2 + gy^2 - 1, gx - gy], gx, gy)
        q, r = reduced(gx^3 + gy, G, gx, gy)
        @test r == 3gy / 2 && iszero(expand(q[1] * G[1] + q[2] * G[2] + r - (gx^3 + gy)))
        q, r = reduced(gx^2 * gy - gr * gx, [gx * gy - gr, gy^2 - 1], gx, gy; order = :grevlex)
        @test iszero(expand(q[1] * (gx * gy - gr) + q[2] * (gy^2 - 1) + r - (gx^2 * gy - gr * gx)))
        @test_throws ArgumentError groebner([gx], gx; order = :revlex)
        @test_throws ArgumentError groebner([gx^2 - 0.5], gx)
        @test_throws ArgumentError groebner([sin(gx) - gy], gx, gy)
        @test_throws ArgumentError groebner([gx < gy], gx, gy)
    end

    @testset "resultants" begin
        @syms ra rb rc rd re rf rg rh ry rz
        for m in (:auto, :modular, :subresultant)
            R(p, q) = resultant(p, q, x; method = m)
            @test R(x^2 - 2, x^2 - 3) == 1
            @test R(x^2 + ry^2 - 1, x - ry) == 2ry^2 - 1
            @test R(ra * x^2 + rb * x + rc, 2ra * x + rb) == 4ra^2 * rc - ra * rb^2
            # res(Π(x - aᵢ), Π(x - cⱼ)) = Π(aᵢ - cⱼ)
            @test iszero(expand(R((x - ra) * (x - rb), (x - rc) * (x - rd)) - (ra - rc) * (ra - rd) * (rb - rc) * (rb - rd)))
            @test R((x - 1) * (x + 2), (x - 1) * (x^2 + 5)) == 0                  # a common root
            @test R(x^2 / 2 - 1 // 3, x / 5 + 1) == Sym(73 // 150)
            @test R(x^3 + ry, x^2 * ry - rz) == ry^5 - rz^3
            @test R(Sym(3), x^2 + 1) == 9 && R(x^2 + 1, Sym(3)) == 9 && R(Sym(2), Sym(3)) == 1
            @test R(ry * x^2 + x + 1, x - ry) == ry^3 + ry + 1                     # lc vanishes at ry = 0
            @test R(x - ry^100, x^3 - 2) == ry^300 - 2
            @test R(Sym(0), x + 1) == 0
            # res(B, A) = (-1)^(mn) res(A, B)
            @test R(x^2 + x + 1, x^3 - ry) == R(x^3 - ry, x^2 + x + 1)          # mn even
            @test R(x + ry, x^3 - 2) == -R(x^3 - 2, x + ry)
        end
        # the methods agree on random dense and sparse inputs
        st = Ref(UInt64(7))                       # a small LCG: no Random dependency
        rnd(r) = (st[] = st[] * 6364136223846793005 + 1442695040888963407; first(r) + Int((st[] >> 33) % length(r)))
        rp(vs, d, k) = sum(rnd(-9:9) * prod(v^rnd(0:d) for v in vs) for _ in 1:k)
        for (vs, d, k) in (([x], 12, 10), ([x, ry], 5, 10), ([x, ry, rz], 3, 8))
            for _ in 1:3
                p, q = rp(vs, d, k), rp(vs, d, k)
                has(p, x) && has(q, x) || continue
                r1 = resultant(p, q, x; method = :modular)
                @test iszero(expand(r1 - resultant(p, q, x; method = :subresultant)))
            end
        end
        @test nterms(Poly(resultant(ra*x^3 + rb*x^2 + rc*x + rd, re*x^3 + rf*x^2 + rg*x + rh, x),
                          ra, rb, rc, rd, re, rf, rg, rh)) == 34                  # the known 34 terms
        # polynomials, coefficients in functions, floats
        P = resultant(Poly(x^2 + ry^2 - 1, x, ry), Poly(x - ry, x, ry), x)
        @test P isa Poly && Sym(P) == 2ry^2 - 1
        @test Sym(resultant(Poly(x^2 - 2, x), Poly(x - 3, x))) == 7
        @test resultant(x^2 - sin(ry), x - cos(ry), x) == cos(ry)^2 - sin(ry)
        @test N(resultant(x^2 - 2.0, x - 1.5, x)) ≈ 0.25
        @test_throws ArgumentError resultant(sin(x), x - 1, x)
        @test_throws ArgumentError resultant(x, x, x; method = :fast)
        # discriminants
        @test discriminant(ra * x^2 + rb * x + rc, x) == rb^2 - 4ra * rc
        @test discriminant(x^3 + ra * x + rb, x) == -4ra^3 - 27rb^2
        @test discriminant((x - 1)^2 * (x + 1), x) == 0
        @test discriminant(3x + 1, x) == 1
        @test Sym(discriminant(Poly(x^2 - 5, x))) == 20
        @test_throws ArgumentError discriminant(Sym(4), x)
    end

    @testset "series" begin
        @syms sa
        S(e, args...) = Sym(series(e, args...))
        @test S(sin(x), x) == x - x^3 / 6 + x^5 / 120
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
        @test_throws ArgumentError series(f(x), x)                         # no derivative rule
        @test_throws ArgumentError series(sin(x), x + 1)
    end

    @testset "subs inside products" begin
        @test subs(x * exp(x), exp(x) => y) == x * y
        @test subs(3x^2 * y, x^2 => z) == 3y * z
    end

    @testset "diff" begin
        @test diff(x^3, x) == 3x^2
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

        # Exponents beyond the uniform packing (16 generators: 4 bits each, max 14). A
        # field too small gets a new layout in the same word (fields sized per generator);
        # past the word, construction and gcd, factor, cancel, ... widen to UInt128 or
        # BigInt, while direct arithmetic (type-stable) raises.
        g = [Sym(Symbol(:v, i)) for i in 1:16]
        h = Poly(g[1] + 1, g...)
        @test h.layout == SymJ._uniform(16, 4)
        @test nterms(h^14) == 15 && (h^14).layout == h.layout
        h15 = h^15                                       # g[1] needs 5 bits: new layout
        @test SymJ._monotype(h15) === UInt64 && h15.layout != h.layout
        @test nterms(h15) == 16 && Sym(h15) == expand((g[1] + 1)^15) && h15 == h^14 * h
        @test hash(h15) == hash(Poly((g[1] + 1)^15, g...))
        hs = Poly(sum(g) + 1, g...)
        @test_throws SymJ._PackOverflow hs^15            # 16 fields of 5 bits: not in a UInt64
        @test occursin("widen", sprint(showerror, SymJ._PackOverflow(15)))
        @test SymJ._monotype(Poly(prod(g)^15, g...)) === UInt128
        @test SymJ._monotype(Poly(prod(g)^(2^40), g...)) === BigInt
        x40 = Poly(g[1]^(2^40) + g[16], g...)           # one high degree among 16: one UInt64
        @test SymJ._monotype(x40) === UInt64 && degree(x40) == 2^40 && Sym(x40) == g[1]^(2^40) + g[16]
        w = convert(SymJ.Poly{SymNumber,UInt128}, h^14)
        @test w == h^14 && hash(w) == hash(h^14) && w * convert(SymJ.Poly{SymNumber,UInt128}, h) == h15
        @test Poly(g[1]^100 * g[2] * g[3], g...) == Poly(g[1]^100 * g[2], g...) * Poly(g[3], g...)
        @test diff(Poly(g[1]^30, g...), g[1]) == Poly(30 * g[1]^29, g...)
        @test gcd(Poly(g[1]^20 - 1, g...), Poly(g[1]^10 - 1, g...)) == Poly(g[1]^10 - 1, g...)
        # division whose remainder outgrows the inputs' degrees (x^k by x + y^5)
        qd, rd = divrem(Poly(g[1]^4, g...), Poly(g[1] + g[2]^5, g...))
        @test rd == Poly(g[2]^20, g...) && qd * Poly(g[1] + g[2]^5, g...) + rd == Poly(g[1]^4, g...)
        big14 = prod(g[3:16])^(2^40)                     # BigInt words
        num, den = Poly(g[1]^4 + big14, g...), Poly(g[1] + g[2]^5, g...)
        qd, rd = divrem(num, den)
        @test SymJ._monotype(num) === BigInt
        @test qd * den + rd == num && rd == Poly(g[2]^20 + big14, g...)
        # more than 64 generators
        g70 = [Sym(Symbol(:u, i)) for i in 1:70]
        P70 = Poly(sum(g70)^2, g70...)
        @test nterms(P70) == 70 * 71 ÷ 2 && Sym(P70) == expand(sum(g70)^2)
        @test gcd(P70, Poly(sum(g70) * g70[1], g70...)) == Poly(sum(g70), g70...)
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

    @testset "calling with keywords" begin
        @parameters pa
        e = x^2 + y
        @test e(; x = 2, y = 3) === 7                    # all given: a Julia number
        @test e(; x = 2.0, y = 3) === 7.0 && (x / 3 + y)(; x = 1, y = 2) === 7 // 3
        @test e(; x = 2) == y + 4                        # something left: substitution
        @test e(; x = y, y = 3) == y^2 + 3               # simultaneous, like subs
        @test (x + 1)(; z = 3, x = 1) === 2              # names not in ex are ignored
        @test (x + 1)(; z = 3) == x + 1 && Sym(5)(; z = 1) === 5
        @test (pa * x^2)(; x = 2, pa = 3) === 12         # parameters by name too
        @test (pa * x^2)(; x = 2) == 4pa
        @test (sin(x) + y)(; x = big(1.0), y = 2) isa BigFloat
        @test (sqrt(Sym(2)) * x)(; x = 1.0) ≈ sqrt(2)
        @test (x^2 ~ 4)(; x = 2) === true && (x + y ≫ 1)(; x = 0) == (y ≫ 1)
        @test Poly(x^2 + y, x, y)(; x = 1, y = 1) === 2
        @test_throws ArgumentError e(x => 1; y = 2)      # not mixed with other forms
        @test_throws ArgumentError e(1; y = 2)
        @test_throws ArgumentError (pa * x + parameter(:x))(; x = 1)   # x names both
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
        u = lambdify(x^2 * y, (x, y))                       # values in one tuple or vector
        @test u((1, 2)) == u([1, 2]) == u(1, 2) == 2
        @test_throws ArgumentError u(1)
        @test_throws ArgumentError u(1, 2, 3)
        @test_throws ArgumentError u((1, 2, 3))
        @test lambdify(x^2)([1 2; 3 4]) == [7 10; 15 22]   # one variable: the value as it is
        @test string(u) == "lambdified function of (x, y)"
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
        # several logarithms of x: combined (rational coefficients), candidates checked
        @test solve(log(x) + log(x - 1) ~ log(Sym(6)), x) == [3]       # not -2 (complex logs)
        @test solve(2log(x) - log(x + 2) ~ 0, x) == [2]
        @test solve(log(x) + log(x + 3) ~ 2log(Sym(2)), x) == [1]
        @test solve(log(x^2) - log(x) ~ log(Sym(5)), x) == [5]
        @test solve(log(x) - log(x - 1) ~ 1, x) == [Sym(ℯ) / (Sym(ℯ) - 1)]
        @test_throws ArgumentError solve(log(x) + sin(x) ~ 1, x)
        @test_throws ArgumentError solve(x * log(x) + log(x - 1) ~ 1, x)
        # N outside the real domain, with an exact integer argument
        @test N(log(Sym(-3))) ≈ log(complex(-3.0)) && N(log(Sym(-2)) + log(Sym(-3)) - log(Sym(6))) ≈ 2π * im
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

        # linear systems: one Dict in a vector; free variables omitted; [] when inconsistent
        holds(eqs, d) = all(eq -> cancel(subs(lhs(eq) - rhs(eq), d)) == 0 ||
                                  abs(N(subs(lhs(eq) - rhs(eq), d))) < 1e-12, eqs)
        @test solve([x + y ~ 3, x - y ~ 1], [x, y]) == [Dict(x => 2, y => 1)]
        sys3 = [2x + 3y - z ~ 1, x - y + 2z ~ 4, 3x + y + z ~ 6]
        @test solve(sys3, [x, y, z]) == [Dict(x => Sym(11 // 5), y => Sym(-1), z => Sym(2 // 5))]
        @test solve([x + y + z ~ 1, x - z ~ 0], [x, y, z]) == [Dict(x => z, y => 1 - 2z, z => z)]
        @test solve([x + y ~ 1, x + y ~ 2], [x, y]) == Dict{Sym,Sym}[]
        @test solve([x + y ~ 1, 2x + 2y ~ 2], [x, y]) == [Dict(x => 1 - y, y => y)]
        @test solve((x + y ~ 3, x - y ~ 1), (x, y)) == [Dict(x => 2, y => 1)]       # tuples
        for (eqs, vs) in (([a * x + y ~ 1, x - y ~ b], [x, y]), ([a * x + b * y ~ 1, b * x - a * y ~ 0], [x, y]),
                          ([x / 2 + y / 3 ~ 1, x ~ y], [x, y]), ([x + y ~ sin(a), x - y ~ cos(a)], [x, y]),
                          ([x + 2y + z ~ a, y - z ~ b, x + z ~ 0], [x, y, z]))
            ds = solve(eqs, vs)
            @test length(ds) == 1 && length(ds[1]) == length(vs) && holds(eqs, ds[1])
        end
        @test solve([a * x + b * y ~ 1, b * x - a * y ~ 0], [x, y]) == [Dict(x => a / (a^2 + b^2), y => b / (a^2 + b^2))]
        d = only(solve([0.5x + y ~ 1, x - y ~ 0.25], [x, y]))
        @test N(d[x]) ≈ 5 / 6 && N(d[y]) ≈ 7 / 12

        # nonlinear systems: through a lex Gröbner basis
        sameset(A, B) = length(A) == length(B) && all(a -> any(b -> keys(a) == keys(b) &&
                                                         all(k -> iszero(expand(a[k] - b[k])), keys(a)), B), A)
        @test sameset(solve([x^2 + y^2 ~ 1, x ~ y], [x, y]),
                      [Dict(x => sqrt(Sym(2)) / 2, y => sqrt(Sym(2)) / 2), Dict(x => -sqrt(Sym(2)) / 2, y => -sqrt(Sym(2)) / 2)])
        @test sameset(solve([x * y ~ 1, x + y ~ 3], [x, y]),
                      [Dict(x => (3 + sqrt(Sym(5))) / 2, y => (3 - sqrt(Sym(5))) / 2),
                       Dict(x => (3 - sqrt(Sym(5))) / 2, y => (3 + sqrt(Sym(5))) / 2)])
        @test length(solve([x^2 ~ 1, y^2 ~ 1], [x, y])) == 4                     # not in shape position
        @test sameset(solve([x^2 + y^2 ~ 4, (x - 1)^2 + y^2 ~ 4], [x, y]),
                      [Dict(x => 1 // 2, y => sqrt(Sym(15)) / 2), Dict(x => 1 // 2, y => -sqrt(Sym(15)) / 2)])
        @test isempty(solve([x^2 + y^2 ~ 1, x^2 + y^2 ~ 4], [x, y]))
        @test sameset(solve([x^2 + y^2 ~ 1], [x, y]), [Dict(x => sqrt(1 - y^2), y => y), Dict(x => -sqrt(1 - y^2), y => y)])   # y free
        @test sameset(solve([1 / x + 1 / y ~ 1, x ~ y], [x, y]), [Dict(x => 2, y => 2)])         # rational
        @test sameset(solve([y ~ x^2 / (x - 1), x * (x - 1) ~ 0], [x, y]), [Dict(x => 0, y => 0)])   # pole dropped
        @test length(solve([x + y + z ~ 0, x * y + y * z + z * x ~ 0, x * y * z ~ 1], [x, y, z])) == 6
        @test length(solve([x^2 + y^2 ~ a^2, x ~ y], [x, y])) == 2                # a parameter
        for (eqs, vs, n) in (([x^2 + x * y + y^2 ~ 3, x^2 - y^2 ~ x - 1], [x, y], 4),
                             ([x^2 + y^2 + z^2 ~ 6, x * y ~ 2, y + z ~ 3], [x, y, z], 4),
                             ([0.5x^2 + y^2 ~ 1.5, x ~ 2y], [x, y], 2),                     # floats
                             ([y ~ x^5 - x - 1, y ~ 0], [x, y], 5),                        # numeric roots
                             ([x + 2y + 2z ~ 1, x^2 + 2y^2 + 2z^2 ~ x, 2x * y + 2y * z ~ y], [x, y, z], 4))
            ds = solve(eqs, vs)
            @test length(ds) == n && all(d -> holds(eqs, d), ds)
        end
        @test_throws ArgumentError solve([sin(x) ~ y, x ~ y^2], [x, y])
        @test solve([1 / x ~ 1], [x]) == [Dict(x => 1)]
        # the call's form fixes the result's: one equation and one unknown, or vectors of both
        @test solve(x ~ 2, x) == [Sym(2)]
        @test solve([x ~ 2], [x]) == [Dict(x => 2)]
        @test_throws ArgumentError solve(x ~ 2, [x])
        @test_throws ArgumentError solve([x ~ 2], x)
        @test_throws ArgumentError solve((x ~ 2,), x)
        @test_throws ArgumentError solve(x + y ~ 2, x, y)
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
        # calling a set: a Bool with every unknown given, else the slice in the rest
        @test C2(x => 1, y => 0.5) === true && C2(x => 1, y => 4) === false
        @test string(C2(y => 1)) == "(-∞, 1)" && C2(y => 1) isa SymJ.SymSet && C2(y => 1)(0)
        T3 = solve([x ≧ 0, y ≧ 0, z ≧ 0, x + 2y + 3z ≦ 6], [x, y, z])
        @test string(T3(z => 1)) == "y ∈ [0, 3/2], x ∈ [0, -2*y + 3]"
        @test string(T3(z => 1, y => 0)) == "[0, 3]" && isempty(T3(x => 7))
        @test T3(x => 1, y => 1, z => 1) === true && T3(x => 2, y => 1, z => 1) === false
        @test isempty(solve([x + y ≪ 1, x + y ≫ 2], [x, y])(x => 0))
        @test_throws ArgumentError C2(z => 1)                                # not an unknown
        S1 = solve(x^2 ≫ 1, x)
        @test S1(2) && !S1(0) && S1(x => 2) && S1(Sym(3) // 2) && S1.var == x
        @test_throws ArgumentError S1(y => 2)
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
            @test !any(m -> SymJ._mono_divides(b.monos[1], m, b.layout, 3), r.monos)
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
        # Zippel (opt-in): the same results, sparse or not, falling back when its
        # interpolation cannot apply (a content in the other variables, a dense gcd)
        w = Sym(:w)
        P4(e) = Poly(e, x, y, z, w)
        for (g, c1, c2) in ((x^10 * y^3 * z + 2y^7 * z^5 * w - x^2 * w^8 + 3y^9 + x * y * z * w,
                             x^9 * z^4 + y^8 * w^2 + 5, y^10 * w^3 - x^6 * w^5 + 7z),
                            ((y + z) * x^7 + y^9 * w^4 + z^6 * w^5 + 1, x^4 * w + y^5 + 2, (w^3 - y) * x^5 + z^8),
                            ((y^2 + 1) * (x^3 + z), x + z * w, x^2 - w),
                            ((x + y + z + w + 1)^4, (x - y + 2z * w + 1)^2 + y, (y * z - x + 3)^2 + w),
                            (big(3)^20 * x * y + big(5)^15 * z^4 * w - 7, x - y * z, y^3 + w))
            G = P4(expand(g))
            A, B = G * P4(c1), G * P4(c2)
            @test gcd(A, B; method = :zippel) == gcd(A, B) == G
        end
        saved = state[]                  # leave the later random tests unchanged
        for _ in 1:20
            a, b, g = randpoly(4, 4), randpoly(4, 3), randpoly(3, 3)
            @test gcd(a * g, b * g; method = :zippel) == gcd(a * g, b * g)
        end
        state[] = saved
        @test gcd(P(x // 2 + y // 3), P((3x + 2y) * z); method = :zippel) == P(x + 2y // 3)
        @test_throws ArgumentError gcd(P(x), P(y); method = :sparse)
        # :auto (the default): dense gcds stay with Brown, sparse or many-variable ones
        # take Zippel; the results agree
        Gd = P4(expand((x + y + z + w + 1)^6))
        Ad, Bd = Gd * P4(x - y + 2z + 1), Gd * P4(y * z - x + 3)
        @test !SymJ._use_zippel(Ad, Bd)
        vs = [Sym(Symbol(:v, i)) for i in 1:12]
        gv = expand(sum(vs)^2 + vs[1] * vs[2])
        Av, Bv = Poly(expand(gv * (vs[1] - vs[3])), vs...), Poly(expand(gv * (vs[2] + vs[4] - 1)), vs...)
        @test SymJ._use_zippel(Av, Bv)
        @test gcd(Av, Bv) == gcd(Av, Bv; method = :zippel) == Poly(gv, vs...)
        @test !SymJ._use_zippel(P(x^2 * y * z - 1), P(x * y * z^2 + 1))      # tiny: dense
        @test occursin(":auto", sprint(showerror, try gcd(P(x), P(y); method = :sparse) catch e; e end))
        # the keyword through the functions built on gcd (the same results)
        zg = x^10 * y^3 * z + 2y^7 * z^5 * w - x^2 * w^8 + 3y^9 + x * y * z * w
        za, zb = expand(zg * (x^9 * z^4 + y^8 * w^2 + 5)), expand(zg * (y^10 * w^3 - x^6 * w^5 + 7z))
        # (the sign follows the order of the generators, which is by hash and so differs
        # between Julia versions)
        @test gcd(za, zb; method = :zippel) == gcd(za, zb) &&
              (iszero(expand(gcd(za, zb) - zg)) || iszero(expand(gcd(za, zb) + zg)))
        @test lcm(za, zb; method = :zippel) == lcm(za, zb)
        @test lcm(P4(za), P4(zb); method = :zippel) == lcm(P4(za), P4(zb))
        @test cancel(za / zb; gcd_method = :zippel) == cancel(za / zb)
        zs = expand(zg^2 * (x * y - z))
        @test factor_squarefree(zs; gcd_method = :zippel) == factor_squarefree(zs)
        @test factor_squarefree(P4(zs); gcd_method = :zippel) == factor_squarefree(P4(zs))
        @test squarefree_part(zs; gcd_method = :zippel) == squarefree_part(zs)
        zf = expand((x * y + z^2 + 1)^2 * (x - y * w))
        @test factor(zf; gcd_method = :zippel) == factor(zf)
        @test factor(P4(zf); gcd_method = :zippel) == factor(P4(zf))
        @test apart(1 / ((x - y) * (x + y * z)), x; gcd_method = :zippel) == apart(1 / ((x - y) * (x + y * z)), x)
        @test apart(1 / (x^2 - 1); gcd_method = :zippel) == apart(1 / (x^2 - 1))
        @test cancel(Sym(3); gcd_method = :zippel) == 3 && factor(6; gcd_method = :zippel) == 6
        @test_throws ArgumentError cancel(za / zb; gcd_method = :sparse)
        @test_throws ArgumentError factor(zf; gcd_method = :sparse)

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
        # beyond the 64-bit packing: the slower path, not a refusal or a silent no-op
        @syms a1 a2 a3 a4 a5 a6 a7 a8
        @test gcd(x^130 * a1 - a2 + a3 + a4 + a5 + a6 + a7 + a8, x^65 + a1) == 1
        wg = x^70 + a1 * a2 + a3^2
        @test gcd(expand(wg * (x - a4 + a5 * a6)), expand(wg * (x + a7 - a8))) == wg
        @test cancel((x^40 - a1^40) / (x - a1) + a2 + a3 + a4 + a5 + a6 + a7 + a8 + y + z) ==
              expand(sum(x^k * a1^(39 - k) for k in 0:39)) + a2 + a3 + a4 + a5 + a6 + a7 + a8 + y + z
        @test factor(expand((x^40 + a1) * (x - a2) * (a3 + a4 + a5 + a6 + a7 + a8 + y + z))) ==
              (x - a2) * (x^40 + a1) * (a3 + a4 + a5 + a6 + a7 + a8 + y + z)

        @test gcd(x^2 - 1, x^2 + 2x + 1) == x + 1
        @test gcd(x^2 * y - y, x * y + y) == x * y + y
        @test lcm(x^2 - 1, x + 1) == x^2 - 1
        @test_throws ArgumentError gcd(sin(x), x)
    end

    @testset "numer_denom and cancel" begin
        # square roots of rationals: cancelled through their conjugates, kept when the
        # denominator would not get simpler
        s2, s3 = sqrt(Sym(2)), sqrt(Sym(3))
        @test cancel((x^2 + 2s2 * x + 2) / (x + s2)) == x + s2
        @test cancel((x^2 - 2) / (x - s2)) == x + s2
        @test cancel((x^2 - 5 - 2s2 * s3) / (x - s2 - s3)) == x + s2 + s3
        @test cancel((x^2 - 3) / ((x + s3) * (x + 2))) == (x - s3) / (x + 2)
        @test iszero(cancel(cancel((x + s2 * y) / (x^2 - 2y^2)) - 1 / (x - s2 * y)))
        @test cancel(1 / (x - s2)) == 1 / (x - s2)                     # not (x + √2)/(x² - 2)
        @test cancel((x + s2) / (x^2 + 1)) == (x + s2) / (x^2 + 1)
        @test cancel((sqrt(x) + 1) / (x - 1)) == (sqrt(x) + 1) / (x - 1)  # sqrt(x) is no number
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

    @testset "solve with assumptions" begin
        @syms a b n
        @test solve(x^2 ~ 4, x; assume = (x ≫ 0,)) == [2] == solve(x^2 ~ 4, x; assume = @assume(x > 0))
        @test solve(x^2 ~ 2, x; assume = @assume(x in Rationals)) == Sym[]
        @test (@assuming x > 0 begin
            solve(x^2 ~ 4, x)
        end) == [2]
        @test @assuming(x > 0, x != 2, begin solve(x^2 ~ 4, x) end) == Sym[]
        @test (@assuming x >= 1 begin
            @assuming x < 3 begin solve(x^2 ~ 9, x) end              # nested: both apply
        end) == Sym[]
        @test solve(x^2 ~ 4, x; assume = Assumptions(x ≶ 2)) == [-2]
        @test solve(x^2 ~ -1, x; assume = (x ∈ ℝ,)) == Sym[]
        @test solve(x^3 ~ 8, x; assume = (x ∈ ℝ,)) == [2]
        @test solve(x^2 - 5x + 6 ~ 0, x; assume = (x ≫ 2,)) == [3]
        @test solve(x^2 ~ a, x; assume = (x ≫ 0, a ≫ 0)) == [sqrt(a)]
        @test solve(x^2 ~ a, x; assume = (x ≫ 0, a ~ 9)) == [3]        # a substituted
        @test solve(x^2 ~ -1, x) == [-IM, IM]                          # none: unchanged
        assuming(x ≫ 0) do
            @test solve(x^2 ~ 4, x) == [2]
            @test solve(x^2 ~ 4, x; assume = (x ≪ 1,)) == Sym[]       # the keyword adds
            @test string(solve(x^2 ≫ 1, x)) == "(1, ∞)"
        end
        # cases: a sign split when a is real, the fact and its negation when not
        cs = solve(x^2 ~ a, x; assume = (x ∈ ℝ, a ∈ ℝ))
        @test cs isa SymJ.Cases && length(cs) == 3
        @test cs[1] == (Any[a ≪ 0] => Sym[]) && cs[2] == (Any[a ~ 0] => [Sym(0)])
        @test cs(a => 4) == [-2, 2] && cs(a => 0) == [0] && cs(a => -1) == Sym[]
        cs = solve(x^2 ~ a, x; assume = (x ≫ 0,))
        @test length(cs) == 2 && cs[1] == (Any[a ≫ 0] => [sqrt(a)])
        @test cs(a => 4) == [2] && cs(a => -1) == Sym[] && cs[a => 4] == cs[1]
        @test length(solve(a * x ~ 1, x; assume = (x ≫ 0,))) == 2
        cs = solve(x^2 ~ 9, x; assume = (x ≫ a,))
        @test [c.second for c in cs] == [[-3, 3], [3], [3], Sym[]]
        @test cs(a => -5) == [-3, 3] && cs(a => 0) == [3] && cs(a => 4) == Sym[]
        cs = solve(x^2 - 3x + 2 ~ 0, x; assume = (n ∈ ℤ, x ≪ n))
        @test cs(n => 1) == Sym[] && cs(n => 2) == [1] && cs(n => 5) == [1, 2]
        # rationality of constants: radicals, π and ℯ
        @test solve(x^2 ~ 2, x; assume = (x ∈ ℚ,)) == Sym[] == solve(x^2 ~ 2, x; assume = (x ∈ ℤ,))
        @test solve(4x^2 ~ 1, x; assume = (x ∈ ℚ,)) == [-1 // 2, 1 // 2]
        @test solve(4x^2 ~ 1, x; assume = (x ∈ ℤ,)) == Sym[]
        @test solve(x^2 - x - 1 ~ 0, x; assume = (x ∈ ℚ,)) == Sym[]
        @test solve(x^3 ~ 2, x; assume = (x ∈ ℚ,)) == Sym[]          # complex roots too
        s2, s3, s6 = sqrt(Sym(2)), sqrt(Sym(3)), sqrt(Sym(6))
        @test SymJ._const_rational(s2 * s3 * s6) === true && SymJ._const_rational(s2 * s3) === false
        @test SymJ._const_rational(2 + 3 * PI) === false && SymJ._const_rational(PI * E) === nothing
        @test SymJ._const_rational(s2 + s3) === nothing                # several radicals: unknown
        @test !(s2 ∈ ℚ) && !(PI ∈ ℤ) && (s2 * s3 * s6 ∈ ℤ)
        cs = solve(2x ~ n, x; assume = (x ∈ ℤ,))
        @test length(cs) == 2 && cs[1].second == [n / 2]
        @test occursin("a > 0  =>  [sqrt(a)]", sprint(show, MIME"text/plain"(), solve(x^2 ~ a, x; assume = (x ≫ 0,))))
        # inequalities: within the unknown's set
        @test string(solve(x^2 ≫ 1, x; assume = (x ≪ 0,))) == "(-∞, -1)"
        # a fact about a product: no recursion (a's sign through the bounds of a*b)
        P = @assume(a * b > 0, a in Reals, b in Reals)
        @test SymJ._signs(a, P) == 0x07 && SymJ._signs(a * b, P) == 0x04
        @test solve(x^2 ~ a * b, x; assume = P) == [-sqrt(a * b), sqrt(a * b)]
        # a product's zero split one factor at a time; linear conditions solved for a symbol
        @test SymJ._solved_for(Relation(:lt, a - b, Sym(0))) in (a ≪ b, b ≫ a)
        @test SymJ._zero_factor(a * b, Assumptions()) in (a, b)
        cs = solve(x^2 ~ 9, x; assume = @assume(x > a - b))
        @test length(cs) == 4 && cs(a => 0, b => 5) == [-3, 3] && cs(a => 5, b => 0) == Sym[]
        @test occursin("a - 3 < b < a + 3", sprint(show, MIME"text/plain"(), cs))
        @test_throws ArgumentError solve(x^2 ~ 4, x; assume = (x ≫ 1, x ≪ 0))   # contradictory
    end

    @testset "limit" begin
        @syms a b
        @test limit(sin(x) / x, x => 0) == 1 && limit((1 + 1 / x)^x, x => Inf) == E
        @test limit(1 / x, x => 0; dir = :+) == Sym(Inf) && limit(-x, x => Inf) == Sym(-Inf)
        @test limit(1 / x, x => 0; dir = :+) isa Sym && limit(x^2, x => Sym(Inf)) == Sym(Inf)
        @test limit(x * exp(1 / x), x => 0; dir = -) == 0
        @test_throws GruntzError limit(sin(x), x => Inf)
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

    @testset "refine and simplify" begin
        @syms a n
        R(e, as...) = refine(e; assume = Assumptions(as...))
        @test R(sqrt(x^2), x ≫ 0) == x && R(sqrt(x^2), x ∈ ℝ) == abs(x) && R(sqrt(x^2)) == sqrt(x^2)
        @test R(sqrt(x^4), x ∈ ℝ) == x^2 && R((x^3)^(1 // 3), x ≫ 0) == x
        @test R((x^3)^(1 // 3), x ∈ ℝ) == (x^3)^(1 // 3)          # a real cube root is not x^(1/3)'s branch
        @test R(abs(x), x ≪ 0) == -x && R(abs(2x * y), x ≫ 0) == 2x * abs(y) && R(sign(x - 1), x ≫ 2) == 1
        @test R(log(exp(x)), x ∈ ℝ) == x && R(log(exp(x))) == log(exp(x)) && R(log(x^2), x ≫ 0) == 2log(x)
        @test R((-1)^(2n), n ∈ ℤ) == 1 && R((-1)^(2n + 1), n ∈ ℤ) == -1
        @test R(sin(n * PI), n ∈ ℤ) == 0 && R(cos(n * PI), n ∈ ℤ) == (-1)^n
        @test R(sin(x + 2n * PI), n ∈ ℤ) == sin(x) && R(cos(x + n * PI), n ∈ ℤ) == (-1)^n * cos(x)
        @test R(floor(n), n ∈ ℤ) == n && R(floor(x + n), n ∈ ℤ) == floor(x) + n
        @test R(sqrt(exp(x)), x ∈ ℝ) == exp(x / 2) && R(x^2 + a, a ~ 3) == x^2 + 3
        @test R(x ≫ 0, x ≫ 1) === true && R(x ≫ 2, x ≫ 1) == (x ≫ 2) && R(x^2 ≧ 0, x ∈ ℝ) === true
        @test (@assuming x > 0 begin refine(sqrt(x^2)) end) == x
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

    @testset "real root isolation" begin
        RR, RI = SymJ.real_roots, SymJ.root_intervals
        @test RR(x^3 - 2x - 5, x) == [2.0945514815423265] && RR(x^2 + 1, x) == Float64[] && RR(Sym(5), x) == Float64[]
        @test RR((x - 1)^2 * (x^2 - 2), x) == [-sqrt(2), 1.0, 1.0, sqrt(2)]
        @test RR(prod(x - k for k in 1:20), x) == Float64.(1:20)                   # Wilkinson
        @test RR(x * (x - 1 // 2) * (x + 3 // 7), x) == [-3 / 7, 0.0, 0.5]
        @test RR(x^5 - x - 1, x; T = BigFloat)[1] ≈ big"1.167303978261418684256045899854842180720560371525489039140082449275651903429519"
        I = RI((x - 1)^2 * (x^2 - 2), x)
        @test length(I) == 3 && I[2] == (1, 1, 2) && all(t -> t[1] < t[2] || t[1] == t[2], I)
        @test all(t -> t[2] - t[1] < 1 // 100, RI(x^2 - 2, x; width = 1 // 100))
        # an exact root next to an interval of the same factor (its endpoint)
        rs = [-11 // 7, 0, 5 // 7, -27 // 7, -2, -23 // 7, 18 // 7]
        @test RR(prod(x - r for r in rs), x) == sort(Float64.(rs))
        @test_throws ArgumentError RR(x^2 - sqrt(Sym(2)), x)
        @test_throws ArgumentError RR(Sym(0), x)
        @test_throws ArgumentError RR(sin(x), x)
    end

    @testset "linear algebra" begin
        @syms a b c d λ
        LA = SymJ.LinearAlgebra
        # det: Bareiss for polynomials, elimination over rational functions otherwise
        @test LA.det([a b; c d]) == a * d - b * c && LA.det([x 1; 1 x]) == x^2 - 1 && LA.det(Sym[1 2; 3 4]) == -2
        @test LA.det(Sym[1 2; 2 4]) == 0 && LA.det(Matrix{Sym}(undef, 0, 0)) == 1
        @test LA.det([1/x 1; 1 x]) == 0 && LA.det([sqrt(Sym(2)) 1; 1 sqrt(Sym(2))]) == 1
        G = [Sym(Symbol("g", i, j)) for i in 1:5, j in 1:5]
        @test length(LA.det(G).args) == 120                                # all 5! terms
        @test LA.det([Sym(1) / (i + j - 1) for i in 1:5, j in 1:5]) == 1 // 266716800000   # Hilbert
        @test LA.det([x 1; 1 x] - λ * LA.I(2)) == expand((x - λ)^2 - 1)  # characteristic polynomial
        # inv: the adjugate for polynomial entries, Gauss-Jordan over rational functions else
        @test inv([a b; c d]) == [d -b; -c a] / (a * d - b * c) && inv(Sym[1 2; 3 4]) == [-2 1; 3//2 -1//2]
        @test inv([x 1; 1 x]) == [x -1; -1 x] / (x^2 - 1) && inv(Sym[2;;]) == Sym[1//2;;]
        @test_throws LA.SingularException inv(Sym[1 2; 2 4])
        @test_throws LA.SingularException inv([1/x 1; 1 x])
        G3 = G[1:3, 1:3]
        @test all(iszero, map(simplify, inv(G3) * G3 - LA.I))
        @test map(simplify, inv([sin(x) cos(x); -cos(x) sin(x)])) == [sin(x) -cos(x); cos(x) sin(x)]
        # lu: a LinearAlgebra.LU, first-nonzero pivoting
        F = LA.lu([a b; c d])
        @test F isa LA.LU && F.L == [1 0; c/a 1] && F.p == [1, 2] && LA.det(F) == a * d - b * c
        @test map(simplify, F.L * F.U - [a b; c d][F.p, :]) == zeros(Sym, 2, 2)
        P = LA.lu(Sym[0 1; 1 0])
        @test P.p == [2, 1] && P.U == Sym[1 0; 0 1]
        @test_throws LA.SingularException LA.lu(Sym[1 2; 2 4])
        @test LA.lu(Sym[1 2; 2 4]; check = false).info == 2
        # \ exactly
        @test [a b; c d] \ [x, y] == [(d * x - b * y) / (a * d - b * c), (a * y - c * x) / (a * d - b * c)]
        @test Sym[1 2; 3 4] \ Sym[5, 6] == [-4, 9 // 2]
        # matrices with I and with expressions as scalars
        @test [x 1; 1 x] - LA.I == [x-1 1; 1 x-1] && 2x * [x 1; 1 x] == [2x^2 2x; 2x 2x^2]
        @test [x 1; 1 x] / x == [1 1/x; 1/x 1]
        @test_throws ArgumentError λ * LA.I
        # rref, rank, nullspace (generic: symbolic pivots nonzero)
        @test SymJ.rref(Sym[1 2 3; 4 5 6; 7 8 9]) == Sym[1 0 -1; 0 1 2; 0 0 0] && LA.rank(Sym[1 2 3; 4 5 6; 7 8 9]) == 2
        @test LA.nullspace(Sym[1 2 3; 4 5 6; 7 8 9]) == Sym[1; -2; 1;;]
        @test LA.rank([a b; 2a 2b]) == 1 && LA.nullspace([a b; 2a 2b]) == [-b/a; 1;;] && SymJ.rref([a 1; 1 1]) == Sym[1 0; 0 1]
        # algebraic constants are numbers, not independent atoms: sqrt(2)^2 - 2 is zero
        s2 = sqrt(Sym(2))
        @test LA.rank([s2 1; 2 s2]) == 1 && LA.det([s2 1; 2 s2]) == 0 && LA.nullspace([s2 1; 2 s2]) == [-s2/2; 1;;]
        @test_throws LA.SingularException inv([s2 1; 2 s2])
        @test inv([s2 1; 1 s2]) == [s2 -1; -1 s2]
        # eigenvalues with multiplicity; eigenvectors; eigen for diagonalizable matrices
        @test LA.eigvals(Sym[2 1; 1 2]) == [1, 3] && LA.eigvals(Sym[1 1; 0 1]) == [1, 1]
        @test LA.eigvals(Sym[1 2 0; 2 1 0; 0 0 3]) == [-1, 3, 3] && Set(LA.eigvals([a b; b a])) == Set([a + b, a - b])
        @test LA.eigvals(Sym[0 -1; 1 0]) == [-IM, IM] && LA.eigvecs(Sym[3 0; 0 3]) == Sym[1 0; 0 1]
        @test size(LA.eigvecs(Sym[1 1; 0 1])) == (2, 1)                     # defective
        @test_throws ArgumentError LA.eigen(Sym[1 1; 0 1])
        for A in (Sym[2 1; 1 2], [a b; b a], Sym[2 0 0; 0 3 4; 0 4 9], [x 1; 1 x], Sym[1 2; 3 4], Sym[0 -1; 1 0])
            F = LA.eigen(A)
            @test F isa LA.Eigen && all(iszero, map(simplify, A * F.vectors - F.vectors * LA.Diagonal(F.values)))
        end
        @test SymJ.charpoly([a b; c a], x) == expand((x - a)^2 - b * c)
        # powers of im with the sign outside: im^3 = 1/im = -im
        @test IM^3 == -IM && 1 / IM == -IM && IM^-1 + IM == 0 && (-1)^(-1 // 3 * one(Sym)) == -(-1)^(2 // 3 * one(Sym))
        @test simplify(1 / (1 + IM)) == (1 - IM) / 2
        @test LA.eigvals([cos(x) -sin(x); sin(x) cos(x)]) == [cos(x) - IM * sin(x), cos(x) + IM * sin(x)]
        # QR: Gram-Schmidt, exact
        for A in (Sym[1 1; 1 0], Sym[3 1; 4 2], Sym[1 2; 3 4; 5 6])
            Q, R = LA.qr(A)
            @test all(iszero, map(simplify, Q * R - A)) && map(simplify, transpose(Q) * Q) == Matrix{Sym}(LA.I(size(A, 2)))
            @test all(i -> iszero(R[i, 1]), 2:size(R, 1))
        end
        @test LA.qr(Sym[3 1; 4 2]).Q == Sym[3//5 -4//5; 4//5 3//5] && LA.qr(Sym[3 1; 4 2]).R == Sym[5 11//5; 0 2//5]
        @test_throws ArgumentError LA.qr(Sym[1 2; 2 4])
        # SVD: exact, checked numerically where simplify can't (nested radicals)
        F = LA.svd(Sym[3 0; 4 5])
        @test F.S == [3sqrt(Sym(5)), sqrt(Sym(5))] && map(simplify, F.U * LA.Diagonal(F.S) * F.Vt) == Sym[3 0; 4 5]
        @test LA.svdvals(Sym[1 1; 0 0]) == [sqrt(Sym(2)), 0] && LA.svdvals(Sym[2 0; 0 -3]) == [3, 2]
        numf(M) = map(e -> ComplexF64(N(e)), M)
        for A in (Sym[1 2; 3 4], Sym[1 im; 0 1], Sym[2 1 0; 1 2 1], Sym[1 1; 1 1; 0 0], Sym[1 0 1; 0 1 1])
            F = LA.svd(A)
            U, S, Vt = numf(F.U), numf(F.S), numf(F.Vt)
            @test isapprox(U * LA.Diagonal(S) * Vt, numf(A); atol = 1e-12) && isapprox(U' * U, LA.I; atol = 1e-12) &&
                  isapprox(Vt * Vt', LA.I; atol = 1e-12) && isapprox(real.(S), LA.svdvals(numf(A)); atol = 1e-12)
        end
        F = LA.svd(Sym[1 1; 1 1; 0 0]; full = true)
        @test size(F.U) == (3, 3) && map(simplify, transpose(F.U) * F.U) == Sym[1 0 0; 0 1 0; 0 0 1]
        @test_throws ArgumentError LA.svd(Sym[1 2 3; 4 5 6; 7 8 10])          # a cubic: refused, not slow
        # symbolic singular values refined under the assumptions
        @syms a b
        @test LA.svdvals(Sym[a 0; 0 b]) == [sqrt(b^2), sqrt(a^2)]
        @test LA.svdvals(Sym[a 0; 0 b]; assume = @assume(a > 0, b > 0)) == [b, a]
        @test assuming(() -> LA.svdvals(Sym[a 0; 0 b]), @assume(a ∈ ℝ, b ∈ ℝ)) == [abs(b), abs(a)]
        @test (@assuming a > 0 begin LA.svd(Sym[a a; a a]).S end) == [2a, 0]
        # Jordan form
        for A in (Sym[1 1; 0 1], Sym[2 1 0; 0 2 0; 0 0 3], Sym[5 4 2 1; 0 1 -1 -1; -1 -1 3 0; 1 1 -1 2], Sym[2 1; 1 2])
            P, J = SymJ.jordan_form(A)
            @test all(iszero, map(simplify, P * J * inv(P) - A))
        end
        @test SymJ.jordan_form(Sym[5 4 2 1; 0 1 -1 -1; -1 -1 3 0; 1 1 -1 2])[2] == Sym[1 0 0 0; 0 2 0 0; 0 0 4 1; 0 0 0 4]
        @test LA.eigvals(Sym[2 0 0; 0 3 4; 0 4 9]) == [1, 2, 11]                 # constants sorted
        # matrix functions through the Jordan form
        @syms t n
        @test exp(t * Sym[0 1; -1 0]) == [cos(t) sin(t); -sin(t) cos(t)]
        @test exp(Sym[1 1; 0 1]) == [E E; 0 E] && exp(t * Sym[2 1; 0 2]) == [exp(2t) t*exp(2t); 0 exp(2t)]
        @test sin(Sym[0 1; 0 0]) == Sym[0 1; 0 0] && cos(Sym[0 1; 0 0]) == Sym[1 0; 0 1]
        Fib = Sym[1 1; 1 0]^n
        @test [round(Int, N(subs(Fib[1, 2], n => k))) for k in 1:10] == [1, 1, 2, 3, 5, 8, 13, 21, 34, 55]
        @test isapprox(N.(subs.(exp(t * Sym[1 2; 3 4]), t => 0.3)), exp(0.3 * [1 2; 3 4]); rtol = 1e-10)
    end

    @testset "summation" begin
        @syms k n r a
        Σ = SymJ.summation
        # polynomials (Faulhaber), presented factored
        @test Σ(k, k, 1, n) == n * (n + 1) / 2 && Σ(k^2, k, 1, n) == n * (n + 1) * (2n + 1) / 6
        @test Σ(k^3, k, 1, n) == n^2 * (n + 1)^2 / 4 && Σ(2k + 1, k, 0, n) == (n + 1)^2 && Σ(Sym(1), k, 1, n) == n
        @test Σ(k^2, k, 1, 10) == 385 && Σ(a, k, 1, n) == a * n
        # Gosper: geometric and polynomial-geometric, telescoping, factorials
        @test Σ(r^k, k, 0, n) == (r^(n + 1) - 1) / (r - 1) && Σ(1 / (k * (k + 1)), k, 1, n) == n / (n + 1)
        @test Σ(k * factorial(k), k, 0, n) == factorial(n + 1) - 1
        @test SymJ._isfn1(Σ(factorial(k), k, 0, n), :Sum) || Σ(factorial(k), k, 0, n).name === :Sum   # no closed form
        # every finite closed form agrees with direct summation
        for f in (k^2, k^3 + 2k, r^k, k * 2^k, 1 / (k * (k + 1)), k * factorial(k), k / ((k + 1) * (k + 2)), 1 / k, binomial(n, k))
            s = Σ(f, k, 1, n)
            v7 = subs(s, n => 7, r => 3)
            @test isapprox(N(v7), N(_add_many_test(subs(f, n => 7, r => 3), k, 7)); rtol = 1e-12)
        end
        # binomial theorem; rational functions: harmonic numbers, combined
        @test Σ(binomial(n, k), k, 0, n) == 2^n && Σ(binomial(n, k) * x^k, k, 0, n) == (1 + x)^n
        @test Σ(1 / k, k, 1, n) == SymJ.harmonic(n) && SymJ.harmonic(4) == 25 // 12
        @test Σ(k / ((k + 1) * (k + 2)), k, 1, n) == SymJ.harmonic(n + 1) + 2 / (n + 2) - 2
        # infinite series: limits, zeta values, convergence by cases
        @test Σ(1 / (k * (k + 1)), k, 1, Inf) == 1 && Σ(1 / k^2, k, 1, Inf) == PI^2 / 6 && Σ(1 / k^4, k, 1, Inf) == PI^4 / 90
        @test Σ(1 / (k^2 - 1), k, 2, Inf) == 3 // 4 && Σ(1 / k, k, 1, Inf) == Sym(Inf)
        @test Σ((Sym(1) / 2)^k, k, 0, Inf) == 2 && Σ((-Sym(1) / 2)^k, k, 0, Inf) == 2 // 3 && Σ(2^k, k, 0, Inf) == Sym(Inf)
        cs = Σ(r^k, k, 0, Inf)
        @test cs isa SymJ.Cases && length(cs) == 3 && cs(r => 1 // 2) == 2
        @test_throws ArgumentError cs(r => 2)
        @test Σ(r^k, k, 0, Inf; assume = @assume(0 < r < 1)) == 1 / (1 - r)
        @test Σ(k * x^k, k, 0, Inf; assume = @assume(-1 < x < 1)) == x / (x - 1)^2
        # hypergeometric series: exp, cos, sin, cosh
        @test Σ(1 / factorial(k), k, 0, Inf) == E && Σ(x^k / factorial(k), k, 0, Inf) == exp(x)
        @test Σ((-1)^k * x^(2k) / factorial(2k), k, 0, Inf) == cos(x) && Σ((-1)^k * x^(2k + 1) / factorial(2k + 1), k, 0, Inf) == sin(x)
        @test Σ(x^(2k) / factorial(2k), k, 0, Inf) == cosh(x)
        @test SymJ._isfn1(Σ(sin(k), k, 1, n), :Sum) || Σ(sin(k), k, 1, n).name === :Sum
        # Zeilberger: definite sums of hypergeometric terms with no antidifference
        @syms m
        @test Σ(binomial(n, k)^2, k, 0, n) == binomial(2n, n)
        @test Σ(k * binomial(n, k), k, 0, n) == n * 2^(n - 1) && Σ(k^2 * binomial(n, k), k, 0, n) == n * (n + 1) * 2^(n - 2)
        @test Σ(binomial(n, k) * binomial(m, n - k), k, 0, n) == binomial(m + n, m)            # Vandermonde
        @test Σ(binomial(n, k) * binomial(m, k), k, 0, n) == binomial(m + n, m)
        @test Σ(binomial(2k, k) * binomial(2n - 2k, n - k), k, 0, n) == 4^n
        @test Σ((-1)^k * binomial(2n, k)^2, k, 0, 2n) == (-1)^n * binomial(2n, n)
        @test Σ(binomial(n + k, k) / 2^k, k, 0, n) == 2^n
        @test Σ(binomial(n, k)^2, k, 1, n) == binomial(2n, n) - 1                             # from 1
        # a second-order recurrence (Franel numbers), and bounds that cut the summand off
        @test Σ(binomial(n, k)^3, k, 0, n).name === :Sum                # Franel: no hypergeometric closed form
        @test Σ(binomial(n, k)^3, k, 0, n)(n => 3) == 56
        # a second-order recurrence with hypergeometric solutions (Hyper): Binet's formula
        fib = Σ(binomial(n - k, k), k, 0, n)
        @test all(m -> simplify(fib(n => m)) == [1, 1, 2, 3, 5, 8, 13, 21][m+1], 0:7)
    end

    @testset "recurrences (Hyper)" begin
        @syms n y() a k
        R, Σ = SymJ.rsolve, SymJ.summation
        φ, ψ = (1 + sqrt(Sym(5))) / 2, (1 - sqrt(Sym(5))) / 2
        @test simplify(R(y(n + 2) ~ y(n + 1) + y(n), y(n), y(0) => 0, y(1) => 1) - (φ^n - ψ^n) / sqrt(Sym(5))) == 0
        @test R(y(n + 2) - 4y(n + 1) + 4y(n) ~ 0, y(n), y(0) => 1, y(1) => 4) == (n + 1) * 2^n   # a double root
        @test R(y(n + 1) ~ (n + 1) * y(n), y(n), y(0) => 1) == factorial(n)
        @test R((n + 1) * y(n + 1) ~ 2 * (2n + 1) * y(n), y(n), y(0) => 1) == binomial(2n, n)
        @test R(y(n) ~ a * y(n - 1), y(n)) == a^n * y(0)                                  # shifted; general
        @test R(y(n + 1) ~ 2y(n), y(n), y(1) => 6) == 3 * 2^n                                 # from y(1)
        # a hypergeometric right side
        @test R(y(n + 1) ~ 2y(n) + 1, y(n)) == 2^n * (y(0) + 1) - 1 && R(y(n + 1) ~ 2y(n) + 1, y(n), y(0) => 0) == 2^n - 1
        @test R(y(n + 1) ~ y(n) + n, y(n), y(0) => 0) == n * (n - 1) / 2
        @test R(y(n + 1) ~ 3y(n) + 2^n, y(n), y(0) => 0) == 3^n - 2^n
        # the internal search: solutions of (E - 2)^2, none for the Franel recurrence
        @test length(SymJ._hyper(Sym[4, -4, 1], n)) == 2
        @test isempty(SymJ._hyper(Sym[-8(n + 1)^2, -(7n^2 + 21n + 16), (n + 2)^2], n))
        @test_throws ArgumentError R((n + 2)^2 * y(n + 2) - (7n^2 + 21n + 16) * y(n + 1) - 8(n + 1)^2 * y(n) ~ 0, y(n), y(0) => 1, y(1) => 2)
        # a right side that is a sum of hypergeometric terms (one class at a time)
        @test R(y(n + 1) - y(n) ~ n * 2^n + 3^n, y(n), y(0) => 0) == 2^n * (n - 2) + 3^n / 2 + 3 // 2
        @test simplify(R(y(n + 2) - y(n) ~ 1 + (-1)^n, y(n), y(0) => 0, y(1) => 0) - (n / 2 + (-1)^n * n / 2)) == 0
        @test simplify(R(y(n + 1) - 2y(n) ~ n + 2^n, y(n), y(0) => 0) - (2^(n - 1) * n + 2^n - n - 1)) == 0
        # algebraic factors: A = n - sqrt(2) dividing n^2 - 2
        @test R(y(n + 2) - (2n + 1) * y(n + 1) + (n^2 - 2) * y(n) ~ 0, y(n), y(0) => 1, y(1) => -sqrt(Sym(2))) ==
              factorial(n - sqrt(Sym(2)) - 1) / factorial(-sqrt(Sym(2)) - 1)
        @test simplify(subs(R(y(n + 1) ~ (n^2 - 2) * y(n), y(n), y(0) => 1), n => 3)) == 4
        # definite sums with bounds other than the summand's own: inhomogeneous recurrences
        @test Σ(binomial(n, k), k, 0, n - 1) == 2^n - 1 && Σ(binomial(n, k)^2, k, 0, n - 1) == binomial(2n, n) - 1
        @test Σ(binomial(n, k) / (k + 1), k, 0, n) == (2 * 2^n - 1) / (n + 1)
        @test Σ(binomial(n, k) * k, k, 1, n - 1) == 2^(n - 1) * n - n
        @test Σ(binomial(2n, k), k, 0, n) == 4^n / 2 + binomial(2n, n) / 2
        @test Σ(binomial(n, k) / binomial(2n, k), k, 0, n) == (2n + 1) / (n + 1)
        @test binomial(n, n + 1) == 0 && binomial(n, n) == 1
        # d'Alembertian solutions: reduction of order, iterated sums
        @test R(y(n + 1) - y(n) ~ 1 / (n + 1), y(n), y(0) => 0) == SymJ.harmonic(n)
        @test R((n + 2) * y(n + 2) - (2n + 3) * y(n + 1) + (n + 1) * y(n) ~ 0, y(n)) == (y(1) - y(0)) * SymJ.harmonic(n) + y(0)
        @test R(y(n + 1) - y(n) ~ 1 / (n + 1)^2, y(n), y(0) => 0) == PI^2 / 6 - SymJ.apply_function(:polygamma, [Sym(1), n + 1])
        let s = R(y(n + 1) - (n + 1) * y(n) ~ 1, y(n), y(0) => 1)          # n! Σ 1/k!
            @test [subs(s, n => m) for m in 0:4] == [1, 2, 5, 16, 65]
        end
        @test R(y(n + 2) - 2y(n + 1) + y(n) ~ 1 / (n + 1), y(n), y(0) => 0, y(1) => 0) == n * SymJ.harmonic(n) - n
        # harmonic numbers summed by parts
        let H = SymJ.harmonic
            @test Σ(H(k), k, 1, n) == (n + 1) * H(n) - n
            @test Σ(k * H(k), k, 1, n) == n * (n + 1) * H(n) / 2 - n * (n - 1) / 4
            @test simplify(Σ(k^2 * H(k), k, 1, n) - (n * (n + 1) * (2n + 1) * H(n) / 6 - n * (n - 1) * (4n + 1) / 36)) == 0
            @test Σ(H(k + 1), k, 0, n) == (n + 2) * H(n + 1) - n - 1
            @test Σ(H(k) / ((k + 1) * (k + 2)), k, 0, n) == n / (n + 1) - H(n) / (n + 2)
            @test Σ(H(k), k, 1, 10) == Sym(55991 // 2520)
            @test SymJ.kind(Σ(2^k * H(k), k, 0, n)) == FN               # Σ 2^(k+1)/(k+1) does not close
        end
        @test Σ((-1)^(k + 1) * binomial(n, k) / k, k, 1, n) == SymJ.harmonic(n)
        @test Σ((-1)^k * binomial(n, k) / (k + 1)^2, k, 0, n) == SymJ.harmonic(n + 1) / (n + 1)
        @test Σ((-1)^k * binomial(4, k) / (k + 1)^2, k, 0, 4) == SymJ.harmonic(Sym(5)) / 5
        @test SymJ.kind(Σ(binomial(n, k)^3, k, 0, n)) == FN                 # Franel: no d'Alembertian solution
        @test_throws ArgumentError R((n + 2)^2 * y(n + 2) - (7n^2 + 21n + 16) * y(n + 1) - 8(n + 1)^2 * y(n) ~ 0, y(n))
        @test_throws ArgumentError R(y(n + 1) ~ y(n)^2, y(n))
        @test_throws ArgumentError R(y(n + 2) ~ y(n + 1) + y(n), y(n), y(1) => 1)
    end

    @testset "differential equations" begin
        @syms y()
        D = SymJ.dsolve
        C1, C2, C3 = Sym(:C1), Sym(:C2), Sym(:C3)
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

    @testset "Laplace transforms" begin
        @syms t s w y() u() v()
        H, δ = SymJ.heaviside, SymJ.dirac
        L(f) = SymJ.laplace_transform(f, t, s)
        IL(F) = SymJ.inverse_laplace_transform(F, s, t)
        C1, C2 = Sym(:C1), Sym(:C2)
        @test L(Sym(1)) == 1 / s && L(t^3) == 6 / s^4 && L(exp(2t)) == 1 / (s - 2) && L(sin(3t)) == 3 / (s^2 + 9)
        @test L(cos(w * t)) == s / (s^2 + w^2) && L(sinh(t)) == 1 / (s^2 - 1) && L(t * sin(t)) == 2s / (s^2 + 1)^2
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
        @test_throws ArgumentError IL(s^2)
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
        @test P(SymJ.summation(a^Sym(:k), Sym(:k), 0, Inf)) == P((a ≫ -1, a ≪ 1) => 1 / (1 - a))
        # refine and limits
        @test refine(P(x ≪ 0 => -x, true => x); assume = @assume(x > 1)) == x
        @test refine(H(x - 1); assume = @assume(x > 2)) == 1
        @test limit(f, x => 0; dir = :+) == 0 && limit(P(x ≪ 1 => x, true => 2 - x), x => 1) == 1
        @test limit(H(x), x => 0; dir = :-) == 0 && limit(f, x => Inf) == Sym(Inf)
        @test_throws GruntzError limit(H(x), x => 0)
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

    @testset "Latexify" begin
        @syms a
        L(v) = repr(MIME"text/latex"(), v)
        @test L(x^2 + 1) == raw"$x^{2} + 1$" && L(sqrt(Sym(2)) * x) == raw"$\sqrt{2} \cdot x$"
        @test L(exp(x)) == raw"$e^{x}$" && L(IM * x) == raw"$i \cdot x$"          # ℯ, im as letters
        @test L(x ~ 2) == raw"$x = 2$" && L(x ≫ 0) == raw"$x > 0$" && L(a ≦ 1) == raw"$a \leq 1$"
        @test L(Sym(2)^100) == raw"$1267650600228229401496703205376$"           # a big literal
        @test L(solve(x^2 ~ 4, x)) == raw"$\left[ -2, \; 2 \right]$" && L(Sym[]) == raw"$\left[ \; \right]$"
        @test startswith(L([x 1; y x]), raw"$\left[")                           # a matrix: an array
        @test L((x, x^2)) == raw"$\left( x, \; x^{2} \right)$"
        @test L(SymNumber(3 // 4)) == raw"$\frac{3}{4}$" && L(Poly(x^2 - 1, x)) == raw"$x^{2} - 1$"
        @test L(series(exp(x), x, 0, 3)) == raw"$1 + x + \frac{x^{2}}{2} + O\left( x^{3} \right)$"
        @test L(solve([x + y ~ 3, x - y ~ 1], [x, y])) ==
              raw"$\left[ \begin{array}{l} \left\{ x = 2, \; y = 1 \right\} \end{array} \right]$"
        @test L(solve(x^2 ≫ 1, x)) == raw"$\left( -\infty, -1 \right) \cup \left( 1, \infty \right)$"
        @test L(solve([x ≧ 0, y ≧ 0, x + y ≦ 1], [x, y])) ==
              raw"$y \in \left[ 0, 1 \right], \; x \in \left[ 0,  - y + 1 \right]$"
        @test L(solve(x^2 ~ a, x; assume = @assume(x > 0))) ==
              raw"$\begin{cases} \left[ \sqrt{a} \right] & a > 0 \\ \left[ \; \right] & \lnot \left( a > 0 \right) \end{cases}$"
        @test occursin(raw"-3 < a < 3", L(solve(x^2 ~ 9, x; assume = @assume(x > a))))
        @test latexify(x / 2) == latexify(:(x / 2))
        # Julia's numbers and mixed tuples are left to their own display
        @test !showable(MIME"text/latex"(), 1.5) && !showable(MIME"text/latex"(), (x, 1.0)) &&
              !showable(MIME"text/latex"(), [1, 2])
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

    @testset "ambiguities" begin
        @test isempty(Test.detect_ambiguities(SymJ; recursive = true))
    end

    include("extensions.jl")
end
