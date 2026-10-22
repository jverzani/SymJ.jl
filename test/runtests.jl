using SymJ
using SymJ: PI, E, IM, ADD, MUL, POW, FN, NUMBER, kind
using TermInterface
using Test

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
        @test f(x) == f(x)
        @test f(x) != g(x)
        @test f(x) != f(y)
        @test f(x) + f(x) == 2f(x)
        @test string(f(x, 2)) == "f(x, 2)"
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
        # no arguments: the value of a constant expression, as a Julia number
        @test_throws ArgumentError u()
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
end
