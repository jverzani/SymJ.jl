# Building, substituting, calling and printing expressions (included by runtests.jl)

@testset "@syms forms" begin
    @syms xs[1:10]
    @test xs isa Vector{Sym} && length(xs) == 10 && string(xs[10]) == "xs₁₀"
    @syms a0[0:3]
    @test string.(a0) == ["a0₀", "a0₁", "a0₂", "a0₃"] && a0[1] == Sym(Symbol("a0₀"))   # indices in the names only
    @syms M[1:2, 1:3]
    @test size(M) == (2, 3) && string(M[2, 3]) == "M₂_₃"
    @syms w[-1:1]
    @test string(w[1]) == "w₋₁"
    @syms p1, q1, h1()
    @test p1 == Sym(:p1) && q1 == Sym(:q1) && h1 isa SymFunction
    @syms fs[1:3]()
    @test fs[2] isa SymFunction && string(fs[2](x)) == "fs₂(x)"
    m = 3
    @syms c3[1:m]                                           # a range from a variable
    @test length(c3) == 3
    @test (@syms p2 q2[1:2]) == (Sym(:p2), [Sym(Symbol("q2₁")), Sym(Symbol("q2₂"))])
    @test_throws ArgumentError @eval @syms bad[1.0:2.0]
    # free_symbols in the natural order of the names
    @syms zs[1:12] ws[-2:2] x2 x10
    @test free_symbols(sum(zs)) == zs
    @test free_symbols(sum(ws)) == ws
    @test free_symbols(x10 + x2 + x) == [x, x2, x10]
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

@testset "subs inside products" begin
    @test subs(x * exp(x), exp(x) => y) == x * y
    @test subs(3x^2 * y, x^2 => z) == 3y * z
end

@testset "subs of compound targets" begin
    @syms a b w
    # powers: an integer multiple of the exponent, and exponents with terms
    @test subs(x^4, x^2 => w) == w^2
    @test subs(1 / x^2, x^2 => w) == 1 / w
    @test subs(x^3, x^2 => w) == x^3                       # 3/2 is not an integer
    @test subs(y * x^(3 // 2), sqrt(x) => w) == w^3 * y
    @test subs(exp(x), exp(-x) => w) == 1 / w
    @test subs(exp(2x), exp(x) => w) == w^2
    @test subs(exp(x + 1), exp(x) => w) == E * w
    @test subs(exp(2a * x), exp(a * x) => w) == w^2
    @test subs(exp(x / 2), exp(x) => w) == exp(x / 2)
    @test subs(exp(x - exp(-x)) - exp(x), exp(-x) => w) == exp(-w) / w - 1 / w
    # sums: every term of the target, the same number of times
    @test subs(a + b + x, a + b => w) == w + x
    @test subs(2a + 2b, a + b => w) == 2w
    @test subs(3 + 2log(x), 1 + log(x) => w) == 2w + 1
    @test subs(2a + b, a + b => w) == 2a + b
    @test subs(sin(a + b + x), a + b => w) == sin(w + x)
    # products: an integer power of the target
    @test subs(a * b * x, a * b => w) == w * x
    @test subs(a * b^2 * x, a * b => w) == b * w * x
    @test subs(y / (a * b), a * b => w) == y / w
    @test subs(x^2 * log(x)^2, x * log(x) => w) == w^2
    ex = (1 + log(x)) * sqrt(1 + (x * log(x))^2)
    @test subs(ex, x * log(x) => w, 1 + log(x) => z) == z * sqrt(1 + w^2)
    @test (exp(2x) + 1)(exp(x) => w) == w^2 + 1                 # calling
end

@testset "subs and calling" begin
    # what a substitution does not reach is kept as it is (not rebuilt)
    let big = sum(sin(k * x + y) * exp(z / k) for k in 1:20) + log(x + 2)
        @test subs(big, sin(z) => 1) === big
        @test subs(big, x => 0) == sum(sin(y) * exp(z / k) for k in 1:20) + log(Sym(2))
        @test subs(x^2 * sin(y), x^2 => z) == z * sin(y)                    # a power factor matched whole
    end
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

@testset "callable / Function" begin
    @test x isa Function
    @test f isa Function
    @test (x .+ [1, 2]) == [x + 1, x + 2]
end

@testset "expand" begin
    @test expand((x + 1)^2) == x^2 + 2x + 1
    @test expand((x + y) * (x - y)) == x^2 - y^2
    @test expand(x * (y + z)) == x * y + x * z
    @test expand((sqrt(Sym(2)) + 1)^2) == 3 + 2sqrt(Sym(2))
    @test length(expand((x + y + z + 1)^10).args) == 285
    @test expand(sin((x + 1)^2)) == sin(x^2 + 2x + 1)
end

@testset "nsimplify" begin
    NS = SymJ.nsimplify
    @test NS(0.75) == Sym(3 // 4)
    @test NS(1 / 3) == Sym(1 // 3)
    @test NS(0.1234) == Sym(617 // 5000)
    @test NS(2.0) == 2
    @test NS(sqrt(2)) == sqrt(Sym(2))
    @test NS(-sqrt(7) / 3) == -sqrt(Sym(7)) / 3
    @test NS((1 + sqrt(5)) / 2) == (1 + sqrt(Sym(5))) / 2
    @test NS(cbrt(2)) == Sym(2)^(1 // 3)
    @test NS(sqrt(2) + sqrt(3)) == sqrt(Sym(2)) + sqrt(Sym(3))         # degree 4, denested
    @test NS(0.5x + 1.4142135623730951) == x / 2 + sqrt(Sym(2))
    @test NS(2pi / 3; constants = [PI]) == 2PI / 3
    @test NS(1 + 2pi; constants = [PI]) == 1 + 2PI
    @test NS(pi / 4 + log(2) / 3; constants = [PI, log(Sym(2))]) == PI / 4 + log(Sym(2)) / 3
    @test NS(3 + 2.718281828459045im; constants = [E]) == 3 + E * im
    @test NS(1.414214; tol = 1e-6) == sqrt(Sym(2))
    @test NS(1.414214) == Sym(1.414214)                                  # too few digits for the default tol
    for v in (log(2), float(pi), exp(1), 0.7853981633974483, 2.288989899226962)
        r = NS(v)                                                        # left alone
        @test SymJ.isnumber(r)
        @test SymJ.isfloat(r.num)
    end
    @test NS(float(pi); constants = [PI]) == PI
    @test NS(x + y) == x + y
    @test_throws ArgumentError NS(0.5; constants = [x])
    @test_throws ArgumentError NS(0.5; tol = 2)
end

@testset "gather, coeff, rationalize" begin
    @syms a b c t w()
    G, C = SymJ.gather, SymJ.coeff
    @test G(a * x^2 + b * x^2 + a * x - b * x + c, x) == (a + b) * x^2 + (a - b) * x + c
    @test G(x * y + a * x * y + x + b * y, [x, y]) == (a + 1) * x * y + b * y + x
    @test G(a * (x + 1)^2 + sin(x) + x * sin(x), x) == a * x^2 + 2a * x + a + sin(x) + x * sin(x)
    @test G(a / x + b / x + c * sqrt(x) + sqrt(x), x) == (a + b) / x + (c + 1) * sqrt(x)
    @test G(a * w(t)^2 + b * w(t)^2 + w(t), w(t)) == (a + b) * w(t)^2 + w(t)
    @test G(a * x^c + b * x^c, x) == (a + b) * x^c
    @test_throws ArgumentError G(x + 1, x + 1)
    @test C(a * x^2 + b * x + c, x, 2) == a && C(a * x^2 + b * x + c, x) == b
    @test C((x + y)^3, x^2) == 3y && C((x + y)^3, x, 0) == y^3
    @test C(a * x + sin(x) + b + 3, x, 0) == b + 3
    @test C(a / x + 2, x, -1) == a && C(a * x^2, x) == 0 && C(3x, x) == 3
    @test rationalize(0.5x + 1.25) == x / 2 + Sym(5 // 4)
    @test rationalize(x^0.5 + sin(0.1x)) == sqrt(x) + sin(x / 10)
    @test rationalize(Sym(3.14159) * x; tol = 1e-3) == Sym(201 // 64) * x
    @test rationalize(Sym(1.5) + 2.5im * x) == Sym(3 // 2) + Sym(5 // 2) * IM * x
    @test rationalize(Sym(1 // 3) * x) == x / 3
    @test rationalize(BigInt, Sym(0.1) * x) == x / 10
end

@testset "printing" begin
    # a float coefficient 1.0 is shown: 1.0*k is not k (SymPy #27198)
    let k = Sym(:k), x = Sym(:x), y = Sym(:y)
        @test string(Sym(1.0) * k) == "1.0*k"
        @test string(-Sym(1.0) * k) == "-1.0*k"
        @test string((1.0x - y)^2 - (x - y)^2) == "(1.0*x - y)^2 - (x - y)^2"
    end
    # NaN prints (its sign is not an integer)
    @test string(Sym(NaN)) == "NaN"
    @test string(Sym(Inf) - Sym(Inf)) == "NaN"
    @test string(x + Sym(NaN)) == "x + NaN"
    # division by zero: SymJ.DivisionByZero, with a message
    @test_throws SymJ.DivisionByZero Sym(1) / Sym(0)
    @test_throws SymJ.DivisionByZero subs(1 / x, x => 0)
    @test_throws SymJ.DivisionByZero Sym(0)^(-1 // 2)
    @test occursin("division by zero", sprint(showerror, SymJ.DivisionByZero("division by zero")))
    @test string(x^2 + 2x + 1) == "x^2 + 2*x + 1"
    @test string(x - y) == "x - y"
    @test string(-x / (2y)) == "-x/(2*y)"
    @test string(x^-2) == "1/x^2"
    @test string(sqrt(Sym(8))) == "2*sqrt(2)"
    # a float exponent prints as a float, not as sqrt (SymPy #29566)
    @test string(x^Sym(0.5)) == "x^0.5"
    @test string(x^Sym(0.5) - sqrt(x)) == "-sqrt(x) + x^0.5"
    @test string(Sym(2)^(1 // 3)) == "2^(1/3)"
    @test string(exp(x)) == "exp(x)"
    @test string(2IM) == "2*im"
    @test string((x + y)^2) == "(x + y)^2"
    @test string(PI) == "π"
    @test sprint(show, MIME"text/plain"(), x + 1) == "x + 1"
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
