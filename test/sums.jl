# Sums and recurrences (included by runtests.jl)

@testset "summation" begin
    # reversed bounds by Karr's convention, as the closed forms: Σ_{a}^{b} = -Σ_{b+1}^{a-1}
    let k = Sym(:k), n = Sym(:n), x = Sym(:x)
        @test SymJ.summation(k, k, 5, 1) == -9
        @test SymJ.summation(k, k, 5, 4) == 0
        @test SymJ.summation(k, k, 5, 0) == subs(SymJ.summation(k, k, n, 0), n => 5)
        # a Sym infinity is a bound as Inf is (SymPy #27893)
        @test SymJ.summation(1 / (x^2 + 1), x, 0, Sym(Inf)) == SymJ.summation(1 / (x^2 + 1), x, 0, Inf)
        @test SymJ.summation(1 / (x^2 + 1), x, Sym(-Inf), 0) == SymJ.summation(1 / (x^2 + 1), x, 0, Inf)
    end
    @syms k n r a
    Σ = SymJ.summation
    # rational summands: shifts both integer and not (the γ of ψ cancels only together)
    @test Σ(1 / ((k + 1) * (2k + 1)), k, 1, Inf) == log(Sym(4)) - 1
    @test Σ(1 / ((k + 1) * (2k + 1)), k, 0, Inf) == log(Sym(4))
    @test abs(N(Σ(k / (k^3 + 1), k, 1, Inf)) - 1.1116439382244) < 1e-6
    # an irreducible denominator through its roots, digamma pairs by reflection (SymPy #29670)
    @test Σ(1 / (k^2 + 1), k, 1, Inf) == PI * coth(PI) / 2 - 1 // 2
    @test Σ(1 / (k^2 + k + 1), k, 0, Inf) == sqrt(Sym(3)) * PI * tanh(sqrt(Sym(3)) * PI / 2) / 3
    @test abs(N(subs(Σ(1 / (k^2 + a^2), k, 1, Inf), a => 2)) - (2π * coth(2π) - 1) / 8) < 1e-12
    # the index is bound: not substituted, and renamed when a value would be captured
    let y = Sym(:y), s = Σ(1 / (k^2 + y), k, 1, n)
        @test string(subs(s, y => k)) == "Sum(1/(k′^2 + k), k′, 1, n)"
        @test subs(s, k => 5) == s
        @test subs(subs(s, y => k), n => 2) == 1 / (k + 1) + 1 / (k + 4)
    end
    # the binomial theorem with a polynomial factor, exact for every n ≥ 0 (SymPy #29731)
    @test Σ(binomial(n, k) * (-1)^(n - k) * k, k, 0, n) == SymJ.piecewise((n ~ 1) => 1, true => 0)
    @test Σ(binomial(n, k) * (-1)^(n - k) * k / factorial(n), k, 0, n) == SymJ.piecewise((n ~ 1) => 1, true => 0)
    @test Σ(binomial(n, k) * (-1)^(n - k), k, 0, n) == Sym(0)^n
    @test Σ(binomial(n, k) * k * 3^k, k, 0, n) == 3n * 4^(n - 1)
    # a constant ratio whose power does not cancel: (-1)^k x^k is (-x)^k, x^(2k) is (x^2)^k
    let x = Sym(:x)
        @test Σ((-1)^k * x^k, k, 1, Inf; assume = (x ≫ 0, x ≪ 1)) == -x / (x + 1)
        @test Σ(x^(2k), k, 1, Inf; assume = (x ≫ 0, x ≪ 1)) == x^2 / (1 - x^2)
        cs = Σ((-1)^k * x^k, k, 1, Inf)
        @test cs isa SymJ.Cases
        @test cs(x => 1 // 2) == -1 // 3
    end
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
    # (a finite sum of sin(k): closed form, SymPy #12107)
    @test N(subs(Σ(sin(k), k, 1, n), n => 7)) ≈ sum(sin(j) for j in 1:7)
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
    @test R(y(n + 1) - y(n) ~ 1 / (n + 1)^2, y(n), y(0) => 0) == SymJ.harmonic(n, 2)
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
        # generalized harmonic numbers
        @test H(Sym(3), 2) == Sym(49 // 36) && H(n, 1) == H(n) && Σ(1 / k^3, k, 1, n) == H(n, 3)
        @test Σ(1 / (k + 1)^2, k, 0, n - 1) == H(n, 2)
        @test Σ(H(k) / k, k, 1, n) == (H(n)^2 + H(n, 2)) / 2
        @test Σ(H(k) / (k + 1), k, 0, n) == (H(n + 1)^2 - H(n + 1, 2)) / 2
        @test Σ(H(k, 2) / k^2, k, 1, n) == (H(n, 2)^2 + H(n, 4)) / 2
        @test Σ(H(k, 2), k, 1, n) == (n + 1) * H(n, 2) - H(n)
        @test simplify(Σ(k * H(k, 2), k, 1, n) - (n * (n + 1) * H(n, 2) / 2 - n / 2 + H(n) / 2)) == 0
        @test Σ(k / (k + 1), k, 0, n - 1) == n - H(n)                 # a polynomial part
        @test subs(Σ(H(k + 1) / k, k, 1, n), n => 6) == Σ(H(k + 1) / k, k, 1, 6)
        @test SymJ.kind(Σ(H(k) / k^2, k, 1, n)) == FN                 # an Euler sum
    end
    @test Σ((-1)^(k + 1) * binomial(n, k) / k, k, 1, n) == SymJ.harmonic(n)
    @test Σ((-1)^k * binomial(n, k) / (k + 1)^2, k, 0, n) == SymJ.harmonic(n + 1) / (n + 1)
    @test Σ((-1)^k * binomial(4, k) / (k + 1)^2, k, 0, 4) == SymJ.harmonic(Sym(5)) / 5
    @test SymJ.kind(Σ(binomial(n, k)^3, k, 0, n)) == FN                 # Franel: no d'Alembertian solution
    # no closed form: the recurrence itself (Franel), evaluated by unrolling
    let s = R((n + 2)^2 * y(n + 2) - (7n^2 + 21n + 16) * y(n + 1) - 8(n + 1)^2 * y(n) ~ 0, y(n), y(0) => 1, y(1) => 2)
        @test SymJ.kind(s) == FN && s.name === :DifferenceRoot
        @test [subs(s, n => m) for m in 0:6] == [sum(binomial(m, j)^3 for j in 0:m) for m in 0:6]
        @test string(s) == "DifferenceRoot((n + 2)^2*y(n + 2) - (7*n^2 + 21*n + 16)*y(n + 1) - 8*(n + 1)^2*y(n) ~ 0, y(0) => 1, y(1) => 2)"
        @test subs(subs(s, n => n + 1), n => 2) == 56
        @test subs(s, n => 200) == sum(binomial(big(200), j)^3 for j in 0:200)
    end
    let s = R(y(n + 2) ~ (n^2 + a) * y(n + 1) + y(n), y(n))              # free initial values, a parameter
        @test expand(subs(s, n => 3) - ((a^2 + a + 1) * y(1) + (a + 1) * y(0))) == 0
    end
    # Liouvillian solutions: interlaced, piecewise in mod(n, m)
    let s = R(y(n + 2) ~ (n + 1) * y(n), y(n), y(0) => 1, y(1) => 1)        # (n - 1)!!
        @test [subs(s, n => m) for m in 0:7] == [1, 1, 1, 2, 3, 8, 15, 48]
        @test SymJ.kind(s) == FN && s.name === :piecewise
    end
    let s = R(y(n + 2) ~ (n + 1) * y(n), y(n))
        @test [subs(s, n => m) for m in 2:5] == [y(0), 2y(1), 3y(0), 8y(1)]
    end
    let s = R(y(n + 3) ~ (n + 1) * y(n), y(n), y(0) => 1, y(1) => 1, y(2) => 1)
        @test [simplify(subs(s, n => m)) for m in 0:8] == [1, 1, 1, 1, 2, 3, 4, 10, 18]
    end
    let s = R(y(n + 2) ~ (n + 1) * y(n) + 1, y(n), y(0) => 0, y(1) => 0)  # sums in the pieces
        @test [subs(s, n => m) for m in 0:6] == [0, 0, 1, 1, 4, 5, 21]
    end
    @test Σ((-1)^k * binomial(n, k)^2, k, 0, n) ==
          SymJ.piecewise((mod(n, 2) ~ 0) => (-1)^(n / 2) * binomial(n, n / 2), true => 0)
    let s = Σ((-1)^k * binomial(n, k)^3, k, 0, n)                         # Dixon
        @test [subs(s, n => m) for m in 0:6] == [Σ((-1)^k * binomial(m, k)^3, k, 0, m) for m in 0:6]
        @test s == SymJ.piecewise((mod(n, 2) ~ 0) => (-1)^(n / 2) * factorial(3n / 2) / factorial(n / 2)^3, true => 0)
    end
    let g = SymJ._present_product(factorial(n - Sym(1 // 3)) * factorial(n - Sym(2 // 3)) / factorial(n))   # Gauss
        @test !has(g, factorial(n - Sym(1 // 3)))
        @test all(m -> subs(g, n => m) == subs(2PI * 3^(-3n - Sym(1 // 2)) * factorial(3n) / factorial(n)^2, n => m), 1:4)
    end
    @test_throws ArgumentError R(y(n + 1) ~ y(n)^2, y(n))
    @test_throws ArgumentError R(y(n + 2) ~ y(n + 1) + y(n), y(n), y(1) => 1)
end
