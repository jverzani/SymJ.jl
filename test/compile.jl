# Compiling expressions (included by runtests.jl)

@testset "horner, cse" begin
    @syms a t1
    H, CSE = SymJ.horner, SymJ.cse
    @test H(2x^3 - x + 5, x) == (2x^2 - 1) * x + 5
    @test H(x^2 * y + x * y^2 + x + 1) == x * (x * y + y^2 + 1) + 1
    @test H(x^5 + 1, x) == x^5 + 1 && H(Sym(7), x) == 7
    @test H(a * x^2 + 2a * x + 3, x) == (a * x + 2a) * x + 3
    let p = 3x^4 - 2x^3 + x - 7, h = H(p, x)
        @test h == ((3x - 2) * x^2 + 1) * x - 7 && expand(h - p) == 0
    end
    @test_throws ArgumentError H(sin(x) + x, x)
    @test CSE((x + y)^2 + sin(x + y)) == ([Sym(:t1) => x + y], Sym(:t1)^2 + sin(Sym(:t1)))
    let (r, e) = CSE([exp(x + y) * sin(x + y), cos(x + y)^2 + exp(x + y)])
        @test length(r) == 2 && r[2].second == exp(r[1].first)       # exp(u) shared as a factor too
        back(u) = subs(subs(u, r[2]), r[1])                              # the replacements undone
        @test e isa Vector && back.(e) == [exp(x + y) * sin(x + y), cos(x + y)^2 + exp(x + y)]
    end
    @test first(CSE((t1 + x)^2 + cos(t1 + x)))[1].first == Sym(:t2)   # names avoid the symbols
    @test CSE(x + y) == (Pair{Sym,Sym}[], x + y)
    # powers of a root of unity reduced first: ω² = 1/ω for ω³ = 1 (SymPy #28051)
    let ω = IM * sqrt(Sym(3)) / 2 - Sym(1) // 2, es = [ω^2 * x, ω^-2 * y, ω * x, ω^-1 * y]
        @test SymJ._unity_powers(ω^2) == 1 / ω
        @test SymJ._unity_powers(ω^-2) == ω
        r, e = CSE(es)
        @test length(r) == 2                                      # ω and 1/ω
        @test all(d -> abs(N(subs(d, x => 1.5, y => 2.5))) < 1e-12,
                 foldl((acc, q) -> subs.(acc, Ref(q)), reverse(r); init = e) .- es)
    end
    let (r, e) = CSE(sqrt(x^2 + y^2) / (x^2 + y^2) + log(x^2 + y^2); prefix = :s)
        @test r == [Sym(:s1) => x^2 + y^2] && e == log(Sym(:s1)) + 1 / sqrt(Sym(:s1))
    end
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
