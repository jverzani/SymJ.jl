# SpecialFunctions, HypergeometricFunctions and PolyLog (dependencies) and the
# RischIntegrate extension.
import SpecialFunctions as SF
import RischIntegrate
import HypergeometricFunctions as HF
import PolyLog as PL

@testset "extensions" begin
    @syms x n a b c
    PI = SymJ.PI

    # ambiguities introduced by the methods SymJ adds to a package's functions (in a
    # submodule for a dependency, the extension for RischIntegrate)
    for (pkg, E) in ((SF, SymJ.SpecialFunctionsSupport), (HF, SymJ.HypergeometricSupport), (PL, SymJ.PolyLogSupport),
                     (RischIntegrate, Base.get_extension(SymJ, :SymJRischIntegrateExt)))
        @test E isa Module
        amb = filter(p -> p[1].module === E || p[2].module === E, Test.detect_ambiguities(pkg, E))
        @test isempty(amb)
    end

    # central difference of f at v
    fd(f, v; h = 1e-6) = (f(v + h) - f(v - h)) / (2h)

    @testset "SpecialFunctions" begin
        # spherical Bessel derivatives: f_ν' = f_{ν-1} - (ν + 1)/x f_ν (SymPy #26958)
        @test diff(SF.sphericalbesselj(2, x), x) == SF.sphericalbesselj(1, x) - 3SF.sphericalbesselj(2, x) / x
        @test lambdify(diff(SF.sphericalbessely(2, x), x), x)(1.3) ≈ (SF.sphericalbessely(2, 1.3 + 1e-6) - SF.sphericalbessely(2, 1.3 - 1e-6)) / 2e-6 rtol = 1e-6
        # (function, arguments before x, point): symbolic value at the point == the package's
        cases = [(SF.gamma, (), 0.4), (SF.gamma, (1.5,), 0.4), (SF.loggamma, (), 0.4), (SF.digamma, (), 0.4),
                 (SF.trigamma, (), 0.4), (SF.invdigamma, (), 0.4), (SF.polygamma, (1,), 0.4), (SF.beta, (1.5,), 0.4),
                 (SF.logbeta, (1.5,), 0.4), (SF.logfactorial, (), 3), (SF.zeta, (), 2.5), (SF.zeta, (2.5,), 0.4),
                 (SF.eta, (), 2.5), (SF.erf, (), 0.4), (SF.erf, (0.1,), 0.4), (SF.erfc, (), 0.4), (SF.erfcx, (), 0.4),
                 (SF.erfi, (), 0.4), (SF.erfinv, (), 0.4), (SF.erfcinv, (), 0.4), (SF.dawson, (), 0.4),
                 (SF.logerf, (0.1,), 0.4), (SF.logerfc, (), 0.4), (SF.logerfcx, (), 0.4), (SF.faddeeva, (), 0.4),
                 (SF.expint, (), 0.4), (SF.expint, (2,), 0.4), (SF.expintx, (), 0.4), (SF.expinti, (), 0.4),
                 (SF.sinint, (), 0.4), (SF.cosint, (), 0.4), (SF.airyai, (), 0.4), (SF.airyaiprime, (), 0.4),
                 (SF.airybi, (), 0.4), (SF.airybiprime, (), 0.4), (SF.airyaix, (), 0.4), (SF.airyaiprimex, (), 0.4),
                 (SF.airybix, (), 0.4), (SF.airybiprimex, (), 0.4), (SF.besselj, (1.5,), 0.4), (SF.bessely, (1.5,), 0.4),
                 (SF.besseli, (1.5,), 0.4), (SF.besselk, (1.5,), 0.4), (SF.besseljx, (1.5,), 0.4),
                 (SF.besselyx, (1.5,), 0.4), (SF.besselix, (1.5,), 0.4), (SF.besselkx, (1.5,), 0.4),
                 (SF.besselj0, (), 0.4), (SF.besselj1, (), 0.4), (SF.bessely0, (), 0.4), (SF.bessely1, (), 0.4),
                 (SF.besselh, (1.5, 1), 0.4), (SF.besselhx, (1.5, 2), 0.4), (SF.hankelh1, (1.5,), 0.4),
                 (SF.hankelh2, (1.5,), 0.4), (SF.hankelh1x, (1.5,), 0.4), (SF.hankelh2x, (1.5,), 0.4),
                 (SF.sphericalbesselj, (1,), 0.4), (SF.sphericalbessely, (1,), 0.4), (SF.jinc, (), 0.4),
                 (SF.ellipk, (), 0.4), (SF.ellipe, (), 0.4)]
        for (f, pre, v) in cases
            ex = f(pre..., x)
            @test kind(ex) == FN
            @test ex(v) ≈ f(pre..., v)                     # calling the expression
            @test lambdify(ex)(v) ≈ f(pre..., v)           # compiled
            v isa AbstractFloat && @test N(f(pre..., Sym(v))) ≈ f(pre..., v) || isnan(f(pre..., v))
        end
        @test SF.sinint(Sym(2.0)) == Sym(SF.sinint(2.0))  # a float argument evaluates
        @test string(SF.besselj(n, x)) == "besselj(n, x)"

        # derivatives against central differences
        for (f, pre, v) in [(SF.gamma, (), 1.3), (SF.gamma, (1.5,), 0.7), (SF.loggamma, (), 1.3), (SF.digamma, (), 1.3),
                            (SF.trigamma, (), 1.3), (SF.invdigamma, (), 0.4), (SF.polygamma, (1,), 1.3),
                            (SF.beta, (1.5,), 0.7), (SF.logbeta, (1.5,), 0.7), (SF.zeta, (2.5,), 0.7),
                            (SF.erf, (), 0.4), (SF.erf, (0.1,), 0.4), (SF.erfc, (), 0.4), (SF.erfcx, (), 0.4),
                            (SF.erfi, (), 0.4), (SF.erfinv, (), 0.4), (SF.erfcinv, (), 0.4), (SF.dawson, (), 0.4),
                            (SF.logerf, (0.1,), 0.4), (SF.logerfc, (), 0.4), (SF.logerfcx, (), 0.4),
                            (SF.expint, (), 0.4), (SF.expint, (2,), 0.4), (SF.expinti, (), 0.4), (SF.sinint, (), 0.4),
                            (SF.cosint, (), 0.4), (SF.airyai, (), 0.4), (SF.airyaiprime, (), 0.4), (SF.airybi, (), 0.4),
                            (SF.airybiprime, (), 0.4), (SF.besselj, (1.5,), 0.4), (SF.bessely, (1.5,), 0.4),
                            (SF.besseli, (1.5,), 0.4), (SF.besselk, (1.5,), 0.4), (SF.hankelh1, (1.5,), 0.4),
                            (SF.besselj0, (), 0.4), (SF.besselj1, (), 0.4), (SF.bessely0, (), 0.4),
                            (SF.bessely1, (), 0.4), (SF.ellipk, (), 0.4), (SF.ellipe, (), 0.4)]
            d = diff(f(pre..., x), x)
            @test d(v) ≈ fd(t -> f(pre..., t), v) rtol = 1e-6
        end
        @test diff(SF.besselj(n, x), x) == (SF.besselj(n - 1, x) - SF.besselj(n + 1, x)) / 2
        @test kind(diff(SF.besselj(n, x), n)) == FN       # no closed form in the order: unevaluated
        @test diff(SF.sinint(x), x) == sin(x) / x

        # exact values and symmetry
        @test SF.gamma(Sym(5)) == 24 && SF.gamma(Sym(1 // 2)) == sqrt(PI) && SF.gamma(Sym(7 // 2)) == 15sqrt(PI) / 8
        # (a pole: an error, as tan(π/2) is; SymPy #23156)
        @test_throws SymJ.DivisionByZero SF.gamma(Sym(-1))
        @test kind(SF.gamma(x)) == FN
        @test SF.zeta(Sym(2)) == PI^2 / 6 && SF.zeta(Sym(0)) == Sym(-1 // 2)
        @test SF.beta(Sym(2), Sym(3)) == Sym(1 // 12)
        @test SF.erf(Sym(0)) == 0 && SF.erfc(Sym(0)) == 1 && SF.erf(-x) == -SF.erf(x) && SF.erfi(-x) == -SF.erfi(x)
        @test SF.sinint(Sym(0)) == 0 && SF.sinint(-x) == -SF.sinint(x)
        @test SF.besselj(0, Sym(0)) == 1 && SF.besselj(2, Sym(0)) == 0 && SF.besselj0(Sym(0)) == 1
        @test SF.ellipk(Sym(0)) == PI / 2 && SF.ellipe(Sym(0)) == PI / 2

        # rewriting through gamma and factorials
        @test SymJ.rewrite(factorial(n) + SF.beta(a, b), SF.gamma) == SF.gamma(n + 1) + SF.gamma(a) * SF.gamma(b) / SF.gamma(a + b)
        @test SymJ.rewrite(SF.gamma(x) + SF.beta(a, b), factorial) == factorial(x - 1) + factorial(a - 1) * factorial(b - 1) / factorial(a + b - 1)

        # gamma simplification: integer shifts, reflection, duplication
        GS = SymJ.simplify_gamma
        @test GS(SF.gamma(x + 1) / SF.gamma(x)) == x && GS(SF.gamma(x) / SF.gamma(x - 2)) == (x - 1) * (x - 2)
        @test GS(SF.gamma(x + 2) / (x * SF.gamma(x))) == x + 1
        @test GS(SF.gamma(x) * SF.gamma(1 - x)) == PI / sin(PI * x) && GS(SF.gamma(x + 1) * SF.gamma(1 - x)) == PI * x / sin(PI * x)
        @test GS(SF.gamma(x + Sym(1 // 2)) * SF.gamma(Sym(1 // 2) - x)) == PI / cos(PI * x)
        @test GS(SF.gamma(x) * SF.gamma(x + Sym(1 // 2))) == 2^(1 - 2x) * sqrt(PI) * SF.gamma(2x)
        @test GS(SF.gamma(x) + SF.gamma(1 - x)) == SF.gamma(x) + SF.gamma(1 - x)          # not a product
    end

    @testset "HypergeometricFunctions" begin
        @test HF._₂F₁(a, b, c, Sym(0)) == 1 && kind(HF._₂F₁(a, b, c, x)) == FN
        @test HF.pFq((a, b), (c,), x) == HF._₂F₁(a, b, c, x) && HF.pFq((a,), (b,), x) == HF._₁F₁(a, b, x)
        @test HF.pFq((a, b, c), (n, 2), x) == HF._₃F₂(a, b, c, n, 2, x)
        for (f, ps, v) in ((HF._₁F₁, (0.5, 1.5), 0.3), (HF._₂F₁, (0.5, 1, 1.5), 0.3), (HF._₃F₂, (0.5, 1, 1.5, 2, 2.5), 0.3))
            @test f(ps..., x)(v) ≈ f(ps..., v)
            @test diff(f(ps..., x), x)(v) ≈ fd(t -> f(ps..., t), v) rtol = 1e-6
        end
        @test kind(diff(HF._₂F₁(a, b, c, x), a)) == FN
        @test HF._₂F₁(a, b, c, x)(a => 0.5, b => 1, c => 1.5)(0.3) ≈ HF._₂F₁(0.5, 1, 1.5, 0.3)
    end

    @testset "PolyLog" begin
        @test PL.li2(Sym(1)) == PI^2 / 6 && PL.li2(Sym(-1)) == -PI^2 / 12 && PL.li3(Sym(0)) == 0
        @test PL.li(n, Sym(0)) == 0 && kind(PL.li(n, x)) == FN
        for f in (PL.li0, PL.li1, PL.li2, PL.li3, PL.li4, PL.li5, PL.li6, PL.reli1, PL.reli2, PL.reli3, PL.reli4)
            @test f(x)(0.3) ≈ f(0.3)
            @test diff(f(x), x)(0.3) ≈ fd(f, 0.3) rtol = 1e-6
        end
        @test PL.li(3, x)(0.3) ≈ PL.li(3, 0.3) && PL.reli(3, x)(0.3) ≈ PL.reli(3, 0.3)
        @test diff(PL.li(n, x), x) == PL.li(n - 1, x) / x
        @test diff(PL.li3(x), x) == PL.li2(x) / x
        # series at 0, and at the branch point z = 1 (not a Taylor series; SymPy #25159)
        @test Sym(series(PL.li2(x), x, 0, 4)) == x + x^2 / 4 + x^3 / 9
        @test Sym(series(PL.li(Sym(3) / 2, exp(-x)), x, 0, 2)) == SF.zeta(Sym(3) / 2) - 2sqrt(PI) * sqrt(x) - SF.zeta(Sym(1) / 2) * x
        @test Sym(series(PL.li1(1 - x), x, 0, 3)) == -log(x)
        @test Sym(series(PL.li(Sym(-1), exp(-x)), x, 0, 3)) == x^-2 - Sym(1) / 12 + x^2 / 240
        @test N(subs(Sym(series(PL.li(Sym(3) / 2, exp(-x)), x, 0, 5)), x => Sym(1) / 10)) ≈ sum(exp(-0.1k) / k^1.5 for k in 1:2000)
        @test N(subs(Sym(series(PL.li3(exp(-x)), x, 0, 5)), x => Sym(1) / 10)) ≈ PL.li3(exp(-0.1))
        # a non-integer order numerically: the series near 0, the expansion near 1
        @test N(PL.li(Sym(1) // 2, Sym(1) // 3)) ≈ sum((1 / 3)^k / sqrt(k) for k in 1:100)
        @test N(PL.li(Sym(1) // 2, Sym(9) // 10)) ≈ sum(0.9^k / sqrt(k) for k in 1:3000)
        @test SF.zeta(Sym(-3)) == Sym(1) / 120
        @test SF.gamma(Sym(-3) / 2) == 4sqrt(PI) / 3
    end

    # (the special functions are in RischIntegrate 1.2 on; Julia 1.10 resolves 1.0)
    isdefined(RischIntegrate, :clausen2) && @testset "RischIntegrate: special functions in answers" begin
        ∫ = RischIntegrate.integrate
        # Clausen functions and the like are expressions: diff and N work on the answer
        F = ∫(x / sin(x), x)
        @test F isa Sym
        @test abs(N((diff(F, x) - x / sin(x))(0.7))) < 1e-12
        @test !occursin("-1 *", string(F))
        G = ∫(x^3 * tan(x), x)
        @test abs(N((diff(G, x) - x^3 * tan(x))(1.3))) < 1e-12
        H = ∫(x / (5 + 4cos(x)), x)
        @test abs(N((diff(H, x) - x / (5 + 4cos(x)))(0.7))) < 1e-12
        @test diff(RischIntegrate.clausen2(x), x) == -log(2 - 2cos(x)) / 2
        @test diff(RischIntegrate.clausen(3, x), x) == -RischIntegrate.clausen2(x)
        @test diff(RischIntegrate.sinhint(x), x) == sinh(x) / x
        @test RischIntegrate.clausen2(x)(1.0) ≈ RischIntegrate.clausen2(1.0)
    end

    @testset "RischIntegrate: piecewise, heaviside, dirac" begin
        P, H, D = SymJ.piecewise, SymJ.heaviside, SymJ.dirac
        ∫ = RischIntegrate.integrate
        @test ∫(H(x), x) == P(x ≪ 0 => 0, true => x)
        @test ∫(x * H(x - 1), x) == P(x ≪ 1 => 0, true => x^2 / 2 - 1 // 2)              # continuous
        @test ∫(H(x) - H(x - 1), x) == P(x ≪ 0 => 0, x ≪ 1 => x, true => 1)
        @test ∫(P(x ≪ 0 => -x, true => x^2), (x, -1, 2)) == 19 // 6
        @test ∫(abs(x - 1), (x, 0, 3)) == 5 // 2 && ∫(abs(x), x) == P(x ≪ 0 => -x^2 / 2, true => x^2 / 2)
        @test ∫(H(x - a), x) == P(x ≪ a => 0, true => x - a)                               # a parameter
        @test ∫(cos(x) * D(x - 1), x) == cos(Sym(1)) * H(x - 1) && ∫(exp(x) * D(2x - 2), x) == E / 2 * H(x - 1)
        @test ∫(x^2 * D(x - 2, 1), x) == 4D(x - 2) - 4H(x - 2)                            # δ′: g(a)δ - g′(a)H (SymPy #17383)
        # δ(u(x)) = Σ δ(x - r)/|u′(r)| over the simple real roots of a polynomial u (SymPy #29001)
        @test ∫(D(x^2 - 1) * x^2, (x, -3, 3)) == 1
        @test ∫(D(x^2 - 1) * x^2, (x, 0, 3)) == 1 // 2
        @test ∫(D(x^3 - x), (x, -2, 2)) == 2
        @test ∫(D(x^2 + 1), (x, -3, 3)) == 0
        @test_throws ArgumentError ∫(D(x^2), (x, -3, 3))                 # a double root
        @test_throws ArgumentError ∫(D(sin(x)), (x, -10, 10))            # not a polynomial
        # float breakpoints (SymPy #29064)
        @test N(∫(P(x ≪ 0.5 => 0, true => 1), (x, 0, 1))) ≈ 0.5
        @test N(∫(P(x ≪ 0.25 => x, x ≪ 0.75 => 1, true => x^2), (x, 0, 1))) ≈ 1 / 32 + 1 / 2 + 37 / 192
        # complex pieces made continuous where the one-sided limit fails (SymPy #28190)
        # (complex constants need a recent RischIntegrate; Julia 1.10 resolves an older one)
        ri_complex = try ∫(exp(IM * x), x); true catch; false end
        ri_complex && let rect = P((-Sym(1) // 2 ≦ x, x ≦ Sym(1) // 2) => 1, true => 0)
            @test ∫(rect * exp(-2IM * x), (x, -Inf, Inf)) == sin(Sym(1))
            @test ∫(rect * exp(-2IM * PI * a * x), (x, -Inf, Inf)) == sin(PI * a) / (PI * a)   # generic a
            @test ∫(rect * exp(-2IM * PI * a * x), (x, -Inf, Inf); assume = (a ~ 0,)) == 1
        end
        # no piece holds past the last condition: not defined there (SymPy #28469)
        @test ∫(P(x ≪ 1 => x), x) == P(x ≪ 1 => x^2 / 2)
        let (a1, a2, a3, t) = Sym.((:a1, :a2, :a3, :t))
            @test ∫(P(x ≪ 1 => a1, x ≪ 2 => a2, x ≪ 3 => a3), (x, 0, t))(t => Sym(5) // 2) == a1 + a2 + a3 / 2
        end
        # a condition of one sign by bounds: |…| ≥ 1 (SymPy #27283)
        @test ∫(abs((3 + cos(3x)^5 / 2) - (1 + sin(2x)^3 / 2)), (x, 0, 2PI)) == 4PI
        # dirac in a product with a sum, abs in the rest (SymPy #27108)
        @test ∫((exp(-abs(x)) - 2D(x)) * exp(-abs(x)), (x, -Inf, Inf)) == -1
        # a periodic integrand: its antiderivative for every period, not one
        let F = ∫(abs(sin(x)), x)
            @test subs(F, x => 3PI) - subs(F, x => 0) == 6
            @test subs(F, x => 0) - subs(F, x => -2PI) == 4
        end
        @test_throws ArgumentError ∫(x * abs(sin(x)), x)              # no finite form: an error
        # x only in x + y: integrated in y as a shift (SymPy #27234)
        let G = ∫(abs(cos(x + y)), y; assume = (x ∈ ℝ, y ∈ ℝ))
            @test subs(G, x => 1, y => 2PI - 1) - subs(G, x => 1, y => -1) == 4
        end
        @test simplify(∫(abs(cos(x + y)), (y, 0, 2PI))) == 4
        # a condition quadratic in x with roots in y: by the discriminant (SymPy #27358)
        @test ∫(P(((x - 5)^2 + (y - 5)^2 ≦ 1) => 1, true => 0), (x, 0, 9), (y, 0, 9)) == PI
        @test ∫(P((x^2 ≧ y) => 1, true => 0), (x, 0, 1), (y, 0, 1)) == 1 // 3
    end

    @testset "RischIntegrate: definite integrals" begin
        ∫ = RischIntegrate.integrate
        Inf_ = Sym(Inf)
        # equal limits; two deltas over the line; a series term by term (SymPy #5098, #5588, #4831)
        let f = SymFunction(:f), a = Sym(:a), x1 = Sym(:x1), x2 = Sym(:x2)
            @test iszero(∫(f(x), (x, a, a)))
            # (sifted by either delta, by the factors' order, which differs between Julia versions)
            @test ∫(x2 * SymJ.dirac(x - x2) * SymJ.dirac(x2 - x1), (x2, -Inf, Inf)) in (x1 * SymJ.dirac(x1 - x), x * SymJ.dirac(x - x1))
        end
        @test Sym(∫(series(sin(x), x, 0, 6), x)) == x^2 / 2 - x^4 / 24 + x^6 / 720
        # the derivative of an answer with root sums (a MethodError; SymPy #5231)
        isdefined(RischIntegrate, :HypTerm) && @test iszero(SymJ.cancel(diff(∫((4x - 5) / (x^5 - x + 1), x), x) - (4x - 5) / (x^5 - x + 1)))
        # a variable that is not a symbol is refused (SymPy #4459); ∫₀^(π/2) cos(x)^(2n) for
        # an integer n > 0 by the limit at π/2, not 0 (SymPy #5804)
        @test_throws ArgumentError ∫(sin(x), sin(x))
        VERSION >= v"1.12" && let n = Sym(:n)                            # (the older RischIntegrate cannot)
            @test iszero(simplify(∫(cos(x)^(2n), (x, 0, PI / 2); assume = (n ∈ ℤ, n ≫ 0)) -
                                  sqrt(PI) * SymJ.SpecialFunctionsSupport.SF.gamma(n + Sym(1) // 2) / (2SymJ.SpecialFunctionsSupport.SF.gamma(n + 1))))
        end
        # Fourier-type integrals by residues (SymPy #6136, #6348); a Gaussian e^(-i t k²)
        # with t > 0 (SymPy #6123)
        let k = Sym(:k), a = Sym(:a), t = Sym(:t)
            @test ∫(exp(SymJ.IM * k) / (k^2 + 1), (k, -Inf, Inf)) == PI / E
            @test ∫(cos(x) / (1 + x^2), (x, -Inf, Inf)) == PI / E
            @test ∫(x * sin(x) / (1 + x^2), (x, -Inf, Inf)) == PI / E
            @test ∫(cos(k * x) / (k^2 + a^2), (k, -Inf, Inf); assume = (a ≫ 0, x ≫ 0)) == PI * exp(-a * x) / a
            @test N(subs(∫(exp(-SymJ.IM * t * k^2), (k, -Inf, Inf); assume = (t ≫ 0,)), t => 2)) ≈ sqrt(π / (2im))
        end
        # pieces whose conditions are free of x (SymPy #6747); breakpoints ±y ordered by
        # y > 0 (SymPy #7158)
        let y = Sym(:y)
            @test ∫(SymJ.piecewise((y ~ 0) => 1, true => x * y), x) == SymJ.piecewise((y ~ 0) => x, true => x^2 * y / 2)
            F = ∫(log(abs(x^2 - y^2)), x; assume = (y ≫ 0,))
            @test N(subs(diff(F, x), x => 1, y => 3)) ≈ log(8)
            @test N(subs(diff(F, x), x => 4, y => 3)) ≈ log(7)
            @test N(subs(∫(log(abs(x^2 - y^2)), (x, 0, 2y); assume = (y ≫ 0,)), y => 3)) ≈ 11.0708601 rtol = 1e-6
        end
        # floats in the limits: a float value (SymPy #7334)
        let t = Sym(:t)
            @test N(∫(∫(SymJ.piecewise(x / t ≪ 0 => x^2, true => 1 + x^2), (t, 0.8, 1)), (x, -1, 1))) ≈ 1 / 3
        end
        # y'² = x³ by its roots (SymPy #7093)
        let u = SymFunction(:u)
            @test SymJ.dsolve(diff(u(x), x)^2 ~ x^3, u(x)) isa AbstractVector
        end
        # u(x) δ(x - a), u unknown, sifted (SymPy #6968); f'g + fg' as (fg)' (SymPy #6958)
        let f = SymFunction(:f), g = SymFunction(:g), a = Sym(:a)
            @test ∫(f(x) * SymJ.dirac(x - a), (x, -Inf, Inf)) == f(a)
            @test ∫(diff(f(x) * g(x), x), x) == f(x) * g(x)
            @test ∫(diff(f(x) * g(x), x), (x, 0, 1)) == f(1) * g(1) - f(0) * g(0)
        end
        # a derivative at a point linear in the variable (SymPy #6407; a stack overflow)
        let f = SymFunction(:f), t = Sym(:t)
            df = subs(diff(f(x), x), x => x - t)
            @test ∫(df, t) == -f(x - t)
            @test ∫(df, (t, 0, x)) == f(x) - f(0)
            @test ∫(diff(f(2t), t), t) == f(2t)
        end
        # x sin(a/x) at 0 by a limit, not substituted (SymPy #6154; a DivisionByZero)
        let t = Sym(:t), a = Sym(:a), w = Sym(:w)
            VERSION >= v"1.12" && @test_throws RischIntegrate.NotElementary ∫(sin(a / t) * cos(w * t) / t, (t, 0, Inf); assume = (a ≫ 0, w ≫ 0))
        end
        # facts that fix a parameter's value in indefinite integrals (SymPy #6985)
        let k = Sym(:k), n = Sym(:n)
            @test ∫(exp(k * x), x; assume = (k ~ 0,)) == x
            @test ∫(x^n, x; assume = (n ~ -1,)) == log(x)
        end
        # a complex antiderivative on part of a long interval is caught (SymPy #7155)
        VERSION >= v"1.12" && @test_throws ArgumentError ∫(x * acos(1 / x), (x, Sym(1) // 10, 100))  # (the older RischIntegrate refuses acos(1/x))
        # floor and ceil integrated by unit intervals (SymPy #8420, #8480)
        @test ∫(1 / floor(x)^2, (x, 1, Inf)) == SymJ.PI^2 / 6
        @test ∫(x - floor(x), (x, 0, Inf)) == Inf_
        @test ∫(x * floor(x), (x, 0, 3)) == Sym(13) // 2
        # |x| in a product's exponent (SymPy #8849)
        VERSION >= v"1.12" && @test ∫(exp(-abs(x)) * cos(2x), (x, -Inf, Inf)) == Sym(2) // 5
        VERSION >= v"1.12" && @test iszero(simplify(∫(exp(-abs(x)) * x, (x, -1, 2)) - (2 / exp(Sym(1)) - 3 / exp(Sym(2)))))
        # the inner integral's cases at y = ±√(-h²) impossible for real y, h > 0 (SymPy #10211)
        # (the older RischIntegrate that Julia 1.10 resolves cannot do the inner integral)
        VERSION >= v"1.12" && let (y, h, w) = (Sym(:y), Sym(:h), Sym(:w))
            r = ∫(1 / sqrt((y - x)^2 + h^2)^3, (x, 0, w), (y, 0, w); assume = (h ≫ 0, w ≫ 0))
            @test iszero(simplify(r - 2(sqrt(h^2 + w^2) - h) / h^2))
        end
        # an answer whose only non-elementary part is a ₂F₁ term: an expression, so limit
        # and the rest take it (SymPy #24879)
        # ₂F₁ antiderivatives where their radicand is positive, limits through the
        # connection formulas (SymPy #24879)
        isdefined(RischIntegrate, :HypTerm) && @test N(∫(1 / sqrt(1 + x^3), (x, -1, Inf_))) ≈ N(SF.gamma(Sym(1) // 3) * sqrt(PI) / SF.gamma(Sym(5) // 6))
        isdefined(RischIntegrate, :HypTerm) && @test_throws ArgumentError ∫(x / sqrt(1 - x^3), (x, 0, 2))
        # Bose–Einstein and Fermi–Dirac integrals (SymPy #24431)
        @test ∫(x^2 / (exp(x) - 1), (x, 0, Inf_)) == 2SF.zeta(Sym(3))
        @test ∫(x^2 / (exp(x) + 1), (x, 0, Inf_)) == 3SF.zeta(Sym(3)) / 2
        # a definite integral refined under the assumptions (n integer: sin(nπ) = 0)
        let n = Sym(:n)
            cs = ∫(cos(n * x) * cos(x), (x, -PI, PI); assume = (n ∈ ℤ, n ≧ 0))
            @test last(cs.cases).second == 0
        end
        # (₂F₁ answers and the larger radicals are in RischIntegrate 1.2 on)
        isdefined(RischIntegrate, :HypTerm) && @test ∫(1 / sqrt(1 + x^3), x) isa Sym
        # simplify of an answer by cases (SymPy #23909 extra)
        let a = Sym(:a), b = Sym(:b), t = Sym(:t)
            @test simplify(∫(exp((b - a) * x - b * t), (x, 0, t))) isa Sym
        end
        # a large antiderivative is not expanded by simplify (#23967)
        isdefined(RischIntegrate, :HypTerm) && let s = Sym(:s), n1 = Sym(:n1)
            f = (sqrt(4s^2 + 1) / 2 + Sym(3) / 2) / (sqrt(4s^2 + 1) / 2 + Sym(1) / 2)^(Sym(3) // 2)
            @test N(subs(∫(f, (s, 0, n1)), n1 => 1)) ≈ sum(N(subs(f, s => (k - Sym(1) // 2) / 2000)) for k in 1:2000) / 2000 rtol = 1e-6
        end
        # exact with the RischIntegrate Julia 1.12 gets; the older one gives approximate constants
        same(v, w) = VERSION >= v"1.12" ? v == w : N(v) ≈ N(w)
        @test ∫(x^2, (x, 0, 3)) == 9
        @test ∫(sin(x), (x, 0, PI)) == 2
        # Gaussian integrals with a parameter (SymPy #28485)
        @test ∫(exp(-a * x^2), (x, -Inf, Inf); assume = (a ≫ 0,)) == sqrt(PI) / sqrt(a)
        @test ∫(x^2 * exp(-a * x^2), (x, 0, Inf); assume = (a ≫ 0,)) == sqrt(PI) / (4a^(3 // 2))
        @test iszero(simplify(∫(exp(-a * x^2 + b * x), (x, -Inf, Inf); assume = (a ≫ 0,)) - sqrt(PI / a) * exp(b^2 / (4a)); assume = (a ≫ 0,)))
        @test ∫(exp(-a * x^2), (x, -Inf, Inf))(a => -1) == Inf
        # erfc(u) = 1 - erf(u) (SymPy #28900, #28894)
        @test simplify(SF.erf(x) + SF.erfc(x)) == 1
        @test simplify(SF.erf(sqrt(Sym(2)) / 2) + SF.erfc(sqrt(Sym(2)) / 2)) == 1
        # no elementary antiderivative: Beta over [0, π/2] (SymPy #30664)
        @test N(∫(1 / sqrt(sin(x)), (x, 0, PI / 2))) ≈ 2.6220575542921198
        @test N(∫(sqrt(tan(x)), (x, 0, PI / 2))) ≈ pi / sqrt(2)
        @test N(∫(sin(x)^(1 // 3) * cos(x)^(2 // 3), (x, 0, PI / 2))) ≈ SF.gamma(2 / 3) * SF.gamma(5 / 6) / sqrt(pi)
        @test ∫(1 / (1 + x^2), (x, 0, Inf)) == PI / 2                   # a limit, not atan(Inf)
        @test ∫(1 / (1 + x^2), (x, -Inf, Inf)) == PI
        @test ∫(exp(-x), (x, 0, Inf)) == 1
        @test ∫(x * exp(-x), (x, 0, Inf)) == 1
        @test same(∫(exp(-x^2), (x, -Inf, Inf)), sqrt(PI))             # erf at ±∞
        @test ∫(sin(x) / x, (x, 1, Inf)) == PI / 2 - SF.sinint(Sym(1))   # sinint at ∞
        @test ∫(cos(x) / x, (x, 1, Inf)) == -SF.cosint(Sym(1))            # cosint at ∞
        @test same(∫(1 / (2 + cos(x)), (x, 0, 2PI)), 2sqrt(Sym(3)) * PI / 3)
        @test ∫(log(x), (x, 0, 1)) == -1                                 # a limit at 0
        @test ∫(1 / sqrt(x), (x, 0, 1)) == 2
        @test ∫(1 / x^2, (x, 1, Inf)) == 1
        @test same(∫(1 / (1 + x^4), (x, 0, Inf)), sqrt(Sym(2)) * PI / 4)
        @test ∫(1 / x, (x, -2, -1)) == log(Sym(1 // 2))
        # divergent: ±∞, or an error
        @test ∫(1 / x^2, (x, -1, 1)) == Inf_                             # was -2
        @test ∫(1 / x^2, (x, 1, -1)) == -Inf_
        @test ∫(1 / x, (x, 1, Inf)) == Inf_
        @test_throws ArgumentError ∫(1 / x, (x, -1, 1))
        @test_throws ArgumentError ∫(1 / x^3, (x, -1, 1))
        @test_throws ArgumentError ∫(sin(x), (x, 0, Inf))
        @test_throws ArgumentError ∫(tan(x), (x, 0, PI))
        # parameters: by cases, or decided by the assumptions
        cs = ∫(exp(-a * x), (x, 0, Inf))
        @test cs isa SymJ.Cases
        @test cs(a => 2) == Sym(1 // 2)
        @test cs(a => -1) == Inf_
        @test cs(a => 0) == Inf_
        @test ∫(exp(-a * x), (x, 0, Inf); assume = @assume(a > 0)) == 1 / a
        @test (@assuming a > 0 begin ∫(exp(-a * x), (x, 0, Inf)) end) == 1 / a
        cs = ∫(x^n, (x, 0, 1))
        @test cs(n => 2) == Sym(1 // 3)
        @test cs(n => -1) == Inf_
        @test cs(n => -2) == Inf_
        if VERSION >= v"1.12"                                            # (the older RischIntegrate cannot)
            cs = ∫(1 / (x^2 + a), (x, 0, Inf))
            @test cs(a => 4) == PI / 4
            @test cs(a => 0) == Inf_
            @test_throws ArgumentError cs(a => -1)                       # a pole at 1
        end
        @test N(∫(1 / (x - a), (x, 0, 1))(a => 2)) ≈ -log(2)
        # iterated: an inner integral by cases on the outer variable keeps the cases its
        # interval allows, or is split where they change
        @syms y m
        pos = @assume(n > 0, m > 0)
        @test ∫(x^n * y^m, (y, 0, x), (x, 0, 1); assume = pos) == 1 / ((m + 1) * (m + n + 2))
        @test ∫(∫(x^n * y^m, (y, 0, x); assume = pos), (x, 0, 1); assume = pos) == 1 / ((m + 1) * (m + n + 2))
        cs = SymJ.Cases{Sym}([Any[x ≪ 0] => -x, Any[x ~ 0] => Sym(0), Any[x ≫ 0] => x])
        @test ∫(cs, (x, -1, 2)) == Sym(5 // 2)
        @test ∫(cs, (x, 2, -1)) == Sym(-5 // 2)
        @test ∫(cs, (x, 1, 2)) == Sym(3 // 2)
        if VERSION >= v"1.12"                                            # (the older RischIntegrate cannot)
            @test same(∫(sqrt(1 - x^2), (x, y, sqrt(1 - y^2)), (y, 0, sqrt(Sym(2)) / 2)), 1 - sqrt(Sym(2)) / 2)   # constants simplified
        end
        # log(u), u < 0 at an end, as log(-u) + iπ, in a power or squared too
        IM = SymJ.IM
        @test ∫(log(x), (x, -1, 0)) == IM * PI - 1
        @test ∫(x * log(x), (x, -1, 0)) == Sym(1 // 4) - IM * PI / 2
        @test ∫(log(x)^2, (x, -1, 0)) == 2 - PI^2 - 2IM * PI
        # (the principal branch: e^{iπm}/(m + 1), in whatever form; at m = 1/3)
        @test isapprox(ComplexF64(N(subs(∫(x^m, (x, -1, 0); assume = @assume(m > 0)), m => Sym(1) // 3))),
                       exp(im * π / 3) / (4 / 3))
        # a variable end: the integral by cases, its limit taking the case that holds (Gabriel's horn)
        if VERSION >= v"1.12"                                            # (the older RischIntegrate cannot)
            @syms M
            horn = ∫(2PI * (1 / x) * sqrt(1 + (-1 / x)^2), (x, 1, M))
            @test horn isa Sym
            @test limit(horn, M => Inf) == Inf_
        end
        # the variable is real: conj(x) is x, as from dot (a line integral, Green's theorem)
        r, ⋅ = [cos(x), sin(x)], SymJ.LinearAlgebra.dot
        @test ∫([-r[2], r[1]] ⋅ diff.(r, x), (x, 0, 2PI)) == 2PI
        @test ∫((r / (r ⋅ r)^(Sym(3) / 2)) ⋅ diff.(r, x), (x, 0, 2PI)) == 0
        @test ∫([x - x^2, x^2 - x^4, x^2 - x^6] ⋅ diff.([x, x^2, x^3], x), (x, 0, 1)) == Sym(3 // 5)
        @test ∫(conj(a) * x, (x, 0, 1)) == a / 2                                          # parameters real, too
        @test ∫(conj(x) * x + real(x) + imag(x), x) == x^3 / 3 + x^2 / 2
        # an integrand that fails is simplified and tried again
        @test ∫(cos(x)^2 / (cos(x)^2 + sin(x)^2), (x, 0, 2PI)) == PI
        @test ∫((cos(a)^2 * x + sin(a)^2 * x) * x - cos(a)^2 * x^2 - sin(a)^2 * x^2, (x, 0, 1)) == 0
        # a radical's nonnegative factors come out, the variable's as |x| (the cone's lateral area)
        @test ∫(sqrt(x^2 * (b^2 + 1)), x) == SymJ.piecewise(x ≪ 0 => -x^2 * sqrt(b^2 + 1) / 2, true => x^2 * sqrt(b^2 + 1) / 2)
        @test ∫(sqrt(x^2 * (b^2 + 1)), (x, 0, a / b), (c, 0, 2PI); assume = @assume(a > 0, b > 0)) == PI * a^2 * sqrt(b^2 + 1) / b^2
        @test ∫(sqrt(4x^2), (x, -1, 2)) == 5
        # an antiderivative that divides by a parameter: generic, the zero only when assumed
        # (SymPy #27280, #27073)
        @syms t
        @test iszero(cancel(∫(exp(c * t), (t, 0, 1)) - (exp(c) - 1) / c))
        @test ∫(exp(c * t), (t, 0, 1); assume = (c ~ 0,)) == 1
        # (complex constants in the integrand: RischIntegrate 1.2 on; not on Julia 1.10)
        (try ∫(exp(IM * x), x); true catch; false end) && @test ∫(exp(im * b * x), (x, -1, 1); assume = (b ~ 0,)) == 2
        @test subs(∫((b - a) * exp(t * c + (1 - t) * a), (t, 0, 1); assume = (c ~ a,)), c => a) == (b - a) * exp(a)
        @test ∫(x / a, (x, 0, 1)) == 1 / (2a)               # the integrand's own denominator: no case
        # drop_removable_singularities = true: a case that is the formula's limit is merged into it
        @test iszero(cancel(∫(exp(c * t), (t, 0, 1); drop_removable_singularities = true) - (exp(c) - 1) / c))
        (try ∫(exp(IM * x), x); true catch; false end) &&
            @test ∫(exp(im * b * x), (x, -1, 1); drop_removable_singularities = true) isa Sym
        @test ∫(x^t, (x, 0, 1); drop_removable_singularities = true) isa SymJ.Cases   # t ≤ -1 is no such case
        # an end at a pole of the antiderivative (tan θ at π/2): a one-sided limit (SymPy #27185)
        @test ∫((8 - (4 - 4cos(x)^2)^(Sym(3) // 2)) / 3, (x, 0, PI / 2)) == 4PI / 3 - Sym(16) // 9
    end

    @testset "RischIntegrate: principal values" begin
        ∫(f, t...; kw...) = RischIntegrate.integrate(f, t...; principal_value = true, kw...)
        same(v, w) = VERSION >= v"1.12" ? v == w : N(v) ≈ N(w)
        @test ∫(1 / x, (x, -1, 1)) == 0
        @test same(∫(1 / x, (x, -1, 2)), log(Sym(2)))
        @test same(∫(1 / (x - 1), (x, 0, 3)), log(Sym(2)))
        @test same(∫(1 / x^3, (x, -1, 2)), Sym(3) / 8)
        @test same(∫(1 / (x^2 - 1), (x, 0, 2)), -log(Sym(3)) / 2)
        @test ∫(x, (x, -Inf, Inf)) == 0                                  # symmetric at ±∞
        @test ∫(x / (1 + x^2), (x, -Inf, Inf)) == 0
        @test ∫(sin(x), (x, -Inf, Inf)) == 0
        @test ∫(tan(x), (x, 0, PI)) == 0
        @test ∫(1 / (1 + x^2), (x, -Inf, Inf)) == PI                     # convergent: unchanged
        @test ∫(1 / x^2, (x, -1, 1)) == Sym(Inf)                          # still divergent
        @test ∫(1 / x, (x, 0, 1)) == Sym(Inf)                             # a singular end is not crossed
        # sin and cos inside an expression: their zeros in every period, not just [0, 2π]
        @test ∫(1 / (sin(x) + 1), (x, -2, 0)) == Sym(Inf)                 # -π/2
        @test ∫(1 / (sin(x) + 1), (x, 4, 5)) == Sym(Inf)                  # 3π/2
        @test ∫(1 / (1 - cos(x)), (x, 5, 7)) == Sym(Inf)                  # 2π
        @test N(∫(1 / (sin(x) + 1), (x, 5, 7))) ≈ 6.450851174313285
        # a periodic condition changes in every period: split at all of its points inside
        @test ∫(abs(sin(2x)), (x, 0, 2PI)) == 4
        @test ∫(abs(sin(3x)), (x, 0, 3PI)) == 6
        @test ∫(abs(cos(3x)), (x, 0, PI)) == 2
        @test ∫(abs(sin(2x)), (x, 2PI, 0)) == -4
        @test ∫(sign(sin(x)), (x, 0, 3PI)) == PI
        # pieces separated by a singular point are not matched there (no Inf constants)
        @test !occursin("Inf", string(RischIntegrate.integrate(1 / abs(sin(x)), x)))
        @test N(∫(1 / abs(sin(x)), (x, 1, 2))) ≈ 1.0476051700585143
        @test ∫(1 / abs(sin(x)), (x, 1, 4)) == Sym(Inf)
        @test ∫(1 / abs(x - 1), (x, 2, 3)) == log(Sym(2))             # the piece is fixed on [2, 3]
        if VERSION >= v"1.12"          # an antiderivative known only numerically (₂F₁, root sums): a float
            @test N(∫(x / sqrt(1 - x^3), (x, 0, 1 // 2))) ≈ 0.12832329160032097
            @test N(∫(1 / (x^5 - x + 1), (x, 0, 1))) ≈ 1.5919853723336097
            @test N(∫(1 / (x^5 - x + 1), (x, 1, 0))) ≈ -1.5919853723336097
            @test_throws ArgumentError ∫(x / sqrt(1 - x^3), (x, 0, 2))     # singular at 1
        end
        @test_throws ArgumentError RischIntegrate.integrate(1 / x, (x, -1, 1))   # only when asked
        # parameters: the cases, with logs of negative arguments made real
        v = ∫(1 / (x - a), (x, 0, 2); assume = @assume(a > 0, a < 2))
        @test isreal(subs(v, a => Sym(1) / 2))
        @test N(subs(v, a => 0.5)) ≈ log(1.5 / 0.5)
        cs = ∫(1 / (x - a), (x, 0, 2))
        @test cs isa SymJ.Cases
        @test N(cs(a => 3)) ≈ log(1 / 3)
        @test cs(a => 0) == Sym(Inf)
    end

    @testset "RischIntegrate: assumptions decide the form" begin
        E = Base.get_extension(SymJ, :SymJRischIntegrateExt)
        if E._has_oracle()                                               # RischIntegrate with a sign oracle
            ∫ = RischIntegrate.integrate
            @test ∫(1 / (x^2 + a), x) isa SymJ.Cases                         # on the sign of a
            F = ∫(1 / (x^2 + a), x; assume = @assume(a < 0))
            @test !SymJ._hasfn(F, :atan)
            @test SymJ._hasfn(F, :log)
            @test isequal((@assuming a < 0 begin ∫(1 / (x^2 + a), x) end), F)
            @test abs(N(subs(diff(F, x) - 1 / (x^2 + a), x => 0.3, a => -2))) < 1e-10
        end
    end

    @testset "RischIntegrate: differential equations" begin
        @syms y()
        C1, C2 = Sym(:C1), Sym(:C2)
        D = SymJ.dsolve
        # an unknown forcing: its integral from the conditions' point (SymPy #27479)
        let t = Sym(:t), P = SymFunction(:P), S = SymFunction(:S)
            r = D(diff(S(t), t) ~ P(t) - S(t), S(t), S(0) => 10)
            @test !occursin("Integral(P(0)", string(r))
            @test SymJ._ode_tidy(subs(r.rhs, t => 0)) == 10
            @test iszero(simplify(diff(r.rhs, t) - (P(t) - r.rhs)))
            r2 = D(diff(y(t), t, 2) ~ P(t), y(t), y(1) => 0, y'(1) => 2)
            @test SymJ._ode_tidy(subs(r2.rhs, t => 1)) == 0
            @test SymJ._ode_tidy(subs(diff(r2.rhs, t), t => 1)) == 2
        end
        @test D(y'(x) ~ 1 / (1 + x^2), y(x)) == (y(x) ~ atan(x) + C1)
        @test D(y''(x) ~ sin(x), y(x)) == (y(x) ~ C1 * x + C2 - sin(x))
        @test D(y'(x) ~ (x^2 + y(x)^2) / (x * y(x)), y(x)) == (y(x) ~ sqrt((C1 + 2log(x)) * x^2))   # homogeneous
        r = D(y''(x) + y(x) ~ tan(x), y(x))                                                   # variation of parameters
        @test iszero(simplify(diff(r.rhs, x, 2) + r.rhs - tan(x)))
        r = D(y'(x) + y(x) ~ x * y(x)^3, y(x))                                                 # Bernoulli
        @test iszero(simplify(diff(r.rhs, x) + r.rhs - x * r.rhs^3))
        @test D(y'(x) ~ exp(x + y(x)), y(x)) == (-exp(-y(x)) - exp(x) ~ C1)                     # separable, implicit
        # one constant, in tan or exp: the equation's solve; a constant in several
        # exponents, e^(k C), absorbed
        @test D(y'(x) ~ (1 + x^2) * (1 + y(x)^2), y(x), y(0) => 1) == (y(x) ~ tan(x^3 / 3 + x + PI / 4))
        @test D(y'(x) ~ y(x) * (1 - y(x) / 100), y(x), y(0) => 20) == (y(x) ~ 100exp(x) / (exp(x) + 4))
        @syms a
        @test D(y'(x) ~ a * y(x) * (1 - y(x)), y(x)) == (y(x) ~ C1 * exp(a * x) / (C1 * exp(a * x) - 1))
        @test D(y'(x) ~ a * y(x) * (1 - y(x)), y(x), y(0) => Sym(1) / 2) == (y(x) ~ exp(a * x) / (exp(a * x) + 1))
    end

    @testset "wrong answers (SymPy review, round 4)" begin
        @syms t A
        # the antiderivative jumps where the integrand is smooth: split there too (#26283)
        # (RischIntegrate 1.2 on; Julia 1.10 resolves 1.0)
        if isdefined(RischIntegrate, :clausen2)
            @test RischIntegrate.integrate(1 / (A - cos(t)), (t, 0, 2PI); assume = (A ≫ 1,)) == 2PI / sqrt(A^2 - 1)
            # (by value: RischIntegrate 1.3 gives it over one denominator)
            @test iszero(simplify(RischIntegrate.integrate(cos(t) / (A - cos(t)), (t, 0, 2PI); assume = (A ≫ 1,)) - (2PI * A / sqrt(A^2 - 1) - 2PI)))
            @test RischIntegrate.integrate(1 / (3 - cos(t)), (t, 0, 2PI)) == PI / sqrt(Sym(2))
        end
        # Bessel functions of non-integer order at 0: x^ν times a series in x² (#25678)
        @test Sym(series(SF.besselj(Sym(1) // 2, x), x, 0, 4)) == expand(sqrt(Sym(2)) / sqrt(PI) * (sqrt(x) - x^(Sym(5) // 2) / 6))
        @test N(subs(Sym(series(SF.bessely(Sym(5) // 3, x), x, 0, 4)), x => Sym(1) // 100)) ≈ SF.bessely(5 / 3, 0.01) rtol = 1e-8
        @test N(subs(Sym(series(SF.besselk(Sym(1) // 3, x), x, 0, 5)), x => Sym(1) // 100)) ≈ SF.besselk(1 / 3, 0.01) rtol = 1e-10
        @test Sym(series(SF.besseli(Sym(1) // 3, 2x), x, 0, 3)) == expand(x^(Sym(1) // 3) * (1 + 3x^2 / 4) / SF.gamma(Sym(4) // 3))
        @test Sym(series(SF.besselj(1, x), x, 0, 5)) == x / 2 - x^3 / 16
        @test_throws ArgumentError series(SF.bessely(1, x), x, 0, 4)       # logarithmic terms
    end

    @testset "missing cases (SymPy review, round 4)" begin
        @syms z a b r R n tf
        F = SymFunction(:F)
        # ₂F₁ antiderivatives at the ends: continuity, Gauss's sum at 1 (#26477)
        if isdefined(RischIntegrate, :clausen2)
            v = RischIntegrate.integrate(z^(a - 1) * (1 - z)^(b - 1), (z, 0, 1); assume = (a ≫ 0, b ≫ 0))
            @test N(subs(v, a => Sym(3) // 2, b => Sym(5) // 2)) ≈ SF.beta(1.5, 2.5)
            v = RischIntegrate.integrate(z^(a - Sym(1) // 2) * sqrt(1 - z), (z, 0, 1); assume = (a ≫ 0,))
            @test N(subs(v, a => 2)) ≈ SF.beta(2.5, 1.5)
        end
        # unknown functions of x: their integrals held, the rest integrated (#25457, #26249, #26116)
        @test RischIntegrate.integrate(F(r) + r, r) == r^2 / 2 + SymJ._mkfn_node(:Integral, Sym[F(r), r])
        @test RischIntegrate.integrate(3F(r) * sin(a), r) == 3sin(a) * SymJ._mkfn_node(:Integral, Sym[F(r), r])
        @test diff(RischIntegrate.integrate(F(r) + r, r), r) == F(r) + r
        @test RischIntegrate.integrate(R * F(R), (R, b, r)) == SymJ._mkfn_node(:Integral, Sym[R * F(R), R, b, r])
    end

    @testset "bigger missing cases (SymPy review, round 4)" begin
        @syms z zR
        # |exp(w)|² for complex w refined to a real Gaussian (#26102)
        v = RischIntegrate.integrate(abs(exp(-IM / 2 * x^2 / (IM * zR + z)))^2, (x, -Inf, Inf); assume = (zR ≫ 0, z ∈ ℝ))
        @test N(subs(v, zR => 2, z => 1)) ≈ sqrt(π) * sqrt(5) / sqrt(2)
    end

    @testset "missing cases (SymPy review, round 6)" begin
        @syms y C tau
        ∫ = RischIntegrate.integrate
        # Γ, ψ at their poles: an error, as tan(π/2) is; Laurent series there (#23156, #23635)
        @test_throws SymJ.DivisionByZero SF.gamma(Sym(0))
        @test_throws SymJ.DivisionByZero SF.digamma(Sym(-2))
        @test Sym(series(SF.gamma(y), y, 0, 1)) == 1 / y + SF.digamma(Sym(1))
        @test Sym(series(SF.polygamma(1, y), y, 0, 1)) == 1 / y^2 + PI^2 / 6
        @test N(subs(Sym(series(SF.gamma(y), y, -1, 3)), y => -1 + Sym(1) // 100)) ≈ SF.gamma(-0.99) rtol = 1e-5
        # Euler's constant as -ψ(1)
        @test Sym(Base.MathConstants.eulergamma) == -SF.digamma(Sym(1))
        # Airy functions at +∞ (#23563)
        v = N(subs(Sym(series(SF.airyai(x), x, Inf, 4)), x => 10))
        @test v ≈ SF.airyai(10.0) rtol = 1e-4
        v = N(subs(Sym(series(SF.airybi(x), x, Inf, 4)), x => 10))
        @test v ≈ SF.airybi(10.0) rtol = 1e-4
        # ∫ dx/|x| real for x < 0 (#22533)
        @test ∫(1 / abs(x), x; assume = (x ∈ ℝ,)) == SymJ.piecewise(x ≪ 0 => -log(-x), true => log(x))
        # an answer by cases integrated case by case (#23610)
        r = ∫(∫(x^a, (x, 0, 1)), a)
        @test r isa SymJ.Cases
        # the end a pole of a real integrand: ±∞, not "no limit" (#23337)
        @test ∫(1 / (x - 3), (x, 1, 3)) == -Inf
        # odd terms over (-a, a) (#23246)
        @test ∫((x^3 * cos(x / 2) + Sym(1) // 2) * sqrt(4 - x^2), (x, -2, 2)) == PI
        # Fresnel, Si and parts (#22487)
        @test ∫(sin(x^2), (x, 0, Inf)) == sqrt(2PI) / 4
        @test N(∫((cos(x^2) - cos(x)) / x^2, (x, -Inf, Inf))) ≈ π - sqrt(2π)
        # breakpoints that move with a parameter: a convolution of boxes (#23761)
        gp = SymJ.piecewise((0 ≦ x, x ≦ 1) => 1, true => 0)
        @test ∫(subs(gp, x => tau) * subs(gp, x => x - tau), (tau, -Inf, Inf)) ==
              SymJ.piecewise(x ≦ 0 => 0, x ≦ 1 => x, x ≦ 2 => 2 - x, true => 0)
        # reduction of order (#23351)
        u = SymFunction(:u)
        sol = SymJ.dsolve(x^2 * (2 - x) * diff(u(x), x, 2) + 2x * diff(u(x), x) - 2u(x) ~ exp(-x) * x^2 * (2 - x)^2, u(x))
        @test sol isa SymJ.Relation
        # 1/(C + cos x) over a period by cases on C, divergent for |C| ≤ 1 (#23613)
        if isdefined(RischIntegrate, :clausen2)
            r = ∫(1 / (C + cos(x)), (x, 0, 2PI))
            @test r isa SymJ.Cases
            @test any(cv -> cv.second == 2PI / sqrt(C^2 - 1), r)
            @test_throws ArgumentError ∫(1 / (Sym(1) // 2 + cos(x)), (x, 0, 2PI))
            @test_throws ArgumentError ∫(1 / (-Sym(1) // 3 + cos(x)), (x, 0, 2PI))
        end
    end

    @testset "wrong answers (SymPy review, round 7)" begin
        @syms a
        # a limit not found is not divergence: the sign of a decided, a > 0 converges (#21054)
        r = RischIntegrate.integrate(sin(x) * exp(-a * x), (x, 0, Inf))
        @test r isa SymJ.Cases
        @test any(cv -> cv.second == 1 / (a^2 + 1), r)
    end

    @testset "missing cases (SymPy review, round 7)" begin
        @syms a b
        # a condition in a second variable: the integral a piecewise function of it, the
        # same in either order (#22023)
        f = SymJ.piecewise((a ≦ 1 // 2, b ≦ 1 // 2) => 1, true => 0)
        g = f + SymJ.piecewise((a ≧ 1 // 2, b ≧ 1 // 2) => 1, true => 0)
        @test simplify(f) isa Sym
        @test subs(RischIntegrate.integrate(f, a), a => 0, b => 0) == 0
        r = RischIntegrate.integrate(f, (a, 0, 1))
        @test r isa SymJ.Cases
        @test any(cv -> cv.second == Sym(1) // 2, r)
        @test RischIntegrate.integrate(RischIntegrate.integrate(g, (a, 0, 1)), (b, 0, 1)) == Sym(1) // 2
        @test RischIntegrate.integrate(RischIntegrate.integrate(g, (b, 0, 1)), (a, 0, 1)) == Sym(1) // 2
        # an impossible case pruned: l ≠ k and l = -k with k, l ≥ 0 (#22066)
        @syms k l L
        r = RischIntegrate.integrate(cos(x * k * PI / L) * cos(x * l * PI / L), (x, -L, L);
                                     assume = (L ≫ 0, k ∈ ℤ, l ∈ ℤ, k ≧ 0, l ≧ 0))
        @test r isa SymJ.Cases
        @test length(r) == 3
    end

    @testset "crashes and hangs (SymPy review, round 8)" begin
        @syms a b
        # a case where the integrand is not defined (b = 0): an error for that case (#19672)
        r = RischIntegrate.integrate(exp(-(x - a)^2 / b), (x, -Inf, Inf))
        @test r isa SymJ.Cases
        @test any(cv -> cv.second isa Sym && iszero(simplify(cv.second - sqrt(b * PI); assume = (b ≫ 0,))), r)   # (√b √π or √(bπ))
        # abs of functions with a pole: no division by zero (#20755)
        @test RischIntegrate.integrate(abs(1 / x), (x, 1, 2)) == log(Sym(2))
        @test RischIntegrate.integrate(abs(atan(1 / x)), (x, -1, 1)) == log(Sym(2)) + PI / 2
    end

    @testset "wrong answers (SymPy review, round 8)" begin
        # |acot x| over (-1, 1): acot jumps at 0, so the abs is kept (#20755)
        @test RischIntegrate.integrate(abs(acot(x)), (x, -1, 1)) == log(Sym(2)) + PI / 2
        @test RischIntegrate.integrate(abs(acot(x)), (x, -1, 0)) == log(Sym(2)) / 2 + PI / 4
        @test RischIntegrate.integrate(abs(acot(x)), x) isa Sym
    end

    @testset "missing cases (SymPy review, round 8)" begin
        @syms t y
        # dirac with pieces: the pieces' value at its point (#19962, #20872)
        @test RischIntegrate.integrate(SymJ.dirac(t - 5) * SymJ.heaviside(t - 1), (t, 0, 10)) == 1
        @test RischIntegrate.integrate(SymJ.piecewise((x ≧ 1 // 50, x ≦ 3 // 100) => 100 * SymJ.dirac(x - 1 // 40), true => 0),
                                       (x, 1 // 50, 3 // 100)) == 100
        # improper integrals: divergence, Γ moments (#20831, #20685, #20043)
        @test RischIntegrate.integrate(exp(x) * exp(log(x)^2), (x, 1, Inf)) == Inf
        @test RischIntegrate.integrate(exp(-x^3), (x, 0, Inf)) == SymJ.apply_function(:gamma, Sym[Sym(4) // 3])
        @test isapprox(Float64(N(RischIntegrate.integrate(t * log(t)^2 * exp(-t), (t, 0, Inf)))), 0.8236806608528794; atol = 1e-12)
        # floats in, floats out
        r = RischIntegrate.integrate(-0.000425170068027211 * x / (0.000212585034013605 * x^2 - 1.0), x)
        @test !SymJ._isexactex(r)
    end

    @testset "missing cases (SymPy review, round 9)" begin
        @syms n
        # the even part over (-∞, ∞) (#19175)
        same(v, w) = VERSION >= v"1.12" ? v == w : N(v) ≈ N(w)            # (the older RischIntegrate: numerically)
        @test same(RischIntegrate.integrate(exp(-x^2) / (1 + exp(x)), (x, -Inf, Inf)), sqrt(PI) / 2)
        # floats in a definite integral: a float value; x^n kept as a power
        @test RischIntegrate.integrate(2.5x, (x, 0, 1)) == Sym(1.25)
        @test !has(RischIntegrate.integrate(x^n * log(x), x), SymJ.E)
        # a nested piecewise solved as one (#18359)
        pdf = SymJ.piecewise((0 ≦ x) & (x ≪ 1) => 1 // 2, (x ≧ 1) & (x ≪ 2) => 0, (x ≧ 2) & (x ≪ 3) => 1 // 2, true => 0)
        F = RischIntegrate.integrate(pdf, (x, 0, x))
        @test solve(SymJ.piecewise(x ≧ 0 => F, true => 0) ~ 1 // 2, x) isa SymJ.SymSet
    end

    @testset "wrong answers (SymPy review, round 10)" begin
        # the boundary terms of a dirac derivative (#17383)
        d1 = SymJ.apply_function(:dirac, Sym[x, Sym(1)])
        @test RischIntegrate.integrate(exp(x) * d1, x) == SymJ.dirac(x) - SymJ.heaviside(x)
        @syms z
        @test RischIntegrate.integrate(exp(x) * d1, (x, -1, z)) == SymJ.dirac(z) - SymJ.heaviside(z)
    end

    @testset "missing cases (SymPy review, round 10)" begin
        @syms t s
        # Frullani integrals; no false divergence from opposite infinite terms (#17058)
        @test RischIntegrate.integrate((atan(4x) - atan(2x)) / x, (x, 0, Inf)) == PI * log(Sym(2)) / 2
        @test RischIntegrate.integrate((exp(-2x) - exp(-5x)) / x, (x, 0, Inf)) == log(Sym(5) / 2)
        # the bound variable of ∫₀ᵗ may be t itself (#17094)
        u = SymFunction(:u)
        @test SymJ.laplace_transform(RischIntegrate.integrate(u(t), (t, 0, t)), t, s) ==
              SymJ.apply_function(:LaplaceTransform, Sym[u(t), t, s]) / s
    end

    @testset "wrong answers (SymPy review, round 11)" begin
        @syms w t
        f = SymFunction(:f)
        # initial conditions at a singular point: an error, not free constants (#15113)
        L = (1 - x^2) * diff(f(x), x, 2) - 2x * diff(f(x), x) + 2f(x) ~ 0
        @test_throws ArgumentError SymJ.dsolve(L, f(x), f(-1) => -1, f(1) => 1)
        @test SymJ.dsolve(L, f(x), f(0) => 0, f(Sym(1) // 2) => Sym(1) // 2) == (f(x) ~ x)
        # a Gaussian with a complex coefficient: decided on its real part (#15292)
        @test N(RischIntegrate.integrate(exp(-(1 - IM) * x^2), (x, 0, Inf))) ≈ sqrt(pi) / (2sqrt(1 - im))
        # a parameter declared complex is not taken real (#16242)
        @test_throws ArgumentError RischIntegrate.integrate(exp(-t * (1 + IM * w)), (t, 0, Inf); assume = (w ∈ ℂ,))
    end

    @testset "wrong answers (SymPy review, round 12)" begin
        @syms t a
        f, g = SymFunction(:f), SymFunction(:g)
        ∫ = RischIntegrate.integrate
        # an antiderivative that jumps at points not located, or a log crossing its cut:
        # no value, not half of it (#13609; RischNotImplemented with some versions)
        @test_throws Exception ∫(exp(2IM * PI * t) / (exp(2IM * PI * t) - 1 // 2), (t, 0, 1))
        # F' ≠ f on a branch cut: no value (#14147)
        @test_throws ArgumentError ∫(asin(x), (x, 0, 9))
        @test_throws ArgumentError ∫(acos(x), (x, 2, 3))
        @test ∫(asin(x), (x, 0, 1)) == PI / 2 - 1
        @test ∫(sqrt(x), (x, -1, 1)) == Sym(2) / 3 + 2IM / 3
        # (log(1 + i) - log(i - 1); -iπ/2 from RischIntegrate 1.3)
        @test ∫(1 / (x + IM), (x, -1, 1)) in (log(1 + IM) - log(IM - 1), -IM * PI / 2)
        # g'(x) is not the unknown's derivative (#13691); ∫ g'(x) dx = g(x)
        @test SymJ.dsolve(diff(f(x), x) ~ a * diff(g(x), x), f(x)) == (f(x) ~ Sym(:C1) + a * g(x))
        @test ∫(diff(g(x), x), (x, 0, 1)) == g(1) - g(0)
        # 1 - erf(15): a cancellation below 256 bits (#14637)
        @test N(1 - SF.erf(Sym(15))) ≈ 7.212994172451207e-100 rtol = 1e-10
    end

    @testset "crashes and hangs (SymPy review, round 12)" begin
        # Si and Ci have no BigFloat methods: N keeps the Float64 value (#14716)
        # (∫ from 1/2 to 1 of log(x + 5)cos(πx))
        Si(u) = SymJ.apply_function(:sinint, Sym[u])
        v = log(Sym(2) // 11) / PI - Si(11PI / 2) / PI + Si(6PI) / PI
        @test N(v) ≈ -0.5604789635702104 rtol = 1e-10
    end

    @testset "missing cases (SymPy review, round 12)" begin
        @syms y z sg
        f = SymFunction(:f)
        ∫ = RischIntegrate.integrate
        # values of an unknown function free of x are constants (#14796)
        @test ∫(f(x) + f(1) + 10, x) == (f(1) + 10) * x + SymJ.apply_function(:Integral, Sym[f(x), x])
        # a Gaussian times cos or sin of a linear argument (#14573)
        @test ∫(cos(x) * exp(-x^2 / (2sg^2)), (x, -Inf, Inf); assume = (sg ≫ 0,)) == sqrt(2PI) * sg * exp(-sg^2 / 2)
        @test ∫(sin(2x + 1) * exp(-x^2), (x, -Inf, Inf)) == sqrt(PI) * sin(Sym(1)) / SymJ.E
        # an iterated integral with the outer variables in their intervals (#14096)
        @test ∫(1 / (1 + x + y + z)^2, (x, 0, 1), (y, 0, 1), (z, 0, 1)) == log(Sym(19683) // 16384)
    end

    @testset "crashes and hangs (SymPy review, round 11)" begin
        @syms a t RC C R1 L Vs V0 I0
        # arithmetic and series on Cases (#16182)
        cs = RischIntegrate.integrate(1 / (a + x), (x, 0, 1))
        @test (1 / cs)(a => 1) == 1 / log(Sym(2))
        # a circuit with symbolic initial values: the constants through symbols for the
        # radicals (#16125)
        V, II = SymFunction(:V), SymFunction(:II)
        r = SymJ.dsolve([diff(V(t), t) ~ -V(t) / RC + II(t) / C, diff(II(t), t) ~ -R1 / L * II(t) - V(t) / L + Vs / L],
                        [V(t), II(t)], V(0) => V0, II(0) => I0)
        vals = (RC => 2, C => 3, R1 => 5, L => 7, Vs => 11, V0 => 1, I0 => 2)
        @test N(subs(subs(r[1].rhs, vals...), t => 0)) ≈ 1
    end

    @testset "missing cases (SymPy review, round 11)" begin
        @syms n
        # x q(x) with q symmetric about the middle; log sin over [0, π/2] (#16396, #16084)
        @test RischIntegrate.integrate(x * sin(x) / (1 + cos(x)^2), (x, 0, PI)) == PI^2 / 4
        @test N(RischIntegrate.integrate(log(sin(x)), (x, 0, PI / 2))) ≈ -pi * log(2) / 2
        # xⁿ e^(-x) log x with no assumption: Cases, divergence in its own case (#15432)
        @test RischIntegrate.integrate(x^n * exp(-x) * log(x), (x, 0, Inf)) isa SymJ.Cases
    end

    @testset "crashes and hangs (SymPy review, round 13)" begin
        @syms x y z c w s L0 r mu
        # a case's fact undefined at another case's value (w = 0 in F3/w): no DivisionByZero
        # (#12568, #11992, #12249)
        @test RischIntegrate.integrate(1 / (1 + x^2 * y^2), (y, 0, z)) isa SymJ.Cases
        if VERSION >= v"1.12"                                            # (the older RischIntegrate cannot)
            @test RischIntegrate.integrate(1 / sqrt(c^2 + (3 - w * s)^2), (s, 0, L0)) isa SymJ.Cases
            @test RischIntegrate.integrate((z - r * mu) / (z^2 + r^2 - 2r * z * mu)^(Sym(3) / 2), (mu, -1, 1)) isa SymJ.Cases
        end
    end

    @testset "missing cases (SymPy review, round 13)" begin
        @syms x a b
        # B(a, b) with a log factor: ∂B/∂a (#12456)
        ψ(w) = SymJ.apply_function(:digamma, Sym[w])
        Γ(w) = SymJ.apply_function(:gamma, Sym[w])
        @test RischIntegrate.integrate(x^(a - 1) * (1 - x)^(b - 1) * log(x), (x, 0, 1); assume = (a ≫ 0, b ≫ 0)) ==
              Γ(a) * Γ(b) / Γ(a + b) * (ψ(a) - ψ(a + b))
        @test N(RischIntegrate.integrate(x^(Sym(1) // 3) * (1 - x)^(Sym(1) // 2) * log(x), (x, 0, 1))) ≈ -0.4528196802991289
    end

    @testset "missing cases (SymPy review, round 14)" begin
        @syms x y a
        Γ(w) = SymJ.apply_function(:gamma, Sym[w])
        # e^(c - a xᵖ) over [0, b]: a lower incomplete gamma (#11845)
        @test N(RischIntegrate.integrate(exp(-x^3), (x, 0, 1))) ≈ 0.8075111821397486
        @test N(RischIntegrate.integrate(exp(-2x^4), (x, 0, 1))) ≈ 0.7490153580905826
        @test N(subs(RischIntegrate.integrate(exp(-y - x^3), (x, 0, 1)), y => 0)) ≈ 0.8075111821397486
        # xᵐ log(1/x)ˢ over [0, 1]: Γ(s + 1)/(m + 1)^(s + 1) (#11876)
        @test RischIntegrate.integrate(log(1 / x)^(a - 1), (x, 0, 1); assume = (a ≫ 0,)) == Γ(a)
        @test RischIntegrate.integrate(sqrt(log(1 / x)), (x, 0, 1)) == sqrt(PI) / 2
        @test RischIntegrate.integrate(x^2 * (-log(x))^3, (x, 0, 1)) == Sym(2) // 27
    end

    @testset "wrong answers (SymPy review, round 15)" begin
        @syms x d k
        # a pole at an end by the facts: ±∞, not F(1) - F(0) with log 0 in it (#9753)
        @test RischIntegrate.integrate(1 / (d - k * x), (x, 0, 1); assume = (d ~ k, k ≫ 0)) == Sym(Inf)
        @test RischIntegrate.integrate(1 / (d - k * x), (x, 0, 1); assume = (d ~ k, k ≪ 0)) == Sym(-Inf)
    end

    @testset "Borwein integrals: a large antiderivative only refined" begin
        # (RischIntegrate's step TX, from its #130, gives the antiderivatives quickly;
        # before it, k = 7 runs for an hour, so the tests are skipped)
        Bs(k) = prod(sin(x / i) / (x / i) for i in 1:2:(2k - 1))
        if isdefined(RischIntegrate, :integrate_trig_products)
            @test RischIntegrate.integrate(Bs(7), (x, 0, Inf)) == PI / 2
            @test RischIntegrate.integrate(Bs(8), (x, 0, Inf)) ==
                  Sym(467807924713440738696537864469 // 935615849440640907310521750000) * PI
        end
    end

    @testset "abs of complex arguments: |exp(2πix) + exp(2πiy)|" begin
        @syms y
        I = RischIntegrate.integrate
        z = exp(2 * PI * im * x) + exp(2 * PI * im * y)
        f = abs(z)
        # |z| is ±z only for real z
        @test SymJ.piecewise(f) == f
        @test SymJ.piecewise(abs(x - y)) != abs(x - y)
        # ∫ over one period of |cos| whatever the shift (0 at y = ±2, negative at 0.7)
        g = I(abs(cos(PI * (x - y))), (x, 0, 1))
        for y0 in (-2, 0.7, 1.2, 2)
            @test N(subs(g, y => y0)) ≈ 2 / π
        end
        @test N(subs(I(abs(sin(x - y)), (x, 0, 1)), y => 1 // 2)) ≈ 2 * (1 - cos(0.5))
        # an antiderivative, for real y, once RischIntegrate's are continuous across the
        # zeros of 1 + cos x (its sign_jump_correction)
        if isdefined(RischIntegrate, :sign_jump_correction)
            F = I(f, x)
            for (x0, y0) in ((0.1, 0.3), (0.4, 0.75), (0.9, -0.2))
                @test abs(N(subs(diff(F, x) - f, x => x0, y => y0))) < 1e-12
            end
            h = I(f, (x, 0, 1 // 2))
            @test N(subs(h, y => 0.8)) ≈ (2 - sin(0.3π) - sin(0.2π)) * 2 / π
        end
        # √(1 + cos x), from RischIntegrate 1.2
        if pkgversion(RischIntegrate) >= v"1.2"
            # (once a stack overflow)
            @test I(abs(1 + exp(2 * PI * im * x)), (x, 0, 1)) == 4 / PI
            # the inner integral is 4/π for every y: one period of 2|cos(π(x - y))| (it was 0)
            @test I(f, (x, 0, 1)) == 4 / PI
            @test I(f, (x, 0, 1), (y, 0, 1)) == 4 / PI
            @test I(f, (y, 0, 1), (x, 0, 1)) == 4 / PI
        end
    end

    @testset "asymptotic series at ±∞: exponential and oscillating factors" begin
        @syms y
        S(f, x0, n) = Sym(series(f, x, x0, n))
        rp = 1 / sqrt(PI)
        same(a, b) = iszero(expand(a - b))
        # erf, erfc: exponentially small terms are kept (SymPy #27109)
        @test same(S(SF.erf(x), Inf, 4), 1 - exp(-x^2) * (rp / x - rp / (2x^3)))
        @test same(S(SF.erfc(x), Inf, 4), exp(-x^2) * (rp / x - rp / (2x^3)))
        @test same(S(SF.erf(x), -Inf, 2), -1 - exp(-x^2) * rp / x)
        @test string(series(SF.erfc(x), x, Inf, 2)) == "1/exp(x^2)*(1/(x*sqrt(π)) + O(1/x^2))"
        # factors cancel: e^(x²) erfc(x) is a power series; so is erfcx
        @test same(S(exp(x^2) * SF.erfc(x), Inf, 4), rp / x - rp / (2x^3))
        @test same(S(SF.erfcx(x), Inf, 4), rp / x - rp / (2x^3))
        @test S(sin(x)^2 + cos(x)^2, Inf, 3) == 1
        # exponential integrals, Si and Ci (SymPy #26257, #26207)
        @test same(S(SF.expinti(x), Inf, 4), exp(x) * (1 / x + 1 / x^2 + 2 / x^3))
        @test same(S(SF.expint(x), Inf, 4), exp(-x) * (1 / x - 1 / x^2 + 2 / x^3))
        @test same(S(SF.sinint(x), Inf, 4), PI / 2 - cos(x) * (1 / x - 2 / x^3) - sin(x) / x^2)
        @test same(S(SF.cosint(x), Inf, 4), sin(x) * (1 / x - 2 / x^3) - cos(x) / x^2)
        # numerically, at moderate x
        for (f, x0, xv, n, tol) in ((SF.sinint(x), -Inf, -20.0, 6, 1e-5), (SF.erf(x), -Inf, -3.0, 6, 1e-6),
                                    (SF.besselj(0, x), Inf, 30.0, 6, 1e-8), (SF.bessely(1, x), Inf, 30.0, 6, 1e-8),
                                    (SF.besseli(1, x), Inf, 30.0, 6, 1e-8), (SF.besselk(1, x), Inf, 30.0, 6, 1e-8),
                                    (SF.expinti(x), -Inf, -30.0, 6, 1e-4), (SF.gamma(x), Inf, 10.0, 5, 1e-7))
            v = N(subs(S(f, x0, n), x => xv))
            @test abs(v / N(subs(f, x => xv)) - 1) < tol
        end
        # modified Bessel functions (SymPy #26210); K of half-integer order terminates
        @test N(subs(S(SF.besselk(Sym(1) // 2, x), Inf, 3), x => 3)) ≈ SF.besselk(0.5, 3.0) rtol = 1e-14
        # Stirling's series (SymPy #26096): Γ(x) ≈ √(2π) x^(x-1/2) e^(-x) (1 + 1/(12x))
        @test same(S(SF.loggamma(x), Inf, 4), (x - Sym(1) // 2) * log(x) - x + log(2PI) / 2 + 1 / (12x) - 1 / (360x^3))
        @test N(subs(S(SF.gamma(x), Inf, 2), x => 20)) ≈ sqrt(2π) * 20.0^19.5 * exp(-20.0) * (1 + 1 / 240) rtol = 1e-12
        # derivatives and sums of series with factors
        @test Sym(diff(series(SF.erf(x), x, Inf, 6), x)) == 2rp * exp(-x^2)
        @test Sym(series(SF.erf(x), x, Inf, 4) + series(SF.erfc(x), x, Inf, 4)) == 1
        # not supported: a reciprocal of a sum with such factors; finite points as before
        @test_throws ArgumentError series(1 / (1 + exp(x)), x, Inf, 3)
        @test_throws ArgumentError series(exp(1 / x), x, 0, 3)
    end

    @testset "missing cases (SymPy review, round 19)" begin
        @syms y
        # a Poly integrated as its expression (#4196)
        @test RischIntegrate.integrate(Poly(x + x^2 * y + y^3, x, y), (y, 0, PI)) ==
              RischIntegrate.integrate(x + x^2 * y + y^3, (y, 0, PI))
        @test RischIntegrate.integrate(Poly(x^2 + 1, x), x) == x^3 / 3 + x
    end

    @testset "wrong answers (SymPy review, round 19)" begin
        @syms y
        # a piecewise ramp's squared difference quotient: the antiderivative at y = 1 (and
        # y = 0) was NaN, a 0·log 0 at the point (#4313)
        u = SymJ.piecewise(x ≪ 0 => 0, x ≪ 1 => x, true => 1)
        g = (u - subs(u, x => y))^2 / (x - y)^2
        F = RischIntegrate.integrate(g, x)
        @test subs(F, y => 1) == RischIntegrate.integrate(subs(g, y => 1), x)
        @test subs(F, y => 0) == RischIntegrate.integrate(subs(g, y => 0), x)
        # a piecewise in the parameters only: a constant for the integrator (it crashed)
        uy = subs(u, x => y)
        G = RischIntegrate.integrate((x - uy)^2 / (x - y)^2, x)
        @test iszero(simplify(subs(diff(G, x) - (x - uy)^2 / (x - y)^2, y => Sym(1) // 2)))
    end

    @testset "wrong answers (SymPy review, round 20)" begin
        # definite integrals along a complex segment: the cuts it crosses and its ends by
        # limits along it (#31017; ∫ 1/z from -i to -1 was 3iπ/2)
        @syms z
        @test RischIntegrate.integrate(1 / z, (z, -im, -1)) == -im * PI / 2
        @test N(RischIntegrate.integrate(log(z), (z, -im, -1))) ≈ 1 + π / 2 + im * (π - 1) rtol = 1e-12
        @test N(RischIntegrate.integrate(log(z) / z^2, (z, -im, -1))) ≈ 1 + π / 2 - im * (π - 1) rtol = 1e-12
        @test RischIntegrate.integrate(1 / z, (z, -1 - im, -1 + im)) == -im * PI / 2
        @test N(RischIntegrate.integrate(sqrt(z), (z, -1 - im, -1 + im))) ≈ 0.4752076 * im rtol = 1e-6
        @test N(RischIntegrate.integrate(log(z + 2), (z, -3 - im, -3 + im))) ≈ 0.26394350735484 * im rtol = 1e-10
        @test N(RischIntegrate.integrate(csc(z), (z, 1 - im, 1 + im))) ≈ 1.8987949537455 * im rtol = 1e-10
        # log(x^i) = i log(x) only for x > e^-π: no value rather than a wrong one (#31017)
        # (the older RischIntegrate of Julia 1.10 refuses the integrand itself)
        @test_throws Union{ArgumentError, RischIntegrate.RischNotImplemented} RischIntegrate.integrate(log(z^im)^2, (z, 0, 1))
    end

    @testset "crashes (SymPy review, round 20)" begin
        @syms a
        # ±√(-a²), a real, are no points of the line: no case √(-a²) = 1, whose limits
        # divided by zero (#2537)
        @test RischIntegrate.integrate(1 / (x^2 + a^2), (x, -1, 1)) == 2atan(1 / a) / a
        @test RischIntegrate.integrate(a / (x^2 + a^2), (x, -1, 1)) == 2atan(1 / a)
        @test RischIntegrate.integrate(1 / (x^2 + a^2), (x, 0, 1)) == atan(1 / a) / a
        # (the older RischIntegrate of Julia 1.10 cannot integrate it: a parameter in the radicand)
        VERSION >= v"1.12" && @test RischIntegrate.integrate((a / sqrt(x^2 + a^2))^3, (x, -1, 1)) == 2a / sqrt(a^2 + 1)
    end

    @testset "no elementary antiderivative: the definite integral unevaluated" begin
        @syms y
        v = RischIntegrate.integrate(exp(x^3), (x, 0, 1))
        @test SymJ._isdefint(v)
        @test N(v) ≈ 1.3419044179774198 rtol = 1e-12
        # ∫₀¹∫₀¹ e^(-x²y²): the inner one in erf, the outer unevaluated (SymPy #5326; the
        # older RischIntegrate of Julia 1.10 cannot do the inner one)
        if VERSION >= v"1.12"
            w = RischIntegrate.integrate(exp(-x^2 * y^2), (x, 0, 1), (y, 0, 1))
            @test SymJ._isdefint(w)
            @test N(w) ≈ 0.9059404763223628 rtol = 1e-10
        end
        @test N(RischIntegrate.integrate(exp(-x^2) * log(x + 2), (x, 0, Inf))) ≈ 0.8231099519447602 rtol = 1e-10
        # integrate(ex): the unevaluated integrals in ex carried out, innermost first; those
        # that cannot be kept
        # (numerically: the older RischIntegrate of Julia 1.10 gives √π as a float)
        @test N(RischIntegrate.integrate(SymJ.Integral(exp(-x^2), (x, 0, Inf)))) ≈ sqrt(π) / 2 rtol = 1e-14
        @test RischIntegrate.integrate(SymJ.Integral(x * y, (x, 0, 1), (y, 0, 2))) == 1
        @test RischIntegrate.integrate(1 + SymJ.Integral(x^2, x)) == 1 + x^3 / 3
        H = SymJ.Integral(exp(x^3), (x, 0, 1))
        @test isequal(RischIntegrate.integrate(H), H)
        @test RischIntegrate.integrate(2H + SymJ.Integral(x, (x, 0, 1))) == 2H + Sym(1) // 2
        @test_throws ArgumentError RischIntegrate.integrate(x^2)
        # not when it diverges (e^(-x) sin(x)/x² at 0), nor with a parameter
        @test_throws RischIntegrate.NotElementary RischIntegrate.integrate(exp(-x) * sin(x) / x^2, (x, 0, Inf))
        @test_throws RischIntegrate.NotElementary RischIntegrate.integrate(y * exp(x^3), (x, 0, 1))
    end

    @testset "indefinite integrals over ℝ: cases on parameter signs" begin
        @syms t
        # (RischIntegrate's sign oracle: the installed one of Julia ≥ 1.12)
        if VERSION >= v"1.12"
            r = RischIntegrate.integrate(1 / (x^2 + a), x)
            @test r isa SymJ.Cases
            @test Set(join(string.(c.first), ", ") for c in r.cases) == Set(["a ≪ 0", "a ≫ 0"])
            @test isequal(only(c.second for c in r.cases if string(only(c.first)) == "a ≫ 0"), atan(x / sqrt(a)) / sqrt(a))
            # each case's form differentiates back, under its condition
            ln = only(c.second for c in r.cases if string(only(c.first)) == "a ≪ 0")
            @test iszero(simplify(diff(ln, x) - 1 / (x^2 + a); assume = @assume(a < 0)))
            # a sign given: that case only
            @test isequal(RischIntegrate.integrate(1 / (x^2 + a), x; assume = @assume(a > 0)), atan(x / sqrt(a)) / sqrt(a))
            # a sign the integrator asks as -1/(4ab): the condition as a*b
            s = RischIntegrate.integrate(1 / (a * x^2 + b), x)
            @test Set(join(string.(c.first), ", ") for c in s.cases) == Set(["a*b ≪ 0", "a*b ≫ 0"])
            # the parameters real: √(a²) is |a| (the older integrator refuses abs)
            @test isequal(RischIntegrate.integrate(sqrt(a^2), x), abs(a) * x)
            # a fact on two intervals (t² > 1) decides its own sign, so each case gets its form
            s = RischIntegrate.integrate(1 / (x^2 + t^2 - 1), x)
            k = findfirst(c -> join(string.(c.first), ", ") == "t^2 - 1 ≫ 0", s.cases)
            @test SymJ._hasfn(s.cases[k].second, :atan)
            # the two signs giving one answer, up to expansion, are one case
            @test RischIntegrate.integrate(sqrt(a^2 - x^2), x) isa Sym
        end
        # a case where the integrand is nowhere real is dropped: √(a x² + b) with a, b < 0
        # (a < 0 implied by -b/a < 0 and b < 0)
        E = Base.get_extension(SymJ, :SymJRischIntegrateExt)
        R = SymJ.FactAnd(Any[a ∈ SymJ.ℝ, b ∈ SymJ.ℝ])
        @test E._nowhere_real(sqrt(a * x^2 + b), x, SymJ._with_fact(SymJ.Assumptions(), SymJ.FactAnd(Any[R, -b / a ≪ 0, b ≪ 0])))
        @test !E._nowhere_real(sqrt(a * x^2 + b), x, SymJ._with_fact(SymJ.Assumptions(), SymJ.FactAnd(Any[R, -b / a ≪ 0, b ≫ 0])))
        @test E._same_answer(-(-sin(x) - x) * a, a * (x + sin(x)))
        # no sign asked: one answer, as before
        @test isequal(RischIntegrate.integrate(x^n, x), x^(n + 1) / (n + 1))
        # t² - 1 > 0 decides the sign of -1/(t² - 1)
        B = SymJ._with_fact(SymJ.Assumptions(), SymJ.Relation(:gt, t^2 - 1, SymJ.ZERO))
        @test SymJ._sign_oracle(B)(-1 / (t^2 - 1)) == -1
    end

    @testset "differentiation under the integral sign" begin
        # (exact on the RischIntegrate of Julia ≥ 1.12; the older one of Julia 1.10 may keep
        # the integral unevaluated, whose value is then the same by quadrature)
        exact = VERSION >= v"1.12"
        # a constant made a parameter: I'(a) = 2π/√(1 - a²), I(0) = 0 (SymPy #30869, #30677)
        v = RischIntegrate.integrate(log((1 + cos(x) / 2) / (1 - cos(x) / 2)) / cos(x), (x, 0, PI))
        @test N(v) ≈ π^2 / 3 rtol = 1e-12
        exact && @test isequal(v, PI^2 / 3)
        w = RischIntegrate.integrate(log(1 + cos(x) / 3), (x, 0, PI))       # π log((3 + 2√2)/6)
        @test N(w) ≈ π * log((3 + 2sqrt(2)) / 6) rtol = 1e-12
        exact && @test !SymJ._isdefint(w)
        u = RischIntegrate.integrate(exp(-x) * sin(2x) / x, (x, 0, Inf))       # atan(2)
        @test N(u) ≈ atan(2) rtol = 1e-12
        exact && @test isequal(u, atan(Sym(2)))
        # I'(a) = 2atan(√(1 - a)/√(1 + a))/√(1 - a²), refined to acos(a)/√(1 - a²) for
        # 0 < a < 1/2, integrates in a: π²/8 - acos(1/2)²/2 = 5π²/72
        h = RischIntegrate.integrate(log(1 + cos(x) / 2) / cos(x), (x, 0, PI / 2))
        @test N(h) ≈ 5π^2 / 72 rtol = 1e-12
        exact && @test isequal(h, 5PI^2 / 72)
        # none of the constants works: the integral unevaluated, as before
        @test SymJ._isdefint(RischIntegrate.integrate(exp(x^3), (x, 0, 1)))
        # with the assumptions' antiderivative a root sum, the generic one (atan's) instead
        @syms a
        exact && @test isequal(RischIntegrate.integrate(2 / (1 - a^2 * cos(x)^2), (x, 0, PI); assume = (a ≫ 0, a ≪ 1)),
                               2PI / sqrt(1 - a^2))
        # with no assumptions, by cases: 2π/√(1 - a²) for |a| < 1 (the limit at acos(1/a)
        # for |a| > 1 raised GMP's DivideError in GruntzLimit < 1.8.1)
        if exact
            r = RischIntegrate.integrate(2 / (1 - a^2 * cos(x)^2), (x, 0, PI))
            @test r isa SymJ.Cases
            @test any(c -> c.second isa Sym && isequal(c.second, 2PI / sqrt(1 - a^2)), r.cases)
            # the conditions on the singular points acos(±1/a) as conditions on a: five
            # cases (a = ±1: ∞; a > 1, a < -1: a pole inside; |a| < 1), no contradictory
            # one (acos(1/a) > π/2 with acos(1/a) < acos(-1/a)) whose limit failed
            @test length(r.cases) == 5
            @test all(c -> c.second isa Sym || occursin("does not converge", sprint(showerror, c.second)), r.cases)
            @test !any(c -> any(f -> occursin("acos", string(f)), c.first), r.cases)
            # each case's condition on a itself: -1/a - 1 ~ 0 as a ~ -1, 1/a + 1 > 0 with
            # a < 0 as a < -1, 1/a² - 1 > 0 as a² < 1
            @test Set(join(string.(c.first), ", ") for c in r.cases) == Set(["a ~ -1", "a ~ 1", "a ≪ -1", "a ≫ 1", "a^2 ≪ 1"])
            # singular points √(-a): their comparisons with the ends as conditions on a
            # (√(-a) ~ 0 is a ~ 0, √(-a) - 1 > 0 is a < -1), and the iπ of log(-u) + iπ
            # cancelling between the ends collected
            s = RischIntegrate.integrate(1 / (x^2 + a), (x, 0, 1))
            @test Set(join(string.(c.first), ", ") for c in s.cases) == Set(["a ~ 0", "a ~ -1", "a ≪ -1", "a ≪ 0, a ≫ -1", "a ≫ 0"])
            v = only(c.second for c in s.cases if join(string.(c.first), ", ") == "a ≪ -1")
            @test !has(v, IM)
            @test N(subs(v, a => -4)) ≈ log(1 / 3) / 4 rtol = 1e-12
        end
    end

    @testset "missing cases (SymPy review, round 20)" begin
        @syms a
        y = SymFunction(:y)
        # a Riccati equation through a solution of its linear form (#31017); numerically:
        # the older RischIntegrate of Julia 1.10 gives √π as a float
        F = SymJ.dsolve(diff(y(x), x) ~ x^2 - 5 - y(x)^2, y(x)).rhs
        @test abs(N(subs(diff(F, x) - (x^2 - 5 - F^2), x => Sym(7) // 10, Sym(:C1) => Sym(1) // 3))) < 1e-9
        # homogeneous with a radical, for x > 0: √(x² + x²v²) = x√(1 + v²) (Kamke 1.113;
        # #31017)
        G = SymJ.dsolve(a * sqrt(x^2 + y(x)^2) + x * diff(y(x), x) - y(x) ~ 0, y(x); assume = (x ≫ 0,)).rhs
        res = a * sqrt(x^2 + G^2) + x * diff(G, x) - G
        @test abs(N(subs(res, x => 2, a => Sym(1) // 3, Sym(:C1) => Sym(1) // 5))) < 1e-12
    end

    @testset "Fresnel integrals" begin
        S, C = SymJ.fresnels, SymJ.fresnelc
        @syms a
        # ∫ sin(ax² + bx + c), ∫ cos(…) in S and C (RischIntegrate's were erf of (-1)^(1/4) multiples)
        F = RischIntegrate.integrate(sin(x^2), x)
        @test N(subs(F, x => 2)) ≈ sqrt(π / 2) * S(sqrt(2 / π) * 2) rtol = 1e-14
        for f in (sin(x^2), cos(x^2), sin(3x^2), cos(x^2 + x + 1), sin(1 - 2x^2), 2x + 5cos(x^2))
            F = RischIntegrate.integrate(f, x)
            @test SymJ._hasfn(F, :fresnels) || SymJ._hasfn(F, :fresnelc)
            @test abs(N(subs(diff(F, x) - f, x => Sym(7) // 10))) < 1e-14
        end
        @test SymJ._hasfn(RischIntegrate.integrate(sin(a * x^2), x; assume = (a ≫ 0,)), :fresnels)
        @test RischIntegrate.integrate(cos(PI * x^2 / 2), (x, 0, 1)) == C(Sym(1))
        @test N(RischIntegrate.integrate(sin(x^2), (x, 0, Inf))) ≈ sqrt(2π) / 4 rtol = 1e-14
        @test N(RischIntegrate.integrate(sin(x^2), (x, -Inf, Inf))) ≈ sqrt(2π) / 2 rtol = 1e-14
        @test latexify(S(x)) == raw"$\mathrm{S}\left( x \right)$"
    end
end
