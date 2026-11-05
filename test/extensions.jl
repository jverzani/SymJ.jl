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
        @test kind(SF.gamma(Sym(-1))) == FN && kind(SF.gamma(x)) == FN
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
        @test ∫(x^2 * D(x - 2, 1), x) == -4H(x - 2)                                       # δ′: -g′(a)
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
            ft = ∫(rect * exp(-2IM * PI * a * x), (x, -Inf, Inf))         # with the case a = 0
            @test ft(a => 0) == 1
            @test last(ft.cases).second == sin(PI * a) / (PI * a)
        end
        # no piece holds past the last condition: not defined there (SymPy #28469)
        @test ∫(P(x ≪ 1 => x), x) == P(x ≪ 1 => x^2 / 2)
        let (a1, a2, a3, t) = symbols("a1 a2 a3 t")
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
        # an answer whose only non-elementary part is a ₂F₁ term: an expression, so limit
        # and the rest take it (SymPy #24879)
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
            @test simplify(∫(exp((b - a) * x - b * t), (x, 0, t))) isa SymJ.Cases
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
        @test ∫(exp(-a * x^2 + b * x), (x, -Inf, Inf); assume = (a ≫ 0,)) == sqrt(PI / a) * exp(b^2 / (4a))
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
        @test ∫(x^m, (x, -1, 0); assume = @assume(m > 0)) == exp(IM * m * PI) / (m + 1)
        # a variable end: the integral by cases, its limit taking the case that holds (Gabriel's horn)
        if VERSION >= v"1.12"                                            # (the older RischIntegrate cannot)
            @syms M
            horn = ∫(2PI * (1 / x) * sqrt(1 + (-1 / x)^2), (x, 1, M))
            @test horn isa SymJ.Cases
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
        # an antiderivative that divides by a parameter: a case where it vanishes (SymPy #27280, #27073)
        @syms t
        cs = ∫(exp(c * t), (t, 0, 1))
        @test cs isa SymJ.Cases
        @test cs(c => 0) == 1
        @test cs(c => 2) == (exp(Sym(2)) - 1) / 2
        # (complex constants in the integrand: RischIntegrate 1.2 on; not on Julia 1.10)
        (try ∫(exp(IM * x), x); true catch; false end) && @test ∫(exp(im * b * x), (x, -1, 1))(b => 0) == 2
        @test ∫((b - a) * exp(t * c + (1 - t) * a), (t, 0, 1))(c => a) == (b - a) * exp(a)
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
            @test ∫(1 / (x^2 + a), x) == atan(x / sqrt(a)) / sqrt(a)         # generic
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
            @test RischIntegrate.integrate(cos(t) / (A - cos(t)), (t, 0, 2PI); assume = (A ≫ 1,)) == 2PI * A / sqrt(A^2 - 1) - 2PI
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
end
