# The package extensions for SpecialFunctions, HypergeometricFunctions, PolyLog and
# RischIntegrate.
import SpecialFunctions as SF
import RischIntegrate
import HypergeometricFunctions as HF
import PolyLog as PL

@testset "extensions" begin
    @syms x n a b c
    PI = SymJ.PI

    # ambiguities introduced by an extension's methods
    for (pkg, ext) in ((SF, :SymJSpecialFunctionsExt), (HF, :SymJHypergeometricFunctionsExt), (PL, :SymJPolyLogExt))
        E = Base.get_extension(SymJ, ext)
        @test E isa Module
        amb = filter(p -> p[1].module === E || p[2].module === E, Test.detect_ambiguities(pkg, E))
        @test isempty(amb)
    end

    # central difference of f at v
    fd(f, v; h = 1e-6) = (f(v + h) - f(v - h)) / (2h)

    @testset "SpecialFunctions" begin
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
    end

    @testset "RischIntegrate: differential equations" begin
        @syms y()
        C1, C2 = Sym(:C1), Sym(:C2)
        D = SymJ.dsolve
        @test D(y'(x) ~ 1 / (1 + x^2), y(x)) == (y(x) ~ atan(x) + C1)
        @test D(y''(x) ~ sin(x), y(x)) == (y(x) ~ C1 * x + C2 - sin(x))
        @test D(y'(x) ~ (x^2 + y(x)^2) / (x * y(x)), y(x)) == (y(x) ~ sqrt((C1 + 2log(x)) * x^2))   # homogeneous
        r = D(y''(x) + y(x) ~ tan(x), y(x))                                                   # variation of parameters
        @test iszero(simplify(diff(r.rhs, x, 2) + r.rhs - tan(x)))
        r = D(y'(x) + y(x) ~ x * y(x)^3, y(x))                                                 # Bernoulli
        @test iszero(simplify(diff(r.rhs, x) + r.rhs - x * r.rhs^3))
        @test D(y'(x) ~ exp(x + y(x)), y(x)) == (-exp(-y(x)) - exp(x) ~ C1)                     # separable, implicit
    end
end
