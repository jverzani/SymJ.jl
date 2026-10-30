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

    @testset "RischIntegrate: definite integrals" begin
        ∫ = RischIntegrate.integrate
        Inf_ = Sym(Inf)
        # exact with the RischIntegrate Julia 1.12 gets; the older one gives approximate constants
        same(v, w) = VERSION >= v"1.12" ? v == w : N(v) ≈ N(w)
        @test ∫(x^2, (x, 0, 3)) == 9
        @test ∫(sin(x), (x, 0, PI)) == 2
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
end
