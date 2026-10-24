# A precompile workload: the operations downstream code reaches first, so their
# compiled code ships with the package instead of being compiled on first use.
# Mixed Sym/number arithmetic covers the number types generic code hands in
# (Int, Rational, Rational{BigInt}, BigInt, Float64, π).

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    nums = Any[2, -1, 1 // 2, big(3) // 4, big(5), 2.5, π, ℯ]
    @compile_workload begin
        x, y = Sym(:x), Sym(:y)
        for c in nums
            c * x; x * c; c + x; x + c; x - c; c - x; x / c; c / x; x^c; c^x
            Sym(c) * x; Sym(c)()
        end
        e = (x + 1)^2 * exp(-x / 2) + sin(x) * cos(y) - log(x) / x + sqrt(Sym(2)) * x^(1 // 3)
        for f in (exp, log, sin, cos, tan, atan, asin, acos, sinh, cosh, tanh, abs, sqrt, sign,
                  floor, ceil, sec, csc, cot, cbrt)
            f(x); f(Sym(2)); f(-x)
        end
        log(Sym(2), x); x^y; inv(x); -x; x == y; hash(x)
        string(e); sprint(show, MIME"text/plain"(), e)
        e(x => 1, y => 2); (x^2 + 1)(3); (x^2 + 1)(3.0); (x^2 + 1)(); N(sqrt(Sym(2)))
        iscall(e); operation(e); arguments(e); sorted_arguments(e)
        for a in arguments(e)
            iscall(a) && maketerm(Sym, operation(a), Any[arguments(a)...], nothing)
        end
        maketerm(Sym, +, Any[x, 1], nothing); maketerm(Sym, *, Any[2, x], nothing)
        maketerm(Sym, exp, Any[2], nothing); maketerm(Sym, ^, Any[x, 1 // 2], nothing)
        expand((x + y + 1)^3); expand(e); diff(e, x); subs(e, x => y + 1)
        cancel((x^2 - 1) / (x + 1)); factor_squarefree(x^3 - x^2 - x + 1)
        factor(x^4 - 1); factor(x^2 * y^2 - 1); apart(1 / (x^2 - 1))
        z = Sym(:z)                                  # multivariate factoring (EEZ), more forms
        factor(x^6 - y^6); factor(expand((x^2 + y * z + 1) * (x * y - z^2 + 3) * (x + y + z)))
        apart(1 / ((x - parameter(:a)) * (x - y)), x); apart((x^3 + 1) / (x^2 - 1), x)
        gcd(x^2 - 1, x^2 + 2x + 1); lcm(x^2 - 1, x + 1)
        solve(exp(x) ~ 2, x); solve(sqrt(x + 1) ~ x - 1, x); solve(x^4 - 10x^2 + 1 ~ 0, x)
        solve(x^2 - 2 ~ 0, x); solve(x^3 - 3x + 1 ~ 0, x); solve(x^2 ≫ 1, x); solve(x * exp(x) ~ 1, x)
        solve([x + y ~ 3, x - y ~ 1], [x, y]); solve([parameter(:a) * x + y ~ 1, x - y ~ 2], [x, y])
        a = parameter(:a); (a * x^2 + x)(3); (a * x^2 + x)(x => 2); free_parameters(a * x)
        p = Poly((x + y + 1)^3, x, y)
        p * p; p + p; p - p; divrem(p, Poly(x + 1, x, y)); gcd(p, Poly((x + y + 1) * (x - 1), x, y))
        Sym(p)
        # limits (GruntzLimit's code for SymJ terms)
        limit(sin(x) / x, x => 0); limit((1 - cos(x)) / x^2, x => 0; dir = :-)
        limit(x * exp(1 / x), x => 0; dir = :-); limit(x^x, x => 0; dir = :+)
        limit(exp(x) * (exp(1 / x - exp(-x)) - exp(1 / x)), x => Inf)
        limit(log(x)^2 / sqrt(x), x => Inf); limit(x * ((1 + 1 / x)^x - ℯ), x => Inf)
        limit(sqrt(x^2 + x) - x, x => Inf); limit((1 + a / x)^x, x => Inf)
        limit(x + sqrt(x^2 + x), x => -Inf); limit((x - Sym(π) / 2) * tan(x), x => Sym(π) / 2)
        limit(x^a, x => Inf); limit(x^a, x => Inf; assume = Assumptions(a ≫ 1))
        # simplification
        simplify((x^2 - 1) / (x - 1)); simplify(sqrt(Sym(2)) * sqrt(Sym(3))); simplify(x^2 + 2x + 1)
        simplify(sin(x)^2 + cos(x)^2); simplify(1 - cos(x)^2); simplify(2sin(x) * cos(x))
        simplify(sin(x) * cos(y) + cos(x) * sin(y)); simplify(tan(x) * cos(x))
        refine(sqrt(x^2); assume = Assumptions(x ≫ 0)); simplify(sqrt(x^2) / x; assume = Assumptions(x ∈ ℝ))
    end
end
