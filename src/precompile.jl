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
        p = Poly((x + y + 1)^3, x, y)
        p * p; p + p; p - p; divrem(p, Poly(x + 1, x, y)); gcd(p, Poly((x + y + 1) * (x - 1), x, y))
        Sym(p)
    end
end
