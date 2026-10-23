# RischIntegrate is generic over term types; this workload runs it on SymJ terms, so the
# compiled code ships with the extension instead of being compiled on the first integral
# (about 4 s otherwise).
module SymJRischIntegrateExt

using SymJ
using RischIntegrate: integrate
using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    @compile_workload begin
        x, a, b, c, d = Sym(:x), Sym(:a), Sym(:b), Sym(:c), Sym(:d)
        for f in (x^3 - 2x + 1, 1 / (x^2 - 1), (2x + 3) / (x^2 + 2x + 5), 1 / (1 + x^3),
                  sqrt(1 + 3x), x / sqrt(4 - x^2), x^2 * (1 + x)^(1 // 3),    # algebraic
                  sqrt(5 - x^2), sqrt(7 - x^2), cos(x) / sqrt(sin(x)^3),
                  sin(x)^2 * cos(x), 1 / (2 + cos(x)), tan(x),                # trigonometric
                  1 / (3 + sin(x) + 2cos(x)), 1 / (7 + 3sin(x) - cos(x)), (sin(x) + 2cos(x)) / (sin(x) - cos(x)),
                  x * exp(2x), exp(x) * sin(x), log(x)^2, x * log(1 + x),     # exp and log
                  sin(x) * exp(cos(x)), exp(atan(x)) / (1 + x^2),
                  atan(x), x * asin(x),                                     # inverse trigonometric
                  1 / (a + x^2), 1 / (b + a * x^2), x^a * log(x),             # parameters
                  x^a * log(b * x), x^(a - 1) * cos(x^a), exp(a * x) * sin(b * x),
                  (a^2 - x^2)^(3 // 2), 1 / (1 + a * sin(x)),
                  (a * sin(x) + cos(x)) / (sin(x) + b * cos(x)), 1 / (a^2 * cos(x)^2 + sin(x)^2),
                  (a * sin(x) + b * cos(x)) / (c * sin(x) + d * cos(x)))
            integrate(f, x)
        end
        # answers with special functions (Ei, ₁F₁): not every term type has them
        for f in (exp(x) / x, exp(a * x) / x, 1 / log(x)^(2 + a))
            try integrate(f, x) catch end
        end
    end
end

end
