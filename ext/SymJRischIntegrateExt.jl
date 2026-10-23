# RischIntegrate is generic over term types; this workload runs it on SymJ terms, so the
# compiled code ships with the extension instead of being compiled on the first integral
# (about 4 s otherwise).
module SymJRischIntegrateExt

using SymJ
using RischIntegrate: integrate
using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    @compile_workload begin
        x, a = Sym(:x), Sym(:a)
        for f in (x^3 - 2x + 1, 1 / (x^2 - 1), (2x + 3) / (x^2 + 2x + 5),     # rational
                  sqrt(1 + 3x), x / sqrt(4 - x^2), x^2 * (1 + x)^(1 // 3),    # algebraic
                  sin(x)^2 * cos(x), 1 / (2 + cos(x)), tan(x),                # trigonometric
                  x * exp(2x), exp(x) * sin(x), log(x)^2, x * log(1 + x),     # exp and log
                  atan(x), x * asin(x),                                     # inverse trigonometric
                  1 / (a + x^2), x^a * log(x))                              # parameters
            integrate(f, x)
        end
    end
end

end
