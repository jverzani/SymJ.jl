# GruntzLimit is generic over term types; this workload runs it on SymJ terms, so the
# compiled code ships with the extension instead of being compiled on the first limit
# (about 3 s otherwise).
module SymJGruntzLimitExt

using SymJ
using GruntzLimit: gruntz_limit
using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    @compile_workload begin
        x, a = Sym(:x), Sym(:a)
        gruntz_limit(sin(x) / x, x, 0)                                 # series at a point
        gruntz_limit((1 - cos(x)) / x^2, x, 0; dir = :-)               # one-sided
        gruntz_limit(x * exp(1 / x), x, 0; dir = :-)
        gruntz_limit(exp(x) * (exp(1 / x - exp(-x)) - exp(1 / x)), x, Inf)   # nested exponentials
        gruntz_limit(log(x)^2 / sqrt(x), x, Inf)
        gruntz_limit(x * ((1 + 1 / x)^x - ℯ), x, Inf)
        gruntz_limit(sqrt(x^2 + x) - x, x, Inf)                        # radicals
        gruntz_limit((1 + a / x)^x, x, Inf)                            # a parameter
        gruntz_limit(x + sqrt(x^2 + x), x, -Inf)
        gruntz_limit((x - Sym(π) / 2) * tan(x), x, Sym(π) / 2)          # a symbolic point
        gruntz_limit(x^x, x, 0)
    end
end

end
