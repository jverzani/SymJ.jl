# RischIntegrate is generic over term types; this workload runs it on SymJ terms, so the
# compiled code ships with the extension instead of being compiled on the first integral
# (about 4 s otherwise).
module SymJRischIntegrateExt

using SymJ
using RischIntegrate: RischIntegrate, integrate
using PrecompileTools: @setup_workload, @compile_workload

# piecewise functions, heaviside and dirac: interval by interval (see SymJ's piecewise.jl);
# abs, sign, max and min, which the integrator does not take, as piecewise functions
# (for real x). The assumptions (assume = …, or an @assuming block) decide the signs of
# parameter expressions the answer's form depends on, through the integrator's sign
# oracle, when the installed RischIntegrate has one; signs they leave open take its
# default.
function RischIntegrate.integrate(ex::Sym, x::Sym; assume = nothing, kw...)
    A = SymJ._assumptions(assume)
    if A !== nothing && _has_oracle() && !haskey(kw, :sign)
        kw = (; sign = SymJ._sign_oracle(A), kw...)
    end
    integ = f -> RischIntegrate.risch_integrate(f, x; kw...)
    SymJ._haspieces(ex, x) && return SymJ._integrate_pieces(ex, x, integ)
    try
        integ(ex)
    catch err
        err isa InterruptException && rethrow()
        pw = SymJ.piecewise(ex)
        SymJ._haspieces(pw, x) || rethrow()
        SymJ._integrate_pieces(pw, x, integ)
    end
end

# does risch_integrate take a sign oracle (RischIntegrate with PR #70)?
const _ORACLE = Ref{Union{Nothing,Bool}}(nothing)
function _has_oracle()
    _ORACLE[] === nothing && (_ORACLE[] = any(m -> :sign in Base.kwarg_decl(m), methods(RischIntegrate.risch_integrate)))
    _ORACLE[]::Bool
end

# definite integrals: split at the integrand's singular points, one-sided limits at them
# and at ±∞, by cases on parameters; principal_value = true takes Cauchy principal values
# at the singular points inside and at ±∞ (see SymJ's definite.jl)
function RischIntegrate.integrate(ex::Sym, xab::Tuple, more::Tuple...; assume = nothing,
                                  principal_value::Bool = false, kw...)
    antider = (f, x, A) -> RischIntegrate.integrate(f, x; assume = A, kw...)
    v = SymJ._definite_integral(ex, xab, assume, antider; principal_value)
    for t in more
        v isa Sym || throw(ArgumentError("integrate: an iterated integral whose inner integral has cases"))
        v = SymJ._definite_integral(v, t, assume, antider; principal_value)
    end
    v
end

# the integrator SymJ uses inside (dsolve)
__init__() = (SymJ._INTEGRATOR[] = (f, x) -> integrate(f, x))

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
            # (an older RischIntegrate, the one Julia 1.10 gets, cannot do all of them)
            try integrate(f, x) catch end
        end
        P, H = SymJ.piecewise(x ≪ 0 => -x, true => x^2), SymJ.heaviside(x - 1)
        try integrate(P, x); integrate(x * H, x); integrate(cos(x) * SymJ.dirac(x - 1), x) catch end
        # answers with special functions (Ei, ₁F₁): not every term type has them
        for f in (exp(x) / x, exp(a * x) / x, 1 / log(x)^(2 + a))
            try integrate(f, x) catch end
        end
    end
end

end
