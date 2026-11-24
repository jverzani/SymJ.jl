# HypergeometricFunctions on expressions: `_₁F₁`, `_₂F₁` and `_₃F₂` with a `Sym`
# argument build unevaluated functions, evaluated numerically by HypergeometricFunctions
# when asked for. `pFq(α, β, z)` maps to the named function of its sizes.
module HypergeometricSupport

import HypergeometricFunctions as HF
using ..SymJ
using ..SymJ: Sym, apply_function, register_function!, _symbolic_signatures, ONE

# name => (p, q): pFq with p numerator and q denominator parameters, then z
const FUNCTIONS = [:_₁F₁ => (1, 1), :_₂F₁ => (2, 1), :_₃F₂ => (3, 2)]

for (f, (p, q)) in FUNCTIONS, sig in _symbolic_signatures(p + q + 1)
    args = [Symbol(:a, i) for i in 1:(p + q + 1)]
    @eval HF.$f($([:($a::$T) for (a, T) in zip(args, sig)]...)) =
        apply_function($(QuoteNode(f)), Sym[$(args...)])
end

# pFq(α, β, z) of the sizes above, one method per size (HypergeometricFunctions
# dispatches on the tuple sizes too)
HF.pFq(α::Tuple{Any}, β::Tuple{Any}, z::Sym) = HF._₁F₁(α..., β..., z)
HF.pFq(α::NTuple{2,Any}, β::Tuple{Any}, z::Sym) = HF._₂F₁(α..., β..., z)
HF.pFq(α::NTuple{3,Any}, β::NTuple{2,Any}, z::Sym) = HF._₃F₂(α..., β..., z)

# ₂F₁ at BigFloat arguments. Near z = 1 and z = ∞ HypergeometricFunctions divides
# 1/Γ(u) - 1/Γ(u + ϵ) by ϵ = c - a - b - m (or b - a - m, m an integer), and at BigFloat a
# parameter like 1/3 leaves ϵ ≈ 1e-77 where it is 0: every digit cancels (₂F₁(1, 1/3;
# 4/3; 3/4) gave 1.4228 for 1.3773). There the series is summed, in z or in z/(z - 1)
# (Pfaff), whichever is smaller; elsewhere, and at other types, HypergeometricFunctions.
function _hyp2f1(a, b, c, z)
    v = _hyp2f1_big(a, b, c, z)
    v === nothing ? HF._₂F₁(a, b, c, z) : v
end

function _hyp2f1_big(a, b, c, z)
    # (SymJ's evaluation keeps integer arguments as Ints)
    big_(t) = real(typeof(t)) === BigFloat
    any(big_, (a, b, c, z)) && all(t -> t isa Integer || big_(t), (a, b, c, z)) || return nothing
    tol = big(2.0)^(8 - precision(BigFloat)) * max(1, abs(a), abs(b), abs(c))
    nearint(t) = abs(t - round(real(t))) <= tol
    (nearint(c - a - b) || nearint(b - a)) || return nothing
    nearint(c) && real(c) < 1 / 2 && return nothing                # a pole in c
    w = z / (z - 1)
    min(abs(z), abs(w)) < 1 - big(2.0)^-20 || return nothing
    p = precision(BigFloat)
    s = setprecision(BigFloat, p + 64) do
        A, B, C, Z, W = big(a), big(b), big(c), big(z), big(w)
        abs(Z) <= abs(W) ? _hyp2f1_series(A, B, C, Z) :
            (v = _hyp2f1_series(A, C - B, C, W); v === nothing ? nothing : (1 - Z)^(-A) * v)
    end
    s === nothing && return nothing
    s isa Complex ? Complex(BigFloat(real(s); precision = p), BigFloat(imag(s); precision = p)) :
        BigFloat(s; precision = p)
end

# Σ (a)ₙ(b)ₙ/((c)ₙ n!) zⁿ for |z| < 1, to the working precision; nothing past 10⁷ terms
function _hyp2f1_series(a, b, c, z)
    s, t = one(z), one(z)
    tol = eps(BigFloat)
    n0 = 2 * (abs(a) + abs(b) + abs(c)) + 4                       # terms may grow first
    for n in 0:10^7
        t *= (a + n) * (b + n) / ((c + n) * (n + 1)) * z
        s += t
        (iszero(t) || n > n0 && abs(t) <= tol * abs(s)) && return s
    end
    nothing
end

function __init__()
    for (f, (p, q)) in FUNCTIONS
        g = getfield(HF, f)
        register_function!(f; julia = f === :_₂F₁ ? _hyp2f1 : g,
            # pFq(α; β; 0) = 1; ₁F₁(a; a; z) = e^z (1/₁F₁(1, 1, x) at 0 went to Gruntz,
            # which has no ₁F₁; SymPy #11465)
            eval = args -> iszero(args[end]) ? ONE : (f === :_₁F₁ && SymJ._eq(args[1], args[2])) ? exp(args[3]) : nothing,
            # d/dz pFq(α; β; z) = Πα / Πβ pFq(α + 1; β + 1; z); no closed form in the parameters
            diff = (args, i) -> i == p + q + 1 ?
                prod(args[1:p]) / prod(args[(p+1):(p+q)]) * g((args[1:(end-1)] .+ 1)..., args[end]) : nothing)
    end
end

end
