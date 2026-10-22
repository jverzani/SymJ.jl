# HypergeometricFunctions on expressions: `_₁F₁`, `_₂F₁` and `_₃F₂` with a `Sym`
# argument build unevaluated functions, evaluated numerically by HypergeometricFunctions
# when asked for. `pFq(α, β, z)` maps to the named function of its sizes.
module SymJHypergeometricFunctionsExt

import HypergeometricFunctions as HF
using SymJ
using SymJ: Sym, apply_function, register_function!, _symbolic_signatures, ONE

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

function __init__()
    for (f, (p, q)) in FUNCTIONS
        g = getfield(HF, f)
        register_function!(f; julia = g,
            eval = args -> iszero(args[end]) ? ONE : nothing,           # pFq(α; β; 0) = 1
            # d/dz pFq(α; β; z) = Πα / Πβ pFq(α + 1; β + 1; z); no closed form in the parameters
            diff = (args, i) -> i == p + q + 1 ?
                prod(args[1:p]) / prod(args[(p+1):(p+q)]) * g((args[1:(end-1)] .+ 1)..., args[end]) : nothing)
    end
end

end
