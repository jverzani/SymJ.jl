# PolyLog on expressions: `li0` ... `li6`, `li(n, z)` and the real versions `reli1` ...
# `reli4`, `reli(n, x)` with a `Sym` argument build unevaluated functions, evaluated
# numerically by PolyLog when asked for.
module PolyLogSupport

import PolyLog as PL
using ..SymJ
using ..SymJ: Sym, apply_function, register_function!, _symbolic_signatures, _unary, ZERO, PI

const ONEARG = (:li0, :li1, :li2, :li3, :li4, :li5, :li6, :reli1, :reli2, :reli3, :reli4)

for f in ONEARG
    @eval PL.$f(z::Sym) = apply_function($(QuoteNode(f)), Sym[z])
end
for f in (:li, :reli), sig in _symbolic_signatures(2)
    @eval PL.$f(n::$(sig[1]), z::$(sig[2])) = apply_function($(QuoteNode(f)), Sym[n, z])
end

# d/dz Li_n(z) = Li_{n-1}(z)/z, with Li_1(z) = -log(1 - z) and Li_0(z) = z/(1 - z);
# the same for the real versions
function __init__()
    at0 = _unary(z -> iszero(z) ? ZERO : nothing)          # Li_n(0) = 0
    for (pre, ks) in (("li", 0:6), ("reli", 1:4)), k in ks
        f = Symbol(pre, k)
        prev = k <= 1 ? nothing : getfield(PL, Symbol(pre, k - 1))
        register_function!(f; julia = getfield(PL, f), eval = at0,
            diff = (args, i) -> k == 0 ? inv((1 - args[1])^2) : k == 1 ? inv(1 - args[1]) :
                                prev(args[1]) / args[1])
    end
    register_function!(:li2; julia = PL.li2,                 # with its special values
        eval = _unary(z -> iszero(z) ? ZERO : isone(z) ? PI^2 / 6 : z == -1 ? -PI^2 / 12 : nothing),
        diff = (args, i) -> PL.li1(args[1]) / args[1])
    for f in (:li, :reli)
        g = getfield(PL, f)
        register_function!(f; julia = g,
            eval = args -> iszero(args[2]) ? ZERO : nothing,
            diff = (args, i) -> i == 2 ? g(args[1] - 1, args[2]) / args[2] : nothing)
    end
end

end
