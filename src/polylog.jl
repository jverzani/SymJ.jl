# PolyLog on expressions: `li0` ... `li6`, `li(n, z)` and the real versions `reli1` ...
# `reli4`, `reli(n, x)` with a `Sym` argument build unevaluated functions, evaluated
# numerically by PolyLog when asked for.
module PolyLogSupport

import PolyLog as PL
using ..SymJ
import SpecialFunctions as SF
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
        register_function!(f; julia = f === :li ? _num_li : (n, z) -> g(_intorder(n), z),
            eval = args -> iszero(args[2]) ? ZERO : nothing,
            diff = (args, i) -> i == 2 ? g(args[1] - 1, args[2]) / args[2] : nothing)
    end
end

# PolyLog's li and reli take an integer order; numerically one may come as a float
_intorder(n) = isinteger(n) ? Int(n) : throw(DomainError(n, "reli: the order must be an integer"))

# Li_s(z) for a real order s, numerically (SymPy #25159): PolyLog for an integer s; else
# Σ zᵏ/kˢ for |z| ≤ 1/2, and near 1, with μ = -log z (|μ| < 2π),
# Γ(1 - s) μ^(s-1) + Σ ζ(s - k) (-μ)^k/k!. Li_s(1) = ζ(s) for s > 1, ∞ otherwise.
function _num_li(s, z)
    isinteger(s) && return PL.li(Int(s), z)
    s isa Real || throw(DomainError(s, "li: a complex order is not supported"))
    s, w = float(s), z isa Real ? float(z) : complex(float(z))
    T = typeof(real(w))
    if abs(w) <= 0.5
        acc, p = zero(w), one(w)
        for k in 1:200
            p *= w
            t = p / T(k)^s
            acc += t
            abs(t) <= eps(T) * abs(acc) && break
        end
        return acc
    end
    isone(w) && return s > 1 ? SF.zeta(s) : T(Inf)
    μ = -log(complex(w))
    abs(μ) < 4 || throw(DomainError(z, "li: a non-integer order with |log z| ≥ 4 is not supported"))
    acc = SF.gamma(1 - s) * μ^(s - 1)
    p = one(μ)
    for k in 0:400
        k > 0 && (p *= -μ / k)
        t = SF.zeta(s - k) * p
        acc += t
        k > 10 && abs(t) <= eps(T) * abs(acc) && break
    end
    # a real z < 1 has a real value (z > 1 is on the branch cut: the principal value)
    w isa Real && w < 1 ? real(acc) : acc
end

end
