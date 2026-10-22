# The Lambert W function (LambertW.jl): lambertw(z) the principal branch W₀,
# lambertw(z, k) branch k. W(z) e^W(z) = z, so it solves (αx + β) e^(γx + δ) = K.

import LambertW: lambertw

register_function!(:lambertw; julia = lambertw,
    eval = function (args)
        z = args[1]
        k = length(args) == 2 ? args[2] : ZERO
        (k.kind == NUMBER && _isint(k.num)) || return nothing
        if iszero(k)
            iszero(z) && return ZERO                        # W(0) = 0
            z === E && return ONE                           # W(e) = 1
        end
        (iszero(k) || k == -1) && _eq(z, -inv(E)) && return MINUS_ONE   # W(-1/e) = -1
        nothing
    end,
    # W'(z) = W / (z (1 + W)); no closed form in the branch
    diff = (args, i) -> i == 1 ? (w = apply_function(:lambertw, args); w / (args[1] * (1 + w))) : nothing)

lambertw(z::Sym) = apply_function(:lambertw, Sym[z])
lambertw(z::Sym, k::Integer) = iszero(k) ? lambertw(z) : apply_function(:lambertw, Sym[z, Sym(k)])
