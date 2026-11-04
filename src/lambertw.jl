# The Lambert W function (LambertW.jl): lambertw(z) the principal branch W₀,
# lambertw(z, k) branch k. W(z) e^W(z) = z, so it solves (αx + β) e^(γx + δ) = K.

import LambertW: lambertw

# the numeric value: LambertW's, or where it fails (a MethodError for some complex z,
# -0.5 + 0im among them) Halley's iteration on w e^w = z, from the series at the branch
# point -1/e or from log z + 2πik - log(log z + 2πik)
function _lambertw_num(z::Number, k::Integer = 0)
    try
        return lambertw(z, k)
    catch err
        err isa Union{MethodError,DomainError} || rethrow()
    end
    c = complex(float(z))
    w = if k == 0 && abs(c + exp(-one(real(c)))) < 0.5
        p = sqrt(2 * (exp(one(real(c))) * c + 1))
        -1 + p - p^2 / 3
    elseif k == 0
        log1p(c)
    else
        L = log(c) + 2π * im * k
        L - log(L)
    end
    for _ in 1:100
        ew = exp(w)
        f = w * ew - c
        d = ew * (w + 1) - (w + 2) * f / (2w + 2)
        δ = f / d
        w -= δ
        abs(δ) <= 4eps(real(typeof(c))) * max(1, abs(w)) && break
    end
    w
end

register_function!(:lambertw; julia = _lambertw_num,
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
