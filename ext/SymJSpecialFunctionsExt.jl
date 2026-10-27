# SpecialFunctions on expressions: `SpecialFunctions.f(args...)` with a `Sym` argument
# builds the unevaluated function `f(args...)`. Its numeric method is SpecialFunctions'
# own `f`, used when asked for (N, calling the expression with numbers, lambdify) or
# when an argument is a float. Exact special values and derivatives are given where
# standard.
module SymJSpecialFunctionsExt

import SpecialFunctions as SF
using SymJ
using SymJ: Sym, apply_function, register_function!, _symbolic_signatures, _unary, _odd,
            PI, ZERO, ONE, NUMBER, isexact, _isint, _toint, _rat, SymNumber

# name => argument counts
const FUNCTIONS = [
    :gamma => (1, 2), :loggamma => (1,), :digamma => (1,), :trigamma => (1,), :invdigamma => (1,),
    :polygamma => (2,), :beta => (2,), :logbeta => (2,), :logfactorial => (1,),
    :zeta => (1, 2), :eta => (1,),
    :erf => (1, 2), :erfc => (1,), :erfcx => (1,), :erfi => (1,), :erfinv => (1,), :erfcinv => (1,),
    :dawson => (1,), :logerf => (2,), :logerfc => (1,), :logerfcx => (1,), :faddeeva => (1,),
    :expint => (1, 2), :expintx => (1, 2), :expinti => (1,), :sinint => (1,), :cosint => (1,),
    :airyai => (1,), :airyaiprime => (1,), :airybi => (1,), :airybiprime => (1,),
    :airyaix => (1,), :airyaiprimex => (1,), :airybix => (1,), :airybiprimex => (1,),
    :besselj => (2,), :bessely => (2,), :besseli => (2,), :besselk => (2,),
    :besseljx => (2,), :besselyx => (2,), :besselix => (2,), :besselkx => (2,),
    :besselj0 => (1,), :besselj1 => (1,), :bessely0 => (1,), :bessely1 => (1,),
    :besselh => (3,), :besselhx => (3,), :hankelh1 => (2,), :hankelh2 => (2,),
    :hankelh1x => (2,), :hankelh2x => (2,), :sphericalbesselj => (2,), :sphericalbessely => (2,),
    :jinc => (1,), :ellipk => (1,), :ellipe => (1,),
]

for (f, arities) in FUNCTIONS, n in arities, sig in _symbolic_signatures(n)
    args = [Symbol(:a, i) for i in 1:n]
    @eval SF.$f($([:($a::$T) for (a, T) in zip(args, sig)]...)) =
        apply_function($(QuoteNode(f)), Sym[$(args...)])
end

## exact values

isint(a::Sym) = a.kind == NUMBER && isexact(a.num) && _isint(a.num)
ishalfint(a::Sym) = a.kind == NUMBER && isexact(a.num) && _isint(a.num + a.num) && !_isint(a.num)
const SQRTPI = sqrt(PI)

# Γ(n) = (n-1)! and Γ(n + 1/2) = (2n)! / (4^n n!) √π
function gamma_value(args)
    length(args) == 1 || return nothing
    a = args[1]
    if isint(a)
        n = _toint(a.num)
        return n >= 1 ? Sym(factorial(big(n - 1))) : nothing
    elseif ishalfint(a) && a.num > SymNumber(0)
        n = _toint(a.num - SymNumber(1 // 2))
        return Sym(factorial(big(2n)) // (big(4)^n * factorial(big(n)))) * SQRTPI
    end
    nothing
end

# ζ(2n) = (-1)^(n+1) B_2n (2π)^2n / (2 (2n)!), for a few n
const ZETA_EVEN = Dict(2 => PI^2 / 6, 4 => PI^4 / 90, 6 => PI^6 / 945, 8 => PI^8 / 9450)
function zeta_value(args)
    length(args) == 1 && isint(args[1]) || return nothing
    n = _toint(args[1].num)
    n == 0 ? Sym(-1 // 2) : get(ZETA_EVEN, n, nothing)
end

# ψ⁽ᵐ⁾(k) = (-1)^(m+1) m! (ζ(m + 1) - Σ_{j<k} 1/j^(m+1)) at a positive integer k, m ≥ 1
function polygamma_value(args)
    m, a = args
    (isint(m) && m.num >= SymNumber(1) && isint(a) && a.num >= SymNumber(1)) || return nothing
    m, k = _toint(m.num), _toint(a.num)
    k <= 1000 || return nothing
    s = sum((big(1) // big(j)^(m + 1) for j in 1:(k - 1)); init = big(0) // 1)
    Sym(isodd(m) ? factorial(big(m)) : -factorial(big(m))) * (SF.zeta(Sym(m + 1)) - Sym(s))
end

at0(v) = _unary(a -> iszero(a) ? v : nothing)
odd_at0(name) = _unary(a -> iszero(a) ? ZERO : _odd(name, a))

function beta_value(args)
    a, b = args
    (isint(a) && isint(b) && a.num > SymNumber(0) && b.num > SymNumber(0)) || return nothing
    m, n = _toint(a.num), _toint(b.num)
    Sym(factorial(big(m - 1)) * factorial(big(n - 1)) // factorial(big(m + n - 1)))
end

# J_n(0): 1 for n = 0, 0 for an integer n ≠ 0
bessel0(args) = iszero(args[2]) && isint(args[1]) ? (iszero(args[1]) ? ONE : ZERO) : nothing

## derivatives (∂ with respect to argument i; nothing where there is no closed form)

d1(f) = (args, i) -> f(args[1])
# for f(ν, x): only the x derivative
dx(f) = (args, i) -> i == 2 ? f(args[1], args[2]) : nothing

function __init__()
    reg(name; kw...) = register_function!(name; julia = getfield(SF, name), kw...)
    E2 = 2 / SQRTPI

    reg(:gamma; eval = gamma_value,
        diff = (args, i) -> length(args) == 1 ? SF.gamma(args[1]) * SF.digamma(args[1]) :
                            i == 2 ? -args[2]^(args[1] - 1) * exp(-args[2]) : nothing)
    reg(:loggamma; eval = _unary(a -> isone(a) || a == 2 ? ZERO : nothing), diff = d1(SF.digamma))
    reg(:digamma; diff = d1(SF.trigamma))
    reg(:trigamma; diff = d1(a -> SF.polygamma(2, a)))
    reg(:invdigamma; diff = d1(a -> inv(SF.trigamma(SF.invdigamma(a)))))
    reg(:polygamma; eval = polygamma_value, diff = dx((m, a) -> SF.polygamma(m + 1, a)))
    reg(:beta; eval = beta_value,
        diff = (args, i) -> SF.beta(args...) * (SF.digamma(args[i]) - SF.digamma(args[1] + args[2])))
    reg(:logbeta; diff = (args, i) -> SF.digamma(args[i]) - SF.digamma(args[1] + args[2]))
    reg(:logfactorial; eval = _unary(a -> isint(a) && a.num >= SymNumber(0) ? log(Sym(factorial(big(_toint(a.num))))) : nothing),
        diff = d1(a -> SF.digamma(a + 1)))
    reg(:zeta; eval = zeta_value,
        diff = (args, i) -> length(args) == 2 && i == 2 ? -args[1] * SF.zeta(args[1] + 1, args[2]) : nothing)
    reg(:eta)

    reg(:erf; eval = args -> length(args) == 1 ? odd_at0(:erf)(args) : nothing,
        diff = (args, i) -> length(args) == 1 ? E2 * exp(-args[1]^2) :
                            (i == 1 ? -1 : 1) * E2 * exp(-args[i]^2))
    reg(:erfc; eval = at0(ONE), diff = d1(a -> -E2 * exp(-a^2)))
    reg(:erfcx; eval = at0(ONE), diff = d1(a -> 2a * SF.erfcx(a) - E2))
    reg(:erfi; eval = odd_at0(:erfi), diff = d1(a -> E2 * exp(a^2)))
    reg(:erfinv; eval = odd_at0(:erfinv), diff = d1(a -> SQRTPI / 2 * exp(SF.erfinv(a)^2)))
    reg(:erfcinv; eval = _unary(a -> isone(a) ? ZERO : nothing), diff = d1(a -> -SQRTPI / 2 * exp(SF.erfcinv(a)^2)))
    reg(:dawson; eval = odd_at0(:dawson), diff = d1(a -> 1 - 2a * SF.dawson(a)))
    reg(:logerf; diff = (args, i) -> (i == 1 ? -1 : 1) * E2 * exp(-args[i]^2) / SF.erf(args[1], args[2]))
    reg(:logerfc; diff = d1(a -> -E2 * exp(-a^2) / SF.erfc(a)))
    reg(:logerfcx; diff = d1(a -> 2a - E2 / SF.erfcx(a)))
    reg(:faddeeva)

    reg(:expint; diff = (args, i) -> length(args) == 1 ? -exp(-args[1]) / args[1] :
                                     i == 2 ? -SF.expint(args[1] - 1, args[2]) : nothing)
    reg(:expintx)
    reg(:expinti; diff = d1(a -> exp(a) / a))
    reg(:sinint; eval = odd_at0(:sinint), diff = d1(a -> sin(a) / a))
    reg(:cosint; diff = d1(a -> cos(a) / a))

    reg(:airyai; diff = d1(SF.airyaiprime))
    reg(:airyaiprime; diff = d1(a -> a * SF.airyai(a)))
    reg(:airybi; diff = d1(SF.airybiprime))
    reg(:airybiprime; diff = d1(a -> a * SF.airybi(a)))
    foreach(reg, (:airyaix, :airyaiprimex, :airybix, :airybiprimex))

    reg(:besselj; eval = bessel0, diff = dx((ν, a) -> (SF.besselj(ν - 1, a) - SF.besselj(ν + 1, a)) / 2))
    reg(:bessely; diff = dx((ν, a) -> (SF.bessely(ν - 1, a) - SF.bessely(ν + 1, a)) / 2))
    reg(:besseli; eval = bessel0, diff = dx((ν, a) -> (SF.besseli(ν - 1, a) + SF.besseli(ν + 1, a)) / 2))
    reg(:besselk; diff = dx((ν, a) -> -(SF.besselk(ν - 1, a) + SF.besselk(ν + 1, a)) / 2))
    reg(:hankelh1; diff = dx((ν, a) -> (SF.hankelh1(ν - 1, a) - SF.hankelh1(ν + 1, a)) / 2))
    reg(:hankelh2; diff = dx((ν, a) -> (SF.hankelh2(ν - 1, a) - SF.hankelh2(ν + 1, a)) / 2))
    reg(:besselj0; eval = at0(ONE), diff = d1(a -> -SF.besselj1(a)))
    reg(:besselj1; eval = odd_at0(:besselj1), diff = d1(a -> SF.besselj0(a) - SF.besselj1(a) / a))
    reg(:bessely0; diff = d1(a -> -SF.bessely1(a)))
    reg(:bessely1; diff = d1(a -> SF.bessely0(a) - SF.bessely1(a) / a))
    foreach(reg, (:besseljx, :besselyx, :besselix, :besselkx, :besselh, :besselhx, :hankelh1x, :hankelh2x,
                  :sphericalbesselj, :sphericalbessely, :jinc))

    reg(:ellipk; eval = at0(PI / 2), diff = d1(m -> (SF.ellipe(m) - (1 - m) * SF.ellipk(m)) / (2m * (1 - m))))
    reg(:ellipe; eval = at0(PI / 2), diff = d1(m -> (SF.ellipe(m) - SF.ellipk(m)) / (2m)))
    nothing
end

end
