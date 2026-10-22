# Partial fraction decomposition over the rationals.
#
# For n/d in lowest terms: a polynomial part and a proper remainder r (n = q*d + r);
# d = c * Π p_i^m_i by `factor`; then r/d = Σ A_i / p_i^m_i with
# A_i = (r/c) * (d/p_i^m_i)^(-1) mod p_i^m_i (extended Euclid over Q, the factors being
# coprime), and A_i expanded in powers of p_i splits A_i / p_i^m_i into
# Σ c_k / p_i^j with deg c_k < deg p_i.

"""
    apart(ex, x)
    apart(ex)

Partial fraction decomposition of a rational function of `x` with exact numeric
coefficients: a polynomial plus terms `c(x)/p(x)^j`, `p` irreducible over the
rationals and `deg c < deg p`. With no `x`, `ex` must have one free symbol.

```julia
apart(1/(x^2 - 1))              # 1/(2*(x - 1)) - 1/(2*(x + 1))
apart((x^3 + 1)/(x^2 - 1))      # x + 1/(x - 1)
```
"""
function apart(ex::Sym, x::Sym)
    x.kind == SYMBOL || throw(ArgumentError("apart: $x is not a symbol"))
    n, d = numer_denom(cancel(ex))
    N, D = Poly{SymNumber}(n, [x]), Poly{SymNumber}(d, [x])
    q, r = divrem(N, D)
    result = Sym(q)
    iszero(r) && return result
    c, facs = factor(D)
    r = _scale(r, inv(c))
    D = _scale(D, inv(c))
    for (p, m) in facs
        P = p^m
        other = divexact(D, P)
        A = rem(r * _qinvmod(other, P), P)
        # A = Σ c_k p^k, so A / p^m = Σ c_k / p^(m - k)
        k = 0
        while !iszero(A)
            A, ck = divrem(A, p)
            iszero(ck) || (result += Sym(ck) / Sym(p)^(m - k))
            k += 1
        end
    end
    result
end

function apart(ex::Sym)
    vars = free_symbols(ex)
    length(vars) == 1 || throw(ArgumentError("apart: give the variable, e.g. apart(ex, x)"))
    apart(ex, only(vars))
end
apart(ex::Number, x...) = Sym(ex)

# a^(-1) mod m in Q[x], for coprime univariate a and m
function _qinvmod(a::Poly{SymNumber}, m::Poly{SymNumber})
    r0, r1 = m, rem(a, m)
    t0, t1 = zero(a), one(a)
    while !iszero(r1)
        q, r = divrem(r0, r1)
        r0, r1 = r1, r
        t0, t1 = t1, t0 - q * t1
    end
    r0.monos[1] == 0 || throw(ArgumentError("not coprime"))
    _scale(t0, inv(r0.coefs[1]))          # r0 is a nonzero constant
end
