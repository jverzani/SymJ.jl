# Partial fraction decomposition over the rationals.
#
# For n/d in lowest terms: a polynomial part and a proper remainder r (n = q*d + r);
# d = c * Π p_i^m_i by `factor`; then r/d = Σ A_i / p_i^m_i with
# A_i = (r/c) * (d/p_i^m_i)^(-1) mod p_i^m_i (extended Euclid over Q, the factors being
# coprime), and A_i expanded in powers of p_i splits A_i / p_i^m_i into
# Σ c_k / p_i^j with deg c_k < deg p_i.

"""
    apart(ex, x; gcd_method = :dense)
    apart(ex; gcd_method = :dense)

Partial fraction decomposition of a rational function of `x` with exact
coefficients: a polynomial plus terms `c(x)/p(x)^j`, `p` irreducible and
`deg c < deg p`. Coefficients may involve other symbols and parameters (and, as in
[`cancel`](@ref), function applications): the decomposition is then over the field
of rational functions in them, `ℚ(a, b, ...)`. With no `x`, `ex` must have one free
variable. `gcd_method = :zippel` uses Zippel's sparse gcd (see [`gcd`](@ref)) in the
initial [`cancel`](@ref) and in factoring the denominator.

```julia
apart(1/(x^2 - 1))              # 1/(2*(x - 1)) - 1/(2*(x + 1))
apart((x^3 + 1)/(x^2 - 1))      # x + 1/(x - 1)
apart(1/((x - a)*(x - b)), x)   # 1/((a - b)*(x - a)) - 1/((a - b)*(x - b))
```
"""
function apart(ex::Sym, x::Sym; gcd_method::Symbol = :dense)
    _isnamed(x) || throw(ArgumentError("apart: $x is not a symbol"))
    n, d = numer_denom(cancel(ex; gcd_method))
    atoms = Set{Sym}()
    (_poly_atoms!(atoms, n) && _poly_atoms!(atoms, d)) ||
        throw(ArgumentError("apart needs exact coefficients"))
    delete!(atoms, x)
    isempty(atoms) || return _apart_params(n, d, x, sort!(collect(atoms); lt = _lt), gcd_method)
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

function apart(ex::Sym; gcd_method::Symbol = :dense)
    vars = free_symbols(ex)
    length(vars) == 1 || throw(ArgumentError("apart: give the variable, e.g. apart(ex, x)"))
    apart(ex, only(vars); gcd_method)
end
apart(ex::Number, x...; gcd_method::Symbol = :dense) = Sym(ex)

# Over ℚ(params): D factored over ℤ in (x, params); by Gauss's lemma its factors
# involving x are its irreducible factors over ℚ(params)[x], the rest a unit there.
# The same steps as above, in Poly{RatFun}.
function _apart_params(n::Sym, d::Sym, x::Sym, params::Vector{Sym}, gcd_method::Symbol = :dense)
    gs = [x; params]
    N, D = Poly{SymNumber}(n, gs), Poly{SymNumber}(d, gs)
    c, facs = factor(D; gcd_method)
    unit = Poly{SymNumber}(Sym(c), gs)
    xfacs = Pair{Poly{SymNumber},Int}[]
    for (p, m) in facs
        degree(p, x) > 0 ? push!(xfacs, p => m) : (unit *= p^m)
    end
    u = _to_ratfun_poly(unit, 1).coefs[1]                    # a constant in x
    Nx = _scale(_to_ratfun_poly(N, 1), inv(u))
    Dx = _to_ratfun_poly(foldl((r, pm) -> r * pm.first^pm.second, xfacs; init = one(N)), 1)
    q, r = divrem(Nx, Dx)
    result = Sym(q)
    iszero(r) && return result
    for (p, m) in xfacs
        px = _to_ratfun_poly(p, 1)
        P = px^m
        A = rem(r * _qinvmod(divexact(Dx, P), P), P)
        k = 0
        while !iszero(A)
            A, ck = divrem(A, px)
            iszero(ck) || (result += _over_common_den(ck) / Sym(p)^(m - k))
            k += 1
        end
    end
    result
end

# c(x) with coefficients in ℚ(params) as (polynomial)/(factored common denominator)
function _over_common_den(c::Poly{RatFun})
    L = foldl((l, r) -> lcm(l, r.den), c.coefs; init = one(c.coefs[1].den))
    num = _add_many(mul(Sym(r.num * divexact(L, _unify(r.den, L)[1])), _monomial(c, m))
                    for (m, r) in zip(c.monos, c.coefs))
    L.monos[1] == 0 ? num / Sym(L) : factor(num) / factor(Sym(L))
end

# a^(-1) mod m in K[x] (K = ℚ or ℚ(params)), for coprime univariate a and m
function _qinvmod(a::Poly{C}, m::Poly{C}) where {C}
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
