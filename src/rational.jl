# Rational expressions: common denominators and cancellation of common factors.

"""
    numer_denom(ex) -> (numerator, denominator)

Write `ex` as a single fraction. Rational coefficients and factors with negative
exponents move to the denominator, and sums are brought over a common denominator.
Nothing is expanded or cancelled; see [`cancel`](@ref).

```julia
numer_denom(x/y + 1)       # (x + y, y)
numer_denom(3x/4)          # (3*x, 4)
```
"""
function numer_denom(ex::Sym)
    k = ex.kind
    if k == NUMBER
        isexact(ex.num) || return (ex, ONE)
        return Sym(SymNumber(_numerator(ex.num))), Sym(SymNumber(_denominator(ex.num)))
    elseif k == ADD
        n, d = numer_denom(Sym(ex.num))
        for i in eachindex(ex.args)
            ni, di = numer_denom(_with_coef(ex.coefs[i], ex.args[i]))
            if _eq(d, di)
                n = n + ni
            else
                n, d = n * di + ni * d, d * di
            end
        end
        return n, d
    elseif k == MUL
        n, d = numer_denom(Sym(ex.num))
        for i in eachindex(ex.args)
            ni, di = numer_denom(_mkpow(ex.args[i], ex.exps[i]))
            n, d = n * ni, d * di
        end
        return n, d
    elseif k == POW
        b, e = ex.args
        if e.kind == NUMBER && _isint(e.num)
            nb, db = numer_denom(b)
            m = _toint(e.num)
            return m > 0 ? (nb^m, db^m) : (db^(-m), nb^(-m))
        end
        _negexp(e) && return (ONE, pow(b, -e))
    end
    ex, ONE
end
numer_denom(x::Number) = numer_denom(Sym(x))

# Collect the generators for viewing ex as a polynomial: symbols, constants, function
# applications and non-integer powers. Returns false for inexact (float) coefficients.
function _poly_atoms!(s::Set{Sym}, ex::Sym)
    k = ex.kind
    k == NUMBER && return isexact(ex.num)
    if k == ADD
        isexact(ex.num) && all(isexact, ex.coefs) || return false
        return all(t -> _poly_atoms!(s, t), ex.args)
    elseif k == MUL
        isexact(ex.num) || return false
        for i in eachindex(ex.args)
            if _nonneg_int(ex.exps[i])
                _poly_atoms!(s, ex.args[i]) || return false
            else
                push!(s, _mkpow(ex.args[i], ex.exps[i]))
            end
        end
        return true
    elseif k == POW && _nonneg_int(ex.args[2])
        return _poly_atoms!(s, ex.args[1])
    end
    push!(s, ex)
    true
end

"""
    cancel(ex; gcd_method = :dense)

Write `ex` as `p/q` with `p` and `q` expanded polynomials that have no common factor.
`q` has a positive leading coefficient and the integer content is shared out, so a
polynomial input comes back expanded.

Symbols, function applications such as `sin(x)`, and non-integer powers such as
`sqrt(x)` are treated as independent variables. This is always valid, but relations
among them are not used: `(x - 1)/(sqrt(x) - 1)` is not reduced to `sqrt(x) + 1`.
Expressions with floating point coefficients are returned unchanged.
`gcd_method = :zippel` uses Zippel's sparse gcd (see [`gcd`](@ref)).

```julia
cancel((x^2 - 1)/(x + 1))                  # x - 1
cancel(1/x + 1/x^2)                        # (x + 1)/x^2
cancel((sin(x)^2 - 1)/(sin(x) + 1))        # sin(x) - 1
```
"""
cancel(ex::Sym; gcd_method::Symbol = :dense) = _with_widening(_cancel, ex; gcd_method)
function _cancel(ex::Sym; gcd_method::Symbol = :dense)
    _check_gcd_method(gcd_method)
    n, d = numer_denom(ex)
    atoms = Set{Sym}()
    (_poly_atoms!(atoms, n) && _poly_atoms!(atoms, d)) || return ex
    gs = sort!(collect(atoms); lt = _lt)
    P, Q = Poly{SymNumber}(n, gs), Poly{SymNumber}(d, gs)
    g = _gcd_poly(P, Q; method = gcd_method)
    P, Q = divexact(P, g), divexact(Q, g)
    if sign(Q.coefs[1]) < 0
        P, Q = -P, -Q
    end
    if Q.monos[1] == 0                                 # constant denominator
        return Sym(_scale(P, inv(Q.coefs[1])))
    end
    Sym(P) / Sym(Q)
end
cancel(x::Number; gcd_method::Symbol = :dense) = Sym(x)
