# Squarefree factorization over the rationals.
#
# A primitive polynomial f is split into its content with respect to its most
# significant generator x_i (a polynomial in the later generators) and its primitive
# part, whose irreducible factors all involve x_i. Yun's algorithm decomposes the
# primitive part using gcds with ∂f/∂x_i; the content is decomposed recursively. The
# two sets of factors are coprime, so factors of equal multiplicity are multiplied.

"""
    factor_squarefree(p::Poly) -> (unit, [f₁ => 1, f₂ => 2, ...])
    factor_squarefree(ex::Sym)

Squarefree factorization: `p == unit * f₁^m₁ * f₂^m₂ * ...` with the `fᵢ` squarefree,
pairwise coprime, of positive leading coefficient and integer content 1, and the
multiplicities `mᵢ` distinct and increasing. Coefficients must be exact.

For an expression, the result is the expression `unit * f₁^m₁ * ...`. Rational
expressions are first reduced with [`cancel`](@ref), and their numerator and
denominator are factored. As in
[`cancel`](@ref), function applications and non-integer powers are treated as
independent variables:

```julia
factor_squarefree(x^3 - x^2 - x + 1)               # (x + 1)*(x - 1)^2
factor_squarefree((x^2 - 1)/(x^2 + 2x + 1)^2)      # (x - 1)/(x + 1)^3
```
"""
factor_squarefree(p::Poly{SymNumber}; gcd_method::Symbol = :dense) =
    _with_widening(_factor_squarefree_poly, p; gcd_method)
function _factor_squarefree_poly(p::Poly{SymNumber}; gcd_method::Symbol = :dense)
    _check_gcd_method(gcd_method)
    iszero(p) && throw(ArgumentError("the zero polynomial has no squarefree factorization"))
    f = _integer_scaled(p)
    c = _int_content(f)
    f = _scale(f, inv(SymNumber(c)))
    sign(f.coefs[1]) < 0 && (f = -f)
    facs = _sqf(f, gcd_method)
    # the unit makes the product exact: p = unit * Π fᵢ^mᵢ
    lcprod = prod((g.coefs[1]^m for (g, m) in facs); init = N1)
    p.coefs[1] / lcprod, facs
end

"""
    squarefree_part(p)

The product of the distinct irreducible factors of `p` (its radical), with positive
leading coefficient and integer content 1. Also for expressions.
"""
squarefree_part(p::Poly{SymNumber}; gcd_method::Symbol = :dense) =
    _with_widening(_squarefree_part_poly, p; gcd_method)
function _squarefree_part_poly(p::Poly{SymNumber}; gcd_method::Symbol = :dense)
    _, facs = _factor_squarefree_poly(p; gcd_method)
    foldl((r, fm) -> r * fm.first, facs; init = one(p))
end

# f primitive with positive lc; returns its nonconstant factors by multiplicity
function _sqf(f::Poly{SymNumber}, method::Symbol = :dense)
    i = _mainvar(f, f)
    i == 0 && return Pair{Poly{SymNumber},Int}[]
    c = _content_mod(f, i, method)
    pp = _isunit(c) ? f : divexact(f, c)
    facs = _yun(pp, i, method)
    _isunit(c) && return facs
    _merge_factors(facs, _sqf(_positive(c), method))
end

# gcd of the coefficients with respect to x_i, with the modular gcd
function _content_mod(p::Poly{SymNumber}, i::Int, method::Symbol = :dense)
    g = zero(p)
    for c in _coeffs_in(p, i)
        g = _gcd_poly(g, c; method)
        _isunit(g) && break
    end
    g
end

function _diff_in(p::Poly{C}, i::Int) where {C}
    n = length(p.gens)
    M = eltype(p.monos)
    unit = one(M) << _shift(p.layout, n, i)
    monos, coefs = M[], C[]
    for (m, c) in zip(p.monos, p.coefs)
        e = _exponent(m, p.layout, n, i)
        e == 0 && continue
        push!(monos, m - unit)
        push!(coefs, c * SymNumber(e))
    end
    Poly{C}(p.gens, p.layout, monos, coefs)
end

# Yun's algorithm for f primitive in x_i (every factor involves x_i)
function _yun(f::Poly{SymNumber}, i::Int, method::Symbol = :dense)
    out = Pair{Poly{SymNumber},Int}[]
    fp = _diff_in(f, i)
    a = _gcd_poly(f, fp; method)
    b, c = divexact(f, a), divexact(fp, a)
    d = c - _diff_in(b, i)
    k = 1
    while b.monos[1] != 0                    # b is not yet a constant
        a = _gcd_poly(b, d; method)
        a.monos[1] == 0 || push!(out, _positive(a) => k)
        b, c = divexact(b, a), divexact(d, a)
        d = c - _diff_in(b, i)
        k += 1
    end
    out
end

function _merge_factors(a, b)
    out = Dict{Int,Poly{SymNumber}}()
    for (f, m) in Iterators.flatten((a, b))
        out[m] = haskey(out, m) ? out[m] * f : f
    end
    [out[m] => m for m in sort!(collect(keys(out)))]
end

## expressions

factor_squarefree(ex::Sym; gcd_method::Symbol = :dense) =
    _with_widening(_factor_squarefree_sym, ex; gcd_method)
function _factor_squarefree_sym(ex::Sym; gcd_method::Symbol = :dense)
    r = _sqf_parts(ex, gcd_method)
    r === nothing && return ex
    (un, fn), (ud, fd) = r
    _sqf_sym(un, fn) / _sqf_sym(ud, fd)
end

squarefree_part(ex::Sym; gcd_method::Symbol = :dense) =
    _with_widening(_squarefree_part_sym, ex; gcd_method)
function _squarefree_part_sym(ex::Sym; gcd_method::Symbol = :dense)
    r = _sqf_parts(ex, gcd_method)
    r === nothing && throw(ArgumentError("squarefree_part needs exact coefficients"))
    (_, fn), _ = r
    _sqf_sym(N1, [fm.first => 1 for fm in fn])
end

# factorizations of the numerator and denominator, or nothing for inexact coefficients
function _sqf_parts(ex::Sym, gcd_method::Symbol = :dense)
    n, d = numer_denom(_cancel(ex; gcd_method))
    atoms = Set{Sym}()
    (_poly_atoms!(atoms, n) && _poly_atoms!(atoms, d)) || return nothing
    gs = sort!(collect(atoms); lt = _lt)
    _factor_squarefree_poly(Poly{SymNumber}(n, gs); gcd_method), _factor_squarefree_poly(Poly{SymNumber}(d, gs); gcd_method)
end

_sqf_sym(u::SymNumber, facs) = foldl((r, fm) -> r * Sym(fm.first)^fm.second, facs; init = Sym(u))
