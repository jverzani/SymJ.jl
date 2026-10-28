# Functional decomposition of univariate polynomials, f = g ∘ h, in characteristic 0 (the
# "tame" case: Kozen & Landau 1989, von zur Gathen 1990). For f monic of degree n = r·s,
# a right component h of degree s, monic with zero constant term, is unique if it exists:
# f = h^r + (terms of degree at most n - s), so the reversed h is the r-th root of the
# reversed f to order s. The left component g is then the h-adic expansion of f, found by
# dividing by h repeatedly; f decomposes with this h when every remainder is a constant.
# Taking the smallest such s gives an indecomposable h, and decomposing g in turn gives a
# complete decomposition. By Ritt's first theorem every complete decomposition has the
# same length, and the same degrees up to order.

"""
    SymJ.decompose(p, x)

A complete functional decomposition of the polynomial `p` in `x`: polynomials
`[g₁, g₂, …, gₖ]` in `x`, each of degree at least 2, with `p = g₁(g₂(⋯gₖ(x)))`, and none
of them decomposable. Each `gᵢ` after the first is monic with zero constant term. `[p]`
when `p` is not decomposable (of prime degree, for instance). The coefficients may be
symbolic. As SymPy's `decompose`. Not exported.

```julia
SymJ.decompose(x^4 + 2x^3 - x - 1, x)         # [x^2 - x - 1, x^2 + x]
SymJ.decompose(expand((x^2 + x)^3 - 2), x)    # [x^3 - 2, x^2 + x]
```

`solve` uses it for irreducible polynomials of degree above 4: when every component has
degree at most 4, the roots are found through the components, exactly.
"""
function decompose(p, x::Sym)
    _isnamed(x) || throw(ArgumentError("decompose: $x is not a symbol"))
    c = _xcoeffs_sym(Sym(p), x)
    length(c) <= 2 && return Sym[Sym(p)]
    [_from_xcoeffs(g, x) for g in _decompose(c)]
end

# the coefficients c[k+1] of x^k in p, a polynomial in x
function _xcoeffs_sym(p::Sym, x::Sym)
    P = try
        Poly{Sym}(expand(p), [x])
    catch
        throw(ArgumentError("decompose: $p is not a polynomial in $x"))
    end
    d = maximum(first, exponent_vectors(P); init = 0)
    c = fill(ZERO, d + 1)
    for (ev, a) in zip(exponent_vectors(P), coefficients(P))
        c[ev[1] + 1] += a
    end
    c
end

_from_xcoeffs(c::Vector{Sym}, x::Sym) = _add_many(c[k + 1] * x^k for k in eachindex(c) .- 1 if !iszero(c[k + 1]))

_dz(a::Sym) = iszero(a) || iszero(cancel(a))
_dnorm(a::Sym) = _hassymbol(a) ? cancel(a) : a

# a complete decomposition of the polynomial with coefficients c, as coefficient vectors
function _decompose(c::Vector{Sym})
    n = length(c) - 1
    for s in 2:(n ÷ 2)
        n % s == 0 || continue
        gh = _decompose_right(c, s)
        gh === nothing && continue
        g, h = gh
        return [_decompose(g); [h]]
    end
    [c]
end

# f = g ∘ h with h of degree s, monic with zero constant term: (g, h), or nothing
function _decompose_right(c::Vector{Sym}, s::Int)
    n = length(c) - 1
    r = n ÷ s
    lc = c[end]
    F = [_dnorm(a / lc) for a in c]                              # monic
    # the reversed F, P(y) = yⁿ F(1/y) = 1 + F[n-1] y + …, to order s; q = P^(1/r)
    P = [F[n + 1 - k] for k in 0:(s - 1)]
    q = Vector{Sym}(undef, s)
    q[1] = ONE
    α = Sym(1 // r)
    for k in 1:(s - 1)                                           # J. C. P. Miller's recurrence
        q[k + 1] = _dnorm(sum(((α + 1) * j - k) * P[j + 1] * q[k - j + 1] for j in 1:k) / k)
    end
    h = [ZERO; reverse(q)]                                       # x^s + q₁ x^(s-1) + … + q_(s-1) x
    # the h-adic expansion: F = Σ gᵢ hⁱ, each gᵢ a constant
    g = Sym[]
    rem = F
    for _ in 1:r
        quo, rm = _divmonic(rem, h)
        all(_dz, rm[2:end]) || return nothing
        push!(g, rm[1])
        rem = quo
    end
    length(rem) == 1 || return nothing
    push!(g, rem[1])
    ([_dnorm(lc * a) for a in g], h)
end

# division of a by the monic b (coefficient vectors, lowest degree first): (quotient, remainder)
function _divmonic(a::Vector{Sym}, b::Vector{Sym})
    m, d = length(a) - 1, length(b) - 1
    m < d && return (Sym[ZERO], copy(a))
    r = copy(a)
    q = fill(ZERO, m - d + 1)
    for k in (m - d):-1:0
        t = r[k + d + 1]
        q[k + 1] = t
        _dz(t) && continue
        for j in 0:d
            r[k + j + 1] = _dnorm(r[k + j + 1] - t * b[j + 1])
        end
    end
    (q, r[1:d])
end
