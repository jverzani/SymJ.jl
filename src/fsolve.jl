# fsolve: the numeric roots of one equation in one unknown, as Maple's fsolve. solve stays
# exact (or an error); fsolve is the numeric counterpart. Polynomials (and rational
# functions, through their numerators) by their exact roots' numeric values (real_roots,
# _numeric_roots); other equations by the Chebyshev-proxy rootfinder of chebzeros.jl.

"""
    SymJ.fsolve(eq, x, (a, b) = (-Inf, Inf); domain = ℝ)

The roots of the equation `eq` (written with `~`) in `x`, numerically, as Maple's
`fsolve`: by default the real roots in the interval `[a, b]` (the whole line when it is
not given; either end may be infinite), sorted, each once, as floats. `eq` must have no
symbol but `x`. Not exported.

With `domain = ℂ` the complex roots too: for a polynomial (or a rational function) all of
them; for another equation the real ones and the complex ones near the real axis that
the rootfinder can vouch for (each polished on the equation itself), so complex roots far
from the axis may be missing. An interval is for real roots only: it is not taken
together with `domain = ℂ`.

Polynomials are solved through their exact roots (correctly rounded real ones). Other
equations are approximated on the interval by Chebyshev interpolants, piece by piece, and
the roots of those checked against the equation: a feature narrower than the samples'
spacing can be missed. The search stops at about 120 roots, those nearest the middle of
the interval (the origin on the whole line), with a warning: `tan(x) ~ x` has infinitely
many.

```julia
@syms x
SymJ.fsolve(cos(x) ~ x, x)                       # [0.7390851332151607]
SymJ.fsolve(log(x) + sin(x) ~ 1, x, (0, 4))      # [1.1099…, 3.3530…]
SymJ.fsolve(x^5 - x + 1 ~ 0, x)                  # [-1.1673039782614187]
SymJ.fsolve(x^5 - x + 1 ~ 0, x; domain = ℂ)      # five roots
```
"""
function fsolve(eq::Relation, x::Sym, interval::Tuple{Any,Any} = (-Inf, Inf); domain = ℝ)
    eq.op === :eq || throw(ArgumentError("fsolve: $eq is not an equation (write lhs ~ rhs)"))
    _isnamed(x) || throw(ArgumentError("fsolve: $x is not a symbol"))
    (domain == ℝ || domain == ℂ) || throw(ArgumentError("fsolve: the domain is ℝ or ℂ, not $domain"))
    given = interval != (-Inf, Inf)
    domain == ℂ && given && throw(ArgumentError("fsolve: an interval is for real roots; give one or domain = ℂ, not both"))
    lo, hi = Float64.(N.(Sym.(interval)))
    lo > hi && ((lo, hi) = (hi, lo))
    ex = eq.lhs - eq.rhs
    others = [s for s in free_symbols(ex) if !_eq(s, x)]
    isempty(others) || throw(ArgumentError("fsolve: $eq has symbols other than $x ($(join(others, ", "))); substitute values for them"))
    _hasparameter(ex) && throw(ArgumentError("fsolve: $eq has parameters; substitute values for them"))
    has(ex, x) || throw(ArgumentError("fsolve: $eq does not involve $x"))
    r = _fsolve_poly(ex, x, lo, hi, domain)
    r === nothing || return r
    _fsolve_cheb(ex, x, lo, hi, domain)
end
fsolve(eq::Relation, x::Sym, a::Real, b::Real; kws...) = fsolve(eq, x, (a, b); kws...)

# a polynomial or a rational function in x: the numerator's roots (not the denominator's),
# real ones in [lo, hi] (domain ℝ) or all of them (ℂ); nothing for another equation
function _fsolve_poly(ex::Sym, x::Sym, lo::Float64, hi::Float64, domain)
    n, d = numer_denom(cancel(ex))
    atoms = Set{Sym}()
    (_poly_atoms!(atoms, n) && _poly_atoms!(atoms, d) && all(a -> _eq(a, x), atoms)) || return nothing
    P = try Poly{Sym}(expand(n), [x]) catch; return nothing end
    deg = degree(P, x)
    deg == 0 && return Sym[]
    cs = Sym[coeff(n, x, k) for k in 0:deg]
    pole(r) = has(d, x) && abs(N(subs(d, x => r))) <= 1e-12 * max(1, abs(N(r)))
    if domain == ℝ && _isexactex(n)
        rs = unique(real_roots(n, x))
        return Sym[Sym(r) for r in rs if lo <= r <= hi && !pole(Sym(r))]
    end
    zs = _numeric_roots(cs)
    out = Sym[]
    for z in zs
        v = N(z)
        if domain == ℝ
            (v isa Real || abs(imag(v)) <= 1e-12 * max(1, abs(v))) || continue
            w = real(v)
            lo <= w <= hi || continue
            z = Sym(w)
        end
        pole(z) || any(o -> _eq(o, z), out) || push!(out, z)
    end
    _sort_roots!(out)
end

# another equation: the Chebyshev-proxy rootfinder on [lo, hi] (real roots), or with the
# trusted complex roots near the real axis (ℂ)
function _fsolve_cheb(ex::Sym, x::Sym, lo::Float64, hi::Float64, domain)
    f = lambdify(ex, x)
    if domain == ℝ
        return Sym[Sym(r) for r in ChebZeros.cheb_zeros(f, lo, hi)]
    end
    r = ChebZeros.cheb_all_zeros(f, lo, hi)
    _sort_roots!([Sym[Sym(v) for v in r.real]; Sym[Sym(real(c.z)) + Sym(imag(c.z)) * IM for c in r.complex]])
end
