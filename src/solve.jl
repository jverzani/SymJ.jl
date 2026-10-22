# Solving equations: CommonSolve's `solve` for a symbolic equation in one variable.
#
# The equation must be a Relation (`x^2 ~ 4`); a bare expression is not taken to mean
# "= 0". lhs - rhs is brought to lowest terms with `cancel`; the numerator's roots are
# the solutions (numerator and denominator then share no factor, so no root is a
# pole). The numerator is factored (other symbols and parameters may appear in the
# coefficients) and each irreducible factor solved:
#   degree 1 and 2: the usual formulas; a xⁿ + b: (-b/a)^(1/n) (-1)^(2k/n);
#   degree 3: Cardano; degree 4: Ferrari (biquadratic ones as quadratics in x²);
#   degree ≥ 5 with numeric coefficients: numeric roots (Aberth's iteration).

import CommonSolve: solve

"""
    solve(eq::Relation, x)

The solutions of the equation `eq` (written with `~`, e.g. `x^2 ~ 4`) for `x`, as a
vector of distinct expressions. The equation must be polynomial or rational in `x`;
its coefficients may involve other symbols and parameters. Solutions are exact through
degree 4 and for equations that factor into such pieces; irreducible factors of higher
degree with numeric coefficients give numeric (floating point) roots. Numeric
solutions are sorted, real ones first. An equation with no solution gives an empty
vector.

```julia
solve(x^2 - 2 ~ 0, x)          # [-sqrt(2), sqrt(2)]
solve(x^2 ~ -1, x)             # [-im, im]
solve(a*x + b ~ 0, x)          # [-b/a]
```
"""
function solve(eq::Relation, x::Sym)
    eq.op === :eq || throw(ArgumentError("solve: inequalities are not supported yet"))
    _isnamed(x) || throw(ArgumentError("solve: $x is not a symbol"))
    n, _ = numer_denom(cancel(eq.lhs - eq.rhs))
    if !has(n, x)
        iszero(n) && throw(ArgumentError("solve: every value of $x satisfies $eq"))
        return Sym[]
    end
    atoms = Set{Sym}()
    exact = _poly_atoms!(atoms, n)
    bad = [a for a in atoms if !_eq(a, x) && has(a, x)]
    isempty(bad) || return _sort_roots!(unique(_solve_through(eq, n, x, bad)))
    roots = exact ? _solve_exact(n, x, atoms) : _solve_inexact(n, x)
    _sort_roots!(unique(roots))
end
solve(eq::Relation, x::Sym, rest...) = throw(ArgumentError("solve: systems are not supported yet"))
solve(ex::Sym, x::Sym) = throw(ArgumentError(
    "solve needs an equation, e.g. solve($ex ~ 0, $x); an expression is not taken to mean \"= 0\""))

## equations in one function of x: f(g(x)) with f polynomial
#
# u stands for g(x); the polynomial in u is solved, then g(x) = r is inverted into an
# equation for g's argument and solved recursively. Inverting a power or a periodic
# function can introduce extraneous roots (sqrt(x) = -1), so candidates are checked
# by substitution. Trig equations give the principal solutions in one period.

function _solve_through(eq::Relation, n::Sym, x::Sym, bad::Vector{Sym})
    common = _common_power(bad)                     # exp(2x), exp(x): powers of one exp
    if common === nothing
        length(bad) == 1 || throw(ArgumentError(
            "solve: $eq involves $x through several functions ($(join(bad, ", "))); not supported yet"))
        g = only(bad)
        u = Sym(Symbol("_u", objectid(g)))
        nu = subs(n, g => u)
    else
        g, ks = common
        u = Sym(Symbol("_u", objectid(g)))
        nu = subs(n, (a => u^k for (a, k) in zip(bad, ks))...)
    end
    if has(nu, x)
        _issqrt(g) && return _solve_radical(eq, nu, u, g, x)
        throw(ArgumentError("solve: $eq involves both $x and $g; not supported (it may need, e.g., Lambert W)"))
    end
    out = Sym[]
    for r in solve(nu ~ 0, u), (arg, val) in _invert(g, r, x)
        append!(out, solve(arg ~ val, x))
    end
    filter(r -> _satisfies(eq, x, r), out)
end

# atoms b^(e_i) with one base b free of x and exponents rational multiples of one
# another: (b^g, k_i) with b^(e_i) = (b^g)^(k_i), k_i integers; else nothing
function _common_power(atoms::Vector{Sym})
    length(atoms) > 1 || return nothing
    all(a -> a.kind == POW, atoms) || return nothing
    b = atoms[1].args[1]
    es = [a.args[2] for a in atoms]
    if !all(a -> _eq(a.args[1], b), atoms)
        # integer bases that are powers of one integer: 4^x = 2^(2x)
        bases = [a.args[1] for a in atoms]
        all(c -> c.kind == NUMBER && _isint(c.num) && c.num > SymNumber(1), bases) || return nothing
        r = _common_integer_base([big(_numerator(c.num)) for c in bases])
        r === nothing && return nothing
        b = Sym(SymNumber(r))
        es = [e * Sym(SymNumber(_ilog(big(_numerator(c.num)), r))) for (e, c) in zip(es, bases)]
    end
    e1 = es[1]
    ratios = [cancel(e / e1) for e in es]
    all(r -> r.kind == NUMBER && isexact(r.num), ratios) || return nothing
    q = reduce(lcm, (big(_denominator(r.num)) for r in ratios); init = big(1))
    g = e1 / Sym(SymNumber(q))
    ks = [Int(_numerator(r.num) * (q ÷ _denominator(r.num))) for r in ratios]
    (pow(b, g), ks)
end

# the least integer r > 1 of which every base is a power, or nothing
function _common_integer_base(bs::Vector{BigInt})
    m = minimum(bs)
    for k in 64:-1:1
        r = _perfect_root(m, k)
        r === nothing && continue
        all(c -> _ilog(c, r) !== nothing, bs) && return r
    end
    nothing
end
function _ilog(c::BigInt, r::BigInt)                 # k with r^k == c, or nothing
    k = 0
    while c > 1 && c % r == 0
        c ÷= r
        k += 1
    end
    c == 1 ? k : nothing
end

_issqrt(g::Sym) = g.kind == POW && g.args[2] == HALF

# an equation in x and R = sqrt(t): reduced with R^2 = t to A + B R = 0, squared to
# A^2 = B^2 t (polynomial), solved, and the extraneous roots of squaring dropped
function _solve_radical(eq::Relation, nu::Sym, u::Sym, R::Sym, x::Sym)
    t = R.args[1]
    atoms = Set{Sym}()
    _poly_atoms!(atoms, nu)
    P = Poly{Sym}(nu, [u])                               # coefficients in x
    A, B = ZERO, ZERO
    for (ev, c) in zip(exponent_vectors(P), coefficients(P))
        k = ev[1]
        iseven(k) ? (A += c * t^(k ÷ 2)) : (B += c * t^((k - 1) ÷ 2))
    end
    cands = solve(expand(A^2 - B^2 * t) ~ 0, x)
    filter(r -> _satisfies(eq, x, r), cands)
end

# g(x) = r as equations arg ~ value for an argument of g
function _invert(g::Sym, r::Sym, x::Sym)
    k = g.kind
    if k == FN && length(g.args) == 1
        t = g.args[1]
        name = g.name
        name === :sin && return [(t, asin(r)), (t, PI - asin(r))]
        name === :cos && return [(t, acos(r)), (t, -acos(r))]
        name === :tan && return [(t, atan(r))]
        name === :log && return [(t, exp(r))]
        name === :asin && return [(t, sin(r))]
        name === :acos && return [(t, cos(r))]
        name === :atan && return [(t, tan(r))]
        name === :sinh && return [(t, asinh(r))]
        name === :cosh && return [(t, acosh(r)), (t, -acosh(r))]
        name === :tanh && return [(t, atanh(r))]
    elseif k == POW
        b, e = g.args
        if !has(b, x)                                       # b^t = r: t = log(r)/log(b)
            iszero(r) && return Tuple{Sym,Sym}[]
            return [(e, b === E ? log(r) : log(r) / log(b))]
        elseif !has(e, x) && e.kind == NUMBER && isexact(e.num)   # t^(p/q) = r
            return [(b, r^inv(e))]
        end
    end
    throw(ArgumentError("solve: cannot invert $g for $x"))
end

# does r satisfy eq? (exactly, or numerically for constants; kept when undecidable)
function _satisfies(eq::Relation, x::Sym, r::Sym)
    d = try
        subs(eq.lhs - eq.rhs, x => r)
    catch err
        err isa DivideError && return false
        rethrow()
    end
    (iszero(d) || iszero(cancel(d))) && return true
    _hassymbol(d) && return true
    v = try N(d) catch; return true end
    abs(v) <= 1e-9 * max(1, abs(N(r)))
end

function _solve_exact(n::Sym, x::Sym, atoms)
    gens = [x; sort!([a for a in atoms if !_eq(a, x)]; lt = _lt)]
    _, facs = factor(Poly{SymNumber}(n, gens))
    roots = Sym[]
    for (f, _) in facs
        degree(f, x) == 0 && continue
        append!(roots, _solve_irreducible(_xcoeffs(f), f))
    end
    roots
end

# coefficients c[k+1] of x^k (x the first generator), as expressions
_xcoeffs(f::Poly{SymNumber}) = [Sym(_coeff_of(f, 1, k)) for k in 0:_deg_in_all(f, 1)]

# numeric coefficients only (floats): numeric roots of the whole numerator
function _solve_inexact(n::Sym, x::Sym)
    c = try
        _xcoeffs(Poly{SymNumber}(n, [x]))
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("solve: inexact coefficients that are not numbers are not supported"))
    end
    _numeric_roots(c)
end

## an irreducible factor, by degree

function _solve_irreducible(c::Vector{Sym}, f)
    d = length(c) - 1
    d == 1 && return [-c[1] / c[2]]
    d == 2 && return _quadratic(c...)
    if all(iszero, c[2:(end-1)])                              # a xⁿ + b
        r = (-c[1] / c[end])^(1 // d)
        return [r * pow(MINUS_ONE, Sym(SymNumber(2k // d))) for k in 0:(d - 1)]
    end
    d == 3 && return _cardano(c...)
    d == 4 && return _ferrari(c...)
    all(ci -> !_hassymbol(ci), c) || throw(ArgumentError(
        "solve: an irreducible factor of degree $d with symbolic coefficients"))
    _numeric_roots(c)
end

function _quadratic(c0, c1, c2)
    s = sqrt(c1^2 - 4 * c2 * c0)
    [(-c1 - s) / (2c2), (-c1 + s) / (2c2)]
end

# is a constant expression zero? (exactly, or numerically when it has no symbols)
function _iszero_value(e::Sym)
    iszero(e) && return true
    _hassymbol(e) && return false
    v = try N(e) catch; return false end
    abs(v) < 1e-12
end

const _OMEGA = (-1 + sqrt(Sym(-3))) / 2                      # a primitive cube root of 1

function _cardano(d, c, b, a)                                # a x³ + b x² + c x + d
    # three distinct real roots and numeric coefficients: the trigonometric form, exact
    # and real (Cardano's radicals would need complex numbers: casus irreducibilis)
    if !any(_hassymbol, (a, b, c, d))
        p = (3a * c - b^2) / (3a^2)
        q = (2b^3 - 9a * b * c + 27a^2 * d) / (27a^3)
        disc = -(4p^3 + 27q^2)
        if N(disc) > 0
            θ = acos(3q / (2p) * sqrt(-3 / p)) / 3
            return [2sqrt(-p / 3) * cos(θ - 2PI * k / 3) - b / (3a) for k in 0:2]
        end
    end
    Δ0 = b^2 - 3a * c
    Δ1 = 2b^3 - 9a * b * c + 27a^2 * d
    s = sqrt(Δ1^2 - 4Δ0^3)
    C = ((Δ1 + s) / 2)^(1 // 3)
    _iszero_value(C) && (C = ((Δ1 - s) / 2)^(1 // 3))
    [-(b + _OMEGA^k * C + Δ0 / (_OMEGA^k * C)) / (3a) for k in 0:2]
end

function _ferrari(e, d, c, b, a)                             # a x⁴ + b x³ + c x² + d x + e
    p = (8a * c - 3b^2) / (8a^2)
    q = (b^3 - 4a * b * c + 8a^2 * d) / (8a^3)
    shift = -b / (4a)
    if _iszero_value(q)                                      # biquadratic in y = x - shift
        r = (-3b^4 + 256a^3 * e - 64a^2 * b * d + 16a * b^2 * c) / (256a^4)
        out = Sym[]
        for y2 in _quadratic(r, p, ONE), sgn in (-1, 1)
            push!(out, shift + sgn * sqrt(y2))
        end
        return out
    end
    Δ0 = c^2 - 3b * d + 12a * e
    Δ1 = 2c^3 - 9b * c * d + 27b^2 * e + 27a * d^2 - 72a * c * e
    Q = ((Δ1 + sqrt(Δ1^2 - 4Δ0^3)) / 2)^(1 // 3)
    _iszero_value(Q) && (Q = ((Δ1 - sqrt(Δ1^2 - 4Δ0^3)) / 2)^(1 // 3))
    S = sqrt(-2p / 3 + (Q + Δ0 / Q) / (3a)) / 2
    for _ in 1:2                                             # S ≠ 0 for some cube root of Q
        _iszero_value(S) || break
        Q *= _OMEGA
        S = sqrt(-2p / 3 + (Q + Δ0 / Q) / (3a)) / 2
    end
    u, v = sqrt(-4S^2 - 2p + q / S) / 2, sqrt(-4S^2 - 2p - q / S) / 2
    [shift - S - u, shift - S + u, shift + S - v, shift + S + v]
end

## numeric roots: Aberth–Ehrlich iteration

function _numeric_roots(c::Vector{Sym})
    a = ComplexF64[_evalc(ci) for ci in c]
    while length(a) > 1 && iszero(a[end])
        pop!(a)
    end
    n = length(a) - 1
    n == 0 && return Sym[]
    a ./= a[end]                                             # monic
    R = 1 + maximum(abs, a[1:(end-1)])                       # Cauchy's root bound
    z = [R * cis(2π * (k + 0.25) / n) for k in 0:(n - 1)]
    p(t) = evalpoly(t, a)
    dp(t) = evalpoly(t, [k * a[k+1] for k in 1:n])
    for _ in 1:500
        δmax = 0.0
        for k in 1:n
            w = p(z[k]) / dp(z[k])
            s = sum(1 / (z[k] - z[j]) for j in 1:n if j != k; init = 0.0im)
            δ = w / (1 - w * s)
            z[k] -= δ
            δmax = max(δmax, abs(δ) / max(1.0, abs(z[k])))
        end
        δmax < 1e-15 && break
    end
    map(z) do t
        abs(imag(t)) <= 1e-12 * max(1.0, abs(t)) ? Sym(real(t)) : Sym(real(t)) + Sym(imag(t)) * IM
    end
end

# numeric (constant) roots sorted: real ones first, by real then imaginary part
function _sort_roots!(roots::Vector{Sym})
    any(_hassymbol, roots) && return roots
    vals = map(r -> (v = try _evalc(r) catch; return nothing end; v), roots)
    any(isnothing, vals) && return roots
    key(v) = (abs(imag(v)) > 1e-10 * max(1, abs(v)), real(v), imag(v))
    roots[sortperm([key(v) for v in vals])]
end
