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
    (_poly_atoms!(atoms, n) && _poly_atoms!(atoms, d) && all(a -> _eq(a, x) || !_hassymbol(a), atoms)) || return nothing
    P = try Poly{Sym}(expand(n), [x]) catch; return nothing end
    deg = degree(P, x)
    deg == 0 && return Sym[]
    cs = Sym[coeff(n, x, k) for k in 0:deg]
    pole(r) = has(d, x) && abs(N(subs(d, x => r))) <= 1e-12 * max(1, abs(N(r)))
    # (correctly rounded real roots for rational coefficients; numeric ones otherwise)
    if domain == ℝ && all(ci -> ci.kind == NUMBER && isexact(ci.num), cs)
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
        zz = z  # (a new binding: a reassigned capture is boxed)
        pole(zz) || any(o -> _eq(o, zz), out) || push!(out, zz)
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

## systems: one solution, as Maple's fsolve for a system

"""
    SymJ.fsolve(eqs, vars, hints = nothing; domain = ℝ)

Solutions of the system of equations `eqs` (a vector, each written with `~`) in the
unknowns `vars`, numerically, as a vector of `Dict`s as `solve` gives for a system. A
polynomial system (polynomials or rational functions in the unknowns) gives all of its
solutions: the real ones (in the ranges given), or with `domain = ℂ` the complex ones
too, through `solve`'s exact elimination. Another system gives one solution, as Maple's
`fsolve` for a system. `eqs` must have no symbols but `vars`, and as many
equations as unknowns. Not exported.

For a system that is not polynomial, `hints` narrows the search, one pair for an unknown: `x => (a, b)` looks for x in
`[a, b]`, `x => x0` starts from `x0`. An unknown without one is looked for in
`[-10, 10]`. The search is Newton's method, damped, with the exact Jacobian, from the
starting values given, or else from points spread over the ranges (the same each call);
with `domain = ℂ` a complex starting value may lead to a complex solution. A system
whose search finds no solution is an error, not an empty vector, since one may exist
elsewhere.

```julia
@syms x y
SymJ.fsolve([x^2 + y^2 ~ 4, exp(x) ~ y], [x, y])                       # one of its two solutions
SymJ.fsolve([x^2 + y^2 ~ 4, exp(x) ~ y], [x, y], [x => (0, 2)])        # x ≈ 0.4...
SymJ.fsolve([sin(x) ~ y, x + y ~ 1], [x, y], [x => 0.5, y => 0.5])
SymJ.fsolve([x^3 + y^3 ~ 1, x^2 + y ~ 3], [x, y])                      # all its real solutions
```
"""
function fsolve(eqs::AbstractVector{<:Relation}, vars::AbstractVector, hints = nothing; domain = ℝ)
    (domain == ℝ || domain == ℂ) || throw(ArgumentError("fsolve: the domain is ℝ or ℂ, not $domain"))
    xs = Sym[Sym(v) for v in vars]
    all(_isnamed, xs) || throw(ArgumentError("fsolve: the unknowns must be symbols"))
    all(e -> e.op === :eq, eqs) || throw(ArgumentError("fsolve: every equation is written lhs ~ rhs"))
    n = length(xs)
    length(eqs) == n || throw(ArgumentError("fsolve: $(length(eqs)) equation$(length(eqs) == 1 ? "" : "s") in $n unknown$(n == 1 ? "" : "s"); give as many equations as unknowns"))
    F = Sym[e.lhs - e.rhs for e in eqs]
    others = [s for s in free_symbols(F) if !any(v -> _eq(v, s), xs)]
    isempty(others) || throw(ArgumentError("fsolve: the equations have symbols other than the unknowns ($(join(others, ", "))); substitute values for them"))
    any(_hasparameter, F) && throw(ArgumentError("fsolve: the equations have parameters; substitute values for them"))
    # ranges and starting values from the hints
    ranges = [(-10.0, 10.0) for _ in 1:n]
    start = Union{Nothing,Number}[nothing for _ in 1:n]
    for h in something(hints, Pair[])
        h isa Pair || throw(ArgumentError("fsolve: a hint is x => (a, b) or x => x0, not $h"))
        i = findfirst(v -> _eq(v, Sym(h.first)), xs)
        i === nothing && throw(ArgumentError("fsolve: $(h.first) is not an unknown"))
        if h.second isa Tuple
            a, b = Float64.(N.(Sym.(h.second)))
            ranges[i] = a <= b ? (a, b) : (b, a)
        else
            v = N(Sym(h.second))
            domain == ℝ && !(v isa Real) && throw(ArgumentError("fsolve: a complex starting value needs domain = ℂ"))
            start[i] = v
        end
    end
    # a polynomial system: all of its solutions, from solve's (Gröbner bases, then each
    # unknown's polynomial), evaluated, those in the domain and the ranges given
    if all(Fi -> _ispoly_system_eq(Fi, xs), F)
        r = _fsolve_poly_system(eqs, xs, ranges, start, domain)
        r === nothing || return r
    end
    f = lambdify(F, xs)
    J = lambdify(Sym[diff(Fi, xj) for Fi in F, xj in xs], xs)
    T = domain == ℂ ? ComplexF64 : Float64
    inrange(z) = domain == ℂ || all(i -> start[i] !== nothing || ranges[i][1] - 1e-9 <= real(z[i]) <= ranges[i][2] + 1e-9, 1:n)
    # starting points: the values given; for the rest, a Halton sequence over the ranges
    function point(k)
        z = Vector{T}(undef, n)
        for i in 1:n
            z[i] = start[i] !== nothing ? T(start[i]) : T(ranges[i][1] + (ranges[i][2] - ranges[i][1]) * _halton(k, _PRIMES[i]))
        end
        z
    end
    nstarts = all(s -> s !== nothing, start) ? 1 : 200
    for k in 1:nstarts
        z = _damped_newton(f, J, point(k), T)
        z === nothing && continue
        inrange(z) || continue
        sol = Dict{Sym,Sym}()
        for i in 1:n
            sol[xs[i]] = T === Float64 ? Sym(real(z[i])) :
                (abs(imag(z[i])) <= 1e-12 * max(1, abs(z[i])) ? Sym(real(z[i])) : Sym(real(z[i])) + Sym(imag(z[i])) * IM)
        end
        return [sol]
    end
    nstarts == 1 && throw(ArgumentError("fsolve: no solution found from the starting values given ($(join(("$(xs[i]) = $(start[i])" for i in 1:n), ", "))); give other ranges or starting values"))
    throw(ArgumentError("fsolve: no solution found in the ranges searched ($(join(("$(xs[i]) ∈ [$(ranges[i][1]), $(ranges[i][2])]" for i in 1:n if start[i] === nothing), ", "))); give other ranges or starting values"))
end

# Fi a polynomial (or a rational function) in the unknowns with numeric coefficients
function _ispoly_system_eq(Fi::Sym, xs::Vector{Sym})
    n, d = numer_denom(cancel(Fi))
    atoms = Set{Sym}()
    # (constant atoms, as 10^(1/20), are coefficients)
    _poly_atoms!(atoms, n) && _poly_atoms!(atoms, d) && all(a -> any(v -> _eq(a, v), xs) || !_hassymbol(a), atoms)
end

# every solution of a polynomial system, numerically (Dicts of floats, sorted), in the
# domain and the ranges given; nothing when solve finds none of its forms (an error, or
# Cases), so Newton's method is tried; an error for infinitely many solutions
function _fsolve_poly_system(eqs, xs::Vector{Sym}, ranges, start, domain)
    # (numeric roots asked for: solve's elimination with its high-degree factors' roots)
    sols = try with(() -> solve(collect(Relation, eqs), xs), _SOLVE_STRICT => false) catch err; err isa ArgumentError || rethrow(); return nothing end
    sols isa AbstractVector || return nothing
    out = Dict{Sym,Sym}[]
    for s in sols
        # (a free unknown: missing, or a value in the unknowns; x = 0 with y free for x y = 0,
        # x/y = 0 was "cannot evaluate free symbol y"; SymPy #4031)
        (all(v -> haskey(s, v), xs) && !any(v -> any(u -> has(s[v], u), xs), xs)) ||
            throw(ArgumentError("fsolve: the system has infinitely many solutions (a free unknown); solve gives them"))
        vals = [N(s[v]) for v in xs]
        any(v -> v isa Number && !isfinite(v), vals) && continue
        if domain == ℝ
            all(v -> v isa Real || abs(imag(v)) <= 1e-10 * max(1, abs(v)), vals) || continue
            vals = [real(v) for v in vals]
            rv = vals
            all(i -> start[i] !== nothing || ranges[i] == (-10.0, 10.0) || ranges[i][1] - 1e-9 <= rv[i] <= ranges[i][2] + 1e-9, eachindex(xs)) || continue
        end
        d = Dict{Sym,Sym}()
        for (v, w) in zip(xs, vals)
            d[v] = w isa Real ? Sym(w) : (abs(imag(w)) <= 1e-12 * max(1, abs(w)) ? Sym(real(w)) : Sym(real(w)) + Sym(imag(w)) * IM)
        end
        any(o -> all(v -> abs(N(o[v] - d[v])) <= 1e-10 * max(1, abs(N(d[v]))), xs), out) || push!(out, d)
    end
    sort!(out; by = d -> [(real(N(d[v])), imag(N(d[v]))) for v in xs])
end

const _PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71)

# the k-th point of the van der Corput sequence in base b (a coordinate of a Halton point)
function _halton(k::Int, b::Int)
    r, f = 0.0, 1.0
    while k > 0
        f /= b
        r += f * (k % b)
        k ÷= b
    end
    r
end

# Newton's method from z, each step halved until |F| decreases (at most 30 halvings), to
# |F| ≤ 1e-12 (relative to the equations' scale at z) with a step near rounding; nothing
# when it stalls, leaves the finite numbers, or the Jacobian is singular
function _damped_newton(f, J, z::Vector{T}, ::Type{T}) where {T}
    # (a complex value at a real point, x^y with x < 0, is no real point of F: as for a
    # DomainError; Float64 of it was an InexactError; SymPy #4027)
    val(w) = try T.(f(w...)) catch err; err isa Union{DomainError,DivideError,InexactError} || rethrow(); nothing end
    v = val(z)
    (v === nothing || !all(isfinite, v)) && return nothing
    scale = max(1.0, maximum(abs, v))
    for _ in 1:100
        r = maximum(abs, v)
        Jz = try T.(J(z...)) catch err; err isa Union{DomainError,DivideError,InexactError} || rethrow(); return nothing end
        all(isfinite, Jz) || return nothing
        δ = try Jz \ v catch err; err isa Union{LinearAlgebra.SingularException,LinearAlgebra.LAPACKException} || rethrow(); return nothing end
        all(isfinite, δ) || return nothing
        t = 1.0
        znew, vnew = z, v
        for _ in 1:30
            znew = z .- t .* δ
            vnew = val(znew)
            vnew !== nothing && all(isfinite, vnew) && maximum(abs, vnew) < r && break
            t /= 2
        end
        (vnew === nothing || !all(isfinite, vnew)) && return r <= 1e-10 * scale ? z : nothing
        # (at the root already, rounding keeps |F| from decreasing: converged)
        maximum(abs, vnew) >= r && return r <= 1e-10 * scale ? z : nothing
        step = maximum(abs, z .- znew)
        z, v = znew, vnew
        maximum(abs, v) <= 1e-12 * scale && step <= 1e-10 * max(1, maximum(abs, z)) && return z
    end
    nothing
end
