# Linear systems: Gauss–Jordan elimination over an exact field, SymNumber for numeric
# coefficients and RatFun (rational functions in the other symbols and parameters)
# otherwise; RatFun's canonical form makes the pivot zero tests exact. With parameters
# the result is the generic solution (pivots such as `a` assumed nonzero), split by
# cases on those pivots by `solve`.

"""
    solve(eqs::AbstractVector{<:Relation}, vars)

Solve a system of equations (written with `~`) for the variables `vars`. The result is a
vector of solutions, each a `Dict` with every one of `vars` as a key; an inconsistent system
gives an empty vector. Both arguments must be vectors (or tuples): one equation in one
unknown is `solve(eq, x)`, which gives a vector of values instead.

- **Linear systems** are solved by Gauss–Jordan elimination. An underdetermined system
  has one solution, with the pivot variables given in terms of the free ones, which map
  to themselves.
- **Polynomial (and rational) systems** are solved through a lex Gröbner basis
  ([`groebner`](@ref)), which is triangular, by back-substitution. The values are exact
  where the equations in one variable along the way can be solved exactly, numeric
  otherwise. With infinitely many solutions, the variables left free map to
  themselves.

Coefficients may involve parameters. Where the solution divides by an expression in
them (a pivot of the elimination, a coefficient of the Gröbner basis, a denominator of
a solution) that could vanish, the result is a `Cases` value, the case where it
vanishes solved again with it substituted; a case that cannot be solved so holds an
error. `assume` is not supported for systems yet.

```julia
solve([a*x + y ~ 1, x - y ~ 2], [x, y])
# a == 0 => [Dict(x => 3, y => 1)], a == -1 => [], a != 0, a != -1 => [Dict(x => 3/(a + 1), y => (1 - 2a)/(a + 1))]
```

```julia
solve([x + y ~ 3, x - y ~ 1], [x, y])          # [Dict(x => 2, y => 1)]
solve([x + y + z ~ 1, x - z ~ 0], [x, y, z])    # [Dict(x => z, y => 1 - 2z, z => z)]
solve([x^2 + y^2 ~ 1, x ~ y], [x, y])          # [Dict(x => -√2/2, y => -√2/2), Dict(x => √2/2, y => √2/2)]
solve([x*y ~ 1, x + y ~ 3], [x, y])
```
"""
function solve(eqs::AbstractVector{<:Relation}, vars::Union{AbstractVector,Tuple}; assume = nothing)
    assume === nothing ||
        throw(ArgumentError("solve: assumptions are not supported for systems yet; solve the system, then refine its solutions"))
    vars = Sym[Sym(v) for v in vars]
    all(_isnamed, vars) || return with(_NOT_PARAMS => unique!(reduce(vcat, map(_named, vars)))) do
        _solve_for_calls((es, vs, _) -> solve(es, vs), collect(Relation, eqs), vars)
    end
    all(e -> e.op === :eq, eqs) || return _solve_ineq_system(eqs, vars)
    _by_cases(Assumptions()) do B
        _solve_system_case(collect(Relation, eqs), vars, B)
    end
end

# the divisions the generic solution made (pivots, Gröbner coefficients), collected
const _DIVISORS = ScopedValue{Union{Nothing,Vector{Sym}}}(nothing)
_divisor!(c) = (D = _DIVISORS[]; D === nothing || push!(D, Sym(c)); nothing)

# the system in the case B: no solutions where a denominator free of the unknowns
# vanishes; the generic solutions where none of the divisors does; else the system again
# with the vanishing one substituted
function _solve_system_case(eqs::Vector{Relation}, vars::Vector{Sym}, B)
    sub = [k => w for (k, w) in B.subst if !any(v -> _eq(v, k), vars)]
    es = try
        isempty(sub) ? eqs : Relation[Relation(:eq, subs(e.lhs, sub...), subs(e.rhs, sub...)) for e in eqs]
    catch err
        _isdivzero(err) || rethrow()
        return Dict{Sym,Sym}[]
    end
    for e in es, f in _param_factors(numer_denom(cancel(e.lhs - e.rhs))[2], vars)
        _decide(Relation(:ne, f, ZERO), B) || return Dict{Sym,Sym}[]
    end
    D = Sym[]
    sols = with(() -> _solve_generic(es, vars), _DIVISORS => D)
    for s in sols, v in vars
        push!(D, numer_denom(s[v])[2])
    end
    ps = _degenerate_subst(D, B, vars)
    ps === nothing && return sols
    _alternatives(ps) do p
        _solve_system_case(Relation[Relation(:eq, subs(e.lhs, p), subs(e.rhs, p)) for e in es], vars, B)
    end
end

# the generic solutions, coefficients with parameters taken to be nonzero (what internal
# callers use)
function _solve_generic(eqs::AbstractVector{<:Relation}, vars::Union{AbstractVector,Tuple})
    vars = Sym[Sym(v) for v in vars]
    all(_isnamed, vars) ||
        return _solve_for_calls((es, vs, _) -> _solve_generic(es, vs), collect(Relation, eqs), vars)
    all(e -> e.op === :eq, eqs) || return _solve_ineq_system(eqs, vars)
    sols = if _islinear(eqs, vars)
        sol = _solve_linear(eqs, vars)
        sol === nothing ? Dict{Sym,Sym}[] : [sol]
    else
        _solve_nonlinear(collect(Relation, eqs), vars)
    end
    # every solution has every unknown as a key; a free one stands for itself
    for s in sols, v in vars
        haskey(s, v) || (s[v] = v)
    end
    filter!(s -> _defined_in(eqs, s, vars), sols)
end

# is every equation defined at the solution s? cancel can remove a factor of a
# denominator, so (x^2 - 1)/(x - 1) = 2 with y = x has no solution, not x = y = 1.
# (Only decided where the expressions that must not vanish become constants: with
# x free, x ≠ 1 is not recorded.)
function _defined_in(eqs, s::Dict{Sym,Sym}, vars::Vector{Sym})
    ps = [v => s[v] for v in vars]
    for eq in eqs
        e = eq.lhs - eq.rhs
        zs = Sym[]
        foreach(v -> _zero_sets!(zs, e, v), vars)
        for b in unique(zs)
            v = try subs(b, ps...) catch err; _isdivzero(err) && return false; rethrow() end
            (iszero(v) || (!_hassymbol(v) && _iszero_value(v))) && return false
        end
    end
    true
end

# every equation of degree at most 1 in the unknowns, with no unknown in a denominator
# or inside a function
function _islinear(eqs, vars::Vector{Sym})
    for eq in eqs
        n, d = numer_denom(cancel(eq.lhs - eq.rhs))
        any(v -> has(d, v), vars) && return false
        atoms = Set{Sym}()
        _poly_atoms!(atoms, _exactfloats(n)) || return false
        any(a -> !any(v -> _eq(a, v), vars) && any(v -> has(a, v), vars), atoms) && return false
        P = Poly{Sym}(n, vars)
        all(ev -> sum(ev) <= 1, exponent_vectors(P)) || return false
    end
    true
end

function _solve_linear(eqs, vars::Vector{Sym})
    A, b = _linear_system(eqs, vars)
    if all(e -> e.kind == NUMBER, Iterators.flatten((A, b)))
        _gauss_jordan(map(e -> e.num, A), map(e -> e.num, b), vars, Sym)
    else
        atoms = Set{Sym}()
        all(e -> _poly_atoms!(atoms, numer_denom(e)[1]) && _poly_atoms!(atoms, numer_denom(e)[2]),
            Iterators.flatten((A, b))) ||
            throw(ArgumentError("solve: inexact coefficients together with symbols are not supported"))
        gens = sort!(collect(atoms); lt = _lt)
        function rf(e)
            n, d = numer_denom(e)
            _ratfun(Poly{SymNumber}(n, gens), Poly{SymNumber}(d, gens))
        end
        _gauss_jordan(map(rf, A), map(rf, b), vars, Sym)
    end
end
solve(eqs::Tuple{Vararg{Relation}}, vars; kws...) = solve(collect(Relation, eqs), vars; kws...)

# A x = b from the equations: each lhs - rhs in lowest terms, no unknown in the
# denominator, of degree at most 1 in the unknowns
function _linear_system(eqs, vars::Vector{Sym})
    m, k = length(eqs), length(vars)
    A = fill(ZERO, m, k)
    b = fill(ZERO, m)
    for (i, eq) in enumerate(eqs)
        n, d = numer_denom(cancel(eq.lhs - eq.rhs))
        any(v -> has(d, v), vars) && throw(ArgumentError("solve: $eq is not linear in $(join(vars, ", "))"))
        P = Poly{Sym}(n, vars)
        for (ev, c) in zip(exponent_vectors(P), coefficients(P))
            s = sum(ev)
            s <= 1 || throw(ArgumentError("solve: $eq is not linear in $(join(vars, ", ")); nonlinear systems are not supported yet"))
            c = c / d
            s == 0 ? (b[i] = -c) : (A[i, findfirst(==(1), ev)] = c)
        end
    end
    A, b
end

# Gauss–Jordan over the field of the entries (SymNumber or RatFun)
function _gauss_jordan(A::Matrix{C}, b::Vector{C}, vars::Vector{Sym}, ::Type{Sym}) where {C}
    m, k = size(A)
    M = hcat(A, b)
    r = 1
    pivots = Tuple{Int,Int}[]
    for col in 1:k
        r > m && break
        p = findfirst(i -> !iszero(M[i, col]), r:m)
        p === nothing && continue
        p += r - 1
        if p != r
            M[r, :], M[p, :] = M[p, :], M[r, :]
        end
        C === SymNumber || _divisor!(M[r, col])
        piv = inv(M[r, col])
        M[r, :] = [e * piv for e in M[r, :]]
        for i in 1:m
            (i == r || iszero(M[i, col])) && continue
            f = M[i, col]
            M[i, :] = [M[i, j] - f * M[r, j] for j in 1:(k + 1)]
        end
        push!(pivots, (r, col))
        r += 1
    end
    if any(i -> !iszero(M[i, k+1]), r:m)                       # 0 = nonzero: inconsistent
        C === SymNumber || foreach(i -> iszero(M[i, k+1]) || _divisor!(M[i, k+1]), r:m)
        return nothing
    end
    pivcols = Set(c for (_, c) in pivots)
    free = [j for j in 1:k if !(j in pivcols)]
    sol = Dict{Sym,Sym}()
    for (i, c) in pivots
        v = Sym(M[i, k+1])
        for j in free
            iszero(M[i, j]) || (v -= Sym(M[i, j]) * vars[j])
        end
        sol[vars[c]] = cancel(v)
    end
    sol
end
