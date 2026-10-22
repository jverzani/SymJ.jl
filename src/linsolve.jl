# Linear systems: Gauss–Jordan elimination over an exact field, SymNumber for numeric
# coefficients and RatFun (rational functions in the other symbols and parameters)
# otherwise; RatFun's canonical form makes the pivot zero tests exact. With parameters
# the result is the generic solution (pivots such as `a` assumed nonzero).

"""
    solve(eqs::AbstractVector{<:Relation}, vars)

Solve a system of linear equations (written with `~`) for the variables `vars`. The
result is a `Dict` from variables to values; for an underdetermined system the pivot
variables are given in terms of the free ones, which do not appear as keys. An
inconsistent system gives `nothing`. Coefficients may involve other symbols and
parameters; the solution is then the generic one (it assumes that the pivots, such as
`a` in `a*x + y ~ 1`, are nonzero).

```julia
solve([x + y ~ 3, x - y ~ 1], [x, y])        # Dict(x => 2, y => 1)
solve([x + y + z ~ 1, x - z ~ 0], [x, y, z])  # Dict(x => z, y => 1 - 2z)
```
"""
function solve(eqs::AbstractVector{<:Relation}, vars::Union{AbstractVector,Tuple})
    vars = Sym[Sym(v) for v in vars]
    all(_isnamed, vars) || throw(ArgumentError("solve: the unknowns must be symbols or parameters"))
    all(e -> e.op === :eq, eqs) || throw(ArgumentError("solve: systems of inequalities are not supported"))
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
solve(eqs::Tuple{Vararg{Relation}}, vars) = solve(collect(Relation, eqs), vars)

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
    any(i -> !iszero(M[i, k+1]), r:m) && return nothing        # 0 = nonzero: inconsistent
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
