# Linear algebra for matrices of expressions: det, inv, lu.
#
# Julia's generic routines pivot by comparing magnitudes (abs(a) > abs(b)), which has no
# meaning for expressions. These eliminate exactly instead, over a field with an exact
# zero test, pivoting on the first nonzero entry:
#   - numbers: SymNumber (exact rationals);
#   - polynomials (det): fraction-free Bareiss elimination over Poly{SymNumber}, every
#     division exact, no gcds;
#   - otherwise: RatFun, rational functions in the atoms of the entries (x, sin(x),
#     sqrt(2), parameters), whose canonical form makes the zero test exact (as for linear
#     systems, see linsolve.jl).
# Results are converted back to Syms in lowest terms.

# An entry of a matrix with algebraic constants (sqrt(33), im): a Sym, normalized after
# each operation (expanded, radicals combined, denominators rationalized), so that the
# canonical form decides zero exactly. As rational functions, such constants would be
# independent atoms, and sqrt(33)^2 - 33 would not be zero.
struct _AlgF
    v::Sym
end
_algnorm(e::Sym) = expand(simplify_radicals(simplify_powers(expand(e))))
Base.:+(a::_AlgF, b::_AlgF) = _AlgF(_algnorm(a.v + b.v))
Base.:-(a::_AlgF, b::_AlgF) = _AlgF(_algnorm(a.v - b.v))
Base.:-(a::_AlgF) = _AlgF(-a.v)
Base.:*(a::_AlgF, b::_AlgF) = _AlgF(_algnorm(a.v * b.v))
Base.inv(a::_AlgF) = _AlgF(_algnorm(cancel(1 / a.v)))
Base.iszero(a::_AlgF) = iszero(a.v) || iszero(cancel(a.v))
Base.zero(::Type{_AlgF}) = _AlgF(ZERO)
Base.one(::Type{_AlgF}) = _AlgF(ONE)
Base.zero(::_AlgF) = _AlgF(ZERO)
Base.one(::_AlgF) = _AlgF(ONE)

# an atom that is an algebraic number: a radical of a number, im
_isalgebraic_atom(a::Sym) = !_hassymbol(a) && a.kind == POW && a.args[1].kind == NUMBER &&
    a.args[2].kind == NUMBER && !isone(_denominator(a.args[2].num))

# the entries over a field: (matrix, back) with back converting an entry to a Sym
function _field_matrix(A::AbstractMatrix)
    S = Sym[Sym(a) for a in A]
    if all(e -> e.kind == NUMBER && isexact(e.num), S)
        return map(e -> e.num, S), Sym
    end
    atoms = Set{Sym}()
    for e in S
        n, d = numer_denom(e)
        (_poly_atoms!(atoms, n) && _poly_atoms!(atoms, d)) ||
            throw(ArgumentError("matrix entries with floating point numbers together with symbols are not supported"))
    end
    any(_isalgebraic_atom, atoms) && return map(e -> _AlgF(_algnorm(e)), S), (a -> cancel(a.v))
    gens = sort!(collect(atoms); lt = _lt)
    isempty(gens) && (gens = Sym[Sym(:_linalg)])         # constants such as sqrt(2) only
    rf(e) = (nd = numer_denom(e); _ratfun(Poly{SymNumber}(nd[1], gens), Poly{SymNumber}(nd[2], gens)))
    map(rf, S), _from_ratfun
end

# a RatFun as a Sym, the denominator's first term positive (d/(a*d - b*c), not
# -d/(-a*d + b*c))
function _from_ratfun(r::RatFun)
    _ispoly(r) && return Sym(r.num)
    n, d = Sym(r.num), Sym(r.den)
    if d.kind == ADD && !isempty(d.coefs) && d.coefs[1] < N0
        n, d = -n, -d
    end
    n / d
end

# all entries polynomials with numeric denominators: the polynomial matrix, else nothing
function _poly_matrix(A::AbstractMatrix)
    S = Sym[Sym(a) for a in A]
    atoms = Set{Sym}()
    for e in S
        _hassymbol(numer_denom(e)[2]) && return nothing      # 1/x is not a polynomial atom
        _poly_atoms!(atoms, e) || return nothing
    end
    isempty(atoms) && return nothing
    any(_isalgebraic_atom, atoms) && return nothing        # sqrt(2)^2 is not an atom squared
    gens = sort!(collect(atoms); lt = _lt)
    map(e -> Poly{SymNumber}(e, gens), S)
end

_square(A) = (n = LinearAlgebra.checksquare(A); n)

# the identity: A - I and A + 2I for matrices of expressions (as for numbers, x + J.λ).
# A symbolic multiple of I cannot be a UniformScaling (whose scale is a Number): λ*I(n)
# is the diagonal matrix instead.
Base.:+(x::Sym, J::LinearAlgebra.UniformScaling) = x + J.λ
Base.:+(J::LinearAlgebra.UniformScaling, x::Sym) = J.λ + x
Base.:-(x::Sym, J::LinearAlgebra.UniformScaling) = x - J.λ
Base.:-(J::LinearAlgebra.UniformScaling, x::Sym) = J.λ - x
Base.:*(x::Sym, ::LinearAlgebra.UniformScaling) = throw(ArgumentError(
    "a symbolic multiple of I is not a UniformScaling (whose scale must be a Number); use $x*I(n), the n×n diagonal matrix"))
Base.:*(J::LinearAlgebra.UniformScaling, x::Sym) = x * J

# an expression times an array, elementwise (for numbers Base does this; a Sym is not a
# Number): 2x*A, A/x, λ*I(n)
Base.:*(x::Sym, A::AbstractArray) = broadcast(*, x, A)
Base.:*(A::AbstractArray, x::Sym) = broadcast(*, A, x)
Base.:/(A::AbstractArray, x::Sym) = broadcast(/, A, x)

"""
    det(A::AbstractMatrix{Sym})

The determinant of a matrix of expressions, exactly: fraction-free (Bareiss) elimination
when the entries are polynomials, Gaussian elimination over rational functions in the
entries' atoms otherwise, pivoting on the first nonzero entry. The result is in lowest
terms.

```julia
det([a b; c d])                  # a*d - b*c
det([x 1; 1 x])                  # x^2 - 1
```
"""
function LinearAlgebra.det(A::AbstractMatrix{Sym})
    n = _square(A)
    n == 0 && return ONE
    # the polynomials built inside, so that a wider exponent packing applies on a retry
    d = _with_widening(B -> (P = _poly_matrix(B); P === nothing ? nothing : Sym(_bareiss_det(P))), A)
    d === nothing || return expand(d)
    M, back = _field_matrix(A)
    back(_elim_det(M))
end

# Bareiss: M[i,j] <- (M[k,k] M[i,j] - M[i,k] M[k,j]) / prev, an exact division
function _bareiss_det(M::Matrix)
    n = size(M, 1)
    M = copy(M)
    sgn = 1
    prev = one(M[1, 1])
    for k in 1:(n - 1)
        if iszero(M[k, k])
            p = findfirst(i -> !iszero(M[i, k]), (k + 1):n)
            p === nothing && return zero(M[1, 1])
            p += k
            M[k, :], M[p, :] = M[p, :], M[k, :]
            sgn = -sgn
        end
        for i in (k + 1):n, j in (k + 1):n
            M[i, j] = divexact(M[k, k] * M[i, j] - M[i, k] * M[k, j], prev)
        end
        prev = M[k, k]
    end
    sgn == 1 ? M[n, n] : -M[n, n]
end

# Gaussian elimination over a field: the product of the pivots, with the swaps' sign
function _elim_det(M::Matrix)
    n = size(M, 1)
    M = copy(M)
    d = one(M[1, 1])
    for k in 1:n
        p = findfirst(i -> !iszero(M[i, k]), k:n)
        p === nothing && return zero(M[1, 1])
        p += k - 1
        if p != k
            M[k, :], M[p, :] = M[p, :], M[k, :]
            d = -d
        end
        d *= M[k, k]
        piv = inv(M[k, k])
        for i in (k + 1):n
            iszero(M[i, k]) && continue
            f = M[i, k] * piv
            for j in (k + 1):n
                M[i, j] -= f * M[k, j]
            end
        end
    end
    d
end

"""
    inv(A::AbstractMatrix{Sym})

The inverse of a matrix of expressions, exactly, by Gauss–Jordan elimination over
rational functions in the entries' atoms (exact rationals for numbers). A singular
matrix throws a `LinearAlgebra.SingularException`.

```julia
inv([a b; c d])                  # [d -b; -c a] / (a*d - b*c)
```
"""
Base.inv(A::StridedMatrix{Sym}) = _inv(A)       # (more specific than LinearAlgebra's)
Base.inv(A::AbstractMatrix{Sym}) = _inv(A)
function _inv(A::AbstractMatrix{Sym})
    n = _square(A)
    if n <= 6                 # polynomial entries: the adjugate, cofactors fraction-free
        R = _with_widening(B -> (P = _poly_matrix(B); P === nothing ? nothing : _adjugate_inv(P)), A)
        R === nothing || return R
    end
    _inv_field(A)
end

# inv(A) = adj(A) / det(A), the cofactors determinants by Bareiss
function _adjugate_inv(P::Matrix)
    n = size(P, 1)
    D = n == 1 ? P[1, 1] : _bareiss_det(P)
    iszero(D) && throw(LinearAlgebra.SingularException(n))
    Ds = Sym(D)
    n == 1 && return Sym[1 / Ds;;]
    R = Matrix{Sym}(undef, n, n)
    for i in 1:n, j in 1:n
        minor = P[setdiff(1:n, j), setdiff(1:n, i)]           # cofactor (j, i): the transpose
        c = n == 2 ? minor[1, 1] : _bareiss_det(minor)
        R[i, j] = _cancel_sign(isodd(i + j) ? -Sym(c) : Sym(c), Ds)
    end
    R
end
_cancel_sign(n::Sym, d::Sym) = (q = cancel(n / d); nd = numer_denom(q); _from_ratfun_sym(nd[1], nd[2]))
function _from_ratfun_sym(n::Sym, d::Sym)
    if d.kind == ADD && !isempty(d.coefs) && d.coefs[1] < N0
        n, d = -n, -d
    end
    n / d
end

function _inv_field(A::AbstractMatrix{Sym})
    n = _square(A)
    M, back = _field_matrix(A)
    T = eltype(M)
    R = Matrix{T}(undef, n, n)
    for i in 1:n, j in 1:n
        R[i, j] = i == j ? one(T) : zero(T)
    end
    for k in 1:n
        p = findfirst(i -> !iszero(M[i, k]), k:n)
        p === nothing && throw(LinearAlgebra.SingularException(k))
        p += k - 1
        if p != k
            M[k, :], M[p, :] = M[p, :], M[k, :]
            R[k, :], R[p, :] = R[p, :], R[k, :]
        end
        piv = inv(M[k, k])
        for j in 1:n                      # (loops: a RatFun is not a broadcasting scalar)
            M[k, j] *= piv
            R[k, j] *= piv
        end
        for i in 1:n
            (i == k || iszero(M[i, k])) && continue
            f = M[i, k]
            for j in 1:n
                M[i, j] -= f * M[k, j]
                R[i, j] -= f * R[k, j]
            end
        end
    end
    Sym[back(r) for r in R]
end

"""
    lu(A::AbstractMatrix{Sym}; check = true)

The LU factorization of a matrix of expressions, `A[p, :] == L * U`, as a
`LinearAlgebra.LU`: computed exactly over rational functions in the entries' atoms,
pivoting on the first nonzero entry of each column (there is no magnitude to compare).
`F.L`, `F.U`, `F.p`, `det(F)` and `F \\ b` work as usual. A singular matrix throws a
`LinearAlgebra.SingularException` unless `check = false`; `NoPivot()` factors without
row exchanges.
"""
function LinearAlgebra.lu(A::AbstractMatrix{Sym}, pivot::Union{LinearAlgebra.RowMaximum,LinearAlgebra.RowNonZero,LinearAlgebra.NoPivot} = LinearAlgebra.RowNonZero();
                          check::Bool = true, allowsingular::Bool = false)
    m, n = size(A)
    M, back = _field_matrix(A)
    ipiv = collect(1:min(m, n))
    info = 0
    for k in 1:min(m, n)
        p = pivot isa LinearAlgebra.NoPivot ? (iszero(M[k, k]) ? nothing : k) :
            (q = findfirst(i -> !iszero(M[i, k]), k:m); q === nothing ? nothing : q + k - 1)
        if p === nothing
            info == 0 && (info = k)
            continue
        end
        ipiv[k] = p
        if p != k
            M[k, :], M[p, :] = M[p, :], M[k, :]
        end
        piv = inv(M[k, k])
        for i in (k + 1):m
            iszero(M[i, k]) && continue
            M[i, k] *= piv
            f = M[i, k]
            for j in (k + 1):n
                M[i, j] -= f * M[k, j]
            end
        end
    end
    check && info > 0 && !allowsingular && throw(LinearAlgebra.SingularException(info))
    LinearAlgebra.LU{Sym,Matrix{Sym},Vector{Int}}(Sym[back(e) for e in M], ipiv, LinearAlgebra.BlasInt(info))
end

# A \ b exactly, through the LU factorization (Julia's generic `\` would pivot by
# comparing magnitudes)
function Base.:\(A::AbstractMatrix{Sym}, b::AbstractVecOrMat)
    size(A, 1) == size(A, 2) || throw(ArgumentError("\\ for a non-square matrix of expressions is not supported"))
    M, back = _field_matrix(hcat(A, b))
    n = size(A, 1)
    k = size(M, 2) - n
    for c in 1:n
        p = findfirst(i -> !iszero(M[i, c]), c:n)
        p === nothing && throw(LinearAlgebra.SingularException(c))
        p += c - 1
        p != c && ((M[c, :], M[p, :]) = (M[p, :], M[c, :]))
        piv = inv(M[c, c])
        for j in 1:size(M, 2)
            M[c, j] *= piv
        end
        for i in 1:n
            (i == c || iszero(M[i, c])) && continue
            f = M[i, c]
            for j in 1:size(M, 2)
                M[i, j] -= f * M[c, j]
            end
        end
    end
    X = Sym[back(M[i, n + j]) for i in 1:n, j in 1:k]
    b isa AbstractVector ? vec(X) : X
end

## row reduction, rank, null space

# reduced row echelon form over the field: (R, pivot columns)
function _rref_field(A::AbstractMatrix)
    M, back = _field_matrix(A)
    m, n = size(M)
    r = 1
    pivots = Int[]
    for c in 1:n
        r > m && break
        p = findfirst(i -> !iszero(M[i, c]), r:m)
        p === nothing && continue
        p += r - 1
        p != r && ((M[r, :], M[p, :]) = (M[p, :], M[r, :]))
        piv = inv(M[r, c])
        for j in 1:n
            M[r, j] *= piv
        end
        for i in 1:m
            (i == r || iszero(M[i, c])) && continue
            f = M[i, c]
            for j in 1:n
                M[i, j] -= f * M[r, j]
            end
        end
        push!(pivots, c)
        r += 1
    end
    Sym[back(e) for e in M], pivots
end

"""
    rref(A::AbstractMatrix)

The reduced row echelon form of a matrix of expressions, exactly. As for linear
systems, a symbolic pivot is taken to be nonzero (the generic case): `rref([a 1; 1 1])`
divides by `a - 1`.

```julia
rref(Sym[1 2 3; 4 5 6; 7 8 9])     # [1 0 -1; 0 1 2; 0 0 0]
```
"""
rref(A::AbstractMatrix) = _rref_field(A)[1]

"""
    rank(A::AbstractMatrix{Sym})

The rank of a matrix of expressions: its number of pivots in exact elimination (for
symbolic entries, the generic rank).
"""
LinearAlgebra.rank(A::AbstractMatrix{Sym}) = length(_rref_field(A)[2])

"""
    nullspace(A::AbstractMatrix{Sym})

A basis of the null space of a matrix of expressions, as the columns of a matrix: one
column for each free variable of the reduced row echelon form (generic, as for
[`rref`](@ref)).
"""
function LinearAlgebra.nullspace(A::AbstractMatrix{Sym})
    R, pivots = _rref_field(A)
    n = size(A, 2)
    free = [j for j in 1:n if !(j in pivots)]
    N = zeros(Sym, n, length(free))
    for (k, f) in enumerate(free)
        N[f, k] = ONE
        for (i, p) in enumerate(pivots)
            N[p, k] = -R[i, f]
        end
    end
    N
end

## eigenvalues and eigenvectors

"""
    SymJ.charpoly(A, λ)

The characteristic polynomial `det(λ*I - A)` of a square matrix of expressions, in `λ`.
Not exported.
"""
function charpoly(A::AbstractMatrix, λ::Sym)
    n = _square(A)
    expand(LinearAlgebra.det(λ * LinearAlgebra.I(n) - Sym[Sym(a) for a in A]))
end

# (eigenvalue, algebraic multiplicity) pairs, from the factored characteristic polynomial
function _eigen_mult(A::AbstractMatrix)
    λ = Sym(Symbol("_eigλ"))
    p = factor(charpoly(A, λ))
    fs = p.kind == MUL ? collect(zip(p.args, p.exps)) : p.kind == POW ? [(p.args[1], p.args[2])] : [(p, ONE)]
    out = Pair{Sym,Int}[]
    for (f, k) in fs
        has(f, λ) || continue
        m = Int(_numerator(k.num))
        for r in solve(f ~ 0, λ)
            i = findfirst(q -> iszero(expand(q.first - r)), out)
            i === nothing ? push!(out, r => m) : (out[i] = r => out[i].second + m)
        end
    end
    out
end

"""
    eigvals(A::AbstractMatrix{Sym})

The eigenvalues of a matrix of expressions, each repeated by its algebraic multiplicity:
the roots of the characteristic polynomial, exact where `solve` finds them exactly.

```julia
eigvals(Sym[2 1; 1 2])       # [1, 3]
eigvals([a b; b a])          # [a - b, a + b]
```
"""
LinearAlgebra.eigvals(A::AbstractMatrix{Sym}) = Sym[v for (v, m) in _eigen_mult(A) for _ in 1:m]

"""
    eigvecs(A::AbstractMatrix{Sym})

Eigenvectors of a matrix of expressions as the columns of a matrix: a basis of each
eigenspace, the null space of `A - v*I`, eigenvalue by eigenvalue (in the order of
`eigvals`). A defective matrix has fewer columns than rows.
"""
function LinearAlgebra.eigvecs(A::AbstractMatrix{Sym})
    n = _square(A)
    cols = Vector{Sym}[]
    for (v, _) in _eigen_mult(A)
        N = LinearAlgebra.nullspace(Sym[Sym(a) for a in A] - v * LinearAlgebra.I(n))
        for k in 1:size(N, 2)
            push!(cols, N[:, k])
        end
    end
    isempty(cols) ? zeros(Sym, n, 0) : reduce(hcat, cols)
end

"""
    eigen(A::AbstractMatrix{Sym})

The eigenvalues and eigenvectors of a diagonalizable matrix of expressions, as a
`LinearAlgebra.Eigen` (`F.values`, `F.vectors`, with `A * F.vectors == F.vectors *
Diagonal(F.values)`). A defective matrix (too few eigenvectors) is an `ArgumentError`;
`eigvals` and `eigvecs` still apply.
"""
function LinearAlgebra.eigen(A::AbstractMatrix{Sym})
    n = _square(A)
    vals, vecs = Sym[], Vector{Sym}[]
    for (v, m) in _eigen_mult(A)
        N = LinearAlgebra.nullspace(Sym[Sym(a) for a in A] - v * LinearAlgebra.I(n))
        size(N, 2) == m || throw(ArgumentError("the matrix is not diagonalizable: the eigenvalue $v has multiplicity $m but $(size(N, 2)) eigenvector$(size(N, 2) == 1 ? "" : "s"); use eigvals and eigvecs"))
        for k in 1:m
            push!(vals, v)
            push!(vecs, N[:, k])
        end
    end
    length(vals) == n || throw(ArgumentError("not all eigenvalues of the matrix were found"))
    LinearAlgebra.Eigen(vals, reduce(hcat, vecs))
end
