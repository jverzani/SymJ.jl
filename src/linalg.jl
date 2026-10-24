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
