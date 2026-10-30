# Linear algebra for matrices of expressions: det, inv, lu, \, rref, rank, nullspace,
# eigvals, eigvecs, eigen, qr, the Jordan form, exp, sin, cos and A^n.
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

# an atom that is algebraic over the others: a radical (sqrt(2), im, sqrt(-t^2)), whose
# powers reduce in the canonical form
_isalgebraic_atom(a::Sym) = a.kind == POW && a.args[2].kind == NUMBER && isexact(a.args[2].num) &&
    !isone(_denominator(a.args[2].num))

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
    SymJ.rref(A::AbstractMatrix)

The reduced row echelon form of a matrix of expressions, exactly. As for linear
systems, a symbolic pivot is taken to be nonzero (the generic case): `rref([a 1; 1 1])`
divides by `a - 1`. Not exported (it is not a LinearAlgebra function).

```julia
SymJ.rref(Sym[1 2 3; 4 5 6; 7 8 9])     # [1 0 -1; 0 1 2; 0 0 0]
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
[`SymJ.rref`](@ref)).
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

## dot and norm, entrywise as for numbers. (The generic methods iterate an entry, and an
## expression iterates as itself.) Symbols are complex: dot(u, v) = Σ conj(uᵢ) vᵢ and
## norm(v) = sqrt(Σ |vᵢ|²); for real symbols refine, under assumptions, removes the conj
## and the abs.

LinearAlgebra.dot(a::Sym, b::Sym) = conj(a) * b
LinearAlgebra.dot(a::Sym, b::Number) = conj(a) * b
LinearAlgebra.dot(a::Number, b::Sym) = conj(a) * b
LinearAlgebra.norm(a::Sym, p::Real = 2) = abs(a)
function LinearAlgebra.norm(v::AbstractArray{Sym}, p::Real = 2)
    isempty(v) && return ZERO
    p == 2 && return sqrt(_add_many(abs2(a) for a in v))
    p == 1 && return _add_many(abs(a) for a in v)
    p == Inf && return reduce(max, abs.(v))
    p == -Inf && return reduce(min, abs.(v))
    p > 0 || throw(ArgumentError("norm: p = $p is not supported for expressions"))
    _add_many(abs(a)^p for a in v)^(1 / Sym(p))
end
LinearAlgebra.normalize(v::AbstractArray{Sym}, p::Real = 2) = v ./ LinearAlgebra.norm(v, p)

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
        for r in _eig_roots(f, λ)
            i = findfirst(q -> iszero(expand(q.first - r)), out)
            i === nothing ? push!(out, r => m) : (out[i] = r => out[i].second + m)
        end
    end
    # constants in a definite order (by real part, then imaginary): the factors' order
    # is not the same on every Julia version
    if all(p -> !_hassymbol(p.first), out)
        key(p) = (v = _evalc(p.first); (real(v), imag(v)))
        sort!(out; by = key)
    end
    out
end

# the roots of a factor of the characteristic polynomial; a quadratic whose
# discriminant is minus a square, -u^2, has the roots (-b ± i u)/2a (rather than with
# sqrt(-u^2); both roots change together, so the pair is the same)
function _eig_roots(f::Sym, λ::Sym)
    P = try Poly{Sym}(expand(f), [λ]) catch; nothing end
    if P !== nothing && degree(P, λ) == 2
        c = Dict(ev[1] => cf for (ev, cf) in zip(exponent_vectors(P), coefficients(P)))
        a2, b1, c0 = get(c, 2, ZERO), get(c, 1, ZERO), get(c, 0, ZERO)
        u = _square_root(expand(4 * a2 * c0 - b1^2))           # -D = u^2
        u === nothing || return Sym[(-b1 - IM * u) / (2a2), (-b1 + IM * u) / (2a2)]
    end
    solve(f ~ 0, λ)
end
# u with e == u^2 (a square of a polynomial times a rational square), else nothing
function _square_root(e::Sym)
    _hassymbol(e) || return nothing
    q = factor(e)
    fs = q.kind == MUL ? collect(zip(q.args, q.exps)) : q.kind == POW ? [(q.args[1], q.args[2])] : [(q, ONE)]
    c = q.kind == MUL ? Sym(q.num) : ONE
    out = ONE
    for (b, k) in fs
        (k.kind == NUMBER && isexact(k.num) && isone(_denominator(k.num)) && iseven(_numerator(k.num))) || return nothing
        out *= b^(k / 2)
    end
    c.kind == NUMBER && c.num > N0 || return nothing
    r = sqrt(c)
    r.kind == NUMBER || return nothing
    r * out
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

## QR

"""
    SymQR

The QR factorization of a matrix of expressions from [`qr`](@ref LinearAlgebra.qr):
`F.Q` with orthonormal columns and `F.R` upper triangular, `A == Q * R`; `Q, R = F`
destructures it. Not exported.
"""
struct SymQR
    Q::Matrix{Sym}
    R::Matrix{Sym}
end
Base.iterate(F::SymQR, s = 1) = s == 1 ? (F.Q, 2) : s == 2 ? (F.R, 3) : nothing
Base.show(io::IO, ::MIME"text/plain", F::SymQR) = (println(io, "SymQR, Q:"); show(io, MIME"text/plain"(), F.Q);
    println(io, "\nR:"); show(io, MIME"text/plain"(), F.R))

# the inner product of columns: symbols taken as real, numbers with im conjugated
_dot(u, v) = _add_many((_has_im(a) ? conj(a) : a) * b for (a, b) in zip(u, v))
_has_im(e::Sym) = has(e, IM)
_qrnorm(e::Sym) = (c = cancel(expand(e)); simplify(c))

"""
    qr(A::AbstractMatrix{Sym})

The (thin) QR factorization of a matrix of expressions with linearly independent
columns, by Gram–Schmidt with exact arithmetic: `Q` has orthonormal columns (their norms
are square roots), `R` is upper triangular with positive diagonal. Symbols are taken to
be real in the inner product. Returns a `SymJ.SymQR` (`F.Q`, `F.R`, or `Q, R = qr(A)`).

```julia
Q, R = qr(Sym[1 1; 1 0])          # Q = [1 1; 1 -1]/sqrt(2), R = [sqrt(2) 1/sqrt(2); 0 1/sqrt(2)]
```
"""
function LinearAlgebra.qr(A::AbstractMatrix{Sym})
    m, n = size(A)
    n <= m || throw(ArgumentError("qr of a matrix with more columns than rows is not supported"))
    Q = Matrix{Sym}(undef, m, n)
    R = zeros(Sym, n, n)
    for k in 1:n
        v = Sym[Sym(a) for a in A[:, k]]
        for j in 1:(k - 1)
            r = _qrnorm(_dot(Q[:, j], v))
            R[j, k] = r
            iszero(r) || (v = Sym[expand(v[i] - r * Q[i, j]) for i in 1:m])
        end
        n2 = _qrnorm(_dot(v, v))
        iszero(n2) && throw(ArgumentError("qr: the columns are linearly dependent (column $k)"))
        nv = simplify(sqrt(n2))
        R[k, k] = nv
        for i in 1:m
            Q[i, k] = simplify(v[i] / nv)
        end
    end
    SymQR(Q, R)
end

## SVD

# Aᴴ: the transpose, entries with i conjugated (symbols taken to be real, as in _dot)
_adjoint(A::AbstractMatrix) = Sym[(e = Sym(A[i, j]); _has_im(e) ? conj(e) : e) for j in axes(A, 2), i in axes(A, 1)]

# orthonormal columns spanning those of B (Gram–Schmidt, as in qr), dependent ones dropped
function _orthonormal(B::AbstractMatrix{Sym}, others::Vector{Vector{Sym}} = Vector{Sym}[])
    out = Vector{Sym}[]
    for k in axes(B, 2)
        v = Sym[Sym(a) for a in B[:, k]]
        for q in Iterators.flatten((others, out))
            r = _qrnorm(_dot(q, v))
            iszero(r) || (v = Sym[expand(v[i] - r * q[i]) for i in eachindex(v)])
        end
        n2 = _qrnorm(_dot(v, v))
        iszero(n2) && continue
        nv = simplify(sqrt(n2))
        push!(out, Sym[simplify(x / nv) for x in v])
    end
    out
end

# (σ, v) for the nonzero singular values, largest first when they are numbers: σ² an
# eigenvalue of AᴴA, v orthonormal eigenvectors (each eigenspace orthonormalized; the
# eigenspaces of the Hermitian AᴴA are orthogonal to each other)
function _svd_parts(A::AbstractMatrix{Sym})
    B = map(_qrnorm, _adjoint(A) * A)
    n = size(B, 1)
    # roots of a cubic or higher (Cardano's, with complex cube roots for real roots) make
    # null spaces that do not finish: refused
    λ = Sym(Symbol("_svdλ"))
    p = factor(charpoly(B, λ))
    for (f, _) in (p.kind == MUL ? collect(zip(p.args, p.exps)) : p.kind == POW ? [(p.args[1], p.args[2])] : [(p, ONE)])
        has(f, λ) || continue
        d = degree(Poly{Sym}(expand(f), [λ]), λ)
        d <= 2 || throw(ArgumentError("svd: the squared singular values are roots of an irreducible polynomial of degree $d; use svd on numbers (e.g. svd(N.(A)))"))
    end
    parts = Tuple{Sym,Vector{Sym}}[]
    for (λ, _) in _eigen_mult(B)
        λ = _qrnorm(λ)
        iszero(λ) && continue
        N = LinearAlgebra.nullspace(B - λ * LinearAlgebra.I(n))
        σ = _hassymbol(λ) ? refine(sqrt(λ)) : simplify(sqrt(λ))   # (refine: the assumptions in force)
        for v in _orthonormal(N)
            push!(parts, (σ, v))
        end
    end
    if all(p -> !_hassymbol(p[1]), parts)
        sort!(parts; by = p -> -real(_evalc(p[1])))
    end
    parts
end

"""
    svd(A::AbstractMatrix{Sym}; full = false)
    svdvals(A::AbstractMatrix{Sym})

The singular value decomposition of a matrix of expressions, exactly, as a
`LinearAlgebra.SVD` (`F.U`, `F.S`, `F.Vt`, or `U, S, V = svd(A)`), with
`A == F.U * Diagonal(F.S) * F.Vt`. The singular values are the square roots of the
eigenvalues of `A'A` (largest first when they are numbers), the columns of `V`
orthonormal eigenvectors of `A'A`, and those of `U` are `A v / σ`. As for `svd` on
numbers, `U` and `V` have `min(m, n)` columns, completed (for zero singular values) by
orthonormal bases of the null spaces; with `full = true`, `U` and `V` are square.
Symbols are taken to be real in `A'A`. A singular value with a symbol, such as
`sqrt(a^2)`, is refined under the assumptions (`assume`, or an [`assuming`](@ref)
block): `a` for `a > 0`, `abs(a)` for real `a`, and left as it is otherwise. The roots of
the characteristic polynomial of `A'A` must be found, and nested square roots are
common, so this is for small matrices.

```julia
F = svd(Sym[3 0; 4 5])           # F.S == [3*sqrt(5), sqrt(5)]
svdvals(Sym[1 1; 0 0])           # [sqrt(2), 0]
svdvals(Sym[a 0; 0 b])           # [sqrt(b^2), sqrt(a^2)]
svdvals(Sym[a 0; 0 b]; assume = @assume(a > 0, b > 0))   # [b, a]
```
"""
function LinearAlgebra.svd(A::AbstractMatrix{Sym}; full::Bool = false, assume = nothing)
    assume === nothing || return assuming(() -> LinearAlgebra.svd(A; full), _assumptions(assume))
    m, n = size(A)
    if m < n                                     # through A' (the smaller of A'A, AA')
        F = LinearAlgebra.svd(_adjoint(A); full)
        return LinearAlgebra.SVD(Matrix(_adjoint(F.Vt)), F.S, Matrix(_adjoint(F.U)))
    end
    A = Sym[Sym(a) for a in A]
    parts = _svd_parts(A)
    r = length(parts)
    vs = [p[2] for p in parts]
    us = [Sym[simplify(x / σ) for x in A * v] for (σ, v) in parts]
    k = full ? m : n
    # completed by the null spaces of A (for V) and of A' (for U)
    append!(vs, _orthonormal(LinearAlgebra.nullspace(A), vs))
    length(us) < k && append!(us, _orthonormal(_complement_basis(us, m), us))
    S = Sym[[p[1] for p in parts]; zeros(Sym, min(m, n) - r)]
    U = reduce(hcat, us[1:k])
    V = reduce(hcat, vs[1:n])
    LinearAlgebra.SVD(U, S, Matrix(_adjoint(V)))
end

# the columns of the identity and a basis of the orthogonal complement of us, enough to
# complete them (the null space of the matrix with rows us')
_complement_basis(us::Vector{Vector{Sym}}, m::Int) =
    isempty(us) ? Matrix{Sym}(LinearAlgebra.I(m)) : LinearAlgebra.nullspace(reduce(vcat, (permutedims(_adjoint(reshape(u, :, 1))[1, :]) for u in us)))

LinearAlgebra.svdvals(A::AbstractMatrix{Sym}; assume = nothing) = LinearAlgebra.svd(A; assume).S

## Jordan form

"""
    SymJ.jordan_form(A)

`(P, J)` with `A == P * J * inv(P)` and `J` the Jordan form of the square matrix of
expressions `A`: a block for each Jordan chain, the eigenvalue on the diagonal and ones
above it. The chains are built from the kernels of `(A - λI)^k`, exactly. Not exported.

```julia
P, J = SymJ.jordan_form(Sym[1 1; 0 1])     # J = [1 1; 0 1]
```
"""
function jordan_form(A::AbstractMatrix)
    n = _square(A)
    S = Sym[Sym(a) for a in A]
    Id = Matrix{Sym}(LinearAlgebra.I(n))
    blocks = Tuple{Sym,Int}[]
    cols = Vector{Sym}[]
    for (λ, m) in _eigen_mult(S)
        B = S - λ * Id
        # kernels of B^k until their dimension is the multiplicity
        Bk, kers = Id, Matrix{Sym}[zeros(Sym, n, 0)]
        while size(kers[end], 2) < m
            Bk = map(_algnorm, Bk * B)
            push!(kers, LinearAlgebra.nullspace(Bk))
            length(kers) > n + 1 && throw(ArgumentError("jordan_form: could not find the generalized eigenvectors of $λ"))
        end
        p = length(kers) - 1
        chosen = [Vector{Sym}[] for _ in 1:p]        # the chain vectors at each level
        for k in p:-1:1
            nblocks = size(kers[k+1], 2) - size(kers[k], 2)       # chains of length ≥ k
            have = length(chosen[k])
            for c in 1:size(kers[k+1], 2)
                have >= nblocks && break
                v = kers[k+1][:, c]
                span = hcat(kers[k], chosen[k]..., v)
                LinearAlgebra.rank(span) > size(kers[k], 2) + length(chosen[k]) || continue
                # a new chain v, Bv, ..., B^(k-1)v, from level k down to 1
                w = v
                chain = Vector{Sym}[]
                for lev in k:-1:1
                    push!(chosen[lev], w)
                    pushfirst!(chain, w)
                    w = Sym[_algnorm(e) for e in B * w]
                end
                append!(cols, chain)
                push!(blocks, (λ, k))
                have += 1
            end
        end
    end
    length(cols) == n || throw(ArgumentError("jordan_form: not all eigenvalues of the matrix were found"))
    J = zeros(Sym, n, n)
    i = 1
    for (λ, k) in blocks
        for j in 0:(k - 1)
            J[i + j, i + j] = λ
            j < k - 1 && (J[i + j, i + j + 1] = ONE)
        end
        i += k
    end
    reduce(hcat, cols), J
end

## matrix functions

# f(A) through the Jordan form: on a block for λ of size k, f(J) has f^(j)(λ)/j! on its
# j-th superdiagonal
function _matfun(f, A::AbstractMatrix{Sym})
    n = _square(A)
    P, J = jordan_form(A)
    t = Sym(Symbol("_matfun_t"))
    F = zeros(Sym, n, n)
    i = 1
    while i <= n
        λ = J[i, i]
        k = 1
        while i + k <= n && J[i + k - 1, i + k] == ONE
            k += 1
        end
        d = f(t)
        for j in 0:(k - 1)
            v = subs(d, t => λ) / factorial(j)
            for r in i:(i + k - 1 - j)
                F[r, r + j] = v
            end
            d = diff(d, t)
        end
        i += k
    end
    map(simplify, P * F * inv(P))
end

"""
    exp(A::AbstractMatrix{Sym}), sin(A), cos(A)

The matrix exponential (sine, cosine) of a square matrix of expressions, exactly,
through its Jordan form. `exp(t*A)` solves `x' = A x`.

```julia
exp(t * Sym[0 1; -1 0])          # [cos(t) sin(t); -sin(t) cos(t)]
```
"""
Base.exp(A::StridedMatrix{Sym}) = _matfun(exp, A)
Base.exp(A::AbstractMatrix{Sym}) = _matfun(exp, A)
Base.sin(A::StridedMatrix{Sym}) = _matfun(sin, A)
Base.cos(A::StridedMatrix{Sym}) = _matfun(cos, A)

"""
    A^n for a matrix of expressions and a symbolic exponent `n`

`A^n` through the Jordan form, for `n` a symbol (taken to be an integer for a defective
matrix): `Sym[1 1; 1 0]^n` gives the Fibonacci numbers.
"""
Base.:^(A::AbstractMatrix{Sym}, n::Sym) = _matfun(λ -> λ^n, A)
Base.:^(A::StridedMatrix{Sym}, n::Sym) = _matfun(λ -> λ^n, A)
