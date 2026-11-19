# Smith and Hermite normal forms of matrices over a Euclidean domain: the integers, or the
# polynomials in one variable with rational coefficients (ℚ[x]). Row and column operations
# by unimodular matrices, the gcd of two entries through the extended gcd.

# ℚ[x], coefficients from the constant term up, no trailing zeros (zero is empty)
struct _QP
    c::Vector{Rational{BigInt}}
    _QP(c::Vector{Rational{BigInt}}) = (while !isempty(c) && iszero(c[end]); pop!(c); end; new(c))
end
_QP(r::Rational{BigInt}) = _QP([r])
Base.zero(::Type{_QP}) = _QP(Rational{BigInt}[])
Base.one(::Type{_QP}) = _QP([one(Rational{BigInt})])
Base.iszero(p::_QP) = isempty(p.c)
Base.:(==)(p::_QP, q::_QP) = p.c == q.c
_deg(p::_QP) = length(p.c) - 1
function Base.:+(p::_QP, q::_QP)
    n = max(length(p.c), length(q.c))
    z = zero(Rational{BigInt})
    _QP(Rational{BigInt}[(i <= length(p.c) ? p.c[i] : z) + (i <= length(q.c) ? q.c[i] : z) for i in 1:n])
end
Base.:-(p::_QP) = _QP(-p.c)
Base.:-(p::_QP, q::_QP) = p + (-q)
function Base.:*(p::_QP, q::_QP)
    (iszero(p) || iszero(q)) && return zero(_QP)
    c = zeros(Rational{BigInt}, length(p.c) + length(q.c) - 1)
    for i in eachindex(p.c), j in eachindex(q.c)
        c[i + j - 1] += p.c[i] * q.c[j]
    end
    _QP(c)
end
function Base.divrem(p::_QP, q::_QP)
    iszero(q) && throw(DivideError())
    r = copy(p.c)
    dq = _deg(q)
    quo = zeros(Rational{BigInt}, max(0, length(r) - dq))
    while length(r) - 1 >= dq && !isempty(r)
        k = length(r) - 1 - dq
        f = r[end] / q.c[end]
        quo[k + 1] = f
        for j in eachindex(q.c)
            r[k + j] -= f * q.c[j]
        end
        while !isempty(r) && iszero(r[end]); pop!(r); end
    end
    _QP(quo), _QP(r)
end

# the domain's operations on BigInt and _QP: the size that the pivot minimizes, the unit
# that normalizes (a positive integer, a monic polynomial), the extended gcd
_size(a::BigInt) = abs(a)
_size(p::_QP) = _deg(p)
_unit(a::BigInt) = a < 0 ? -one(BigInt) : one(BigInt)
_unit(p::_QP) = iszero(p) ? one(_QP) : _QP(inv(p.c[end]))      # the inverse of the leading coefficient
function _gcdx(a::BigInt, b::BigInt)
    g, s, t = gcdx(a, b)
    g, s, t
end
function _gcdx(a::_QP, b::_QP)
    r0, r1, s0, s1, t0, t1 = a, b, one(_QP), zero(_QP), zero(_QP), one(_QP)
    while !iszero(r1)
        q, r = divrem(r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, s0 - q * s1
        t0, t1 = t1, t0 - q * t1
    end
    u = _unit(r0)
    r0 * u, s0 * u, t0 * u
end
_divides(a, b) = iszero(divrem(b, a)[2])        # a | b
# the reduced remainder for Hermite entries above a pivot: 0 ≤ r < p for integers, deg r < deg p
_hnf_rem(a::BigInt, p::BigInt) = mod(a, p)
_hnf_rem(a::_QP, p::_QP) = divrem(a, p)[2]
_hnf_quo(a::BigInt, p::BigInt) = fld(a, p)
_hnf_quo(a::_QP, p::_QP) = divrem(a, p)[1]

# the matrix's entries as integers, or as polynomials in its one variable over ℚ: (M, x),
# x nothing for integers
function _nf_matrix(A::AbstractMatrix)
    S = Sym.(A)
    syms = Set{Sym}()
    for e in S
        union!(syms, free_symbols(e))
    end
    length(syms) <= 1 || throw(ArgumentError("normal form: the entries are integers or polynomials in one variable, not in $(join(sort!(collect(syms); lt = _lt), ", "))"))
    if isempty(syms) && all(e -> e.kind == NUMBER && isexact(e.num) && _isint(e.num), S)
        return map(e -> BigInt(_rat(e.num)), S), nothing
    end
    x = isempty(syms) ? Sym(:x) : only(syms)
    M = map(S) do e
        P = try Poly{SymNumber}(expand(e), [x]) catch err
            err isa ArgumentError || rethrow()
            throw(ArgumentError("normal form: $e is not a polynomial in $x"))
        end
        d = degree(P, x)
        c = zeros(Rational{BigInt}, max(d, -1) + 1)
        for (ev, cf) in zip(exponent_vectors(P), coefficients(P))
            (isexact(cf)) || throw(ArgumentError("normal form: $e has a float coefficient"))
            c[ev[1] + 1] = _rat(cf)
        end
        _QP(c)
    end
    M, x
end
_nf_sym(a::BigInt, x) = Sym(a)
_nf_sym(p::_QP, x) = _add_many(Sym(SymNumber(c)) * x^(i - 1) for (i, c) in enumerate(p.c))
_nf_out(M, x) = map(a -> _nf_sym(a, x), M)

# row and column operations on A, mirrored on U (rows) and V (columns)
function _rowcombine!(A, U, i, j, a, b, c, d)          # [row i; row j] ← [a b; c d] [row i; row j]
    for M in (A, U)
        for k in axes(M, 2)
            ri, rj = M[i, k], M[j, k]
            M[i, k], M[j, k] = a * ri + b * rj, c * ri + d * rj
        end
    end
end
function _rowscale!(A, U, i, u)                       # row i times the unit u
    for M in (A, U), k in axes(M, 2)
        M[i, k] = u * M[i, k]
    end
end
function _colcombine!(A, V, i, j, a, b, c, d)          # [col i, col j] ← [col i, col j] [a c; b d]
    for M in (A, V)
        for k in axes(M, 1)
            ci, cj = M[k, i], M[k, j]
            M[k, i], M[k, j] = a * ci + b * cj, c * ci + d * cj
        end
    end
end

function _smith!(A::Matrix{T}) where {T}
    m, n = size(A)
    U = T[i == j ? one(T) : zero(T) for i in 1:m, j in 1:m]
    V = T[i == j ? one(T) : zero(T) for i in 1:n, j in 1:n]
    for t in 1:min(m, n)
        while true
            # the smallest nonzero entry of the lower right block, moved to (t, t)
            best = nothing
            for i in t:m, j in t:n
                iszero(A[i, j]) && continue
                (best === nothing || _size(A[i, j]) < _size(A[best...])) && (best = (i, j))
            end
            best === nothing && return A, U, V
            i, j = best
            i != t && _rowcombine!(A, U, t, i, zero(T), one(T), one(T), zero(T))
            j != t && _colcombine!(A, V, t, j, zero(T), one(T), one(T), zero(T))
            # the column below the pivot, then the row right of it, cleared through gcds
            for i in (t + 1):m
                iszero(A[i, t]) && continue
                a, b = A[t, t], A[i, t]
                if _divides(a, b)
                    q = divrem(b, a)[1]
                    _rowcombine!(A, U, t, i, one(T), zero(T), -q, one(T))
                else
                    g, s, r = _gcdx(a, b)
                    _rowcombine!(A, U, t, i, s, r, -divrem(b, g)[1], divrem(a, g)[1])
                end
            end
            for j in (t + 1):n
                iszero(A[t, j]) && continue
                a, b = A[t, t], A[t, j]
                if _divides(a, b)
                    q = divrem(b, a)[1]
                    _colcombine!(A, V, t, j, one(T), zero(T), -q, one(T))
                else
                    g, s, r = _gcdx(a, b)
                    _colcombine!(A, V, t, j, s, r, -divrem(b, g)[1], divrem(a, g)[1])
                end
            end
            # (a column step can refill the column below the pivot)
            any(i -> !iszero(A[i, t]), (t + 1):m) && continue
            # the pivot divides the rest of the block, else a row of it is added and again
            k = findfirst(ij -> !_divides(A[t, t], A[ij[1], ij[2]]), [(i, j) for i in (t + 1):m for j in (t + 1):n])
            if k !== nothing
                i = [(i, j) for i in (t + 1):m for j in (t + 1):n][k][1]
                _rowcombine!(A, U, t, i, one(T), one(T), zero(T), one(T))
                continue
            end
            _rowscale!(A, U, t, _unit(A[t, t]))
            break
        end
    end
    A, U, V
end

function _hermite!(A::Matrix{T}) where {T}
    m, n = size(A)
    U = T[i == j ? one(T) : zero(T) for i in 1:m, j in 1:m]
    r = 1                                              # the next pivot row
    for j in 1:n
        r > m && break
        # the gcd of column j's entries from row r down into row r
        for i in (r + 1):m
            iszero(A[i, j]) && continue
            if iszero(A[r, j])
                _rowcombine!(A, U, r, i, zero(T), one(T), one(T), zero(T))
                continue
            end
            a, b = A[r, j], A[i, j]
            g, s, t = _gcdx(a, b)
            _rowcombine!(A, U, r, i, s, t, -divrem(b, g)[1], divrem(a, g)[1])
        end
        iszero(A[r, j]) && continue                     # no pivot in this column
        _rowscale!(A, U, r, _unit(A[r, j]))
        # the entries above the pivot reduced by it
        for i in 1:(r - 1)
            q = _hnf_quo(A[i, j], A[r, j])
            iszero(q) || _rowcombine!(A, U, i, r, one(T), -q, zero(T), one(T))
        end
        r += 1
    end
    A, U
end

"""
    SymJ.smith_normal_form(A)
    SymJ.smith_normal_decomp(A)

The Smith normal form of a matrix whose entries are integers, or polynomials in one variable
with rational coefficients: the diagonal matrix `D` (of `A`'s size) whose diagonal entries
`d₁, d₂, …` each divide the next, nonnegative integers or monic polynomials, with
`U A V = D` for invertible `U`, `V` over the integers (or the polynomials). `smith_normal_decomp`
gives `(D, U, V)`. The `dᵢ` are `A`'s invariant factors: `d₁ ⋯ dₖ` is the gcd of the `k × k`
minors. As SymPy's `smith_normal_form` and `smith_normal_decomp`. Not exported.

```julia
SymJ.smith_normal_form(Sym[2 4 4; -6 6 12; 10 -4 -16])   # diagonal 2, 6, 12
D, U, V = SymJ.smith_normal_decomp(Sym[x 1; 0 x])         # diagonal 1, x^2
```
"""
smith_normal_form(A::AbstractMatrix) = smith_normal_decomp(A)[1]
function smith_normal_decomp(A::AbstractMatrix)
    M, x = _nf_matrix(A)
    D, U, V = _smith!(copy(M))
    _nf_out(D, x), _nf_out(U, x), _nf_out(V, x)
end

"""
    SymJ.hermite_normal_form(A)
    SymJ.hermite_normal_decomp(A)

The Hermite normal form of a matrix whose entries are integers, or polynomials in one variable
with rational coefficients: the upper triangular (row echelon) `H = U A`, `U` invertible over
the integers (or the polynomials), each pivot positive (monic), the entries above a pivot
reduced by it (in `[0, p)` for integers, of lower degree for polynomials), and zero rows last.
`hermite_normal_decomp` gives `(H, U)`. This is the row form, as Mathematica's
`HermiteDecomposition`; SymPy's `hermite_normal_form` is the column form, `transpose` of this
form of `transpose(A)` (with its zero columns dropped). Not exported.

```julia
SymJ.hermite_normal_form(Sym[2 3 6 2; 5 6 1 6; 8 3 1 1])
```
"""
hermite_normal_form(A::AbstractMatrix) = hermite_normal_decomp(A)[1]
function hermite_normal_decomp(A::AbstractMatrix)
    M, x = _nf_matrix(A)
    H, U = _hermite!(copy(M))
    _nf_out(H, x), _nf_out(U, x)
end
