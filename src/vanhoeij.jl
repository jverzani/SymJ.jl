# van Hoeij's recombination for Zassenhaus: the subsets of the lifted factors that
# give true factors, found by lattice reduction instead of trying subsets (whose
# number grows exponentially with the number of modular factors).
#
# For a monic lifted factor g_i mod M = p^a, the traces are t_ij = ℓ^j P_j(g_i) mod M
# (P_j the j-th power sum of its roots, by Newton's identities; ℓ = lc(f)). For a true
# factor G = ℓ_G Π_{i ∈ S} g_i, Σ_{i∈S} t_ij ≡ ℓ^j P_j(G) (mod M), an integer of size at
# most T_j = n (|ℓ| R)^j (R a root bound of f). In the lattice spanned by the rows
# [e_i | t_i1/D_1 ... t_iN/D_N] and [0 | M/D_j at column j] (D_j ≈ T_j / r, the columns
# scaled down so true factors give short vectors), the vectors (e_S, ...) of the true
# factors have norm at most B (below). After LLL, a basis vector whose Gram–Schmidt
# norm exceeds B cannot be in a vector of norm ≤ B: dropping those leaves a span
# containing the true-factor vectors; once it contains nothing else, its reduced
# echelon form is their 0/1 partition. Otherwise more traces are added.

const _VANHOEIJ_MIN_FACTORS = 9        # fewer modular factors: subsets are cheaper

# LLL reduction (δ = 3/4) with exact rational Gram–Schmidt; returns the basis and the
# squared Gram–Schmidt norms
function _lll(b::Vector{Vector{BigInt}})
    n = length(b)
    Q = Rational{BigInt}
    dot(u, v) = sum(u[i] * v[i] for i in eachindex(u))
    μ = zeros(Q, n, n)
    Bn = zeros(Q, n)
    bs = Vector{Vector{Q}}(undef, n)
    for i in 1:n
        v = Q.(b[i])
        for j in 1:(i - 1)
            μ[i, j] = dot(Q.(b[i]), bs[j]) / Bn[j]
            v .-= μ[i, j] .* bs[j]
        end
        bs[i] = v
        Bn[i] = dot(v, v)
    end
    function reduce!(k, l)
        abs(μ[k, l]) <= 1 // 2 && return
        q = round(BigInt, μ[k, l])
        b[k] .-= q .* b[l]
        μ[k, l] -= q
        for i in 1:(l - 1)
            μ[k, i] -= q * μ[l, i]
        end
    end
    k = 2
    while k <= n
        reduce!(k, k - 1)
        if Bn[k] < (3 // 4 - μ[k, k-1]^2) * Bn[k-1]
            m = μ[k, k-1]
            B = Bn[k] + m^2 * Bn[k-1]
            μ[k, k-1] = m * Bn[k-1] / B
            Bn[k] = Bn[k-1] * Bn[k] / B
            Bn[k-1] = B
            b[k], b[k-1] = b[k-1], b[k]
            for j in 1:(k - 2)
                μ[k, j], μ[k-1, j] = μ[k-1, j], μ[k, j]
            end
            for i in (k + 1):n
                t = μ[i, k]
                μ[i, k] = μ[i, k-1] - m * t
                μ[i, k-1] = t + μ[k, k-1] * μ[i, k]
            end
            k = max(2, k - 1)
        else
            for l in (k - 2):-1:1
                reduce!(k, l)
            end
            k += 1
        end
    end
    b, Bn
end

# power sums P_1..P_N of the roots of a monic g (low to high), mod M
function _power_sums(g::Vector{BigInt}, N::Int, M::BigInt)
    d = length(g) - 1
    c(i) = g[i+1]                                   # coefficient of x^i
    P = zeros(BigInt, N)
    for j in 1:N
        s = big(0)
        for i in 1:min(j - 1, d)
            s += c(d - i) * P[j-i]
        end
        j <= d && (s += j * c(d - j))
        P[j] = mod(-s, M)
    end
    P
end

# reduced row echelon form over ℚ; the nonzero rows
function _rref(A::Matrix{Rational{BigInt}})
    A = copy(A)
    m, n = size(A)
    r = 1
    for c in 1:n
        r > m && break
        piv = findfirst(i -> !iszero(A[i, c]), r:m)
        piv === nothing && continue
        piv += r - 1
        A[r, :], A[piv, :] = A[piv, :], A[r, :]
        A[r, :] ./= A[r, c]
        for i in 1:m
            i != r && !iszero(A[i, c]) && (A[i, :] .-= A[i, c] .* A[r, :])
        end
        r += 1
    end
    A[1:(r-1), :]
end

# The true factors of f (squarefree, primitive) from its monic modular factors (mod p)
# or nothing (the caller falls back to subset search).
function _vanhoeij(f::Vector{BigInt}, modfacs::Vector{Vector{Int}}, p::Int)
    n, r = length(f) - 1, length(modfacs)
    ℓ = f[end]
    # |roots| ≤ R: Fujiwara's bound 2 max (|f_{n-i}|/|ℓ|)^(1/i), much tighter than
    # Cauchy's when the coefficients are large and the roots small
    R = 2 * maximum(i -> _iroot(cld(abs(f[n-i+1]), abs(ℓ)), i) + 1, 1:n)
    T(j) = n * (abs(ℓ) * R)^j
    N = 1
    while N <= 2n
        Ds = [max(big(1), cld(T(j), r)) for j in 1:N]
        # precision: M well beyond D_j 2^(r+N) so that wrong combinations stand out
        need = maximum(Ds) * big(2)^(r + N + 20)
        k, M = 1, big(p)
        while M <= need
            k += 1
            M *= p
        end
        g = _hensel_lift(_zmod(f, M), modfacs, p, k)
        tr = [_power_sums(gi, N, M) for gi in g]
        sym(v) = (w = mod(v, M); w > M >> 1 ? w - M : w)
        rows = Vector{Vector{BigInt}}()
        for i in 1:r
            row = zeros(BigInt, r + N)
            row[i] = 1
            for j in 1:N
                row[r+j] = div(sym(mod(ℓ^j * tr[i][j], M)), Ds[j], RoundNearest)
            end
            push!(rows, row)
        end
        for j in 1:N
            row = zeros(BigInt, r + N)
            row[r+j] = div(M, Ds[j], RoundNearest)
            push!(rows, row)
        end
        b, Bn = _lll(rows)
        # a true factor's vector: 0/1 part ≤ r, each trace entry ≤ T_j/D_j + r + 1
        Bsq = r + sum((T(j) ÷ Ds[j] + r + 1)^2 for j in 1:N)
        s = length(b)
        while s > 0 && Bn[s] > Bsq
            s -= 1
        end
        s == 0 && return nothing
        s == 1 && return [f]                           # only f itself: irreducible
        W = _rref(Rational{BigInt}[b[i][j] for i in 1:s, j in 1:r])
        parts = [findall(!iszero, W[i, :]) for i in 1:size(W, 1)]
        if all(i -> all(v -> v == 1, W[i, parts[i]]), 1:size(W, 1)) &&
           sort!(vcat(parts...)) == collect(1:r)
            # a partition that looks right may still be too coarse a span: the
            # candidates are tested, and if one fails, more traces are added
            out = Vector{BigInt}[]
            rest = f
            for S in parts
                cand = BigInt[ℓ]
                for i in S
                    cand = _zmod(_zmul(cand, g[i]), M)
                end
                h = _zprimitive(_zsym(cand, M))
                q = _zdivexact(rest, h)
                q === nothing && break
                push!(out, h)
                rest = q
            end
            length(out) == length(parts) && length(rest) == 1 && return out
        end
        N *= 2
    end
    nothing
end
