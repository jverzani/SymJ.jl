# Recognizing floats as exact numbers (SymPy's nsimplify): a float is matched by a
# rational, a rational combination of given constants, or an algebraic number of low
# degree written in radicals. The relations are integer relations, found by PSLQ, among
# x, 1 and the constants, or among the powers 1, x, …, x^d. For any x, among the
# (2M + 1)^n integer vectors b of length n with entries at most M, some b·x is as small
# as about 1/(2M + 1)^n relative to the terms; a relation is accepted only when it holds
# to within tol and (2M + 1)^n tol is small, so it is unlikely to hold by chance, and
# random floats are left alone.

"""
    SymJ.nsimplify(ex; constants = (), tol = 1e-14, maxcoeff = 1000, maxdegree = 4)

`ex` with its floats replaced by exact numbers that match them, as SymPy's `nsimplify`.
Each float is recognized, in this order, as

- a rational with small numerator and denominator (`0.75 = 3/4`);
- a rational combination of the `constants` (exact numbers such as `pi`, `ℯ`, `sqrt(Sym(2))`,
  `log(Sym(2))`) with small coefficients (`7.283185307179586 = 1 + 2π` with `constants = [pi]`);
- a root of an integer polynomial of degree at most `maxdegree` with small coefficients,
  given in radicals (or as `solve` gives it) (`1.4142135623730951 = √2`,
  `1.618033988749895 = (1 + √5)/2`).

`tol` is the relative error allowed; a float typed to fewer digits needs a larger
`tol` (`SymJ.nsimplify(1.414214; tol = 1e-6)` is `√2`). A relation (with coefficients at
most `maxcoeff`) is accepted only when its coefficients are small enough that it is
unlikely to hold to within `tol` by chance, so a float that is not recognized is left as
it is. To write every float as
a rational, use [`rationalize`](@ref). Not exported.

```julia
SymJ.nsimplify(0.5x + 1.4142135623730951)              # x/2 + √2
SymJ.nsimplify(2.0943951023931953; constants = [pi])   # 2π/3
SymJ.nsimplify(1.2599210498948732)                     # 2^(1/3)
SymJ.nsimplify(3.0 + 2.718281828459045im; constants = [ℯ])   # 3 + ℯ*im
```
"""
function nsimplify(ex; constants = (), tol::Real = 1e-14, maxcoeff::Integer = 1000, maxdegree::Integer = 4)
    cs = Tuple{Sym,Float64}[]
    for c in constants
        c = Sym(c)
        v = _hassymbol(c) ? nothing : N(c)
        v isa Real || throw(ArgumentError("nsimplify: the constant $c has no real value"))
        push!(cs, (c, Float64(v)))
    end
    rtol = Float64(tol)
    0 < rtol < 1 || throw(ArgumentError("nsimplify: tol is between 0 and 1, not $tol"))
    _mapfloats(f -> something(_identify(f, cs, rtol, maxcoeff, maxdegree), Sym(f)), Sym(ex))
end

# a relation b unlikely to hold to within rtol by chance: few vectors are its size
_unlikely(b, rtol) = (2.0 * maximum(abs, b) + 1)^length(b) * rtol <= 1e-3

function _identify(f::Float64, cs, rtol, maxcoeff, maxdegree)
    iszero(f) && return ZERO
    q = rationalize(BigInt, f; tol = rtol * abs(f))
    _unlikely((numerator(q), denominator(q)), rtol) && return Sym(SymNumber(q))
    setprecision(BigFloat, 128) do
        x = BigFloat(f)
        if !isempty(cs)
            # b₁x + b₂ + Σ bᵢ cᵢ = 0
            v = [x; one(BigFloat); BigFloat[c[2] for c in cs]]
            b = _relation(v, rtol, maxcoeff)
            if b !== nothing && !iszero(b[1])
                r = -(b[2] + sum(b[i + 2] * c[1] for (i, c) in enumerate(cs))) / b[1]
                return r
            end
        end
        for d in 2:maxdegree
            b = _relation(BigFloat[x^k for k in 0:d], rtol, maxcoeff)
            (b === nothing || iszero(b[end])) && continue
            r = _algebraic_root(b, f, rtol)
            r === nothing || return r
        end
        nothing
    end
end

# an integer relation of v with small coefficients, checked against v
function _relation(v::Vector{BigFloat}, rtol, maxcoeff)
    b = _pslq(v, BigFloat(rtol) * maxcoeff, maxcoeff)
    b === nothing && return nothing
    (maximum(abs, b) <= maxcoeff && _unlikely(b, rtol)) || return nothing
    abs(sum(b .* v)) <= rtol * sum(abs.(b) .* abs.(v)) || return nothing
    b
end

# the root of Σ bₖ zᵏ⁻¹ nearest to f, exactly
function _algebraic_root(b::Vector{BigInt}, f::Float64, rtol)
    d = length(b) - 1
    if d == 2
        c, bb, a = b
        disc = bb^2 - 4a * c
        disc < 0 && return nothing
        rs = (Sym(-bb) + sqrt(Sym(disc))) / Sym(2a), (Sym(-bb) - sqrt(Sym(disc))) / Sym(2a)
    else
        z = Sym(:_nsimplify_z)
        p = sum(Sym(b[k + 1]) * z^k for k in 0:d)
        rs = try solve(p ~ 0, z) catch; return nothing end
        rs isa AbstractVector || return nothing
    end
    best, err = nothing, Inf
    for r in rs
        v = try N(r) catch; continue end
        v isa Real || continue
        e = abs(v - f)
        e < err && ((best, err) = (r, e))
    end
    err <= 10rtol * abs(f) || return nothing
    d == 2 ? best : simplify(best)              # solve's radicals denested: √(5 + 2√6) = √2 + √3
end

# PSLQ (Ferguson, Bailey and Arno): an integer vector b, not zero, with b·x below tol
# (x scaled to norm 1), or nothing when every relation has norm above maxcoeff
function _pslq(x::Vector{BigFloat}, tol::BigFloat, maxcoeff::Integer; maxiter::Integer = 100 * length(x)^2)
    n = length(x)
    s = [sqrt(sum(x[j]^2 for j in k:n)) for k in 1:n]
    y = x ./ s[1]
    s ./= s[1]
    for i in 1:n                                    # a zero entry is a relation itself
        abs(y[i]) < tol && return BigInt[j == i for j in 1:n]
    end
    γ = sqrt(BigFloat(4) / 3)
    H = zeros(BigFloat, n, n - 1)
    for j in 1:(n - 1), i in j:n
        H[i, j] = i == j ? s[i + 1] / s[i] : -y[i] * y[j] / (s[j] * s[j + 1])
    end
    A = BigInt[i == j for i in 1:n, j in 1:n]
    B = copy(A)
    function reduce!(i, j)
        iszero(H[j, j]) && return
        t = round(BigInt, H[i, j] / H[j, j])
        iszero(t) && return
        y[j] += t * y[i]
        for k in 1:j
            H[i, k] -= t * H[j, k]
        end
        for k in 1:n
            A[i, k] -= t * A[j, k]
            B[k, j] += t * B[k, i]
        end
    end
    for i in 2:n, j in (i - 1):-1:1
        reduce!(i, j)
    end
    for _ in 1:maxiter
        m = argmax(i -> γ^i * abs(H[i, i]), 1:(n - 1))
        y[m], y[m + 1] = y[m + 1], y[m]
        for k in 1:n
            A[m, k], A[m + 1, k] = A[m + 1, k], A[m, k]
            B[k, m], B[k, m + 1] = B[k, m + 1], B[k, m]
        end
        for k in 1:(n - 1)
            H[m, k], H[m + 1, k] = H[m + 1, k], H[m, k]
        end
        if m <= n - 2
            t0 = hypot(H[m, m], H[m, m + 1])
            t1, t2 = H[m, m] / t0, H[m, m + 1] / t0
            for i in m:n
                t3, t4 = H[i, m], H[i, m + 1]
                H[i, m], H[i, m + 1] = t1 * t3 + t2 * t4, -t2 * t3 + t1 * t4
            end
        end
        for i in (m + 1):n, j in min(i - 1, m + 1):-1:1
            reduce!(i, j)
        end
        j = argmin(abs.(y))
        abs(y[j]) < tol && return B[:, j]
        M = maximum(abs(H[i, i]) for i in 1:(n - 1))
        (iszero(M) || inv(M) > maxcoeff) && return nothing
    end
    nothing
end
