# Real root isolation for polynomials with rational coefficients: SymJ.root_intervals
# and SymJ.real_roots.
#
# The square-free decomposition gives the multiplicities; each square-free factor, with
# integer coefficients, is isolated by Descartes' rule of signs with bisection
# (Vincent–Collins–Akritas), in exact integer arithmetic. The positive roots are mapped
# to (0, 1) by x -> 2ᵏ x (2ᵏ above a root bound), the negative ones through p(-x). For
# (0, 1) the count is the sign variations of (x + 1)ⁿ p(1/(x + 1)); a count of 0 or 1 is
# final, otherwise the interval is halved: the left half's polynomial is 2ⁿ p(x/2), the
# right half's that shifted by 1, so each node is built from its parent with shifts by
# bits and additions. An interval is (k, d) for [k/2ᵈ, (k + 1)/2ᵈ]. Intervals are refined
# by bisection with exact signs (p(a/b) bⁿ, an integer).

"""
    SymJ.root_intervals(p, x; width = nothing)

Isolating intervals for the real roots of the polynomial `p` in `x` (rational
coefficients): a sorted vector of `(lo, hi, m)`, each interval holding exactly one root
of multiplicity `m`, the intervals disjoint, their endpoints exact rationals
(`lo == hi` for a rational root found exactly). With `width`, each interval is refined to
be narrower than it. Not exported.

```julia
@syms x
SymJ.root_intervals(x^3 - 2x - 5, x)            # [(2, 4, 1)]  (one real root)
SymJ.root_intervals((x - 1)^2 * (x^2 - 2), x)   # [(-2, -1, 1), (1, 1, 2), (1, 2, 1)]-like
```
"""
function root_intervals(p, x::Sym; width = nothing)
    facs = _rr_factors(Sym(p), x)
    out = Tuple{Rational{BigInt},Rational{BigInt},Int,Vector{BigInt}}[]
    for (c, m) in facs
        for (lo, hi) in _rr_isolate(c)
            push!(out, (lo, hi, m, c))
        end
    end
    sort!(out; by = t -> t[1])
    # intervals of different factors may overlap (an open interval holds its root
    # strictly inside, so a shared endpoint is no overlap): refine until apart
    i = 1
    while i < length(out)
        a, b = out[i], out[i + 1]
        if b[1] < a[2]
            a[1] == a[2] || (out[i] = (_rr_refine(a[4], a[1], a[2], (a[2] - a[1]) / 2)..., a[3], a[4]))
            b[1] == b[2] || (out[i + 1] = (_rr_refine(b[4], b[1], b[2], (b[2] - b[1]) / 2)..., b[3], b[4]))
            sort!(out; by = t -> t[1])
            i = max(i - 1, 1)
            continue
        end
        i += 1
    end
    w = width === nothing ? nothing : Rational{BigInt}(width isa AbstractFloat ? rationalize(BigInt, width) : width)
    [(w === nothing ? (lo, hi, m) : (_rr_refine(c, lo, hi, w)..., m)) for (lo, hi, m, c) in out]
end

"""
    SymJ.real_roots(p, x; T = Float64)

The real roots of the polynomial `p` in `x` (rational coefficients) as numbers of type
`T` (`Float64`, or `BigFloat` at the current precision), in increasing order, each
repeated by its multiplicity. Each root is the common rounding of the endpoints of an
interval refined until they round alike, so it is correctly rounded. Not exported.

```julia
SymJ.real_roots(x^3 - 2x - 5, x)               # [2.0945514815423265]
SymJ.real_roots((x - 1)^2 * (x^2 - 2), x)      # [-1.4142135623730951, 1.0, 1.0, 1.4142135623730951]
```
"""
function real_roots(p, x::Sym; T::Type{<:AbstractFloat} = Float64)
    out = T[]
    for (c, m) in _rr_factors(Sym(p), x)
        for (lo, hi) in _rr_isolate(c)
            while lo != hi && T(lo) != T(hi)
                lo, hi = _rr_refine(c, lo, hi, (hi - lo) / 2)
            end
            append!(out, fill(T(lo), m))
        end
    end
    sort!(out)
end

# the square-free factors of p, as dense integer coefficient vectors (constant first),
# with their multiplicities
function _rr_factors(p::Sym, x::Sym)
    _isnamed(x) || throw(ArgumentError("root_intervals: $x is not a symbol"))
    P = try Poly{SymNumber}(expand(p), [x]) catch
        throw(ArgumentError("root_intervals: $p is not a polynomial in $x with rational coefficients"))
    end
    all(c -> isexact(c), coefficients(P)) || throw(ArgumentError("root_intervals: the coefficients must be exact rationals"))
    iszero(P) && throw(ArgumentError("root_intervals: the zero polynomial has every number as a root"))
    degree(P, x) == 0 && return Tuple{Vector{BigInt},Int}[]
    _, facs = _factor_squarefree_poly(P)
    [(_rr_dense(f), m) for (f, m) in facs if degree(f, f.gens[1]) > 0]
end

function _rr_dense(f::Poly)
    n = maximum(ev[1] for ev in exponent_vectors(f))
    c = zeros(BigInt, n + 1)
    for (ev, a) in zip(exponent_vectors(f), coefficients(f))
        c[ev[1] + 1] = BigInt(a)
    end
    c
end

# isolating intervals of the real roots of a square-free integer polynomial c (constant
# first): (lo, hi) with one root strictly inside, or lo == hi for an exact root
function _rr_isolate(c::Vector{BigInt})
    out = Tuple{Rational{BigInt},Rational{BigInt}}[]
    if iszero(c[1])                           # a root at 0 (once: square-free)
        push!(out, (big(0) // 1, big(0) // 1))
        c = c[2:end]
    end
    n = length(c) - 1
    n == 0 && return out
    # 2ᵏ above the Cauchy bound 1 + max|cᵢ/cₙ|
    B = 1 + cld(maximum(abs, @view c[1:n]), abs(c[end]))
    k = ndigits(B; base = 2)
    buf = Vector{BigInt}(undef, n + 1)
    for s in (1, -1)
        q = BigInt[(s == -1 && isodd(i) ? -c[i + 1] : c[i + 1]) << (k * i) for i in 0:n]   # p(±2ᵏ x)
        found = Tuple{BigInt,Int,Bool}[]          # (k, d, exact)
        _rr_bisect!(found, q, big(0), 0, buf)
        for (j, d, exact) in found
            lo = (j // big(2)^d) * big(2)^k
            hi = exact ? lo : ((j + 1) // big(2)^d) * big(2)^k
            push!(out, s == 1 ? (lo, hi) : (-hi, -lo))
        end
    end
    sort!(out)
end

signvar(c) = (v = 0; s = 0; for a in c; iszero(a) && continue; t = sign(a); s != 0 && t != s && (v += 1); s = t; end; v)

# c(x) -> c(x + 1), in place, additions only
function _rr_shift1!(c::Vector{BigInt})
    n = length(c)
    @inbounds for i in 1:(n - 1), j in (n - 1):-1:i
        Base.GMP.MPZ.add!(c[j], c[j + 1])
    end
    c
end

# the Descartes count for (0, 1): sign variations of (x + 1)ⁿ q(1/(x + 1))
function _rr_count(q::Vector{BigInt}, buf::Vector{BigInt})
    n = length(q)
    for i in 1:n
        buf[i] = Base.GMP.MPZ.set(q[n + 1 - i])
    end
    signvar(_rr_shift1!(buf))
end

function _rr_bisect!(found, q::Vector{BigInt}, j::BigInt, d::Int, buf::Vector{BigInt})
    length(buf) == length(q) || resize!(buf, length(q))
    v = _rr_count(q, buf)
    v == 0 && return
    v == 1 && (push!(found, (j, d, false)); return)
    n = length(q) - 1
    L = BigInt[q[i + 1] << (n - i) for i in 0:n]            # 2ⁿ q(x/2): [0, 1/2]
    R = _rr_shift1!(BigInt[Base.GMP.MPZ.set(a) for a in L])  # its shift by 1: [1/2, 1]
    if iszero(R[1])                                          # a root at the midpoint
        push!(found, (2j + 1, d + 1, true))
        R = R[2:end]
    end
    _rr_bisect!(found, L, 2j, d + 1, buf)
    _rr_bisect!(found, R, 2j + 1, d + 1, buf)
end

# the sign of c at a rational point: c(a/b) bⁿ, an integer
function _rr_sign(c::Vector{BigInt}, r::Rational{BigInt})
    a, b = numerator(r), denominator(r)
    n = length(c) - 1
    acc = c[end]
    bp = big(1)
    for i in (n - 1):-1:0
        bp *= b
        acc = acc * a + c[i + 1] * bp
    end
    sign(acc)
end

# c(x)/(b x - a) for a root a/b of c, integer coefficients (exact: Gauss's lemma)
function _rr_deflate(c::Vector{BigInt}, r::Rational{BigInt})
    a, b = numerator(r), denominator(r)
    n = length(c) - 1
    q = Vector{BigInt}(undef, n)
    rem = c[end]
    for i in n:-1:1                       # c = (b x - a) q: from the top coefficient down
        q[i] = div(rem, b)
        rem = c[i] + a * q[i]
    end
    q
end

# (lo, hi) bisected until narrower than w (one root inside; an exact root ends it). An
# endpoint that is itself a root of c (an exact root found next to the interval) is
# divided out first, so the signs at the endpoints decide the side.
function _rr_refine(c::Vector{BigInt}, lo::Rational{BigInt}, hi::Rational{BigInt}, w)
    lo == hi && return (lo, hi)
    iszero(_rr_sign(c, lo)) && (c = _rr_deflate(c, lo))
    iszero(_rr_sign(c, hi)) && (c = _rr_deflate(c, hi))
    slo = _rr_sign(c, lo)
    while hi - lo >= w
        m = (lo + hi) / 2
        sm = _rr_sign(c, m)
        sm == 0 && return (m, m)
        sm == slo ? (lo = m) : (hi = m)
    end
    (lo, hi)
end
