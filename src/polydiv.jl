# Multivariate division with respect to lex order.
#
# The dividend is held in a `_MonoTable` with a max-heap of its monomials. Each step
# pops the largest monomial m. If LM(q) divides m, the term t = m / LM(q) goes to the
# quotient and -t*(q - LT(q)) is added to the table; every monomial it touches is smaller
# than m, so a popped monomial is never revisited. Otherwise the term goes to the
# remainder. Both come out in decreasing order, so no final sort is needed.

# max-heap of packed monomials
function _heappush!(h::Vector{M}, x::M) where {M}
    push!(h, x)
    i = length(h)
    @inbounds while i > 1
        j = i >> 1
        h[j] >= x && break
        h[i] = h[j]
        i = j
    end
    @inbounds h[i] = x
    h
end

function _heappop!(h::Vector{M}) where {M}
    @inbounds top = h[1]
    x = pop!(h)
    n = length(h)
    n == 0 && return top
    i = 1
    @inbounds while true
        l = 2i
        l > n && break
        r = l + 1
        c = (r <= n && h[r] > h[l]) ? r : l
        h[c] <= x && break
        h[i] = h[c]
        i = c
    end
    @inbounds h[i] = x
    top
end

# does monomial a divide monomial m?
@inline function _mono_divides(a::M, m::M, L::_Layout, n::Int) where {M}
    for i in 1:n
        s, mask = _shift(L, n, i), _fieldmask(L, i)
        (a >> s) & mask > (m >> s) & mask && return false
    end
    true
end

"""
    divrem(p::Poly, q::Poly)

Quotient and remainder of multivariate division of `p` by `q` in lex order (generators
most significant first): `p == quo*q + rem`, and no term of `rem` is divisible by the
leading monomial of `q`. For univariate polynomials this is ordinary long division.
Coefficients are divided exactly, so integer polynomials may give rational results.
"""
function Base.divrem(p::Poly{C}, q::Poly{C}) where {C}
    p, q = _unify(p, q)
    iszero(q) && throw(DivideError())
    _divrem_fit(p, q)
end

# The remainder's degrees can exceed the inputs' (x^k by x + y^5): a layout with room
# for twice the exponent that did not fit, in the same word, and again; if the word
# cannot hold it, `_PackOverflow` for the caller
function _divrem_fit(p::Poly{C,M}, q::Poly{C,M}) where {C,M}
    while true
        try
            return _divrem(p, q)
        catch err
            (err isa _DivOverflow) || rethrow()
            degs = max.(_maxdegs(p), _maxdegs(q))
            degs[err.i] = max(degs[err.i], 2 * err.e)
            L = _layout_in(M, degs)
            L === nothing && throw(_PackOverflow(err.e))
            p, q = _repack(p, p.gens, L, M), _repack(q, q.gens, L, M)
        end
    end
end

# an exponent of generator i overflowing during division
struct _DivOverflow <: Exception
    e::Int
    i::Int
end

function _divrem(p::Poly{C,M}, q::Poly{C,M}) where {C,M}
    gs, L, n = p.gens, p.layout, length(p.gens)
    lm, lc = q.monos[1], q.coefs[1]
    qdeg = _maxdegs(q)
    qm, qc = q.monos, q.coefs
    work = _MonoTable{C}(2 * (length(p.monos) + length(qm)), M)
    heap = M[]
    for (m, c) in zip(p.monos, p.coefs)
        _addto!(work, m, c)
        _heappush!(heap, m)
    end
    quom, quoc, remm, remc = M[], C[], M[], C[]
    while !isempty(heap)
        m = _heappop!(heap)
        c = _lookup(work, m)
        iszero(c) && continue
        if !_mono_divides(lm, m, L, n)
            push!(remm, m)
            push!(remc, c)
            continue
        end
        t, ct = m - lm, c / lc
        for i in 1:n                      # t * (terms of q) must not overflow a field
            e = _exponent(t, L, n, i) + qdeg[i]
            e > _maxexp(L, i) && throw(_DivOverflow(e, i))
        end
        push!(quom, t)
        push!(quoc, ct)
        @inbounds for j in 2:length(qm)
            mj = t + qm[j]
            _addto!(work, mj, -(ct * qc[j])) && _heappush!(heap, mj)
        end
    end
    Poly{C}(gs, L, quom, quoc), Poly{C}(gs, L, remm, remc)
end

Base.div(p::Poly, q::Poly) = divrem(p, q)[1]
Base.rem(p::Poly, q::Poly) = divrem(p, q)[2]

"""
    divides(p::Poly, q::Poly) -> (Bool, quotient)

Whether `q` divides `p` exactly, and the quotient when it does.
"""
function divides(p::Poly, q::Poly)
    d, r = divrem(p, q)
    iszero(r), d
end

"""
    divexact(p::Poly, q::Poly)

The quotient `p / q`, which must be exact (throws an `ArgumentError` otherwise).
"""
function divexact(p::Poly, q::Poly)
    ok, d = divides(p, q)
    ok || throw(ArgumentError("division is not exact"))
    d
end
