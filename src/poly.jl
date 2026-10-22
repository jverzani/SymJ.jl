# Sparse multivariate polynomials with packed exponents.
#
# A monomial's exponent vector is packed into one UInt64: `bits = 64 ÷ ngens` bits per
# generator, the first generator in the most significant field. Then
#   * monomial multiplication is integer addition (no carries, guarded by degree checks),
#   * lex order (gens[1] > gens[2] > ...) is integer order,
# so terms are kept as parallel vectors sorted by decreasing packed monomial and sums are
# linear merges. Products accumulate in a hash table keyed by the packed monomial.

"""
    Poly{C}

Sparse multivariate polynomial in the generators `gens(p)` with coefficients of type `C`:
`SymNumber` for numeric coefficients, or `Sym` for symbolic ones. Terms are stored in
lex order with `gens(p)[1]` most significant.

    Poly(ex, x, y, ...)

Convert the expression `ex` to a polynomial in `x, y, ...`. Coefficients are `SymNumber`
when they are all numeric, and `Sym` otherwise (so `Poly(a*x^2 + b, x)` has coefficients
`a` and `b`). Use `Poly{C}(ex, x, ...)` to choose the coefficient type. `Sym(p)` converts
back. Each generator may have degree up to `2^(64 ÷ ngens) - 2`; exceeding it throws an
`OverflowError`.
"""
struct Poly{C} <: AbstractSymbolic
    gens::Vector{Sym}
    bits::Int
    monos::Vector{UInt64}   # strictly decreasing
    coefs::Vector{C}        # nonzero
end

_polybits(n::Int) = n == 0 ? 64 : (n <= 64 ? 64 ÷ n : throw(ArgumentError("at most 64 generators are supported")))
_mask(bits::Int) = bits == 64 ? typemax(UInt64) : (UInt64(1) << bits) - 1
# largest allowed exponent; an all-ones field is never used, so the packed monomial
# typemax(UInt64) is free to serve as the empty-slot marker in `_MonoTable`
_maxexp(bits::Int) = Int(min(_mask(bits) - 1, UInt64(typemax(Int))))
@inline _shift(bits::Int, n::Int, i::Int) = bits * (n - i)
@inline _exponent(m::UInt64, bits::Int, n::Int, i::Int) = Int((m >> _shift(bits, n, i)) & _mask(bits))

_zeropoly(::Type{C}, gens, bits) where {C} = Poly{C}(gens, bits, UInt64[], C[])
_constpoly(::Type{C}, gens, bits, c) where {C} =
    iszero(c) ? _zeropoly(C, gens, bits) : Poly{C}(gens, bits, UInt64[0], C[c])
_genpoly(::Type{C}, gens, bits, i::Int) where {C} =
    Poly{C}(gens, bits, UInt64[UInt64(1) << _shift(bits, length(gens), i)], C[one(C)])

## accessors

"""
    gens(p::Poly)

The generators (variables) of `p`, most significant first.
"""
gens(p::Poly) = p.gens

"""
    nterms(p::Poly)

Number of nonzero terms.
"""
nterms(p::Poly) = length(p.monos)

"""
    coefficients(p::Poly)

Coefficients of `p` in lex order (aligned with [`exponent_vectors`](@ref)).
"""
coefficients(p::Poly) = p.coefs

"""
    exponent_vectors(p::Poly)

Exponent vector of each term, in lex order.
"""
exponent_vectors(p::Poly) =
    [[_exponent(m, p.bits, length(p.gens), i) for i in eachindex(p.gens)] for m in p.monos]

function _genindex(p::Poly, x::Sym)
    i = findfirst(g -> _eq(g, x), p.gens)
    i === nothing && throw(ArgumentError("$x is not a generator of the polynomial"))
    i
end

"""
    degree(p::Poly)
    degree(p::Poly, x)

Total degree, or the degree in the generator `x`. The zero polynomial has degree `-1`.
"""
function degree(p::Poly)
    n = length(p.gens)
    maximum((sum(_exponent(m, p.bits, n, i) for i in 1:n; init = 0) for m in p.monos); init = -1)
end
function degree(p::Poly, x::Sym)
    i, n = _genindex(p, x), length(p.gens)
    maximum((_exponent(m, p.bits, n, i) for m in p.monos); init = -1)
end

function _maxdegs(p::Poly)
    n = length(p.gens)
    d = zeros(Int, n)
    for m in p.monos, i in 1:n
        d[i] = max(d[i], _exponent(m, p.bits, n, i))
    end
    d
end

# would multiplying monomials of p and q overflow an exponent field?
function _check_degrees(p::Poly, q::Poly)
    emax = _maxexp(p.bits)
    for (a, b) in zip(_maxdegs(p), _maxdegs(q))
        a + b > emax && throw(OverflowError("exponent exceeds $emax; too many generators for this degree"))
    end
end

Base.iszero(p::Poly) = isempty(p.monos)
Base.isone(p::Poly) = length(p.monos) == 1 && p.monos[1] == 0 && isone(p.coefs[1])
Base.zero(p::Poly{C}) where {C} = _zeropoly(C, p.gens, p.bits)
Base.one(p::Poly{C}) where {C} = _constpoly(C, p.gens, p.bits, one(C))

Base.:(==)(p::Poly, q::Poly) = _samegens(p, q) && p.monos == q.monos && p.coefs == q.coefs
Base.hash(p::Poly, h::UInt) = hash(p.coefs, hash(p.monos, hash(p.gens, h)))

## generators: arithmetic requires a common ring; differing gens are unified

_samegens(p::Poly, q::Poly) = p.gens === q.gens || (length(p.gens) == length(q.gens) && all(_eq.(p.gens, q.gens)))

function _unify(p::Poly{C}, q::Poly{C}) where {C}
    _samegens(p, q) && return (p, q)
    gs = copy(p.gens)
    for g in q.gens
        any(h -> _eq(g, h), gs) || push!(gs, g)
    end
    bits = _polybits(length(gs))
    (_repack(p, gs, bits), _repack(q, gs, bits))
end

function _repack(p::Poly{C}, gs::Vector{Sym}, bits::Int) where {C}
    n, np = length(gs), length(p.gens)
    pos = [findfirst(g -> _eq(g, v), gs) for v in p.gens]
    emax = _maxexp(bits)
    monos = similar(p.monos)
    for (k, m) in enumerate(p.monos)
        r = UInt64(0)
        for i in 1:np
            e = _exponent(m, p.bits, np, i)
            e > emax && throw(OverflowError("exponent $e does not fit in $bits bits"))
            r |= UInt64(e) << _shift(bits, n, pos[i])
        end
        monos[k] = r
    end
    perm = sortperm(monos; rev = true)
    Poly{C}(gs, bits, monos[perm], p.coefs[perm])
end

## arithmetic

function Base.:+(p::Poly{C}, q::Poly{C}) where {C}
    p, q = _unify(p, q)
    a, b, ac, bc = p.monos, q.monos, p.coefs, q.coefs
    na, nb = length(a), length(b)
    monos = Vector{UInt64}(undef, na + nb)
    coefs = Vector{C}(undef, na + nb)
    i, j, k = 1, 1, 0
    @inbounds while i <= na && j <= nb
        if a[i] > b[j]
            k += 1; monos[k] = a[i]; coefs[k] = ac[i]; i += 1
        elseif a[i] < b[j]
            k += 1; monos[k] = b[j]; coefs[k] = bc[j]; j += 1
        else
            s = ac[i] + bc[j]
            if !iszero(s)
                k += 1; monos[k] = a[i]; coefs[k] = s
            end
            i += 1; j += 1
        end
    end
    @inbounds while i <= na
        k += 1; monos[k] = a[i]; coefs[k] = ac[i]; i += 1
    end
    @inbounds while j <= nb
        k += 1; monos[k] = b[j]; coefs[k] = bc[j]; j += 1
    end
    Poly{C}(p.gens, p.bits, resize!(monos, k), resize!(coefs, k))
end

Base.:-(p::Poly{C}) where {C} = Poly{C}(p.gens, p.bits, p.monos, C[-c for c in p.coefs])
Base.:-(p::Poly{C}, q::Poly{C}) where {C} = p + (-q)

function _scale(p::Poly{C}, c::C) where {C}
    iszero(c) && return zero(p)
    isone(c) && return p
    cs = C[x * c for x in p.coefs]
    keep = findall(!iszero, cs)    # inexact or zero-divisor coefficients may vanish
    Poly{C}(p.gens, p.bits, p.monos[keep], cs[keep])
end

# p * c*m for a single term
function _mul_term(p::Poly{C}, m::UInt64, c::C) where {C}
    cs = C[x * c for x in p.coefs]
    keep = findall(!iszero, cs)
    Poly{C}(p.gens, p.bits, p.monos[keep] .+ m, cs[keep])
end

function Base.:*(p::Poly{C}, q::Poly{C}) where {C}
    p, q = _unify(p, q)
    (iszero(p) || iszero(q)) && return zero(p)
    _check_degrees(p, q)
    length(q.monos) == 1 && return _mul_term(p, q.monos[1], q.coefs[1])
    length(p.monos) == 1 && return _mul_term(q, p.monos[1], p.coefs[1])
    if C === SymNumber
        r = _mul_smallint(p, q)
        r === nothing || return r
    end
    t = _MonoTable{C}(length(p.monos) + length(q.monos))
    _mul_into!(t, p.monos, p.coefs, q.monos, q.coefs)
    monos, coefs = _contents(t)
    _sorted_poly(p.gens, p.bits, monos, coefs)
end

# Integer coefficients that fit in Int: accumulate in Int128 when the bound
# max|a| * max|b| * min(#p, #q) guarantees no overflow.
_issmallint(c::SymNumber) = c.tag == T_SMALL && c.d == 1

function _mul_smallint(p::Poly{SymNumber}, q::Poly{SymNumber})
    (all(_issmallint, p.coefs) && all(_issmallint, q.coefs)) || return nothing
    bound = Float64(maximum(c -> abs(c.n), p.coefs)) * Float64(maximum(c -> abs(c.n), q.coefs)) *
            min(length(p.coefs), length(q.coefs))
    bound < 2.0^125 || return nothing
    t = _MonoTable{Int128}(length(p.monos) + length(q.monos))
    _mul_into!(t, p.monos, Int128[c.n for c in p.coefs], q.monos, Int128[c.n for c in q.coefs])
    monos, acc = _contents(t)
    _sorted_poly(p.gens, p.bits, monos, SymNumber[SymNumber(a) for a in acc])
end

# Open-addressing hash table from packed monomial to accumulated coefficient.
const _EMPTY = typemax(UInt64)

mutable struct _MonoTable{T}
    keys::Vector{UInt64}
    vals::Vector{T}
    count::Int
end
function _MonoTable{T}(n::Int) where {T}
    cap = nextpow(2, max(16, 2n))
    _MonoTable{T}(fill(_EMPTY, cap), Vector{T}(undef, cap), 0)
end

@inline _slot(m::UInt64, mask::Int) = Int(((m * 0x9e3779b97f4a7c15) >> 32) % UInt64 & (mask % UInt64)) + 1

# add c to the entry for m; returns true if m was not present
@inline function _addto!(t::_MonoTable{T}, m::UInt64, c::T) where {T}
    keys = t.keys
    mask = length(keys) - 1
    i = _slot(m, mask)
    @inbounds while true
        k = keys[i]
        if k == m
            t.vals[i] += c
            return false
        elseif k == _EMPTY
            keys[i] = m
            t.vals[i] = c
            t.count += 1
            2 * t.count > length(keys) && _grow!(t)
            return true
        end
        i = (i & mask) + 1
    end
end

# the entry for m, which must be present
@inline function _lookup(t::_MonoTable, m::UInt64)
    keys = t.keys
    mask = length(keys) - 1
    i = _slot(m, mask)
    @inbounds while keys[i] != m
        i = (i & mask) + 1
    end
    @inbounds t.vals[i]
end

function _grow!(t::_MonoTable{T}) where {T}
    oldk, oldv = t.keys, t.vals
    t.keys = fill(_EMPTY, 2 * length(oldk))
    t.vals = Vector{T}(undef, 2 * length(oldk))
    t.count = 0
    @inbounds for i in eachindex(oldk)
        oldk[i] == _EMPTY || _addto!(t, oldk[i], oldv[i])
    end
end

function _mul_into!(t::_MonoTable{T}, pm, pc::Vector{T}, qm, qc::Vector{T}) where {T}
    @inbounds for i in eachindex(pm)
        mi, ci = pm[i], pc[i]
        for j in eachindex(qm)
            _addto!(t, mi + qm[j], ci * qc[j])
        end
    end
    t
end

function _contents(t::_MonoTable{T}) where {T}
    monos = Vector{UInt64}(undef, t.count)
    vals = Vector{T}(undef, t.count)
    k = 0
    @inbounds for i in eachindex(t.keys)
        if t.keys[i] != _EMPTY
            k += 1
            monos[k] = t.keys[i]
            vals[k] = t.vals[i]
        end
    end
    monos, vals
end

function _sorted_poly(gs, bits, monos::Vector{UInt64}, coefs::Vector{C}) where {C}
    perm = sortperm(monos; rev = true)
    keep = filter(k -> !iszero(coefs[k]), perm)
    Poly{C}(gs, bits, monos[keep], coefs[keep])
end

function Base.:^(p::Poly{C}, k::Integer) where {C}
    k < 0 && throw(DomainError(k, "negative powers of polynomials are not polynomials"))
    k == 0 && return one(p)
    k == 1 && return p
    emax = _maxexp(p.bits)
    degs = _maxdegs(p)
    any(d -> d * k > emax, degs) && throw(OverflowError("exponent exceeds $emax"))
    t = length(p.monos)
    t == 1 && return Poly{C}(p.gens, p.bits, [p.monos[1] * UInt64(k)], [p.coefs[1]^k])
    # The multinomial theorem computes each composition of k once, which is the result
    # size when few compositions share a monomial. Many collide when the base is dense
    # (e.g. 1 + x + ... + x^9), and then repeated multiplication by the small base is
    # cheaper: compare the number of compositions with the degree bound on the size.
    ncomp = binomial(big(k + t - 1), t - 1)
    sizebound = prod(big(k * d + 1) for d in degs)
    ncomp <= 4 * sizebound && return _multinomial_pow(p, Int(k))
    r = p
    for _ in 2:k
        r = r * p
    end
    r
end

# (c_1 m_1 + ... + c_t m_t)^k = Σ k!/(a_1!...a_t!) Π c_i^a_i m_i^a_i. Packed monomials
# make m_i^a_i a multiplication and the product a sum (degrees checked by the caller).
# Integer coefficients are computed in Int128 when (Σ|c_i|)^k bounds every coefficient
# below 2^125, else in BigInt, and converted once per term: SymNumber's big values are
# rationals, and normalizing them on every multiply-add dominated the cost.
function _multinomial_pow(p::Poly{C}, k::Int) where {C}
    if C === SymNumber && all(_isint, p.coefs)
        T = sum(c -> abs(big(_numerator(c))), p.coefs)^k < big(2)^125 ? Int128 : BigInt
        monos, coefs = _multinomial_terms(p.monos, T[T(_numerator(c)) for c in p.coefs], k, one(T))
        return _sorted_poly(p.gens, p.bits, monos, SymNumber[SymNumber(c) for c in coefs])
    end
    monos, coefs = _multinomial_terms(p.monos, p.coefs, k, N1)
    _sorted_poly(p.gens, p.bits, monos, coefs)
end

function _multinomial_terms(ms::Vector{UInt64}, cs::Vector{T}, k::Int, b1) where {T}
    cpow = [T[a == 0 ? one(T) : c^a for a in 0:k] for c in cs]   # exact 1, even for floats
    table = _MonoTable{T}(min(binomial(k + length(ms) - 1, length(ms) - 1), 1 << 20))
    _multinomial!(table, ms, cpow, 1, k, one(T), UInt64(0), b1)
    _contents(table)
end

# binomial(rem, a) from binomial(rem, a - 1)
_binstep(b::Integer, rem::Int, a::Int) = div(b * (rem - a + 1), a)
_binstep(b::SymNumber, rem::Int, a::Int) = b * SymNumber(rem - a + 1) / SymNumber(a)
_bincoef(::Type{T}, b::Integer) where {T<:Integer} = T(b)
_bincoef(::Type{C}, b::SymNumber) where {C} = _tocoef(C, b)

function _multinomial!(table::_MonoTable{T}, ms, cpow, i::Int, rem::Int, coef::T, mono::UInt64, b1) where {T}
    if rem == 0
        _addto!(table, mono, coef)
    elseif i == length(ms)
        _addto!(table, mono + UInt64(rem) * ms[i], coef * cpow[i][rem+1])
    else
        b = b1                                    # binomial(rem, a), updated in the loop
        for a in 0:rem
            a > 0 && (b = _binstep(b, rem, a))
            ca = coef * _bincoef(T, b) * cpow[i][a+1]
            _multinomial!(table, ms, cpow, i + 1, rem - a, ca, mono + UInt64(a) * ms[i], b1)
        end
    end
    table
end

# scalars and expressions lift into the polynomial's ring
_lift(p::Poly{C}, x::Number) where {C} = _constpoly(C, p.gens, p.bits, _tocoef(C, Sym(x)))
_lift(p::Poly{C}, x::Sym) where {C} = Poly{C}(x, p.gens)
for op in (:+, :-, :*)
    @eval Base.$op(p::Poly, x::Union{Number,Sym}) = $op(p, _lift(p, x))
    @eval Base.$op(x::Union{Number,Sym}, p::Poly) = $op(_lift(p, x), p)
end

"""
    diff(p::Poly, x)

Derivative of `p` with respect to the generator `x`.
"""
function Base.diff(p::Poly{C}, x::Sym) where {C}
    i = findfirst(g -> _eq(g, x), p.gens)
    i === nothing && return zero(p)
    n = length(p.gens)
    unit = UInt64(1) << _shift(p.bits, n, i)
    monos, coefs = UInt64[], C[]
    for (m, c) in zip(p.monos, p.coefs)
        e = _exponent(m, p.bits, n, i)
        e == 0 && continue
        push!(monos, m - unit)     # order is preserved: same unit removed from each
        push!(coefs, c * _tocoef(C, Sym(e)))
    end
    Poly{C}(p.gens, p.bits, monos, coefs)
end

## conversion from expressions

_tocoef(::Type{SymNumber}, x::Sym) =
    x.kind == NUMBER ? x.num : throw(ArgumentError("$x is not a numeric coefficient; use Poly{Sym}"))
_tocoef(::Type{Sym}, x::Sym) = x
_tocoef(::Type{SymNumber}, x::SymNumber) = x
_tocoef(::Type{Sym}, x::SymNumber) = Sym(x)

_nonneg_int(e::Sym) = e.kind == NUMBER && _isint(e.num) && sign(e.num) >= 0

function Poly(ex, gs::Sym...)
    ex = Sym(ex)
    g = collect(gs)
    _has_numeric_coefs(ex, g) ? Poly{SymNumber}(ex, g) : Poly{Sym}(ex, g)
end
Poly{C}(ex, gs::Sym...) where {C} = Poly{C}(Sym(ex), collect(gs))

function Poly{C}(ex::Sym, gs::Vector{Sym}) where {C}
    all(g -> g.kind != NUMBER, gs) || throw(ArgumentError("generators must not be numbers"))
    _topoly(C, ex, gs, _polybits(length(gs)), Dict{Sym,Int}(g => i for (i, g) in enumerate(gs)))
end

function _has_numeric_coefs(ex::Sym, gs::Vector{Sym})
    any(g -> _eq(g, ex), gs) && return true
    k = ex.kind
    k == NUMBER && return true
    k == ADD && return all(t -> _has_numeric_coefs(t, gs), ex.args)
    k == MUL && return all(i -> _nonneg_int(ex.exps[i]) && _has_numeric_coefs(ex.args[i], gs), eachindex(ex.args))
    k == POW && return _nonneg_int(ex.args[2]) && _has_numeric_coefs(ex.args[1], gs)
    false
end

function _topoly(::Type{C}, ex::Sym, gs, bits, index) where {C}
    i = get(index, ex, 0)
    i > 0 && return _genpoly(C, gs, bits, i)
    k = ex.kind
    k == NUMBER && return _constpoly(C, gs, bits, _tocoef(C, ex))
    if k == ADD
        monos, coefs = UInt64[], C[]
        iszero(ex.num) || (push!(monos, 0); push!(coefs, _tocoef(C, ex.num)))
        parts = Poly{C}[]
        for i in eachindex(ex.args)
            m = _monocode(ex.args[i], index, bits, length(gs))
            if m === nothing
                push!(parts, _scale(_topoly(C, ex.args[i], gs, bits, index), _tocoef(C, ex.coefs[i])))
            else
                push!(monos, m)
                push!(coefs, _tocoef(C, ex.coefs[i]))
            end
        end
        push!(parts, Poly{C}(gs, bits, monos, coefs))   # unsorted; `_sum_polys` sorts
        return _sum_polys(C, gs, bits, parts)
    elseif k == MUL
        r = _constpoly(C, gs, bits, _tocoef(C, ex.num))
        for i in eachindex(ex.args)
            r = r * _topoly(C, _mkpow(ex.args[i], ex.exps[i]), gs, bits, index)
        end
        return r
    elseif k == POW && _nonneg_int(ex.args[2])
        return _topoly(C, ex.args[1], gs, bits, index)^_toint(ex.args[2].num)
    end
    any(g -> has(ex, g), gs) && throw(ArgumentError("$ex is not polynomial in $(join(gs, ", "))"))
    _constpoly(C, gs, bits, _tocoef(C, ex))
end

# packed monomial of a product of generator powers, or `nothing`
function _monocode(ex::Sym, index, bits::Int, n::Int)
    i = get(index, ex, 0)
    i > 0 && return UInt64(1) << _shift(bits, n, i)
    k = ex.kind
    if k == POW
        i = get(index, ex.args[1], 0)
        (i > 0 && _nonneg_int(ex.args[2])) || return nothing
        e = _toint(ex.args[2].num)
        e > _maxexp(bits) && throw(OverflowError("exponent $e does not fit in $bits bits"))
        return UInt64(e) << _shift(bits, n, i)
    elseif k == MUL && _isN1(ex.num)
        m = UInt64(0)
        for j in eachindex(ex.args)
            mj = _monocode(_mkpow(ex.args[j], ex.exps[j]), index, bits, n)
            mj === nothing && return nothing
            m += mj               # distinct bases: distinct fields, no carries
        end
        return m
    end
    nothing
end

# sum of many polynomials: concatenate, sort once, combine runs
function _sum_polys(::Type{C}, gs, bits, ps) where {C}
    monos = reduce(vcat, (p.monos for p in ps); init = UInt64[])
    coefs = reduce(vcat, (p.coefs for p in ps); init = C[])
    perm = sortperm(monos; rev = true)
    om, oc = UInt64[], C[]
    for k in perm
        if !isempty(om) && om[end] == monos[k]
            oc[end] += coefs[k]
        else
            push!(om, monos[k])
            push!(oc, coefs[k])
        end
    end
    keep = findall(!iszero, oc)
    Poly{C}(gs, bits, om[keep], oc[keep])
end

## conversion to expressions

function Sym(p::Poly{C}) where {C}
    n = length(p.gens)
    if C === SymNumber && all(_isnamed, p.gens)
        return _sym_fast(p)
    end
    _add_many(mul(Sym(c), _monomial(p, m)) for (m, c) in zip(p.monos, p.coefs))
end

function _monomial(p::Poly, m::UInt64)
    n = length(p.gens)
    r = ONE
    for i in 1:n
        e = _exponent(m, p.bits, n, i)
        e > 0 && (r = mul(r, pow(p.gens[i], Sym(e))))
    end
    r
end

# Monomials in symbol generators are built directly in canonical form: bases visited in
# `_cmp` order are already sorted, and distinct monomials give distinct keys.
function _sym_fast(p::Poly{SymNumber})
    n = length(p.gens)
    order = sortperm(p.gens; lt = _lt)
    keys = Sym[]
    coefs = SymNumber[]
    c0 = N0
    for (m, c) in zip(p.monos, p.coefs)
        if m == 0
            c0 = c
            continue
        end
        bases, exps = Sym[], Sym[]
        for i in order
            e = _exponent(m, p.bits, n, i)
            e > 0 || continue
            push!(bases, p.gens[i])
            push!(exps, Sym(e))
        end
        push!(keys, _mkmul(N1, bases, exps))
        push!(coefs, c)
    end
    perm = sortperm(keys; lt = _lt)
    _mkadd(c0, keys[perm], coefs[perm])
end

Base.convert(::Type{Sym}, p::Poly) = Sym(p)
free_symbols(p::Poly) = free_symbols(Sym(p))
free_parameters(p::Poly) = free_parameters(Sym(p))

# called like the expression it represents (see calling a `Sym`)
(p::Poly)(args...) = Sym(p)(args...)

function Base.show(io::IO, p::Poly)
    print(io, "Poly(")
    if iszero(p)
        print(io, "0")
    end
    for (k, (m, c)) in enumerate(zip(p.monos, p.coefs))
        s = sprint(_show, mul(Sym(c), _monomial(p, m)), P_ADD)
        if k == 1
            print(io, s)
        elseif startswith(s, "-")
            print(io, " - ", chop(s; head = 1, tail = 0))
        else
            print(io, " + ", s)
        end
    end
    isempty(p.gens) || print(io, ", ", join(p.gens, ", "))
    print(io, ")")
end
Base.show(io::IO, ::MIME"text/plain", p::Poly) = show(io, p)
