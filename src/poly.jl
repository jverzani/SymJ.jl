# Sparse multivariate polynomials with packed exponents.
#
# A monomial's exponent vector is packed into one integer, `bits` bits per generator,
# the first generator in the most significant field. Then
#   * monomial multiplication is integer addition (no carries, guarded by degree checks),
#   * lex order (gens[1] > gens[2] > ...) is integer order,
# so terms are kept as parallel vectors sorted by decreasing packed monomial and sums are
# linear merges. Products accumulate in a hash table keyed by the packed monomial.
#
# The word is a UInt64 with `bits = 64 ÷ ngens` (the fast case, nearly always), else a
# UInt128 (`128 ÷ ngens`), else a BigInt (any width): slower, but no degree or number
# of generators is refused. The word type M is a type parameter, and arithmetic never
# changes it, so that a computation stays type-stable. An exponent that does not fit
# throws `_PackOverflow`; the entry points (`Poly(ex, ...)`, gcd, factor, cancel, ...)
# catch it and rerun the computation with wider words (`_with_widening`). BigInt words
# simply grow their fields.

# Field widths per generator (gens[1] in the most significant field), with the shifts
# and masks derived from them
struct _Layout
    widths::Vector{Int}
    shifts::Vector{Int}
    masks::Vector{UInt64}
end
function _Layout(widths::Vector{Int})
    n = length(widths)
    shifts = Vector{Int}(undef, n)
    s = 0
    for i in n:-1:1
        shifts[i] = s
        s += widths[i]
    end
    _Layout(widths, shifts, UInt64[_mask(w) for w in widths])
end
Base.:(==)(a::_Layout, b::_Layout) = a === b || a.widths == b.widths
_total(L::_Layout) = sum(L.widths; init = 0)

"""
    Poly{C}

Sparse multivariate polynomial in the generators `gens(p)` with coefficients of type `C`:
`SymNumber` for numeric coefficients, or `Sym` for symbolic ones. Terms are stored in
lex order with `gens(p)[1]` most significant.

    Poly(ex, x, y, ...)

Convert the expression `ex` to a polynomial in `x, y, ...`. Coefficients are `SymNumber`
when they are all numeric, and `Sym` otherwise (so `Poly(a*x^2 + b, x)` has coefficients
`a` and `b`). Use `Poly{C}(ex, x, ...)` to choose the coefficient type. `Sym(p)` converts
back. Degrees and the number of generators are unlimited; polynomials whose exponents
fit in 64 bits (`64 ÷ ngens` bits each) are the fast case, and larger ones are packed
into wider words. Arithmetic on polynomials keeps their packing: `p * q` or `p^k` whose
exponents do not fit throws an error, while `Poly`, `gcd`, `factor` and the functions
on expressions widen and carry on.
"""
struct Poly{C,M<:Integer} <: AbstractSymbolic
    gens::Vector{Sym}
    layout::_Layout
    monos::Vector{M}        # strictly decreasing
    coefs::Vector{C}        # nonzero
end
Poly{C}(gens::Vector{Sym}, L::_Layout, monos::Vector{M}, coefs::Vector) where {C,M<:Integer} =
    Poly{C,M}(gens, L, monos, coefs)

# Storing into a container of polynomials of another word repacks if the exponents fit
function Base.convert(::Type{Poly{C,M}}, p::Poly{C}) where {C,M}
    p isa Poly{C,M} && return p
    L = _layout_in(M, _maxdegs(p))
    L === nothing && throw(_PackOverflow(maximum(_maxdegs(p); init = 0)))
    _repack(p, p.gens, L, M)
end

# an exponent that does not fit the packing in use: caught by `_with_widening`
struct _PackOverflow <: Exception
    e::Int
end
Base.showerror(io::IO, err::_PackOverflow) = print(io, "OverflowError: exponent $(err.e) ",
    "does not fit the packed exponents of these polynomials; operations such as gcd, ",
    "factor and the functions on expressions widen the packing automatically")

_mask(bits::Int) = bits == 64 ? typemax(UInt64) : (UInt64(1) << bits) - 1
# largest allowed exponent; an all-ones field is never used, so the packed monomial
# typemax(UInt64) (UInt128) is free to serve as the empty-slot marker in `_MonoTable`
_maxexp(bits::Int) = Int(min(_mask(bits) - 1, UInt64(typemax(Int))))
_maxexp(L::_Layout, i::Int) = _maxexp(L.widths[i])
@inline _shift(L::_Layout, n::Int, i::Int) = @inbounds L.shifts[i]
@inline _fieldmask(L::_Layout, i::Int) = @inbounds L.masks[i]
@inline _exponent(m::Integer, L::_Layout, n::Int, i::Int) =
    @inbounds Int((m >> L.shifts[i]) & L.masks[i])
_fits(L::_Layout, degs) = all(i -> degs[i] <= _maxexp(L, i), eachindex(degs))

# the width of a field holding exponent e (≤ 62 bits, so exponents stay Ints)
function _bitsfor(e::Integer)
    e < 2^61 || throw(OverflowError("exponent $e is too large"))
    64 - leading_zeros(UInt64(e + 1))
end

# the word holding a layout; the uniform width of a fixed-size word
_monotype(L::_Layout) = (t = _total(L); t <= 64 ? UInt64 : t <= 128 ? UInt128 : BigInt)
_monotype(p::Poly{C,M}) where {C,M} = M
_capacity(::Type{UInt64}) = 64
_capacity(::Type{UInt128}) = 128
_fullbits(::Type{M}, n::Int) where {M} = n == 0 ? 64 : n == 1 && M === UInt64 ? 64 : min(_capacity(M) ÷ n, 62)

# uniform layouts: the full-width ones for up to 64 generators precomputed (no lookup
# on the common path), others built on demand
const _UNIFORM64 = [_Layout(fill(_fullbits(UInt64, n), n)) for n in 0:64]
const _UNIFORM128 = [_Layout(fill(_fullbits(UInt128, n), n)) for n in 0:64]
function _uniform(n::Int, w::Int)
    n <= 64 && w == _fullbits(UInt64, n) && return @inbounds _UNIFORM64[n+1]
    n <= 64 && w == _fullbits(UInt128, n) && return @inbounds _UNIFORM128[n+1]
    _Layout(fill(w, n))
end

# Within `_with_widening` a floor (`_minbits`) on every field forces wider words.
const _MINBITS = :SymJ_poly_minbits
_minbits() = get(task_local_storage(), _MINBITS, 0)::Int

# The layout for exponents up to degs (one per generator), in the narrowest word at
# least `M0`: the uniform one (the widest equal fields) whenever it holds them, which
# is nearly always and the same layout for all such polynomials; otherwise fields sized
# to the degrees, the spare bits shared out as headroom (say one generator of high
# degree among many linear ones, still in one UInt64).
function _layout_for(degs::AbstractVector{<:Integer}, ::Type{M0} = UInt64) where {M0}
    n = length(degs)
    n == 0 && return _uniform(0, 64)
    fl = _minbits()
    # the common case first, without allocating: the uniform UInt64 layout holds degs
    if fl == 0 && M0 === UInt64
        u = _fullbits(UInt64, n)
        emax = _maxexp(u)
        u >= 1 && all(d -> d <= emax, degs) && return _uniform(n, u)
    end
    need = Int[max(_bitsfor(d), fl, 1) for d in degs]
    for M in (UInt64, UInt128)
        (M0 === UInt128 && M === UInt64 || M0 === BigInt) && continue
        u = _fullbits(M, n)
        u >= 1 && all(<=(u), need) && return _uniform(n, u)
        cap = _capacity(M)
        tot = sum(need)
        tot <= cap || continue
        spare = cap - tot
        w = [min(need[i] + spare ÷ n + (i <= spare % n), 62) for i in 1:n]
        return _Layout(w)
    end
    _Layout(Int[min(max(2 * k, 2), 62) for k in need])          # BigInt: room to grow
end
# a layout for degs within the word M, or nothing if M cannot hold them
function _layout_in(::Type{M}, degs) where {M}
    L = _layout_for(degs, M)
    (M === BigInt || _total(L) <= _capacity(M)) ? L : nothing
end
_polylayout(n::Int) = _layout_for(zeros(Int, n))

# Run f(args...), and if a packed exponent overflows, again with wider words (the
# polynomial arguments repacked): the slow path, instead of failing.
function _with_widening(f, args...; kws...)
    try
        return f(args...; kws...)
    catch err
        err isa _PackOverflow || rethrow()
        _minbits() >= 62 && throw(OverflowError("exponent $(err.e) is too large"))
        floor = max(_minbits() + 1, _bitsfor(err.e))
        return task_local_storage(_MINBITS, floor) do
            _with_widening(f, map(_rewiden, args)...; kws...)
        end
    end
end
_rewiden(x) = x
function _rewiden(p::Poly)
    L = _layout_for(_maxdegs(p), _monotype(p))
    L == p.layout ? p : _repack(p, p.gens, L, _monotype(L) === BigInt ? BigInt : promote_type(_monotype(L), _monotype(p)))
end

_zeropoly(::Type{C}, gens, L::_Layout, ::Type{M} = _monotype(L)) where {C,M} =
    Poly{C,M}(gens, L, M[], C[])
_constpoly(::Type{C}, gens, L::_Layout, c, ::Type{M} = _monotype(L)) where {C,M} =
    iszero(c) ? _zeropoly(C, gens, L, M) : Poly{C,M}(gens, L, M[zero(M)], C[c])
_genpoly(::Type{C}, gens, L::_Layout, i::Int, ::Type{M} = _monotype(L)) where {C,M} =
    Poly{C,M}(gens, L, M[one(M) << _shift(L, length(gens), i)], C[one(C)])

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
    [[_exponent(m, p.layout, length(p.gens), i) for i in eachindex(p.gens)] for m in p.monos]

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
    maximum((sum(_exponent(m, p.layout, n, i) for i in 1:n; init = 0) for m in p.monos); init = -1)
end
function degree(p::Poly, x::Sym)
    i, n = _genindex(p, x), length(p.gens)
    maximum((_exponent(m, p.layout, n, i) for m in p.monos); init = -1)
end

function _maxdegs(p::Poly)
    n = length(p.gens)
    d = zeros(Int, n)
    for m in p.monos, i in 1:n
        d[i] = max(d[i], _exponent(m, p.layout, n, i))
    end
    d
end

# Room for exponents up to degs (for a product, the sums of the degrees): a new layout
# in the same word if it holds them (type-stable), else `_PackOverflow` (BigInt words
# always hold them)
function _fit_product(p::Poly{C,M}, q::Poly{C,M}) where {C,M}
    degs = _maxdegs(p) .+ _maxdegs(q)
    _fits(p.layout, degs) && return p, q
    L = _layout_in(M, degs)
    L === nothing && throw(_PackOverflow(maximum(degs)))
    _repack(p, p.gens, L, M), _repack(q, q.gens, L, M)
end
function _fit_degrees(p::Poly{C,M}, degs) where {C,M}
    _fits(p.layout, degs) && return p
    L = _layout_in(M, degs)
    L === nothing && throw(_PackOverflow(min(maximum(degs), typemax(Int))))
    _repack(p, p.gens, L, M)
end

Base.iszero(p::Poly) = isempty(p.monos)
Base.isone(p::Poly) = length(p.monos) == 1 && p.monos[1] == 0 && isone(p.coefs[1])
Base.zero(p::Poly{C,M}) where {C,M} = _zeropoly(C, p.gens, p.layout, M)
Base.one(p::Poly{C,M}) where {C,M} = _constpoly(C, p.gens, p.layout, one(C), M)

function Base.:(==)(p::Poly, q::Poly)
    _samegens(p, q) || return false
    p.layout == q.layout && return p.monos == q.monos && p.coefs == q.coefs
    length(p.monos) == length(q.monos) && p.coefs == q.coefs &&       # different packings
        exponent_vectors(p) == exponent_vectors(q)
end
# independent of the packing (equal polynomials may be packed differently)
Base.hash(p::Poly, h::UInt) = hash(p.coefs, hash(exponent_vectors(p), hash(p.gens, h)))

## generators: arithmetic requires a common ring; differing gens are unified

_samegens(p::Poly, q::Poly) = p.gens === q.gens || (length(p.gens) == length(q.gens) && all(_eq.(p.gens, q.gens)))

# common generators and layout, in the same word (type-stable): the union of the
# generators must fit it, else `_PackOverflow`
function _unify(p::Poly{C,M}, q::Poly{C,M}) where {C,M}
    if _samegens(p, q)
        p.layout == q.layout && return (p, q)
        degs = max.(_maxdegs(p), _maxdegs(q))
        L = _fits(p.layout, degs) ? p.layout : _fits(q.layout, degs) ? q.layout : _layout_in(M, degs)
        L === nothing && throw(_PackOverflow(maximum(degs)))
        return (_repack(p, p.gens, L, M), _repack(q, p.gens, L, M))
    end
    gs = copy(p.gens)
    for g in q.gens
        any(h -> _eq(g, h), gs) || push!(gs, g)
    end
    degs = zeros(Int, length(gs))
    for (r, d) in ((p, _maxdegs(p)), (q, _maxdegs(q))), (k, v) in enumerate(r.gens)
        j = findfirst(g -> _eq(g, v), gs)
        degs[j] = max(degs[j], d[k])
    end
    L = _layout_in(M, degs)
    L === nothing && throw(_PackOverflow(maximum(degs)))
    (_repack(p, gs, L, M), _repack(q, gs, L, M))
end
# different words: the narrower one widened (not type-stable; only at the boundaries)
# (explicit M1, M2: with `p::Poly{C}, q::Poly{C}` this method would be chosen over the
# one above even for equal words)
function _unify(p::Poly{C,M1}, q::Poly{C,M2}) where {C,M1,M2}
    rank(M) = M === UInt64 ? 1 : M === UInt128 ? 2 : 3
    if rank(M1) < rank(M2)
        q2, p2 = _unify(q, convert(Poly{C,M2}, p))
        return p2, q2
    end
    _unify(p, convert(Poly{C,M1}, q))
end

_repack(p::Poly{C}, gs::Vector{Sym}, L::_Layout) where {C} = _repack(p, gs, L, _monotype(L))
function _repack(p::Poly{C}, gs::Vector{Sym}, L::_Layout, ::Type{M}) where {C,M}
    n, np = length(gs), length(p.gens)
    p isa Poly{C,M} && p.layout == L && p.gens === gs && return p
    pos = [findfirst(g -> _eq(g, v), gs) for v in p.gens]
    monos = Vector{M}(undef, length(p.monos))
    for (k, m) in enumerate(p.monos)
        r = zero(M)
        for i in 1:np
            e = _exponent(m, p.layout, np, i)
            e > _maxexp(L, pos[i]) && throw(_PackOverflow(e))
            r |= M(e) << _shift(L, n, pos[i])
        end
        monos[k] = r
    end
    perm = sortperm(monos; rev = true)
    Poly{C,M}(gs, L, monos[perm], p.coefs[perm])
end

## arithmetic

function Base.:+(p::Poly{C}, q::Poly{C}) where {C}
    p, q = _unify(p, q)
    a, b, ac, bc = p.monos, q.monos, p.coefs, q.coefs
    na, nb = length(a), length(b)
    monos = similar(a, na + nb)
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
    Poly{C}(p.gens, p.layout, resize!(monos, k), resize!(coefs, k))
end

Base.:-(p::Poly{C}) where {C} = Poly{C}(p.gens, p.layout, p.monos, C[-c for c in p.coefs])
Base.:-(p::Poly{C}, q::Poly{C}) where {C} = p + (-q)

function _scale(p::Poly{C}, c::C) where {C}
    iszero(c) && return zero(p)
    isone(c) && return p
    cs = C[x * c for x in p.coefs]
    keep = findall(!iszero, cs)    # inexact or zero-divisor coefficients may vanish
    Poly{C}(p.gens, p.layout, p.monos[keep], cs[keep])
end

# p * c*m for a single term
function _mul_term(p::Poly{C}, m::Integer, c::C) where {C}
    cs = C[x * c for x in p.coefs]
    keep = findall(!iszero, cs)
    Poly{C}(p.gens, p.layout, p.monos[keep] .+ m, cs[keep])
end

# p * x_j^e, widened if the exponent needs it (safe where a packed unit is not)
function _mul_genpow(p::Poly{C}, j::Int, e::Int) where {C}
    (iszero(p) || e == 0) && return p
    degs = _maxdegs(p)
    degs[j] += e
    p = _fit_degrees(p, degs)
    _mul_term(p, _monotype(p)(e) << _shift(p.layout, length(p.gens), j), one(C))
end

function Base.:*(p::Poly{C}, q::Poly{C}) where {C}
    p, q = _unify(p, q)
    (iszero(p) || iszero(q)) && return zero(p)
    p, q = _fit_product(p, q)
    length(q.monos) == 1 && return _mul_term(p, q.monos[1], q.coefs[1])
    length(p.monos) == 1 && return _mul_term(q, p.monos[1], p.coefs[1])
    if C === SymNumber
        r = _mul_smallint(p, q)
        r === nothing || return r
    end
    t = _MonoTable{C}(length(p.monos) + length(q.monos), _monotype(p))
    _mul_into!(t, p.monos, p.coefs, q.monos, q.coefs)
    monos, coefs = _contents(t)
    _sorted_poly(p.gens, p.layout, monos, coefs)
end

# Integer coefficients that fit in Int: accumulate in Int128 when the bound
# max|a| * max|b| * min(#p, #q) guarantees no overflow.
_issmallint(c::SymNumber) = c.tag == T_SMALL && c.d == 1

function _mul_smallint(p::Poly{SymNumber}, q::Poly{SymNumber})
    (all(_issmallint, p.coefs) && all(_issmallint, q.coefs)) || return nothing
    bound = Float64(maximum(c -> abs(c.n), p.coefs)) * Float64(maximum(c -> abs(c.n), q.coefs)) *
            min(length(p.coefs), length(q.coefs))
    bound < 2.0^125 || return nothing
    t = _MonoTable{Int128}(length(p.monos) + length(q.monos), _monotype(p))
    _mul_into!(t, p.monos, Int128[c.n for c in p.coefs], q.monos, Int128[c.n for c in q.coefs])
    monos, acc = _contents(t)
    _sorted_poly(p.gens, p.layout, monos, SymNumber[SymNumber(a) for a in acc])
end

# Open-addressing hash table from packed monomial to accumulated coefficient.
const _EMPTY = typemax(UInt64)
_emptykey(::Type{UInt64}) = typemax(UInt64)
_emptykey(::Type{UInt128}) = typemax(UInt128)
_emptykey(::Type{BigInt}) = big(-1)               # packed monomials are nonnegative

mutable struct _MonoTable{M,T}
    keys::Vector{M}
    vals::Vector{T}
    count::Int
end
function _MonoTable{T}(n::Int, ::Type{M} = UInt64) where {T,M}
    cap = nextpow(2, max(16, 2n))
    _MonoTable{M,T}(fill(_emptykey(M), cap), Vector{T}(undef, cap), 0)
end

@inline _slot(m::UInt64, mask::Int) = Int(((m * 0x9e3779b97f4a7c15) >> 32) % UInt64 & (mask % UInt64)) + 1
@inline _slot(m::UInt128, mask::Int) = _slot((m % UInt64) ⊻ ((m >> 64) % UInt64) * 0xbf58476d1ce4e5b9, mask)
@inline _slot(m::BigInt, mask::Int) = Int(hash(m) & (mask % UInt)) + 1

# add c to the entry for m; returns true if m was not present
@inline function _addto!(t::_MonoTable{M,T}, m::M, c::T) where {M,T}
    keys = t.keys
    mask = length(keys) - 1
    i = _slot(m, mask)
    @inbounds while true
        k = keys[i]
        if k == m
            t.vals[i] += c
            return false
        elseif k == _emptykey(M)
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
@inline function _lookup(t::_MonoTable{M}, m::M) where {M}
    keys = t.keys
    mask = length(keys) - 1
    i = _slot(m, mask)
    @inbounds while keys[i] != m
        i = (i & mask) + 1
    end
    @inbounds t.vals[i]
end

function _grow!(t::_MonoTable{M,T}) where {M,T}
    oldk, oldv = t.keys, t.vals
    t.keys = fill(_emptykey(M), 2 * length(oldk))
    t.vals = Vector{T}(undef, 2 * length(oldk))
    t.count = 0
    @inbounds for i in eachindex(oldk)
        oldk[i] == _emptykey(M) || _addto!(t, oldk[i], oldv[i])
    end
end

function _mul_into!(t::_MonoTable{M,T}, pm, pc::Vector{T}, qm, qc::Vector{T}) where {M,T}
    @inbounds for i in eachindex(pm)
        mi, ci = pm[i], pc[i]
        for j in eachindex(qm)
            _addto!(t, mi + qm[j], ci * qc[j])
        end
    end
    t
end

function _contents(t::_MonoTable{M,T}) where {M,T}
    monos = Vector{M}(undef, t.count)
    vals = Vector{T}(undef, t.count)
    k = 0
    @inbounds for i in eachindex(t.keys)
        if t.keys[i] != _emptykey(M)
            k += 1
            monos[k] = t.keys[i]
            vals[k] = t.vals[i]
        end
    end
    monos, vals
end

function _sorted_poly(gs, L, monos::Vector{<:Integer}, coefs::Vector{C}) where {C}
    perm = sortperm(monos; rev = true)
    keep = filter(k -> !iszero(coefs[k]), perm)
    Poly{C}(gs, L, monos[keep], coefs[keep])
end

function Base.:^(p::Poly{C}, k::Integer) where {C}
    k < 0 && throw(DomainError(k, "negative powers of polynomials are not polynomials"))
    k == 0 && return one(p)
    k == 1 && return p
    p = _fit_degrees(p, [Int(min(d * big(k), big(typemax(Int)))) for d in _maxdegs(p)])
    degs = _maxdegs(p)
    t = length(p.monos)
    t == 1 && return Poly{C}(p.gens, p.layout, [p.monos[1] * _monotype(p)(k)], [p.coefs[1]^k])
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
        return _sorted_poly(p.gens, p.layout, monos, SymNumber[SymNumber(c) for c in coefs])
    end
    monos, coefs = _multinomial_terms(p.monos, p.coefs, k, N1)
    _sorted_poly(p.gens, p.layout, monos, coefs)
end

function _multinomial_terms(ms::Vector{M}, cs::Vector{T}, k::Int, b1) where {M,T}
    cpow = [T[a == 0 ? one(T) : c^a for a in 0:k] for c in cs]   # exact 1, even for floats
    table = _MonoTable{T}(min(binomial(k + length(ms) - 1, length(ms) - 1), 1 << 20), M)
    _multinomial!(table, ms, cpow, 1, k, one(T), zero(M), b1)
    _contents(table)
end

# binomial(rem, a) from binomial(rem, a - 1)
_binstep(b::Integer, rem::Int, a::Int) = div(b * (rem - a + 1), a)
_binstep(b::SymNumber, rem::Int, a::Int) = b * SymNumber(rem - a + 1) / SymNumber(a)
_bincoef(::Type{T}, b::Integer) where {T<:Integer} = T(b)
_bincoef(::Type{C}, b::SymNumber) where {C} = _tocoef(C, b)

function _multinomial!(table::_MonoTable{M,T}, ms, cpow, i::Int, rem::Int, coef::T, mono::M, b1) where {M,T}
    if rem == 0
        _addto!(table, mono, coef)
    elseif i == length(ms)
        _addto!(table, mono + M(rem) * ms[i], coef * cpow[i][rem+1])
    else
        b = b1                                    # binomial(rem, a), updated in the loop
        for a in 0:rem
            a > 0 && (b = _binstep(b, rem, a))
            ca = coef * _bincoef(T, b) * cpow[i][a+1]
            _multinomial!(table, ms, cpow, i + 1, rem - a, ca, mono + M(a) * ms[i], b1)
        end
    end
    table
end

# scalars and expressions lift into the polynomial's ring
_lift(p::Poly{C,M}, x::Number) where {C,M} = _constpoly(C, p.gens, p.layout, _tocoef(C, Sym(x)), M)
_lift(p::Poly{C,M}, x::Sym) where {C,M} =
    _topoly(C, M, x, p.gens, p.layout, Dict{Sym,Int}(g => i for (i, g) in enumerate(p.gens)))
for op in (:+, :-, :*)
    @eval Base.$op(p::Poly, x::Union{Number,Sym}) = $op(p, _lift(p, x))
    @eval Base.$op(x::Union{Number,Sym}, p::Poly) = $op(_lift(p, x), p)
end
# division by a constant (by a polynomial in the generators: divrem or divexact)
function Base.:/(p::Poly, x::Union{Number,Sym})
    c = Sym(x)
    any(g -> has(c, g), p.gens) &&
        throw(ArgumentError("divide by a polynomial in the generators with divrem or divexact"))
    p * (ONE / c)
end
Base.://(p::Poly, x::Union{Number,Sym}) = p / x

"""
    diff(p::Poly, x)

Derivative of `p` with respect to the generator `x`.
"""
function Base.diff(p::Poly{C,M}, x::Sym) where {C,M}
    i = findfirst(g -> _eq(g, x), p.gens)
    i === nothing && return zero(p)
    n = length(p.gens)
    unit = one(M) << _shift(p.layout, n, i)
    monos, coefs = M[], C[]
    for (m, c) in zip(p.monos, p.coefs)
        e = _exponent(m, p.layout, n, i)
        e == 0 && continue
        push!(monos, m - unit)     # order is preserved: same unit removed from each
        push!(coefs, c * _tocoef(C, Sym(e)))
    end
    Poly{C}(p.gens, p.layout, monos, coefs)
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
    _checked_gens(_has_numeric_coefs(ex, g) ? Poly{SymNumber}(ex, g) : Poly{Sym}(ex, g))
end
Poly{C}(ex, gs::Sym...) where {C} = _checked_gens(Poly{C}(Sym(ex), collect(gs)))

# a generator that is an expression (1/x, sin(x)) must not have its symbols in a
# coefficient: x + 1/x is not a polynomial in 1/x (coefficients [1, x])
function _checked_gens(p::Poly)
    for g in p.gens
        g.kind == SYMBOL && continue
        ss = free_symbols(g)
        for c in coefficients(p)
            c isa Sym && any(s -> has(c, s), ss) &&
                throw(ArgumentError("not a polynomial in $g: the coefficient $c depends on it"))
        end
    end
    p
end

function Poly{C}(ex::Sym, gs::Vector{Sym}) where {C}
    all(g -> g.kind != NUMBER, gs) || throw(ArgumentError("generators must not be numbers"))
    index = Dict{Sym,Int}(g => i for (i, g) in enumerate(gs))
    n = length(gs)
    # the uniform layout nearly always holds ex; if not, a layout from bounds on its
    # degrees (products and powers inside find room for themselves)
    L = _polylayout(n)
    try
        return _topoly(C, _monotype(L), ex, gs, L, index)
    catch err
        err isa _PackOverflow || rethrow()
    end
    L = _layout_for(_degbound(ex, index, n))
    _topoly(C, _monotype(L), ex, gs, L, index)
end

# an upper bound on the degree of ex in each generator (exact for an expanded sum)
function _degbound(ex::Sym, index, n::Int)
    d = zeros(Int, n)
    _degbound!(d, zeros(Int, n), ex, index)
    d
end
# d .= max.(d, the degrees of ex); t is scratch
function _degbound!(d, t, ex::Sym, index)
    if ex.kind == ADD
        for a in ex.args
            if _monodeg!(t, a, index)
                d .= max.(d, t)
            else
                d .= max.(d, _degbound(a, index, length(d)))
            end
        end
    elseif _monodeg!(t, ex, index)
        d .= max.(d, t)
    else                                                # a product: degrees add
        s = zeros(Int, length(d))
        factors = ex.kind == MUL ? [_mkpow(ex.args[i], ex.exps[i]) for i in eachindex(ex.args)] : [ex]
        for f in factors
            if f.kind == POW && _nonneg_int(f.args[2])
                s .+= _toint(f.args[2].num) .* _degbound(f.args[1], index, length(d))
            elseif _monodeg!(t, f, index)
                s .+= t
            elseif f.kind == ADD
                s .+= _degbound(f, index, length(d))
            end                         # else not polynomial: `_topoly` reports it
        end
        d .= max.(d, s)
    end
    d
end
# t = the exponents of a product of generator powers or a constant (true), else false
function _monodeg!(t, ex::Sym, index)
    fill!(t, 0)
    i = get(index, ex, 0)
    i > 0 && (t[i] = 1; return true)
    k = ex.kind
    if k == POW
        i = get(index, ex.args[1], 0)
        i > 0 && _nonneg_int(ex.args[2]) && (t[i] = _toint(ex.args[2].num); return true)
        return !any(g -> has(ex, g), keys(index))
    elseif k == MUL
        for j in eachindex(ex.args)
            i = get(index, ex.args[j], 0)
            (i > 0 && _nonneg_int(ex.exps[j])) || return !any(g -> has(ex, g), keys(index))
            t[i] += _toint(ex.exps[j].num)
        end
        return true
    elseif k == ADD
        return false
    end
    true                                                # a number or other constant
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

function _topoly(::Type{C}, ::Type{M}, ex::Sym, gs, L::_Layout, index) where {C,M}
    i = get(index, ex, 0)
    i > 0 && return _genpoly(C, gs, L, i, M)
    k = ex.kind
    k == NUMBER && return _constpoly(C, gs, L, _tocoef(C, ex), M)
    if k == ADD
        monos, coefs = M[], C[]
        iszero(ex.num) || (push!(monos, 0); push!(coefs, _tocoef(C, ex.num)))
        parts = Poly{C,M}[]
        for i in eachindex(ex.args)
            m = _monocode(M, ex.args[i], index, L, length(gs))
            if m === nothing
                push!(parts, _scale(_topoly(C, M, ex.args[i], gs, L, index), _tocoef(C, ex.coefs[i])))
            else
                push!(monos, m)
                push!(coefs, _tocoef(C, ex.coefs[i]))
            end
        end
        push!(parts, Poly{C,M}(gs, L, monos, coefs))   # unsorted; `_sum_polys` sorts
        return _sum_polys(C, M, gs, L, parts)
    elseif k == MUL
        r = _constpoly(C, gs, L, _tocoef(C, ex.num), M)
        for i in eachindex(ex.args)
            r = r * _topoly(C, M, _mkpow(ex.args[i], ex.exps[i]), gs, L, index)
        end
        return r
    elseif k == POW && _nonneg_int(ex.args[2])
        return _topoly(C, M, ex.args[1], gs, L, index)^_toint(ex.args[2].num)
    end
    any(g -> has(ex, g), gs) && throw(ArgumentError("$ex is not polynomial in $(join(gs, ", "))"))
    _constpoly(C, gs, L, _tocoef(C, ex), M)
end

# packed monomial of a product of generator powers, or `nothing`
function _monocode(::Type{M}, ex::Sym, index, L::_Layout, n::Int) where {M}
    i = get(index, ex, 0)
    i > 0 && return one(M) << _shift(L, n, i)
    k = ex.kind
    if k == POW
        i = get(index, ex.args[1], 0)
        (i > 0 && _nonneg_int(ex.args[2])) || return nothing
        e = _toint(ex.args[2].num)
        e > _maxexp(L, i) && throw(_PackOverflow(e))
        return M(e) << _shift(L, n, i)
    elseif k == MUL && _isN1(ex.num)
        m = zero(M)
        for j in eachindex(ex.args)
            mj = _monocode(M, _mkpow(ex.args[j], ex.exps[j]), index, L, n)
            mj === nothing && return nothing
            m += mj               # distinct bases: distinct fields, no carries
        end
        return m
    end
    nothing
end

# sum of many polynomials: concatenate, sort once, combine runs
function _sum_polys(::Type{C}, ::Type{M}, gs, L, ps) where {C,M}
    # parts built in other layouts (none, given the degree bounds): bring all to one
    if any(p -> !(p.layout == L), ps)
        degs = reduce((a, p) -> max.(a, _maxdegs(p)), ps; init = zeros(Int, length(gs)))
        L = _layout_in(M, degs)
        L === nothing && throw(_PackOverflow(maximum(degs)))
        ps = [_repack(p, gs, L, M) for p in ps]
    end
    monos = reduce(vcat, (p.monos for p in ps); init = M[])
    coefs = reduce(vcat, (p.coefs for p in ps); init = C[])
    perm = sortperm(monos; rev = true)
    om, oc = M[], C[]
    for k in perm
        if !isempty(om) && om[end] == monos[k]
            oc[end] += coefs[k]
        else
            push!(om, monos[k])
            push!(oc, coefs[k])
        end
    end
    keep = findall(!iszero, oc)
    Poly{C}(gs, L, om[keep], oc[keep])
end

## conversion to expressions

function Sym(p::Poly{C}) where {C}
    n = length(p.gens)
    if C === SymNumber && all(_isnamed, p.gens)
        return _sym_fast(p)
    end
    _add_many(mul(Sym(c), _monomial(p, m)) for (m, c) in zip(p.monos, p.coefs))
end

function _monomial(p::Poly, m::Integer)
    n = length(p.gens)
    r = ONE
    for i in 1:n
        e = _exponent(m, p.layout, n, i)
        e > 0 && (r = mul(r, pow(p.gens[i], Sym(e))))
    end
    r
end

# Monomials in symbol generators are built directly in canonical form: bases visited in
# `_cmp` order are already sorted, and distinct monomials give distinct keys.
function _sym_fast(p::Poly{SymNumber})
    n = length(p.gens)
    order = sortperm(p.gens; lt = _lt)
    shifts = [_shift(p.layout, n, i) for i in order]
    masks = [_fieldmask(p.layout, i) for i in order]
    gens = p.gens[order]
    keys = Vector{Sym}(undef, length(p.monos))
    coefs = Vector{SymNumber}(undef, length(p.monos))
    c0 = N0
    k = 0
    for (m, c) in zip(p.monos, p.coefs)
        if m == 0
            c0 = c
            continue
        end
        cnt = 0                                   # exact sizes: no growing by push!
        for (s, mask) in zip(shifts, masks)
            (m >> s) & mask == 0 || (cnt += 1)
        end
        bases = Vector{Sym}(undef, cnt)
        exps = Vector{Sym}(undef, cnt)
        j = 0
        for (g, s, mask) in zip(gens, shifts, masks)
            e = Int((m >> s) & mask)
            e > 0 || continue
            j += 1
            bases[j] = g
            exps[j] = Sym(SymNumber(e))
        end
        k += 1
        keys[k] = _mkmul(N1, bases, exps)
        coefs[k] = c
    end
    resize!(keys, k)
    resize!(coefs, k)
    perm = _canonical_perm(keys)
    _mkadd(c0, permute!(keys, perm), permute!(coefs, perm))
end

Base.convert(::Type{Sym}, p::Poly) = Sym(p)
free_symbols(p::Poly) = free_symbols(Sym(p))
free_parameters(p::Poly) = free_parameters(Sym(p))

# called like the expression it represents (see calling a `Sym`)
(p::Poly)(args...; kws...) = Sym(p)(args...; kws...)

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
