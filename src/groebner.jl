# Gröbner bases: Buchberger's algorithm with the Gebauer–Möller criteria and the sugar
# strategy, and the FGLM change of order.
#
# Polynomials are kept in their own representation, `_GP`: per term, the exponents
# packed into a word (`m`, one field per variable plus a guard bit above each field)
# and an order key (`k`), a word whose integer order is the monomial order. Both are
# linear in the exponents, so multiplying by a monomial adds to both:
#   lex      k = m
#   grlex    k = [deg | e1 | ... | en]
#   grevlex  k = [deg | e1+...+e(n-1) | ... | e1+e2 | e1]
# (grevlex: at equal degree, a smaller last exponent wins, which is a larger prefix
# sum). Divisibility is one subtraction against the guard bits, and an exponent or a
# degree that outgrows its field sets a guard bit: `_PackOverflow`, and the entry point
# reruns with wider words (UInt64, UInt128, then BigInt).
#
# Coefficients are integers (`SymNumber`), or, when the polynomials involve other
# symbols, integer polynomials in them (`Poly{SymNumber}`): a basis over Q(params).
# Reduction is fraction-free, with primitive parts taken after each reduction.
#
# `order = :lex` with finitely many solutions computes a grevlex basis, usually far
# cheaper, and converts it with FGLM (Faugère, Gianni, Lazard, Mora), which enumerates
# monomials in lex order and finds the linear relations among their normal forms.

struct _GCtx{W}
    n::Int
    b::Int                  # value bits per field (fields are b + 1 bits wide)
    order::Symbol
    G::W                    # the guard bits of the n exponent fields
    KG::W                   # the guard bit of the key's degree field (or G for lex)
end

function _GCtx(::Type{W}, n::Int, order::Symbol) where {W}
    b = W === BigInt ? 62 : _capacity(W) ÷ (n + 1) - 1
    b >= 1 || throw(_PackOverflow(1))
    G = zero(W)
    for i in 1:n
        G |= one(W) << (_gshift(n, b, i) + b)
    end
    KG = order === :lex ? G : one(W) << (n * (b + 1) + b)
    _GCtx{W}(n, b, order, G, KG)
end

@inline _gshift(n::Int, b::Int, i::Int) = (n - i) * (b + 1)          # field i of 1..n
@inline _gexp(c::_GCtx{W}, m::W, i::Int) where {W} = Int((m >> _gshift(c.n, c.b, i)) & ((one(W) << c.b) - one(W)))

function _gmono(c::_GCtx{W}, e) where {W}
    m = zero(W)
    for i in 1:c.n
        e[i] < (1 << c.b) || throw(_PackOverflow(e[i]))
        m |= W(e[i]) << _gshift(c.n, c.b, i)
    end
    m
end

function _gkey(c::_GCtx{W}, m::W) where {W}
    c.order === :lex && return m
    n, b = c.n, c.b
    d = sum(i -> _gexp(c, m, i), 1:n; init = 0)
    d < (1 << b) || throw(_PackOverflow(d))
    k = W(d) << (n * (b + 1))
    if c.order === :grlex
        return k | m
    end
    s = 0                                   # grevlex: field j holds e1 + ... + e(n-j)
    for i in 1:(n-1)
        s += _gexp(c, m, i)
        k |= W(s) << _gshift(n, b, n - i)
    end
    k
end

_gdeg(c::_GCtx, m) = sum(i -> _gexp(c, m, i), 1:c.n; init = 0)
@inline _gdivides(c::_GCtx, a, b) = ((b | c.G) - a) & c.G == c.G            # a | b
function _glcm(c::_GCtx{W}, a::W, b::W) where {W}
    m = zero(W)
    for i in 1:c.n
        m |= W(max(_gexp(c, a, i), _gexp(c, b, i))) << _gshift(c.n, c.b, i)
    end
    m
end

# terms: keys strictly decreasing, coefficients nonzero
struct _GP{C,W}
    k::Vector{W}
    m::Vector{W}
    c::Vector{C}
end
Base.iszero(p::_GP) = isempty(p.k)
Base.length(p::_GP) = length(p.k)

## coefficients: integers, integer polynomials in the parameters, or (for division)
## elements of a field: rationals, or expressions in the parameters kept canceled

_cgcd(a::SymNumber, b::SymNumber) = SymNumber(gcd(_numerator(a), _numerator(b)))
_cgcd(a::Poly, b::Poly) = gcd(a, b)
_cdiv(a::SymNumber, b::SymNumber) = a / b
_cdiv(a::Poly, b::Poly) = divexact(a, b)
_cunit(a::SymNumber) = abs(a) == 1
_cunit(a::Poly) = _isunit(a)
_csign(a::SymNumber) = sign(a)
_csign(a::Poly) = sign(a.coefs[1])
_cnorm(a) = a
_cnorm(a::Sym) = cancel(a)

# a p - b t g, with t the monomial (tk, tm)
function _gsub(a::C, p::_GP{C,W}, b::C, tk::W, tm::W, g::_GP{C,W}, ctx::_GCtx{W}) where {C,W}
    np, ng = length(p.k), length(g.k)
    K, M, Cs = W[], W[], C[]
    sizehint!(K, np + ng); sizehint!(M, np + ng); sizehint!(Cs, np + ng)
    sa = !isone(a)
    i = j = 1
    @inbounds while i <= np || j <= ng
        if j <= ng
            gk, gm = g.k[j] + tk, g.m[j] + tm
            (gm & ctx.G != 0 || gk & ctx.KG != 0) && throw(_PackOverflow(1 << ctx.b))
        end
        if j > ng || (i <= np && p.k[i] > gk)
            push!(K, p.k[i]); push!(M, p.m[i]); push!(Cs, sa ? _cnorm(a * p.c[i]) : p.c[i])
            i += 1
        elseif i > np || gk > p.k[i]
            push!(K, gk); push!(M, gm); push!(Cs, _cnorm(-(b * g.c[j])))
            j += 1
        else
            v = _cnorm((sa ? a * p.c[i] : p.c[i]) - b * g.c[j])
            if !iszero(v)
                push!(K, gk); push!(M, gm); push!(Cs, v)
            end
            i += 1; j += 1
        end
    end
    _GP{C,W}(K, M, Cs)
end

# p divided by the gcd of its coefficients, leading coefficient positive; and the gcd
function _gprimitive(p::_GP{C,W}) where {C,W}
    iszero(p) && return p, nothing
    g = p.c[1]
    for c in p.c
        _cunit(g) && break
        g = _cgcd(g, c)
    end
    _csign(p.c[1]) * _csign(g) < 0 && (g = -g)
    isone(g) && return p, g
    _GP{C,W}(p.k, p.m, C[_cdiv(c, g) for c in p.c]), g
end

# the normal form of p by the polynomials polys[G] (fraction-free): r = mult * p mod G,
# with every term of r irreducible (only the leading term when full = false)
function _gnf(p::_GP{C,W}, polys, G, ctx::_GCtx{W}; full::Bool = true) where {C,W}
    iszero(p) && return p, nothing
    mult = one(p.c[1])
    i = 1
    while i <= length(p.k)
        mi = p.m[i]
        j = 0
        for g in G
            _gdivides(ctx, polys[g].m[1], mi) && (j = g; break)
        end
        if j == 0
            full || break
            i += 1
            continue
        end
        q = polys[j]
        c, d = p.c[i], q.c[1]
        gg = _cgcd(c, d)
        a, b = _cdiv(d, gg), _cdiv(c, gg)
        p = _gsub(a, p, b, p.k[i] - q.k[1], mi - q.m[1], q, ctx)
        mult = mult * a
    end
    p, mult
end

struct _GPair{W}
    i::Int
    j::Int
    lm::W                   # the lcm of the leading monomials, and its key
    lk::W
    sugar::Int
end

# Gebauer–Möller: add polys[h] to the basis G and its pairs to B (Becker–Weispfenning,
# UPDATE)
function _gupdate!(G::Vector{Int}, B::Vector{_GPair{W}}, h::Int, polys, sug, ctx::_GCtx{W}) where {W}
    mh = polys[h].m[1]
    lcmof(g) = _glcm(ctx, mh, polys[g].m[1])
    coprime(g) = lcmof(g) == mh + polys[g].m[1]
    C = copy(G)
    D = Int[]
    while !isempty(C)
        g1 = popfirst!(C)
        l1 = lcmof(g1)
        if coprime(g1) || !(any(g2 -> _gdivides(ctx, lcmof(g2), l1), C) || any(g2 -> _gdivides(ctx, lcmof(g2), l1), D))
            push!(D, g1)
        end
    end
    filter!(B) do pr
        !(_gdivides(ctx, mh, pr.lm) && _glcm(ctx, polys[pr.i].m[1], mh) != pr.lm &&
          _glcm(ctx, mh, polys[pr.j].m[1]) != pr.lm)
    end
    for g in D
        coprime(g) && continue
        l = lcmof(g)
        s = max(sug[h] - _gdeg(ctx, mh), sug[g] - _gdeg(ctx, polys[g].m[1])) + _gdeg(ctx, l)
        push!(B, _GPair{W}(g, h, l, _gkey(ctx, l), s))
    end
    filter!(g -> !_gdivides(ctx, mh, polys[g].m[1]), G)
    push!(G, h)
    nothing
end

function _gspoly(f::_GP{C,W}, g::_GP{C,W}, pr::_GPair{W}, ctx) where {C,W}
    gg = _cgcd(f.c[1], g.c[1])
    a, b = _cdiv(g.c[1], gg), _cdiv(f.c[1], gg)
    tf, tg = pr.lm - f.m[1], pr.lm - g.m[1]
    af = _GP{C,W}(f.k .+ (pr.lk - f.k[1]), f.m .+ tf, C[a * c for c in f.c])
    _gsub(one(a), af, b, pr.lk - g.k[1], tg, g, ctx)
end

# a reduced Gröbner basis, primitive, ordered by decreasing leading monomial
function _buchberger(F::Vector{_GP{C,W}}, ctx::_GCtx{W}) where {C,W}
    polys = _GP{C,W}[]
    sug = Int[]
    G = Int[]
    B = _GPair{W}[]
    for f in sort(F; by = f -> f.k[1])
        h, _ = _gnf(f, polys, G, ctx)
        iszero(h) && continue
        h, _ = _gprimitive(h)
        push!(polys, h)
        push!(sug, maximum(m -> _gdeg(ctx, m), h.m))
        _gupdate!(G, B, length(polys), polys, sug, ctx)
    end
    while !isempty(B)
        best = 1
        for t in 2:length(B)
            (B[t].sugar, B[t].lk) < (B[best].sugar, B[best].lk) && (best = t)
        end
        pr = B[best]
        deleteat!(B, best)
        s = _gspoly(polys[pr.i], polys[pr.j], pr, ctx)
        h, _ = _gnf(s, polys, G, ctx)
        iszero(h) && continue
        h, _ = _gprimitive(h)
        push!(polys, h)
        push!(sug, pr.sugar)
        if h.k[1] == _gkey(ctx, zero(W))              # a constant: the whole ring
            return [h]
        end
        _gupdate!(G, B, length(polys), polys, sug, ctx)
    end
    _greduce(polys[G], ctx)
end

# minimal, then each tail reduced by the others
function _greduce(F::Vector{_GP{C,W}}, ctx) where {C,W}
    F = sort(F; by = f -> f.k[1])
    keep = _GP{C,W}[]
    for f in F
        any(g -> _gdivides(ctx, g.m[1], f.m[1]), keep) || push!(keep, f)
    end
    out = similar(keep)
    for i in eachindex(keep)
        others = [j for j in eachindex(keep) if j != i]
        r, _ = _gnf(keep[i], keep, others, ctx)
        out[i] = first(_gprimitive(r))
    end
    sort!(out; by = f -> f.k[1], rev = true)
end

## FGLM: a zero-dimensional basis to lex

_gzerodim(G, ctx) = all(i -> any(g -> (m = g.m[1]; _gexp(ctx, m, i) > 0 &&
    m == typeof(m)(_gexp(ctx, m, i)) << _gshift(ctx.n, ctx.b, i)), G), 1:ctx.n)

_q(c::SymNumber) = Rational{BigInt}(big(_numerator(c)), big(_denominator(c)))

# FGLM with multiplication matrices: the grevlex standard monomials s_1, ..., s_D are a
# basis of the quotient; column s of M_i is the normal form of x_i s. A lex candidate
# u = x_i t, t already standard in lex, has the normal form M_i nf(t): a matrix-vector
# product instead of a reduction, which for the high-degree lex monomials is costly.
function _fglm(G::Vector{_GP{SymNumber,W}}, ctx::_GCtx{W}, lex::_GCtx{W}) where {W}
    n = ctx.n
    Q = Rational{BigInt}
    unit(i) = one(W) << _gshift(n, ctx.b, i)
    # the standard monomials of G
    std = W[zero(W)]
    index = Dict{W,Int}(zero(W) => 1)
    k = 1
    while k <= length(std)
        for i in 1:n
            u = std[k] + unit(i)
            u & ctx.G != 0 && throw(_PackOverflow(1 << ctx.b))
            haskey(index, u) && continue
            any(g -> _gdivides(ctx, g.m[1], u), G) && continue
            push!(std, u)
            index[u] = length(std)
        end
        k += 1
    end
    D = length(std)
    idx = collect(eachindex(G))
    function nfvec(m::W)
        v = zeros(Q, D)
        haskey(index, m) && (v[index[m]] = 1; return v)
        r, mult = _gnf(_GP{SymNumber,W}(W[_gkey(ctx, m)], W[m], SymNumber[SymNumber(1)]), G, idx, ctx)
        iszero(r) && return v
        μ = _q(mult)
        for t in eachindex(r.m)
            v[index[r.m[t]]] = _q(r.c[t]) / μ
        end
        v
    end
    cols = [Vector{Union{Nothing,Vector{Q}}}(nothing, D) for _ in 1:n]
    col(i, s) = (c = cols[i][s]; c === nothing ? (cols[i][s] = nfvec(std[s] + unit(i))) : c)
    function mulvar(i, v)
        w = zeros(Q, D)
        for s in 1:D
            iszero(v[s]) && continue
            c = col(i, s)
            for q in 1:D
                iszero(c[q]) || (w[q] += v[s] * c[q])
            end
        end
        w
    end

    stair = W[]                                 # lex standard monomials, increasing
    vecs = Vector{Q}[]                          # their normal forms
    rows = Tuple{Int,Vector{Q},Vector{Q}}[]     # echelon: pivot, vector, combination of vecs
    out = _GP{SymNumber,W}[]
    cand = Tuple{W,Int,Int}[(zero(W), 0, 0)]    # monomial, from stair[t] times x_i
    seen = Set{W}([zero(W)])
    while !isempty(cand)
        c = argmin(first, cand)                 # lex key = the packed monomial
        deleteat!(cand, findfirst(==(c), cand))
        u, t, i = c
        any(g -> _gdivides(lex, g.m[1], u), out) && continue
        v0 = t == 0 ? nfvec(u) : mulvar(i, vecs[t])
        v = copy(v0)
        α = zeros(Q, length(stair) + 1)
        α[end] = 1
        for (piv, rv, rα) in rows
            iszero(v[piv]) && continue
            f = v[piv] / rv[piv]
            for q in 1:D
                iszero(rv[q]) || (v[q] -= f * rv[q])
            end
            for q in eachindex(rα)
                α[q] -= f * rα[q]
            end
        end
        if all(iszero, v)
            # u + Σ α_q stair_q = 0
            ms = [u; stair]
            cs = [α[end]; α[1:end-1]]
            terms = sort([(ms[q], cs[q]) for q in findall(!iszero, cs)]; by = first, rev = true)
            den = reduce(lcm, (denominator(c) for (_, c) in terms); init = big(1))
            push!(out, _GP{SymNumber,W}(W[mm for (mm, _) in terms], W[mm for (mm, _) in terms],
                                        SymNumber[SymNumber(numerator(c * den)) for (_, c) in terms]))
        else
            push!(stair, u)
            push!(vecs, v0)
            push!(rows, (findfirst(!iszero, v), v, α))
            for i2 in 1:n
                w = u + unit(i2)
                w & ctx.G != 0 && throw(_PackOverflow(1 << ctx.b))
                w in seen || (push!(seen, w); push!(cand, (w, length(stair), i2)))
            end
        end
    end
    sort!([first(_gprimitive(f)) for f in out]; by = f -> f.k[1], rev = true)
end

## conversion

# the polynomials in vars with coefficients in params (integer, or integer polynomials)
function _gconvert(exs::Vector{Sym}, vars::Vector{Sym}, params::Vector{Sym}, ctx::_GCtx{W}) where {W}
    nv = length(vars)
    gs = [vars; params]
    if isempty(params)
        F = _GP{SymNumber,W}[]
        for ex in exs
            P = Poly{SymNumber}(ex, gs)
            iszero(P) && continue
            l = reduce(lcm, (big(_denominator(c)) for c in P.coefs); init = big(1))
            ms = W[_gmono(ctx, [_exponent(m, P.layout, nv, i) for i in 1:nv]) for m in P.monos]
            cs = SymNumber[c * SymNumber(l) for c in P.coefs]
            push!(F, _gpsorted(ms, cs, ctx))
        end
        return F
    end
    F = _GP{Poly{SymNumber},W}[]
    np = length(params)
    for ex in exs
        P = Poly{SymNumber}(ex, gs)
        iszero(P) && continue
        l = reduce(lcm, (big(_denominator(c)) for c in P.coefs); init = big(1))
        groups = Dict{W,Tuple{Vector{Vector{Int}},Vector{SymNumber}}}()
        for (m, c) in zip(P.monos, P.coefs)
            e = [_exponent(m, P.layout, nv + np, i) for i in 1:(nv+np)]
            gm = _gmono(ctx, e[1:nv])
            es, cs = get!(() -> (Vector{Int}[], SymNumber[]), groups, gm)
            push!(es, e[(nv+1):end])
            push!(cs, c * SymNumber(l))
        end
        ms = collect(keys(groups))
        cs = Poly{SymNumber}[_gparampoly(params, groups[m]...) for m in ms]
        push!(F, _gpsorted(ms, cs, ctx))
    end
    F
end

function _gparampoly(params::Vector{Sym}, es::Vector{Vector{Int}}, cs::Vector{SymNumber})
    np = length(params)
    L = _layout_for([maximum(e[i] for e in es) for i in 1:np])
    M = _monotype(L)
    ms = M[reduce(|, (M(e[i]) << _shift(L, np, i) for i in 1:np); init = zero(M)) for e in es]
    _sorted_poly(params, L, ms, cs)
end

function _gpsorted(ms::Vector{W}, cs::Vector{C}, ctx::_GCtx{W}) where {C,W}
    ks = W[_gkey(ctx, m) for m in ms]
    perm = sortperm(ks; rev = true)
    _GP{C,W}(ks[perm], ms[perm], cs[perm])
end

_gsym(c::SymNumber) = Sym(c)
_gsym(c::Poly) = Sym(c)
_gsym(c::Sym) = c
_gmonosym(ctx, m, vars) = prod((vars[i]^_gexp(ctx, m, i) for i in 1:ctx.n); init = ONE)

# monic, as an expression
const _GROEBNER_ALGEBRAIC_MSG =
    "groebner: a leading coefficient vanishes for the algebraic numbers in the coefficients, " *
    "which are taken as parameters; replace one (sqrt(3)) by a variable r, last, and add r^2 - 3"
function _gtosym(p::_GP, ctx, vars; monic::Bool = true)
    lc = _gsym(p.c[1])
    monic && iszero(lc) && throw(ArgumentError(_GROEBNER_ALGEBRAIC_MSG))
    _add_many(_cnorm(monic ? _gsym(c) / lc : _gsym(c)) * _gmonosym(ctx, m, vars) for (m, c) in zip(p.m, p.c))
end

# expressions, variables (given or all the symbols) and parameters
function _gsetup(F, vars::Vector{Sym}, who::String)
    exs = Sym[]
    for f in F
        f isa Relation && (f.op === :eq || throw(ArgumentError("$who: equations only, not $f")); f = f.lhs - f.rhs)
        n, d = numer_denom(Sym(f))
        any(v -> has(d, v), vars) && throw(ArgumentError("$who: $f is not a polynomial in $(join(vars, ", "))"))
        push!(exs, n)
    end
    atoms = Set{Sym}()
    for ex in exs
        _poly_atoms!(atoms, ex) || throw(ArgumentError("$who: floating point coefficients are not supported; use exact numbers"))
    end
    if isempty(vars)
        vars = _sort_by_name!([a for a in atoms if _isnamed(a)])
        isempty(vars) && return exs, vars, Sym[]
    end
    all(_isnamed, vars) || throw(ArgumentError("$who: the variables must be symbols"))
    params = Sym[]
    for a in atoms
        any(v -> _eq(a, v), vars) && continue
        any(v -> has(a, v), vars) && throw(ArgumentError("$who: not a polynomial in $(join(vars, ", ")): $a"))
        push!(params, a)
    end
    sort!(params; lt = _lt)
    exs, vars, params
end

const _GORDERS = (:lex, :grlex, :grevlex)

# run f(W) with the narrowest word that holds the computation
function _gwiden(f)
    for W in (UInt64, UInt128, BigInt)
        try
            return f(W)
        catch err
            err isa _PackOverflow || rethrow()
            W === BigInt && rethrow()
        end
    end
end

"""
    groebner(polys, vars...; order = :lex)

The reduced Gröbner basis of the ideal generated by `polys` (expressions, or equations
`lhs ~ rhs`, taken as `lhs - rhs`), in the variables `vars` (most significant first; all
the symbols, in alphabetical order, if none are given), as a vector of expressions with
leading coefficient 1, by decreasing leading monomial.

`order` is `:lex`, `:grlex` or `:grevlex`. A lex basis is triangular: its last
polynomials involve only the last variables, which is how it eliminates variables and
solves systems. When the system has finitely many solutions, a lex basis is computed as
a grevlex one, usually much cheaper, then converted (FGLM).

Other symbols are parameters: the basis is over rational functions in them, right for
generic values of the parameters, which may fail at special values.

The basis is `[1]` exactly when the polynomials have no common zero, over the complex
numbers.

```julia
groebner([x^2 + y^2 - 1, x - y], x, y)       # [x - y, y^2 - 1/2]
groebner([x*y - 1, x^2 + y^2 - 4], x, y; order = :grevlex)
groebner([x^2 + y^2 - r^2, x - y], x, y)     # [x - y, y^2 - r^2/2]
```

See also [`reduced`](@ref).
"""
function groebner(F::AbstractVector, vars::Sym...; order::Symbol = :lex)
    order in _GORDERS || throw(ArgumentError("groebner: the order must be :lex, :grlex or :grevlex, not :$order"))
    exs, vs, params = _gsetup(F, Sym[vars...], "groebner")
    exs = filter(!iszero, exs)
    (isempty(exs) || isempty(vs)) && return isempty(exs) ? Sym[] : Sym[ONE]
    _groebner(exs, vs, params, order)
end
groebner(F::Tuple, vars::Sym...; kws...) = groebner(collect(F), vars...; kws...)

# fglm = false computes a lex basis directly (for testing the conversion)
function _groebner(exs::Vector{Sym}, vs::Vector{Sym}, params::Vector{Sym}, order::Symbol; fglm::Bool = true)
    _gwiden() do W
        ctx = _GCtx(W, length(vs), order)
        if fglm && order === :lex && isempty(params)
            ctxr = _GCtx(W, length(vs), :grevlex)
            G = _buchberger(_gconvert(exs, vs, params, ctxr), ctxr)
            if length(G) == 1 && G[1].m[1] == zero(W)
                return Sym[ONE]
            end
            _gzerodim(G, ctxr) && return Sym[_gtosym(g, ctx, vs) for g in _fglm(G, ctxr, ctx)]
        end
        G = _buchberger(_gconvert(exs, vs, params, ctx), ctx)
        Sym[_gtosym(g, ctx, vs) for g in G]
    end
end

"""
    reduced(f, G, vars...; order = :lex) -> (quotients, remainder)

Divide `f` by the polynomials `G` in the variables `vars` and the monomial `order`:
`f = Σ qᵢ Gᵢ + r` with no term of `r` divisible by a leading term of a `Gᵢ`. When `G` is
a Gröbner basis (from [`groebner`](@ref), with the same variables and order), the
remainder is the normal form of `f`, zero exactly when `f` is in the ideal.

```julia
G = groebner([x^2 + y^2 - 1, x - y], x, y)
reduced(x^3 + y, G, x, y)          # ([x^2 + x*y + y^2 + ..., ...], 3y/2) …
```
"""
function reduced(f, G::AbstractVector, vars::Sym...; order::Symbol = :lex)
    order in _GORDERS || throw(ArgumentError("reduced: the order must be :lex, :grlex or :grevlex, not :$order"))
    exs = Sym[g isa Relation ? g.lhs - g.rhs : Sym(g) for g in [f; collect(G)]]
    _, vs, params = _gsetup(exs, Sym[vars...], "reduced")
    isempty(vs) && throw(ArgumentError("reduced: no variables"))
    _gwiden() do W
        ctx = _GCtx(W, length(vs), order)
        # field coefficients: rationals, or expressions in the parameters kept canceled
        P = [_gfield(e, vs, ctx, isempty(params)) for e in exs]
        p, gs = P[1], P[2:end]
        C = eltype(p.c)
        qs = [_GP{C,W}(W[], W[], C[]) for _ in gs]
        rk, rm, rc = W[], W[], C[]
        while !iszero(p)
            j = findfirst(g -> !iszero(g) && _gdivides(ctx, g.m[1], p.m[1]), gs)
            if j === nothing
                push!(rk, p.k[1]); push!(rm, p.m[1]); push!(rc, p.c[1])
                p = _GP{C,W}(p.k[2:end], p.m[2:end], p.c[2:end])
                continue
            end
            g = gs[j]
            c = _cnorm(p.c[1] / g.c[1])
            tk, tm = p.k[1] - g.k[1], p.m[1] - g.m[1]
            push!(qs[j].k, tk); push!(qs[j].m, tm); push!(qs[j].c, c)
            p = _gsub(one(c), p, c, tk, tm, g, ctx)
        end
        tos(q) = iszero(q) ? ZERO : _gtosym(q, ctx, vs; monic = false)
        Sym[tos(q) for q in qs], tos(_GP{C,W}(rk, rm, rc))
    end
end

# e as a polynomial in vs over a field: rational coefficients, or expressions in the
# parameters
function _gfield(e::Sym, vs::Vector{Sym}, ctx::_GCtx{W}, numeric::Bool) where {W}
    nv = length(vs)
    if numeric
        P = Poly{SymNumber}(e, vs)
        return _gpsorted(W[_gmono(ctx, [_exponent(m, P.layout, nv, i) for i in 1:nv]) for m in P.monos], copy(P.coefs), ctx)
    end
    P = Poly{Sym}(e, vs)
    cs = Sym[cancel(c) for c in P.coefs]
    keep = findall(!iszero, cs)
    _gpsorted(W[_gmono(ctx, [_exponent(m, P.layout, nv, i) for i in 1:nv]) for m in P.monos[keep]], cs[keep], ctx)
end
