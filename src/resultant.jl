# Resultants and discriminants.
#
# res_x(A, B) is the determinant of the Sylvester matrix of A and B as polynomials in x
# of degrees m = deg_x A and n = deg_x B: zero exactly when A and B have a common root
# (or both leading coefficients vanish), and a polynomial in the other variables, so it
# eliminates x from A = B = 0.
#
# Two methods, over Z (rational inputs are scaled to integer ones):
#
#   * modular: the resultant modulo word-size primes, each by evaluating the other
#     variables on a grid of points, Euclid's algorithm in Z_p[x] at every point, and
#     dense interpolation; the images combined by Chinese remaindering up to a bound on
#     the coefficients. Evaluation and reduction mod p commute with the determinant as
#     long as the degrees m and n are kept formal, which `_ures_formal` does, so there
#     are no unlucky primes or points.
#   * subresultant: the subresultant PRS (Cohen, Algorithm 3.3.7) on sparse polynomials,
#     for many variables of low degree, where the dense grid is far larger than the
#     resultant (two quartics with symbolic coefficients: 5^10 points, 219 terms).

const _RES_METHODS = (:auto, :modular, :subresultant)
_check_res_method(m::Symbol) = m in _RES_METHODS ? m :
    throw(ArgumentError("the resultant method must be :auto, :modular or :subresultant, not :$m"))

## over Z_p: univariate, with formal degrees

# res of a and b (coefficients mod p, constant first, trimmed) as polynomials of the
# formal degrees m and n, whose leading coefficients may vanish:
#   res_{m,n}(a, b) = (-1)^(mn) res_{n,m}(b, a)
#   res_{m,n}(a, b) = (-1)^(mn) lc(b)^(m - deg r) res_{n,deg r}(b, r),  r = a mod b,
# the second when deg b = n
function _ures_formal(a::Vector{Int}, m::Int, b::Vector{Int}, n::Int, p::Int)
    acc = 1
    while true
        n == 0 && return _mulmod(acc, powermod(isempty(b) ? 0 : b[1], m, p), p)
        m == 0 && return _mulmod(acc, powermod(isempty(a) ? 0 : a[1], n, p), p)
        if _udeg(b) < n
            _udeg(a) < m && return 0                    # a zero first column
            isodd(m * n) && (acc = _submod(0, acc, p))
            a, m, b, n = b, n, a, m
        end
        r = _urem!(copy(a), b, p)
        isempty(r) && return 0
        dr = _udeg(r)
        f = powermod(b[end], m - dr, p)
        isodd(m * n) && (f = _submod(0, f, p))
        acc = _mulmod(acc, f, p)
        a, m, b, n = b, n, r, dr
    end
end

# coefficients (constant first) of the polynomial of degree ≤ D with the values v at
# 0, 1, ..., D: Newton's divided differences, then Horner into monomial form
function _uinterp!(v::AbstractVector{Int}, invs::Vector{Int}, p::Int)
    D = length(v) - 1
    c = collect(v)
    for j in 1:D, k in D:-1:j
        c[k+1] = _mulmod(_submod(c[k+1], c[k], p), invs[j], p)
    end
    q = zeros(Int, D + 1)
    q[1] = c[D+1]
    for k in (D-1):-1:0                      # q = q * (t - k) + c_k
        for i in (D-k+1):-1:2
            q[i] = _submod(q[i-1], _mulmod(k, q[i], p), p)
        end
        q[1] = _addmod(_submod(0, _mulmod(k, q[1], p), p), c[k+1], p)
    end
    v .= q
end

## over Z: modular

# the terms of an integer polynomial as x-degrees, exponents in the other variables
# (columns) and coefficients
function _res_terms(P::Poly{SymNumber}, i::Int, vs::Vector{Int})
    n = length(P.gens)
    ks = [_exponent(m, P.layout, n, i) for m in P.monos]
    E = [_exponent(m, P.layout, n, v) for v in vs, m in P.monos]
    cs = BigInt[_numerator(c) for c in P.coefs]
    ks, E, cs
end

# Σ_k (Σ |c| over the terms in x^k)^2: |a_k(y)| ≤ ||a_k||_1 on the unit torus, so by
# Hadamard |res(y)| ≤ SA^(n/2) SB^(m/2) there, which bounds the coefficients of res
function _res_normsq(ks::Vector{Int}, cs::Vector{BigInt}, m::Int)
    s = zeros(BigInt, m + 1)
    for (k, c) in zip(ks, cs)
        s[k+1] += abs(c)
    end
    sum(t -> t^2, s)
end

function _res_modular(A::Poly{SymNumber}, B::Poly{SymNumber}, i::Int, m::Int, n::Int)
    gs, ng = A.gens, length(A.gens)
    dA, dB = _maxdegs(A), _maxdegs(B)
    vs = [v for v in 1:ng if v != i && (dA[v] > 0 || dB[v] > 0)]
    D = [n * dA[v] + m * dB[v] for v in vs]           # deg_v res ≤ n deg_v A + m deg_v B
    dims = Tuple(D .+ 1)
    ka, Ea, ca = _res_terms(A, i, vs)
    kb, Eb, cb = _res_terms(B, i, vs)
    target = 4 * _res_normsq(ka, ca, m)^n * _res_normsq(kb, cb, n)^m   # want M^2 > target

    H = zeros(BigInt, dims)
    M = big(1)
    p = _MODPRIME_START + 1
    while M^2 <= target
        p = _prevprime(p)
        V = _res_image(ka, Ea, ca, m, kb, Eb, cb, n, D, dims, p)
        if isone(M)
            H .= _symmetric.(V, p)
        else
            w = invmod(Int(mod(M, p)), p)
            Mp = M * p
            for j in eachindex(H)
                h = H[j]
                x = h + M * _mulmod(_submod(V[j], Int(mod(h, p)), p), w, p)
                H[j] = _symmetric(mod(x, Mp), Mp)
            end
        end
        M *= p
    end

    # the result, in A's generators
    degs = zeros(Int, ng)
    for I in CartesianIndices(H)
        iszero(H[I]) && continue
        for (j, v) in enumerate(vs)
            degs[v] = max(degs[v], I[j] - 1)
        end
    end
    L = _layout_for(degs)
    W = _monotype(L)
    monos, coefs = W[], SymNumber[]
    for I in CartesianIndices(H)
        iszero(H[I]) && continue
        mono = zero(W)
        for (j, v) in enumerate(vs)
            mono |= W(I[j] - 1) << _shift(L, ng, v)
        end
        push!(monos, mono)
        push!(coefs, SymNumber(H[I]))
    end
    _sorted_poly(gs, L, monos, coefs)
end

# the resultant mod p on the grid of points, interpolated: its coefficients mod p
function _res_image(ka, Ea, ca, m, kb, Eb, cb, n, D, dims, p)
    nv = length(D)
    ra, rb = Int[mod(c, p) for c in ca], Int[mod(c, p) for c in cb]
    # pw[j][e+1, α+1] = α^e
    pw = map(1:nv) do j
        emax = max(maximum(Ea[j, :]; init = 0), maximum(Eb[j, :]; init = 0))
        T = ones(Int, emax + 1, D[j] + 1)
        for α in 0:D[j], e in 1:emax
            T[e+1, α+1] = _mulmod(T[e, α+1], α, p)
        end
        T
    end
    V = zeros(Int, dims)
    a, b = zeros(Int, m + 1), zeros(Int, n + 1)
    for I in CartesianIndices(V)
        _res_eval!(a, ka, Ea, ra, pw, I, p)
        _res_eval!(b, kb, Eb, rb, pw, I, p)
        V[I] = _ures_formal(_utrim!(copy(a)), m, _utrim!(copy(b)), n, p)
    end
    # interpolate along each axis in turn
    for j in 1:nv
        invs = Int[invmod(t, p) for t in 1:D[j]]
        pre, post = prod(dims[1:(j-1)]; init = 1), prod(dims[(j+1):end]; init = 1)
        R = reshape(V, pre, dims[j], post)
        for s in 1:pre, t in 1:post
            _uinterp!(view(R, s, :, t), invs, p)
        end
    end
    V
end

function _res_eval!(a::Vector{Int}, ks, E, cs, pw, I, p)
    fill!(a, 0)
    @inbounds for t in eachindex(ks)
        v = cs[t]
        for j in eachindex(pw)
            v = _mulmod(v, pw[j][E[j, t]+1, I[j]], p)
        end
        a[ks[t]+1] = _addmod(a[ks[t]+1], v, p)
    end
    a
end

## over Z: subresultant PRS, with x the first generator

# lc(b)^(δ+1) a mod b in x_1, δ = deg a - deg b
function _prem_full(a::Poly{C}, b::Poly{C}) where {C}
    db, lcb = _deg_in(b, 1), _lc_in(b, 1)
    δ = _deg_in(a, 1) - db
    W = eltype(b.monos)
    unit = one(W) << _shift(b.layout, length(b.gens), 1)
    r, k = a, 0
    while !iszero(r) && (dr = _deg_in(r, 1)) >= db
        r = lcb * r - _lc_in(r, 1) * _mul_term(b, unit * W(dr - db), one(C))
        k += 1
    end
    k <= δ ? r * lcb^(δ + 1 - k) : r
end

function _res_subresultant(A::Poly{SymNumber}, B::Poly{SymNumber})
    a, b = _content_in(A, 1), _content_in(B, 1)
    A, B = divexact(A, a), divexact(B, b)
    t = a^_deg_in(B, 1) * b^_deg_in(A, 1)
    s = 1
    if _deg_in(A, 1) < _deg_in(B, 1)
        A, B = B, A
        isodd(_deg_in(A, 1)) && isodd(_deg_in(B, 1)) && (s = -1)
    end
    g = h = one(A)
    while true
        dA, dB = _deg_in(A, 1), _deg_in(B, 1)
        δ = dA - dB
        isodd(dA) && isodd(dB) && (s = -s)
        R = _prem_full(A, B)
        iszero(R) && return zero(A)
        A, B = B, divexact(R, g * h^δ)
        g = _lc_in(A, 1)
        δ == 0 || (h = divexact(g^δ, h^(δ - 1)))
        _deg_in(B, 1) > 0 || break
    end
    dA = _deg_in(A, 1)
    h = divexact(_lc_in(B, 1)^dA, h^(dA - 1))
    s == 1 ? t * h : -(t * h)
end

## entry points

# p scaled to integer coefficients, and the factor
function _int_scaled(p::Poly{SymNumber})
    all(_isint, p.coefs) && return p, big(1)
    l = reduce(lcm, (big(_denominator(c)) for c in p.coefs); init = big(1))
    _scale(p, SymNumber(l)), l
end

# :auto takes the subresultant PRS when the grid of evaluation points is large against
# the inputs: sparse ones in many variables (symbolic coefficients), where the PRS was
# measured 100-500x faster; dense ones in few variables go the other way, by as much
function _res_use_subres(A::Poly, B::Poly, i::Int, m::Int, n::Int)
    dA, dB = _maxdegs(A), _maxdegs(B)
    G = prod((Float64(n * dA[v] + m * dB[v] + 1) for v in eachindex(dA) if v != i); init = 1.0)
    G > 16.0 * length(A.monos) * length(B.monos)
end

"""
    resultant(p, q, x; method = :auto)

The resultant of `p` and `q` with respect to `x`: the determinant of their Sylvester
matrix as polynomials in `x`. It is zero exactly when `p` and `q` have a common root in
`x` (or both leading coefficients vanish), and it is a polynomial in the other symbols,
so it eliminates `x` from the system `p = q = 0`.

`p` and `q` are polynomials in `x` (expressions, or [`Poly`](@ref)s, which give a
`Poly`) whose coefficients may involve other symbols or functions of them. `method` is
`:modular` (evaluation and interpolation modulo primes), `:subresultant` (the
subresultant PRS) or `:auto`, which chooses.

```julia
resultant(x^2 - 2, x^2 - 3, x)              # 1
resultant(x^2 + y^2 - 1, x - y, x)          # 2y^2 - 1
resultant(a*x^2 + b*x + c, 2a*x + b, x)     # -a*(b^2 - 4a*c)
```
"""
function resultant(p::Poly{SymNumber}, q::Poly{SymNumber}, x::Sym; method::Symbol = :auto)
    _check_res_method(method)
    _with_widening(_resultant_poly, p, q, x; method)
end
function resultant(p::Poly{Sym}, q::Poly, x::Sym; kws...)
    _poly_like(resultant(Sym(p), Sym(q), x; kws...), _unify_gens(p, q))
end
resultant(p::Poly{SymNumber}, q::Poly{Sym}, x::Sym; kws...) =
    _poly_like(resultant(Sym(p), Sym(q), x; kws...), _unify_gens(p, q))
function resultant(p::Poly, q::Poly; kws...)
    gs = _unify_gens(p, q)
    length(gs) == 1 || throw(ArgumentError("resultant: give the variable to eliminate"))
    resultant(p, q, only(gs); kws...)
end

_unify_gens(p::Poly, q::Poly) = (gs = copy(p.gens); for g in q.gens; any(h -> _eq(g, h), gs) || push!(gs, g); end; gs)
_poly_like(ex::Sym, gs::Vector{Sym}) = Poly{Sym}(ex, gs)

function _resultant_poly(p::Poly{SymNumber}, q::Poly{SymNumber}, x::Sym; method::Symbol = :auto)
    p, q = _unify(p, q)
    i = _genindex(p, x)
    (iszero(p) || iszero(q)) && return zero(p)
    m, n = _maxdegs(p)[i], _maxdegs(q)[i]
    (m == 0 || n == 0) && return n == 0 ? q^m : p^n
    (any(isfloat, p.coefs) || any(isfloat, q.coefs)) && return _res_float(p, q, i, m, n)
    A, la = _int_scaled(p)
    B, lb = _int_scaled(q)
    R = if method === :subresultant || (method === :auto && _res_use_subres(A, B, i, m, n))
        _res_subres_in(A, B, i)
    else
        _res_modular(A, B, i, m, n)
    end
    d = la^n * lb^m
    R = isone(d) ? R : _scale(R, SymNumber(1 // d))
    _unify(R, p)[1]
end

# the subresultant PRS on A and B repacked with x first, the result repacked back
function _res_subres_in(A::Poly{SymNumber}, B::Poly{SymNumber}, i::Int)
    gs = A.gens
    gx = [gs[i]; gs[[j for j in eachindex(gs) if j != i]]]
    perm = [i; [j for j in eachindex(gs) if j != i]]
    W = _monotype(A)
    L = _layout_in(W, max.(_maxdegs(A), _maxdegs(B))[perm])
    L === nothing && throw(_PackOverflow(maximum(_maxdegs(A))))
    R = _res_subresultant(_repack(A, gx, L, W), _repack(B, gx, L, W))
    dR = _maxdegs(R)
    degs = zeros(Int, length(gs))
    degs[perm] .= dR
    L2 = _layout_for(degs)
    _repack(R, gs, L2)
end

# floating point coefficients, one variable: the Sylvester determinant by elimination
# with partial pivoting
function _res_float(p::Poly{SymNumber}, q::Poly{SymNumber}, i::Int, m::Int, n::Int)
    count(>(0), max.(_maxdegs(p), _maxdegs(q))) == 1 ||
        throw(ArgumentError("resultant: floating point coefficients are supported for polynomials in one variable only"))
    num(c) = (v = N(Sym(c)); v isa Integer || v isa Rational ? float(v) : v)
    pc = zeros(Sym, m + 1)
    qc = zeros(Sym, n + 1)
    ng = length(p.gens)
    for (mo, c) in zip(p.monos, p.coefs); pc[_exponent(mo, p.layout, ng, i)+1] = Sym(c); end
    for (mo, c) in zip(q.monos, q.coefs); qc[_exponent(mo, q.layout, ng, i)+1] = Sym(c); end
    vals = [num.(pc); num.(qc)]
    T = promote_type(map(typeof, vals)...)
    S = zeros(T, m + n, m + n)
    for r in 1:n, k in 0:m
        S[r, r+m-k] = num(pc[k+1])
    end
    for r in 1:m, k in 0:n
        S[n+r, r+n-k] = num(qc[k+1])
    end
    d = one(T)
    N_ = m + n
    for c in 1:N_
        piv = argmax(abs.(S[c:N_, c])) + c - 1
        iszero(S[piv, c]) && return zero(p)
        if piv != c
            S[[c, piv], :] = S[[piv, c], :]
            d = -d
        end
        d *= S[c, c]
        for r in (c+1):N_
            f = S[r, c] / S[c, c]
            S[r, c:N_] .-= f .* S[c, c:N_]
        end
    end
    _constpoly(SymNumber, p.gens, p.layout, _tocoef(SymNumber, Sym(d)), _monotype(p))
end

function resultant(a, b, x::Sym; method::Symbol = :auto)
    _check_res_method(method)
    _isnamed(x) || throw(ArgumentError("resultant: $x is not a symbol"))
    a, b = Sym(a), Sym(b)
    atoms = Set{Sym}()
    exact = _poly_atoms!(atoms, a) && _poly_atoms!(atoms, b)
    if !exact                                # floating point coefficients
        atoms = Set{Sym}()
        _float_atoms!(atoms, a)
        _float_atoms!(atoms, b)
    end
    for t in atoms
        !_eq(t, x) && has(t, x) && throw(ArgumentError("resultant: the arguments must be polynomials in $x, not involve $t"))
    end
    push!(atoms, x)
    gs = sort!(collect(atoms); lt = _lt)
    P = Poly{SymNumber}(a, gs)
    Q = Poly{SymNumber}(b, gs)
    Sym(resultant(P, Q, x; method))
end

# the atoms of a polynomial expression with floating point coefficients
function _float_atoms!(s::Set{Sym}, ex::Sym)
    k = ex.kind
    k == NUMBER && return s
    if k == ADD
        foreach(t -> _float_atoms!(s, t), ex.args)
    elseif k == MUL
        for i in eachindex(ex.args)
            _nonneg_int(ex.exps[i]) ? _float_atoms!(s, ex.args[i]) : push!(s, _mkpow(ex.args[i], ex.exps[i]))
        end
    elseif k == POW && _nonneg_int(ex.args[2])
        _float_atoms!(s, ex.args[1])
    else
        push!(s, ex)
    end
    s
end

"""
    discriminant(p, x)

The discriminant of `p` in `x`, `(-1)^(n(n-1)/2) resultant(p, p', x) / lc(p)` for `p` of
degree `n ≥ 1` in `x`: the product of `lc(p)^2 (rᵢ - rⱼ)^2` over pairs of roots, zero
exactly when `p` has a repeated root.

```julia
discriminant(a*x^2 + b*x + c, x)    # b^2 - 4a*c
discriminant(x^3 + p*x + q, x)      # -4p^3 - 27q^2
```
"""
function discriminant(p::Poly, x::Sym; kws...)
    n = degree(p, x)
    n >= 1 || throw(ArgumentError("discriminant: degree in $x below 1"))
    r = resultant(p, diff(p, x), x; kws...)
    i = _genindex(p, x)
    lc = _lc_x(p, i)
    d = divexact(_unify(r, lc)...)
    isodd(n * (n - 1) ÷ 2) ? -d : d
end
discriminant(p::Poly; kws...) = length(p.gens) == 1 ? discriminant(p, only(p.gens); kws...) :
    throw(ArgumentError("discriminant: give the variable"))
function discriminant(ex, x::Sym; kws...)
    ex = Sym(ex)
    P = Poly{Sym}(ex, [x])
    n = degree(P, x)
    n >= 1 || throw(ArgumentError("discriminant: degree in $x below 1"))
    lc = P.coefs[1]
    r = resultant(ex, diff(ex, x), x; kws...)
    d = cancel(r / lc)
    isodd(n * (n - 1) ÷ 2) ? -d : d
end

# the coefficient of the top power of x_i, as a polynomial in the same generators
function _lc_x(p::Poly{C}, i::Int) where {C}
    ng = length(p.gens)
    dmax = _maxdegs(p)[i]
    clear = ~(eltype(p.monos)(_fieldmask(p.layout, i)) << _shift(p.layout, ng, i))
    keep = [k for k in eachindex(p.monos) if _exponent(p.monos[k], p.layout, ng, i) == dmax]
    _sorted_poly(p.gens, p.layout, p.monos[keep] .& clear, p.coefs[keep])
end
