# radsimp: denominators free of radicals

"""
    SymJ.radsimp(ex)

`ex` as a quotient whose denominator has no radicals (and no `im`, which is `(-1)^(1/2)`):
the numerator and denominator are multiplied by a cofactor of the denominator, one group
of radicals of one base at a time. Radicals inside function arguments are left alone, and
so is an expression without radicals in its denominator. Unlike [`simplify`](@ref), the
result is not chosen for its size.

```julia
SymJ.radsimp(1 / (1 + sqrt(Sym(2))))                  # sqrt(2) - 1
SymJ.radsimp(1 / (2^(Sym(1)//3) + 2^(Sym(1)//4)))     # no radicals below
SymJ.radsimp((1 + sqrt(Sym(3))*im) / (1 + (2 + sqrt(Sym(3)))*im))
SymJ.radsimp(1 / (sqrt(x) + 1))                       # (1 - √x)/(1 - x)
```
"""
radsimp(ex) = _radsimp(Sym(ex), _RADSIMP_MAXL)
# (maxL: the largest group of radicals of one base, b^(1/L) for L ≤ maxL, eliminated)
function _radsimp(e::Sym, maxL::Int)
    n, d = numer_denom(e)
    d = expand(d)
    _rad_groups(d) === nothing && return e
    for _ in 1:12                                      # a group of radicals at a time
        grp = _rad_groups(d)
        grp === nothing && break
        r = _rad_eliminate(d, grp..., maxL)
        r === nothing && return e                      # (a group too large: as given)
        g, d = r
        n = _rad_reduce(expand(n * g), grp...)
    end
    _rad_groups(d) === nothing || return e
    _rad_normalize(n, d)
end

# n/d with d's sign and integer content into n: d a number, n/d expanded; else the gcd of
# the integer contents divided out and d's leading coefficient made positive
function _rad_normalize(n::Sym, d::Sym)
    d.kind == NUMBER && return expand(n / d)
    cn, cd = _int_content(n), _int_content(d)
    g = cn === nothing || cd === nothing ? N1 : gcd(cn, cd)
    s = _leading_sign(d)
    k = Sym(SymNumber(s * g))
    isone(k) ? n / d : expand(n / k) / expand(d / k)
end
_int_content(e::Sym) = begin
    cs = e.kind == ADD ? [e.coefs..., (iszero(e.num) ? () : (e.num,))...] :
         e.kind == MUL ? [e.num] : e.kind == NUMBER ? [e.num] : [N1]
    all(c -> isexact(c) && _isint(c), cs) ? reduce(gcd, (big(_numerator(c)) for c in cs)) : nothing
end
_leading_sign(e::Sym) = (c = e.kind == ADD ? e.coefs[1] : e.kind == MUL || e.kind == NUMBER ? e.num : N1;
                         (try c < N0 catch; false end) ? -1 : 1)

# the radicals b^(p/q) of e outside function arguments (and im): one base b whose radicals
# lie inside no other radical's base, with L the lcm of their q; nothing if none
function _rad_groups(e::Sym)
    rads = Tuple{Sym,Sym}[]                            # (base, exponent)
    function find(t::Sym)
        t.kind == FN && return
        if t.kind == POW
            b, p = t.args
            if p.kind == NUMBER && isexact(p.num) && !isinteger(p) && b !== E
                push!(rads, (b, p))
            else
                find(b)
            end
        elseif t.kind == MUL
            for (b, p) in zip(t.args, t.exps)
                if p.kind == NUMBER && isexact(p.num) && !isinteger(p) && b !== E
                    push!(rads, (b, p))
                elseif p.kind == NUMBER
                    find(b)
                end
            end
        elseif t.kind == ADD
            foreach(find, t.args)
        end
    end
    find(e)
    isempty(rads) && return nothing
    bases = Sym[]
    for (b, _) in rads
        any(c -> _eq(c, b), bases) || push!(bases, b)
    end
    # outermost: a base that occurs in no other base
    inner(b) = any(c -> !_eq(c, b) && has(c, b), bases)
    k = findfirst(b -> !inner(b), bases)
    k === nothing && return nothing
    b = bases[k]
    L = reduce(lcm, (_denominator(p.num) for (c, p) in rads if _eq(c, b)))
    (b, Int(L))
end

# the radicals of base b in e as powers of t = b^(1/L): b^(p/q) = b^⌊k/L⌋ t^(k mod L), k = pL/q
function _rad_sub(e::Sym, b::Sym, L::Int, t::Sym)
    e.kind == FN && return e
    pw(base, p) = begin
        k = p.num * L
        kk = _toint(k)
        (fld(kk, L) == 0 ? ONE : base^Sym(fld(kk, L))) * t^Sym(mod(kk, L))
    end
    isrel(base, p) = _eq(base, b) && p.kind == NUMBER && isexact(p.num) && !isinteger(p) && isinteger(SymNumber(L) * p.num)
    if e.kind == POW
        base, p = e.args
        isrel(base, p) && return pw(base, p)
        return p.kind == NUMBER && isinteger(p) ? _rad_sub(base, b, L, t)^p : e
    elseif e.kind == MUL
        return Sym(e.num) * prod((isrel(base, p) ? pw(base, p) :
                                  p.kind == NUMBER && isinteger(p) ? _rad_sub(base, b, L, t)^p : base^p
                                  for (base, p) in zip(e.args, e.exps)); init = ONE)
    elseif e.kind == ADD
        return _rebuild_node(e, Sym[_rad_sub(a, b, L, t) for a in e.args])
    end
    e
end

const _RADSIMP_MAXL = 16

# e with its powers of t = b^(1/L) reduced mod t^L - b: b^(7/2) = b³√b, so that a sum's
# terms in one radical combine; e as given if it is no polynomial in t
function _rad_reduce(e::Sym, b::Sym, L::Int)
    t = Sym(Symbol("_radsimp_t"))
    E = expand(_rad_sub(e, b, L, t))
    has(E, t) || return E
    P = try Poly{Sym}(E, [t]) catch err; err isa InterruptException && rethrow(); return e end
    out = ZERO
    for k in 0:degree(P, t)
        ck = coeff(E, t, k)
        (iszero(ck) || has(ck, t)) && (has(ck, t) ? (return e) : continue)
        out += expand(ck * b^Sym(fld(k, L))) * b^(Sym(mod(k, L)) // Sym(L))
    end
    expand(out)
end

# d with the radicals of base b gone: (g, N) with d g = N, N free of them; d as a
# polynomial in t = b^(1/L) mod t^L - b, its multiplication matrix M (column j: d t^(j-1)),
# g the first column of adj(M) and N = det(M) (d g = M adj(M) e₁ = N e₁)
function _rad_eliminate(d::Sym, b::Sym, L::Int, maxL::Int)
    t = Sym(Symbol("_radsimp_t"))
    D = expand(_rad_sub(d, b, L, t))
    P = try Poly{Sym}(D, [t]) catch err; err isa InterruptException && rethrow(); return nothing end
    deg = degree(P, t)
    cs = [expand(coeff(D, t, k)) for k in 0:deg]
    any(c -> has(c, t), cs) && return nothing
    c = fill(ZERO, L)                                    # reduced mod t^L - b
    for (k, ck) in enumerate(cs)
        iszero(ck) && continue
        m = k - 1
        c[mod(m, L) + 1] = expand(c[mod(m, L) + 1] + ck * b^Sym(fld(m, L)))
    end
    # only powers of t^s present: t^s = b^(s/L) generates the same field with L/s
    s = reduce(gcd, (k - 1 for k in 1:L if !iszero(c[k])); init = L)
    s == L && return (ONE, c[1])                         # no radical after all
    if s > 1
        c = c[1:s:L]
        L ÷= s
    end
    L > maxL && return nothing
    tb = b^(Sym(1) // Sym(L))
    if L == 2                                            # the conjugate
        g = c[1] - c[2] * tb
        return (g, expand(c[1]^2 - c[2]^2 * b))
    end
    M = Matrix{Sym}(undef, L, L)
    for j in 1:L, i in 1:L
        m = i - j                                        # t^(i-1) in t^(j-1) d: c_{i-j}, wrapped
        M[i, j] = m >= 0 ? c[m + 1] : expand(b * c[m + L + 1])
    end
    adj = Sym[(isodd(i) ? 1 : -1) * LinearAlgebra.det(M[2:L, [k for k in 1:L if k != i]]) for i in 1:L]
    Nm = expand(sum(M[1, i] * adj[i] for i in 1:L))
    iszero(Nm) && return nothing
    g = _add_many(expand(adj[i]) * tb^Sym(i - 1) for i in 1:L)
    (g, Nm)
end
