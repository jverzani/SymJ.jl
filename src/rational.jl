# Rational expressions: common denominators and cancellation of common factors.

"""
    numer_denom(ex) -> (numerator, denominator)

Write `ex` as a single fraction. Rational coefficients and factors with negative
exponents move to the denominator, and sums are brought over a common denominator.
Nothing is expanded or cancelled; see [`cancel`](@ref).

```julia
numer_denom(x/y + 1)       # (x + y, y)
numer_denom(3x/4)          # (3*x, 4)
```
"""
function numer_denom(ex::Sym)
    k = ex.kind
    if k == NUMBER
        isexact(ex.num) || return (ex, ONE)
        return Sym(SymNumber(_numerator(ex.num))), Sym(SymNumber(_denominator(ex.num)))
    elseif k == ADD
        n, d = numer_denom(Sym(ex.num))
        for i in eachindex(ex.args)
            ni, di = numer_denom(_with_coef(ex.coefs[i], ex.args[i]))
            if _eq(d, di)
                n = n + ni
            else
                n, d = n * di + ni * d, d * di
            end
        end
        return n, d
    elseif k == MUL
        n, d = numer_denom(Sym(ex.num))
        for i in eachindex(ex.args)
            ni, di = numer_denom(_mkpow(ex.args[i], ex.exps[i]))
            n, d = n * ni, d * di
        end
        return n, d
    elseif k == POW
        b, e = ex.args
        if e.kind == NUMBER && _isint(e.num)
            nb, db = numer_denom(b)
            m = _toint(e.num)
            return m > 0 ? (nb^m, db^m) : (db^(-m), nb^(-m))
        end
        _negexp(e) && return (ONE, pow(b, -e))
    end
    ex, ONE
end
numer_denom(x::Number) = numer_denom(Sym(x))

# Collect the generators for viewing ex as a polynomial: symbols, constants, function
# applications and non-integer powers. Returns false for inexact (float) coefficients.
function _poly_atoms!(s::Set{Sym}, ex::Sym)
    k = ex.kind
    k == NUMBER && return isexact(ex.num)
    if k == ADD
        isexact(ex.num) && all(isexact, ex.coefs) || return false
        return all(t -> _poly_atoms!(s, t), ex.args)
    elseif k == MUL
        isexact(ex.num) || return false
        for i in eachindex(ex.args)
            if _nonneg_int(ex.exps[i])
                _poly_atoms!(s, ex.args[i]) || return false
            else
                push!(s, _mkpow(ex.args[i], ex.exps[i]))
            end
        end
        return true
    elseif k == POW && _nonneg_int(ex.args[2])
        return _poly_atoms!(s, ex.args[1])
    end
    push!(s, ex)
    true
end

"""
    cancel(ex; gcd_method = :auto)

Write `ex` as `p/q` with `p` and `q` expanded polynomials that have no common factor.
`q` has a positive leading coefficient and the integer content is shared out, so a
polynomial input comes back expanded.

Symbols, function applications such as `sin(x)`, and non-integer powers such as
`sqrt(x)` are treated as independent variables. This is always valid, but relations
among them are not used: `(x - 1)/(sqrt(x) - 1)` is not reduced to `sqrt(x) + 1`. The
exception is square roots of rational numbers, whose squares are known: a common factor
through them is cancelled, as in `cancel((x^2 + 2sqrt(2)x + 2)/(x + sqrt(2))) ==
x + sqrt(2)`.
Expressions with floating point coefficients are returned unchanged.
`gcd_method` (`:auto`, `:dense` or `:zippel`; see [`gcd`](@ref)) chooses the gcd method.

```julia
cancel((x^2 - 1)/(x + 1))                  # x - 1
cancel(1/x + 1/x^2)                        # (x + 1)/x^2
cancel((sin(x)^2 - 1)/(sin(x) + 1))        # sin(x) - 1
```
"""
cancel(ex::Sym; gcd_method::Symbol = :auto) = _with_widening(_cancel, ex; gcd_method)
function _cancel(ex::Sym; gcd_method::Symbol = :auto)
    _check_gcd_method(gcd_method)
    n, d = numer_denom(ex)
    atoms = Set{Sym}()
    (_poly_atoms!(atoms, n) && _poly_atoms!(atoms, d)) || return ex
    # square roots of rationals last: the leading coefficient, and so the sign
    # normalization, is then decided by the other atoms
    gs = sort!(collect(atoms); lt = (a, b) -> _isrationalsqrt(a) == _isrationalsqrt(b) ? _lt(a, b) : _isrationalsqrt(b))
    P, Q = Poly{SymNumber}(n, gs), Poly{SymNumber}(d, gs)
    g = _gcd_poly(P, Q; method = gcd_method)
    P, Q = divexact(P, g), divexact(Q, g)
    r = _cancel_radicals(P, Q, gs, gcd_method)
    r === nothing || ((P, Q) = r)
    if sign(Q.coefs[1]) < 0
        P, Q = -P, -Q
    end
    if Q.monos[1] == 0                                 # constant denominator
        return Sym(_scale(P, inv(Q.coefs[1])))
    end
    Sym(P) / Sym(Q)
end
cancel(x::Number; gcd_method::Symbol = :auto) = Sym(x)

# Square roots of rational numbers are generators like any other atom, so the gcd above
# does not use sqrt(2)^2 = 2: (x^2 + 2sqrt(2)x + 2)/(x + sqrt(2)) stays as it is. So when
# the denominator Q involves such roots s_i, both sides are multiplied by Q's other
# conjugates (each s_i -> ±s_i), which makes the denominator rational; the gcd of that
# denominator and every coefficient of the numerator (grouped by its s-monomial) is a
# common factor. When only the numerator P involves them, the same with the roles
# swapped: (x + sqrt(2)y)/(x^2 - 2y^2) -> 1/(x - sqrt(2)y). The result is kept only when
# the denominator's degree goes down, so 1/(x - sqrt(2)) is not turned into
# (x + sqrt(2))/(x^2 - 2). nothing if there is no gain.
function _cancel_radicals(P::Poly{SymNumber}, Q::Poly{SymNumber}, gs::Vector{Sym}, gcd_method::Symbol)
    rad = [i for i in eachindex(gs) if _isrationalsqrt(gs[i])]
    isempty(rad) && return nothing
    inQ = [i for i in rad if _maxdegs(Q)[i] > 0]
    inP = [i for i in rad if _maxdegs(P)[i] > 0]
    if !isempty(inQ)
        r = _conj_cancel(Q, P, gs, inQ, gcd_method)       # Q rationalized
        (r === nothing || r === :zero) && return nothing
        Q2, P2 = r
    elseif !isempty(inP)
        r = _conj_cancel(P, Q, gs, inP, gcd_method)       # P rationalized
        r === :zero && return zero(P), one(Q)             # P was zero: sqrt(2)^2 - 2, unexpanded
        r === nothing && return nothing
        P2, Q2 = r
    else
        return nothing
    end
    _degree_without(Q2, rad) < _degree_without(Q, rad) ? (P2, Q2) : nothing
end

# total degree, not counting the generators gs[rad]
_degree_without(p::Poly, rad) =
    maximum((sum(e for (i, e) in enumerate(ev) if !(i in rad); init = 0) for ev in exponent_vectors(p)); init = -1)

# R and S times R's other conjugates in the roots gs[rad] (so R becomes free of them),
# then divided by the gcd of the new R and every s-monomial coefficient of the new S
function _conj_cancel(R::Poly{SymNumber}, S::Poly{SymNumber}, gs::Vector{Sym}, rad::Vector{Int}, gcd_method::Symbol)
    length(rad) > 3 && return nothing                      # 2^k - 1 conjugates
    Rs = Sym(R)
    conjs = Sym[]
    for mask in 1:(2^length(rad) - 1)
        push!(conjs, subs(Rs, (gs[rad[j]] => -gs[rad[j]] for j in eachindex(rad) if (mask >> (j - 1)) & 1 == 1)...))
    end
    c = prod(conjs)
    rn, sn = expand(Rs * c), expand(Sym(S) * c)
    iszero(rn) && return :zero
    any(i -> has(rn, gs[i]), rad) && return nothing
    # the products may leave the generators (sqrt(r^2)^2 = r^2): no gain to be had then
    atoms = Set{Sym}()
    (_poly_atoms!(atoms, rn) && _poly_atoms!(atoms, sn) && all(a -> any(g -> _eq(a, g), gs), atoms)) || return nothing
    Rn, Sn = Poly{SymNumber}(rn, gs), Poly{SymNumber}(sn, gs)
    n = length(gs)
    groups = Dict{Vector{Int},Vector{Tuple{Vector{Int},SymNumber}}}()
    for (ev, cf) in zip(exponent_vectors(Sn), coefficients(Sn))
        push!(get!(groups, ev[rad], Tuple{Vector{Int},SymNumber}[]), (ev, cf))
    end
    g = Rn
    for terms in values(groups)
        A = Poly{SymNumber}(_add_many(Sym(cf) * prod((gs[i]^ev[i] for i in 1:n if !(i in rad) && ev[i] > 0); init = ONE)
                                       for (ev, cf) in terms), gs)
        g = _gcd_poly(g, A; method = gcd_method)
        g.monos[1] == 0 && return nothing
    end
    divexact(Rn, g), divexact(Sn, g)
end
_isrationalsqrt(a::Sym) = a.kind == POW && a.args[2] == HALF && a.args[1].kind == NUMBER &&
    isexact(a.args[1].num) && a.args[1].num > N0
