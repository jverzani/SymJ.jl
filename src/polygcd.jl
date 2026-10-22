# Polynomial GCD over the integers (rational inputs are scaled to integer ones).
#
# Recursive primitive PRS: a polynomial is treated as univariate in its most significant
# active generator x_i, with coefficients that are polynomials in the remaining
# generators (the same `Poly` type, with field i cleared). Because generators before i
# are absent, lex order groups the terms by their power of x_i, so coefficients are
# contiguous runs. Then
#     gcd(a, b) = gcd(cont(a), cont(b)) * pp(PRS(pp(a), pp(b)))
# where cont is the gcd of the coefficients (recursion on fewer generators), pp is the
# primitive part, and each pseudo-remainder is made primitive to control coefficient
# growth.

# most significant generator present in p or q (0 if both are constants)
function _mainvar(p::Poly, q::Poly)
    all = reduce(|, p.monos; init = UInt64(0)) | reduce(|, q.monos; init = UInt64(0))
    n = length(p.gens)
    for i in 1:n
        (all >> _shift(p.bits, n, i)) & _mask(p.bits) != 0 && return i
    end
    0
end

_deg_in(p::Poly, i::Int) = isempty(p.monos) ? -1 : _exponent(p.monos[1], p.bits, length(p.gens), i)

# coefficients of p as a polynomial in x_i (generators before i absent), highest first
function _coeffs_in(p::Poly{C}, i::Int) where {C}
    n = length(p.gens)
    clear = ~(_mask(p.bits) << _shift(p.bits, n, i))
    out = Poly{C}[]
    k = 1
    while k <= length(p.monos)
        e = _exponent(p.monos[k], p.bits, n, i)
        j = k
        while j < length(p.monos) && _exponent(p.monos[j+1], p.bits, n, i) == e
            j += 1
        end
        push!(out, Poly{C}(p.gens, p.bits, p.monos[k:j] .& clear, p.coefs[k:j]))
        k = j + 1
    end
    out
end

_lc_in(p::Poly{C}, i::Int) where {C} = first(_coeffs_in(p, i))

_positive(p::Poly) = !iszero(p) && sign(p.coefs[1]) < 0 ? -p : p
_isunit(p::Poly) = length(p.monos) == 1 && p.monos[1] == 0 && abs(p.coefs[1]) == 1

function _content_in(p::Poly, i::Int)
    g = zero(p)
    for c in _coeffs_in(p, i)
        g = _gcd(g, c)
        _isunit(g) && break
    end
    g
end

# pseudo-remainder of a by b in x_i, up to a factor that is a power of lc(b)
function _prem(a::Poly{C}, b::Poly{C}, i::Int) where {C}
    db, lcb = _deg_in(b, i), _lc_in(b, i)
    unit = UInt64(1) << _shift(b.bits, length(b.gens), i)
    r = a
    while !iszero(r) && (dr = _deg_in(r, i)) >= db
        r = lcb * r - _lc_in(r, i) * _mul_term(b, unit * UInt64(dr - db), one(C))
    end
    r
end

function _primpart_in(p::Poly, i::Int)
    c = _content_in(p, i)
    _positive(_isunit(c) ? p : divexact(p, c))
end

# gcd of integer-coefficient polynomials over the same generators, positive lc
function _gcd(a::Poly{C}, b::Poly{C}) where {C}
    iszero(a) && return _positive(b)
    iszero(b) && return _positive(a)
    i = _mainvar(a, b)
    if i == 0
        g = gcd(_numerator(a.coefs[1]), _numerator(b.coefs[1]))
        return _constpoly(C, a.gens, a.bits, SymNumber(g))
    end
    _deg_in(a, i) == 0 && return _gcd(a, _content_in(b, i))   # a is free of x_i
    _deg_in(b, i) == 0 && return _gcd(_content_in(a, i), b)
    ca, cb = _content_in(a, i), _content_in(b, i)
    c = _gcd(ca, cb)
    pa = _isunit(ca) ? a : divexact(a, ca)
    pb = _isunit(cb) ? b : divexact(b, cb)
    _deg_in(pa, i) < _deg_in(pb, i) && ((pa, pb) = (pb, pa))
    while true
        r = _prem(pa, pb, i)
        iszero(r) && break
        if _deg_in(r, i) == 0                  # gcd of the primitive parts is 1
            pb = one(pb)
            break
        end
        pa, pb = pb, _primpart_in(r, i)
    end
    _positive(c * _primpart_in(pb, i))
end

function _integer_scaled(p::Poly{SymNumber})
    any(isfloat, p.coefs) && throw(ArgumentError("gcd is not defined for floating point coefficients"))
    all(_isint, p.coefs) && return p
    l = reduce(lcm, (_denominator(c) for c in p.coefs); init = 1)
    _scale(p, SymNumber(l))
end

"""
    gcd(p::Poly, q::Poly; method = :dense)

Greatest common divisor of polynomials with numeric (exact) coefficients. For integer
inputs the result is the integer gcd, including the gcd of the integer contents, with a
positive leading coefficient. If either input has non-integer rational coefficients,
the result is monic.

The default `method = :dense` is Brown's modular algorithm. `method = :zippel` adds
Zippel's sparse interpolation: much faster when the gcd has few terms relative to its
degrees in three or more variables, slower when the gcd is dense. The result is the
same either way (it is verified by division).
"""
_check_gcd_method(m::Symbol) = m in (:dense, :zippel) ? m :
    throw(ArgumentError("the gcd method must be :dense or :zippel, not :$m"))

function Base.gcd(p::Poly{SymNumber}, q::Poly{SymNumber}; method::Symbol = :dense)
    _check_gcd_method(method)
    p, q = _unify(p, q)
    g = _gcd_modular(_integer_scaled(p), _integer_scaled(q); sparse = method === :zippel)
    rational = !all(_isint, p.coefs) || !all(_isint, q.coefs)
    rational && !iszero(g) ? _scale(g, inv(g.coefs[1])) : g
end

"""
    lcm(p::Poly, q::Poly; method = :dense)

Least common multiple, normalized like [`gcd`](@ref), whose `method` it takes.
"""
function Base.lcm(p::Poly{SymNumber}, q::Poly{SymNumber}; method::Symbol = :dense)
    p, q = _unify(p, q)
    (iszero(p) || iszero(q)) && return zero(p)
    _positive(divexact(p * q, gcd(p, q; method)))
end

"""
    gcd(a::Sym, b::Sym; method = :dense)
    lcm(a::Sym, b::Sym; method = :dense)

GCD and LCM of expressions that are polynomials with numeric coefficients in their
symbols, computed through [`Poly`](@ref); `method = :zippel` as for polynomials.
"""
Base.gcd(a::Sym, b::Sym; method::Symbol = :dense) = Sym(gcd(_poly_pair(a, b)...; method))
Base.lcm(a::Sym, b::Sym; method::Symbol = :dense) = Sym(lcm(_poly_pair(a, b)...; method))

function _poly_pair(a::Sym, b::Sym)
    gs = Set{Sym}()
    (_numpoly_symbols!(gs, a) && _numpoly_symbols!(gs, b)) ||
        throw(ArgumentError("gcd requires polynomials with numeric coefficients"))
    g = sort!(collect(gs); lt = _lt)
    Poly{SymNumber}(a, g), Poly{SymNumber}(b, g)
end
