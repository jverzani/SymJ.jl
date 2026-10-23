# Exact rational functions over the integers, ℚ(a, b, ...): a coefficient type for
# polynomials whose coefficients involve parameters (partial fractions in x over
# ℚ(a, b); later, linear systems).
#
# Canonical form, so that equality and zero tests are exact: num/den with
#   * den = 1 and num with rational coefficients, when the reduced denominator is a
#     constant;
#   * otherwise num and den integer polynomials with gcd 1 (the gcd includes the
#     integer contents) and den of positive leading coefficient.
# Values with different generator lists are combined through Poly's unification, so
# constants (no generators) mix freely with the rest.

# numerator and denominator in the 64-bit packing (concrete fields keep the arithmetic
# fast); parameters of degree beyond it are not supported
struct RatFun
    num::Poly{SymNumber,UInt64}
    den::Poly{SymNumber,UInt64}
    RatFun(n::Poly{SymNumber}, d::Poly{SymNumber}) = new(_rf64(n), _rf64(d))
end
_rf64(p::Poly{SymNumber,UInt64}) = p
function _rf64(p::Poly{SymNumber})
    try
        convert(Poly{SymNumber,UInt64}, p)
    catch err
        err isa _PackOverflow || rethrow()
        throw(ArgumentError("rational functions of parameters are limited to exponents " *
                            "that fit 64 bits ($(64 ÷ max(length(p.gens), 1)) bits each here)"))
    end
end

const _NOGENS = Sym[]
_constpoly0(c::SymNumber) = _constpoly(SymNumber, _NOGENS, _uniform(0, 64), c, UInt64)

function RatFun(n::Poly{SymNumber}, d::Poly{SymNumber}, ::Val{:normalize})
    iszero(d) && throw(DivideError())
    iszero(n) && return RatFun(_constpoly0(N0), _constpoly0(N1))
    n, d = _unify(n, d)
    if d.monos[1] == 0                                   # constant denominator
        return RatFun(_scale(n, inv(d.coefs[1])), one(d))
    end
    l = reduce(lcm, (_denominator(c) for c in Iterators.flatten((n.coefs, d.coefs))); init = 1)
    isone(l) || (n = _scale(n, SymNumber(l)); d = _scale(d, SymNumber(l)))
    g = _gcd_poly(n, d)
    if !(g.monos[1] == 0 && isone(g.coefs[1]))
        n, d = divexact(n, g), divexact(d, g)
    end
    sign(d.coefs[1]) < 0 && ((n, d) = (-n, -d))
    d.monos[1] == 0 && return RatFun(_scale(n, inv(d.coefs[1])), one(d))
    RatFun(n, d)
end
_ratfun(n, d) = RatFun(n, d, Val(:normalize))
RatFun(p::Poly{SymNumber}) = RatFun(p, one(p))
RatFun(c::SymNumber) = RatFun(_constpoly0(c), _constpoly0(N1))

_ispoly(r::RatFun) = r.den.monos[1] == 0

Base.zero(::Type{RatFun}) = RatFun(N0)
Base.one(::Type{RatFun}) = RatFun(N1)
Base.zero(::RatFun) = zero(RatFun)
Base.one(::RatFun) = one(RatFun)
Base.iszero(r::RatFun) = iszero(r.num)
Base.isone(r::RatFun) = _ispoly(r) && isone(r.num)

function Base.:+(a::RatFun, b::RatFun)
    _ispoly(a) && _ispoly(b) && return RatFun(a.num + b.num)
    _ratfun(a.num * b.den + b.num * a.den, a.den * b.den)
end
Base.:-(a::RatFun) = RatFun(-a.num, a.den)
Base.:-(a::RatFun, b::RatFun) = a + (-b)
function Base.:*(a::RatFun, b::RatFun)
    _ispoly(a) && _ispoly(b) && return RatFun(a.num * b.num)
    _ratfun(a.num * b.num, a.den * b.den)
end
Base.inv(a::RatFun) = _ratfun(a.den, a.num)
Base.:/(a::RatFun, b::RatFun) = a * inv(b)
Base.:^(a::RatFun, k::Integer) = k >= 0 ? RatFun(a.num^k, a.den^k) : inv(a)^(-k)   # reduced stays reduced
Base.:(==)(a::RatFun, b::RatFun) = iszero(a - b)
Base.hash(a::RatFun, h::UInt) = hash(_ispoly(a), h)       # == is by difference: hash coarsely

_tocoef(::Type{RatFun}, x::SymNumber) = RatFun(x)
_tocoef(::Type{RatFun}, x::Sym) = x.kind == NUMBER ? RatFun(x.num) :
    throw(ArgumentError("$x is not a constant rational function coefficient"))

Sym(r::RatFun) = _ispoly(r) ? Sym(r.num) : Sym(r.num) / Sym(r.den)
Base.show(io::IO, r::RatFun) = show(io, Sym(r))

# a multivariate Poly over [gens[i]; params] as a polynomial in gens[i] with
# coefficients in ℚ(params)
function _to_ratfun_poly(P::Poly{SymNumber}, i::Int)
    n = length(P.gens)
    params = P.gens[[j for j in 1:n if j != i]]
    pos = [j for j in 1:n if j != i]
    # the parameters' fields: wide enough for P's degrees in them
    PL = _layout_for(_maxdegs(P)[pos])
    M = _monotype(PL)
    groups = Dict{Int,Tuple{Vector{M},Vector{SymNumber}}}()
    for (m, c) in zip(P.monos, P.coefs)
        e = _exponent(m, P.layout, n, i)
        code = zero(M)
        for (k, j) in enumerate(pos)
            code |= M(_exponent(m, P.layout, n, j)) << _shift(PL, length(params), k)
        end
        ms, cs = get!(() -> (M[], SymNumber[]), groups, e)
        push!(ms, code)
        push!(cs, c)
    end
    es = sort!(collect(keys(groups)); rev = true)
    monos = UInt64[UInt64(e) for e in es]
    coefs = map(es) do e
        ms, cs = groups[e]
        perm = sortperm(ms; rev = true)
        RatFun(Poly{SymNumber}(params, PL, ms[perm], cs[perm]))
    end
    Poly{RatFun}([P.gens[i]], _uniform(1, 64), monos, coefs)
end
