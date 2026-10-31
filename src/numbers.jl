# SymNumber: the scalar used for numeric leaves and for coefficients/constant terms.
#
# Exact rationals with an unboxed `Int` fast path and a `Rational{BigInt}` fallback on
# overflow, plus `Float64` for inexact values. A tag field (rather than an abstract type)
# keeps every container concretely typed and every operation a few predictable branches.

const T_SMALL = 0x00   # n//d with d > 0, gcd(n, d) == 1, n != typemin(Int)
const T_BIG   = 0x01   # Rational{BigInt} that does not fit in T_SMALL
const T_FLOAT = 0x02   # Float64 stored in `f`

"""
    SymJ.DivisionByZero(msg)

Thrown when exact arithmetic divides by zero: `Sym(1)/0`, `Sym(0)/0`, `0^(-1)`, or a
substitution that makes a denominator zero (`subs(1/x, x => 0)`). Not exported.
"""
struct DivisionByZero <: Exception
    msg::String
end
Base.showerror(io::IO, e::DivisionByZero) = print(io, "DivisionByZero: ", e.msg)

# a division by zero, from SymJ's arithmetic or Julia's (what callers that probe for a
# pole catch)
_isdivzero(err) = err isa DivisionByZero || err isa DivideError

"""
    SymNumber <: Real

Numeric scalar used inside symbolic expressions. It is exact (an integer or rational of
arbitrary size) unless constructed from a floating point value. Small values never
allocate; arithmetic promotes to `BigInt` on overflow and demotes back when results fit.
"""
struct SymNumber <: Real
    tag::UInt8
    n::Int
    d::Int
    f::Float64
    big::Union{Nothing,Rational{BigInt}}
end

@inline _small(n::Int, d::Int) = SymNumber(T_SMALL, n, d, 0.0, nothing)

function _frombig(r::Rational{BigInt})
    n, d = numerator(r), denominator(r)
    if typemin(Int) < n <= typemax(Int) && d <= typemax(Int)
        return _small(Int(n), Int(d))
    end
    SymNumber(T_BIG, 0, 1, 0.0, r)
end

SymNumber(x::SymNumber) = x
SymNumber(x::Bool) = _small(Int(x), 1)
function SymNumber(x::Integer)
    typemin(Int) < x <= typemax(Int) ? _small(Int(x), 1) : _frombig(Rational{BigInt}(x))
end
function SymNumber(x::Rational{Int})
    n = numerator(x)
    n == typemin(Int) ? _frombig(Rational{BigInt}(x)) : _small(n, denominator(x))
end
SymNumber(x::Rational) = _frombig(Rational{BigInt}(x))
SymNumber(x::AbstractFloat) = SymNumber(T_FLOAT, 0, 1, Float64(x), nothing)

const N0 = SymNumber(0)
const N1 = SymNumber(1)
const NM1 = SymNumber(-1)

isfloat(x::SymNumber) = x.tag == T_FLOAT
isexact(x::SymNumber) = x.tag != T_FLOAT

# exact 1 / exact 0 (a float 1.0 coefficient is kept, a float 0.0 is not)
@inline _isN1(x::SymNumber) = x.tag == T_SMALL && x.n == 1 && x.d == 1
@inline _isint(x::SymNumber) =
    x.tag == T_SMALL ? x.d == 1 : x.tag == T_BIG && isone(denominator(x.big::Rational{BigInt}))

@inline function _rat(x::SymNumber)
    x.tag == T_SMALL && return Base.unsafe_rational(BigInt(x.n), BigInt(x.d))
    x.tag == T_BIG && return x.big::Rational{BigInt}
    Rational{BigInt}(x.f)
end

_numerator(x::SymNumber) = x.tag == T_SMALL ? x.n : numerator(_rat(x))
_denominator(x::SymNumber) = x.tag == T_SMALL ? x.d : denominator(_rat(x))

# integer exponent as an Int (exponents beyond Int are not supported)
function _toint(x::SymNumber)
    _isint(x) || throw(InexactError(:_toint, Int, x))
    x.tag == T_SMALL ? x.n : Int(numerator(x.big::Rational{BigInt}))
end

function Base.Float64(x::SymNumber)
    x.tag == T_FLOAT && return x.f
    x.tag == T_SMALL && return x.n / x.d
    Float64(x.big::Rational{BigInt})
end
Base.float(x::SymNumber) = Float64(x)
Base.Rational{BigInt}(x::SymNumber) = _rat(x)
Base.BigInt(x::SymNumber) = (_isint(x) || throw(InexactError(:BigInt, BigInt, x)); numerator(_rat(x)))
Base.Int(x::SymNumber) = _toint(x)

Base.iszero(x::SymNumber) = x.tag == T_SMALL ? x.n == 0 : (x.tag == T_FLOAT && iszero(x.f))
Base.isone(x::SymNumber) = x.tag == T_SMALL ? (x.n == 1 && x.d == 1) : (x.tag == T_FLOAT && isone(x.f))
Base.isinteger(x::SymNumber) = x.tag == T_FLOAT ? isinteger(x.f) : _isint(x)
Base.sign(x::SymNumber) =
    x.tag == T_SMALL ? sign(x.n) : x.tag == T_FLOAT ? Int(sign(x.f)) : sign(numerator(x.big::Rational{BigInt}))
Base.signbit(x::SymNumber) = x.tag == T_FLOAT ? signbit(x.f) : sign(x) < 0
Base.abs(x::SymNumber) = signbit(x) ? -x : x
Base.zero(::Type{SymNumber}) = N0
Base.one(::Type{SymNumber}) = N1

## Arithmetic

@inline function _add_small(n1::Int, d1::Int, n2::Int, d2::Int)
    if d1 == 1 && d2 == 1
        s, o = Base.add_with_overflow(n1, n2)
        (o || s == typemin(Int)) && return nothing
        return _small(s, 1)
    end
    g = gcd(d1, d2)
    a1, a2 = d2 ÷ g, d1 ÷ g
    x, o1 = Base.mul_with_overflow(n1, a1)
    y, o2 = Base.mul_with_overflow(n2, a2)
    s, o3 = Base.add_with_overflow(x, y)
    d, o4 = Base.mul_with_overflow(d1, a1)
    (o1 | o2 | o3 | o4 || s == typemin(Int)) && return nothing
    g2 = gcd(s, d)
    _small(s ÷ g2, d ÷ g2)
end

@inline function _mul_small(n1::Int, d1::Int, n2::Int, d2::Int)
    if d1 == 1 && d2 == 1
        p, o = Base.mul_with_overflow(n1, n2)
        (o || p == typemin(Int)) && return nothing
        return _small(p, 1)
    end
    g1, g2 = gcd(n1, d2), gcd(n2, d1)
    p, o1 = Base.mul_with_overflow(n1 ÷ g1, n2 ÷ g2)
    q, o2 = Base.mul_with_overflow(d1 ÷ g2, d2 ÷ g1)
    (o1 | o2 || p == typemin(Int)) && return nothing
    _small(p, q)
end

function Base.:+(a::SymNumber, b::SymNumber)
    if a.tag == T_SMALL && b.tag == T_SMALL
        r = _add_small(a.n, a.d, b.n, b.d)
        r === nothing || return r
    elseif a.tag == T_FLOAT || b.tag == T_FLOAT
        return SymNumber(Float64(a) + Float64(b))
    end
    _frombig(_rat(a) + _rat(b))
end

function Base.:*(a::SymNumber, b::SymNumber)
    if a.tag == T_SMALL && b.tag == T_SMALL
        r = _mul_small(a.n, a.d, b.n, b.d)
        r === nothing || return r
    elseif a.tag == T_FLOAT || b.tag == T_FLOAT
        return SymNumber(Float64(a) * Float64(b))
    end
    _frombig(_rat(a) * _rat(b))
end

function Base.:-(a::SymNumber)
    a.tag == T_SMALL && return _small(-a.n, a.d)
    a.tag == T_FLOAT && return SymNumber(-a.f)
    _frombig(-(a.big::Rational{BigInt}))
end
Base.:-(a::SymNumber, b::SymNumber) = a + (-b)

function Base.inv(a::SymNumber)
    if a.tag == T_SMALL
        a.n == 0 && throw(DivisionByZero("division by zero"))
        return a.n > 0 ? _small(a.d, a.n) : _small(-a.d, -a.n)
    end
    a.tag == T_FLOAT && return SymNumber(inv(a.f))
    _frombig(inv(a.big::Rational{BigInt}))
end
Base.:/(a::SymNumber, b::SymNumber) = a * inv(b)
Base.://(a::SymNumber, b::SymNumber) = a / b

@inline function _checked_pow(x::Int, k::Int)
    r, o = 1, false
    while true
        if isodd(k)
            r, o1 = Base.mul_with_overflow(r, x)
            o |= o1
        end
        k >>= 1
        k == 0 && break
        x, o2 = Base.mul_with_overflow(x, x)
        o |= o2
    end
    r, o
end

function Base.:^(a::SymNumber, k::Integer)
    a.tag == T_FLOAT && return SymNumber(a.f^k)
    k == 0 && return N1
    if k < 0
        iszero(a) && throw(DivisionByZero("division by zero"))
        return inv(a)^(-big(k))
    end
    if a.tag == T_SMALL && k <= typemax(Int)
        n, o1 = _checked_pow(a.n, Int(k))
        d, o2 = _checked_pow(a.d, Int(k))
        (o1 | o2 || n == typemin(Int)) || return _small(n, d)
    end
    _frombig(_rat(a)^k)
end

## Comparison and hashing (`==` is numeric; `isequal`/`hash` are structural so that
## exact 1 and float 1.0 are distinct keys)

function Base.:(==)(a::SymNumber, b::SymNumber)
    a.tag == T_SMALL && b.tag == T_SMALL && return a.n == b.n && a.d == b.d
    (a.tag == T_FLOAT || b.tag == T_FLOAT) && return Float64(a) == Float64(b)
    _rat(a) == _rat(b)
end

function Base.:<(a::SymNumber, b::SymNumber)
    a.tag == T_SMALL && b.tag == T_SMALL && return widemul(a.n, b.d) < widemul(b.n, a.d)
    (a.tag == T_FLOAT || b.tag == T_FLOAT) && return Float64(a) < Float64(b)
    _rat(a) < _rat(b)
end
Base.:<=(a::SymNumber, b::SymNumber) = a < b || a == b
Base.isless(a::SymNumber, b::SymNumber) = a < b

# structural three-way comparison used for canonical ordering
function _numcmp(a::SymNumber, b::SymNumber)
    a.tag != b.tag && return a.tag < b.tag ? -1 : 1
    if a.tag == T_SMALL
        c = cmp(a.n, b.n)
        return c != 0 ? c : cmp(a.d, b.d)
    end
    a.tag == T_BIG && return cmp(a.big::Rational{BigInt}, b.big::Rational{BigInt})
    cmp(reinterpret(UInt64, a.f), reinterpret(UInt64, b.f))
end

Base.isequal(a::SymNumber, b::SymNumber) = _numcmp(a, b) == 0

function Base.hash(x::SymNumber, h::UInt)
    x.tag == T_SMALL && return x.d == 1 ? hash(x.n, h) : hash(x.n, hash(x.d, h))
    x.tag == T_FLOAT && return hash(x.f, h + 0x01)
    hash(x.big::Rational{BigInt}, h)
end

## Conversion / promotion

Base.convert(::Type{SymNumber}, x::SymNumber) = x
Base.convert(::Type{SymNumber}, x::Real) = SymNumber(x)
Base.promote_rule(::Type{SymNumber}, ::Type{<:Integer}) = SymNumber
Base.promote_rule(::Type{SymNumber}, ::Type{<:Rational}) = SymNumber
Base.promote_rule(::Type{SymNumber}, ::Type{<:AbstractFloat}) = SymNumber

function Base.show(io::IO, x::SymNumber)
    if x.tag == T_FLOAT
        print(io, x.f)
    elseif _isint(x)
        print(io, _numerator(x))
    else
        print(io, _numerator(x), "//", _denominator(x))
    end
end

## Integer roots, used to canonicalize n^(p/q)

const _SMALL_PRIMES = let
    s = trues(1000)
    s[1] = false
    for i in 2:31, j in (i * i):i:1000
        s[i] && (s[j] = false)
    end
    findall(s)
end

# largest r with r^k <= m, m > 0
function _iroot(m::Integer, k::Int)
    k == 2 && return isqrt(m)
    mb = BigInt(m)
    r = BigInt(floor(BigFloat(mb)^(1 / BigFloat(k))))
    while (r + 1)^k <= mb
        r += 1
    end
    while r^k > mb
        r -= 1
    end
    oftype(m, r)
end

function _perfect_root(m::Integer, k::Int)
    r = _iroot(m, k)
    big(r)^k == m ? r : nothing
end

# n = a^q * c with c free of q-th powers (exact for prime factors < 1000, best effort beyond)
function _extract_power(n::Integer, q::Int)
    a, c, m = one(n), one(n), n
    for p in _SMALL_PRIMES
        p * p > m && break
        if m % p == 0
            k = 0
            while m % p == 0
                m ÷= p
                k += 1
            end
            a *= oftype(n, p)^(k ÷ q)
            c *= oftype(n, p)^(k % q)
        end
    end
    if m > 1
        r = _perfect_root(m, q)
        r === nothing ? (c *= m) : (a *= r)
    end
    a, c
end
