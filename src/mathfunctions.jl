# More of Base's mathematical functions. Those with a natural symbolic form are
# rewritten into it (sind(x) = sin(π*x/180)); the rest are unevaluated functions with
# rules, as in functions.jl. Functions that need to know whether a symbol is real
# (real, imag, conj, angle) wait for assumptions.

const _DEG = PI / 180

## rewrites

for (f, g) in ((:sind, :sin), (:cosd, :cos), (:tand, :tan), (:secd, :sec), (:cscd, :csc), (:cotd, :cot))
    @eval Base.$f(x::Sym) = $g(_DEG * x)
end
for (f, g) in ((:asind, :asin), (:acosd, :acos), (:atand, :atan), (:asecd, :asec), (:acscd, :acsc), (:acotd, :acot))
    @eval Base.$f(x::Sym) = $g(x) / _DEG
end
Base.sincosd(x::Sym) = (sind(x), cosd(x))
Base.tanpi(x::Sym) = tan(PI * x)
Base.sincospi(x::Sym) = (sinpi(x), cospi(x))
Base.deg2rad(x::Sym) = _DEG * x
Base.rad2deg(x::Sym) = x / _DEG
Base.fourthroot(x::Sym) = x^(1 // 4)
Base.abs2(x::Sym) = abs(x)^2
Base.cis(x::Sym) = exp(IM * x)
Base.cispi(x::Sym) = exp(IM * PI * x)
Base.mod2pi(x::Sym) = mod(x, 2PI)
Base.numerator(x::Sym) = numer_denom(x)[1]
Base.denominator(x::Sym) = numer_denom(x)[2]

const _SymOrNum = Union{Sym,Number}
for (S, T) in ((Sym, Sym), (Sym, Number), (Number, Sym))
    @eval begin
        Base.hypot(x::$S, y::$T) = sqrt(x^2 + y^2)
        Base.:\(x::$S, y::$T) = y / x
        Base.div(x::$S, y::$T) = trunc(x / y)
        Base.fld(x::$S, y::$T) = floor(x / y)
        Base.cld(x::$S, y::$T) = ceil(x / y)
        Base.mod(x::$S, y::$T) = apply_function(:mod, Sym[Sym(x), Sym(y)])
        Base.rem(x::$S, y::$T) = apply_function(:rem, Sym[Sym(x), Sym(y)])
        Base.max(x::$S, y::$T) = apply_function(:max, Sym[Sym(x), Sym(y)])
        Base.min(x::$S, y::$T) = apply_function(:min, Sym[Sym(x), Sym(y)])
        Base.minmax(x::$S, y::$T) = (min(x, y), max(x, y))
    end
end
# muladd and fma for every mix of Sym and Number with at least one Sym (so that
# generic code such as evalpoly works)
for S in (Sym, Number), T in (Sym, Number), U in (Sym, Number)
    S === T === U === Number && continue
    @eval Base.muladd(x::$S, y::$T, z::$U) = x * y + z
    @eval Base.fma(x::$S, y::$T, z::$U) = x * y + z
end
Base.clamp(x::Sym, lo, hi) = min(max(x, lo), hi)

# hypot of one argument, and of three or more (a Sym among the first three)
Base.hypot(x::Sym) = abs(x)
for T in _symbolic_signatures(3)
    @eval Base.hypot(x::$(T[1]), y::$(T[2]), z::$(T[3]), rest::_SymOrNum...) =
        sqrt(_add_many(Sym(a)^2 for a in (x, y, z, rest...)))
end

## atan(y, x): the angle of the point (x, y), in (-π, π], for real x and y. Internally
## the function atan2, so that nothing taking atan's one argument meets it; printed as
## atan. Evaluated for constants of known sign; refine places it under assumptions.

# atan(y/x) in the quadrant of (x, y), from the signs (-1, 0, 1) of y and x
function _atan2_value(y::Sym, x::Sym, sy::Int, sx::Int)
    sx > 0 && return atan(y / x)
    sx < 0 && return atan(y / x) + (sy < 0 ? -PI : PI)
    sy == 0 ? ZERO : sy * PI / 2
end

register_function!(:atan2; julia = atan,
    eval = function (args)
        length(args) == 2 || throw(ArgumentError("atan takes one or two arguments"))
        y, x = args
        (_hassymbol(y) || _hassymbol(x)) && return nothing
        sy, sx = _csign(y), _csign(x)
        sy === nothing || sx === nothing ? nothing : _atan2_value(y, x, sy, sx)
    end,
    diff = function (args, i)
        y, x = args
        (i == 1 ? x : -y) / (x^2 + y^2)
    end)

for (S, T) in ((Sym, Sym), (Sym, Number), (Number, Sym))
    @eval begin
        Base.atan(y::$S, x::$T) = apply_function(:atan2, Sym[Sym(y), Sym(x)])
        Base.atand(y::$S, x::$T) = atan(y, x) / _DEG
    end
end

## floor, ceil, trunc and round: with a rounding mode, to `digits` after the point or
## `sigdigits` significant digits (exactly: round(π; digits = 2) is 157/50), and to an
## integer type

function _round(x::Sym, r::RoundingMode)
    r == RoundNearest && return apply_function(:round, Sym[x])
    r == RoundDown && return apply_function(:floor, Sym[x])
    r == RoundUp && return apply_function(:ceil, Sym[x])
    r == RoundToZero && return apply_function(:trunc, Sym[x])
    r == RoundNearestTiesUp && return floor(x + ONE / 2)
    if r == RoundNearestTiesAway
        x.kind == NUMBER && isexact(x.num) && return Sym(SymNumber(round(BigInt, _rat(x.num), RoundNearestTiesAway)))
        return sign(x) * floor(abs(x) + ONE / 2)
    end
    throw(ArgumentError("round: the rounding mode $r is not supported for expressions"))
end

function Base.round(x::Sym, r::RoundingMode = RoundNearest; digits::Union{Nothing,Integer} = nothing,
                    sigdigits::Union{Nothing,Integer} = nothing, base::Integer = 10)
    digits === nothing && sigdigits === nothing && return _round(x, r)
    digits === nothing || sigdigits === nothing ||
        throw(ArgumentError("round: give digits or sigdigits, not both"))
    if sigdigits !== nothing
        _hassymbol(x) && throw(ArgumentError("round: sigdigits needs a constant, not $x"))
        iszero(x) && return x
        v = N(x)
        v isa Real || throw(ArgumentError("round: $x is not real"))
        digits = sigdigits - 1 - floor(Int, log(base, abs(v)))
    end
    s = Sym(base)^digits
    _round(x * s, r) / s
end
for (f, r) in ((:floor, RoundDown), (:ceil, RoundUp), (:trunc, RoundToZero))
    @eval begin
        Base.$f(x::Sym; kw...) = round(x, $r; kw...)
        Base.$f(::Type{T}, x::Sym) where {T<:Integer} = T(round(x, $r))
    end
end
Base.round(::Type{T}, x::Sym, r::RoundingMode = RoundNearest) where {T<:Integer} = T(round(x, r))

## conversions of numbers: Int(Sym(3)), Rational(Sym(7)/2), BigFloat(PI); an
## InexactError when the value does not fit the type (a symbol, an irrational for an
## integer or a rational, a complex value for a float)

function (::Type{T})(x::Sym) where {T<:Integer}
    x.kind == NUMBER || throw(InexactError(nameof(T), T, x))
    n = x.num
    isfloat(n) && return T(n.f)
    _isint(n) || throw(InexactError(nameof(T), T, x))
    T(_numerator(n))
end
# an Int when the value fits, else a BigInt (as Integer(3.0) is an Int)
function Base.Integer(x::Sym)
    v = BigInt(x)
    typemin(Int) <= v <= typemax(Int) ? Int(v) : v
end

function (::Type{T})(x::Sym) where {T<:Rational}
    x.kind == NUMBER || throw(InexactError(:Rational, T, x))
    n = x.num
    isfloat(n) && return T(n.f)
    r = _rat(n)
    T === Rational || return T(r)
    # Rational{Int} when the value fits
    typemin(Int) <= numerator(r) <= typemax(Int) && denominator(r) <= typemax(Int) ? Rational{Int}(r) : r
end

# Float64 is in manipulation.jl; BigFloat evaluates at the current precision
function (::Type{T})(x::Sym) where {T<:AbstractFloat}
    v = T === BigFloat ? (_hassymbol(x) ? N(x) : evaluate(x, big(0.0))) : N(x)
    v isa Real || throw(InexactError(nameof(T), T, x))
    T(v)
end

## predicates: true when known to hold, from the expression or the facts of an
## enclosing @assuming block, false otherwise (as isinteger: isreal(x) is false for a
## symbol, which is complex unless assumed otherwise). A real expression is finite.

_scope_assumptions() = (A = _assumptions(nothing); A === nothing ? Assumptions() : A)

# a float infinity or NaN among e's numbers
function _hasnonfinite(e::Sym)
    nf(n) = isfloat(n) && !isfinite(n.f)
    k = e.kind
    k == NUMBER && return nf(e.num)
    (k == ADD || k == MUL) && (nf(e.num) || (k == ADD && any(nf, e.coefs))) && return true
    any(_hasnonfinite, e.args) || (k == MUL && any(_hasnonfinite, e.exps))
end

# the numeric value of a constant, or nothing
_const_value(x::Sym) = _hassymbol(x) ? nothing : (try _evalc(x) catch; nothing end)

function Base.isreal(x::Sym)
    x.kind == NUMBER && return true
    _isreal(x, _scope_assumptions()) === true
end
function Base.isfinite(x::Sym)
    x.kind == NUMBER && return !isfloat(x.num) || isfinite(x.num.f)
    _hasnonfinite(x) && return false
    v = _const_value(x)
    v === nothing || return isfinite(v)
    _isreal(x, _scope_assumptions()) === true
end
Base.isinf(x::Sym) = x.kind == NUMBER ? isfloat(x.num) && isinf(x.num.f) :
    (v = _const_value(x); v !== nothing && isinf(v))
Base.isnan(x::Sym) = x.kind == NUMBER ? isfloat(x.num) && isnan(x.num.f) :
    (v = _const_value(x); v !== nothing && isnan(v))

## sinc(x) = sin(πx)/(πx) with sinc(0) = 1, and cosc = sinc'

register_function!(:sinc; julia = sinc,
    eval = _unary(function (a)
        iszero(a) && return ONE
        n = sin(PI * a)
        _hasfn(n, :sin) || return n / (PI * a)
        _even(:sinc, a)
    end),
    diff = (args, i) -> cosc(args[1]))
register_function!(:cosc; julia = cosc,
    eval = _unary(function (a)
        iszero(a) && return ZERO
        c, s = cos(PI * a), sin(PI * a)
        (_hasfn(c, :cos) || _hasfn(s, :sin)) || return c / a - s / (PI * a^2)
        _odd(:cosc, a)
    end),
    diff = (args, i) -> -PI^2 * sinc(args[1]) - 2 * cosc(args[1]) / args[1])

Base.sinc(x::Sym) = apply_function(:sinc, Sym[x])
Base.cosc(x::Sym) = apply_function(:cosc, Sym[x])

## mod and rem: evaluated for constants (through floor and trunc), else unevaluated

for (f, g) in ((:mod, floor), (:rem, trunc))
    register_function!(f; julia = getfield(Base, f),
        eval = function (args)
            length(args) == 2 || throw(ArgumentError("$f takes two arguments"))
            a, b = args
            (!_hassymbol(a) && !_hassymbol(b) && !iszero(b)) || return nothing
            q = g(a / b)
            q.kind == NUMBER ? a - b * q : nothing
        end)
end

## max and min: n-ary, nested calls flattened, constants reduced to the extreme one,
## arguments in canonical order so that max(x, y) == max(y, x)

function _extreme_rule(name::Symbol, better)
    function (args)
        flat = Sym[]
        for a in args
            a.kind == FN && a.name === name ? append!(flat, a.args) : push!(flat, a)
        end
        best, bestv = nothing, 0.0
        rest = Sym[]
        for a in flat
            if !_hassymbol(a)
                v = a.kind == NUMBER ? a.num : (try N(a) catch; nothing end)
                (v isa Real || v isa SymNumber) || (push!(rest, a); continue)
                if best === nothing || better(v, bestv)
                    best, bestv = a, v
                end
            else
                any(r -> _eq(r, a), rest) || push!(rest, a)
            end
        end
        best === nothing || push!(rest, best)
        length(rest) == 1 && return rest[1]
        sort!(rest; lt = _lt)
        length(rest) == length(args) && all(_eq.(rest, args)) && return nothing
        _mkfn_node(name, rest)
    end
end
register_function!(:max; julia = max, eval = _extreme_rule(:max, (v, w) -> v > w))
register_function!(:min; julia = min, eval = _extreme_rule(:min, (v, w) -> v < w))

## factorial and binomial

register_function!(:factorial; julia = factorial,
    eval = _unary(function (a)
        (a.kind == NUMBER && _isint(a.num) && 0 <= sign(a.num)) || return nothing
        n = _toint(a.num)
        n <= 100_000 ? Sym(SymNumber(factorial(big(n)))) : nothing
    end))
Base.factorial(x::Sym) = apply_function(:factorial, Sym[x])

# binomial(n, k) for an integer k >= 0 is the polynomial n(n-1)...(n-k+1)/k!
register_function!(:binomial; julia = binomial,
    eval = function (args)
        length(args) == 2 || throw(ArgumentError("binomial takes two arguments"))
        n, k = args
        # binomial(n, n + m) = 0 for a positive integer m (1/Γ(1 - m) = 0), binomial(n, n) = 1
        d = expand(k - n)
        if d.kind == NUMBER && _isint(d.num) && (_hassymbol(n) || _hassymbol(k))
            d.num > N0 && return ZERO
            iszero(d.num) && return ONE
        end
        (k.kind == NUMBER && _isint(k.num)) || return nothing
        kk = _toint(k.num)
        kk < 0 && return ZERO
        n.kind == NUMBER && _isint(n.num) && return Sym(SymNumber(binomial(big(_toint(n.num)), kk)))
        kk <= 64 || return nothing
        r = ONE
        for i in 0:(kk - 1)
            r = r * (n - i)
        end
        r / Sym(SymNumber(factorial(big(kk))))
    end)
for (S, T) in ((Sym, Sym), (Sym, Integer), (Integer, Sym))
    @eval Base.binomial(n::$S, k::$T) = apply_function(:binomial, Sym[Sym(n), Sym(k)])
end
