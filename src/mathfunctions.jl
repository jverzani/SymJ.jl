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
