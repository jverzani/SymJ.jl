# Unevaluated function applications (kind FN) and the rules that simplify them.
#
# A function is identified by its name. Behaviour is looked up in `FUNCTION_RULES`:
#   julia: a Julia function for numeric evaluation (and TermInterface.operation)
#   eval:  args -> Union{Nothing,Sym}; return a simplified value or `nothing`
#   diff:  (args, i) -> partial derivative with respect to the i-th argument
# Names with no rule (user functions) stay unevaluated.

"""
    SymJ.SymFunction(name)

A symbolic function. Calling it builds an unevaluated application:

```julia
f = SymFunction(:f)     # or `@syms f()`
f(x, 2)                 # f(x, 2)
```

A function is its name: one named like a built-in function (`SymFunction(:gamma)`,
`@syms sin()`) is that function, with its values and derivatives (`SymFunction(:gamma)(1)`
is 1).
"""
struct SymFunction <: AbstractSymbolic
    name::Symbol
end
(f::SymFunction)(args...) = apply_function(f.name, Sym[Sym(a) for a in args])
Base.show(io::IO, f::SymFunction) = print(io, f.name)
Base.show(io::IO, ::MIME"text/plain", f::SymFunction) = show(io, f)

# f' (f'', …): the derivative of a function of one argument; f'(x) is diff(f(x), x) and
# f'(0) the derivative at 0 (for initial conditions, `y'(0) => 1`)
struct SymDerivative <: AbstractSymbolic
    f::SymFunction
    k::Int
end
Base.adjoint(f::SymFunction) = SymDerivative(f, 1)
Base.adjoint(d::SymDerivative) = SymDerivative(d.f, d.k + 1)
function (d::SymDerivative)(a)
    a = Sym(a)
    t = a.kind == SYMBOL ? a : Sym(:x)
    e = d.f(t)
    for _ in 1:d.k
        e = _mkfn_node(:Derivative, Sym[e, t])
    end
    a.kind == SYMBOL ? e : apply_function(:Derivative, Sym[e.args[1], t, a])
end
Base.show(io::IO, d::SymDerivative) = print(io, d.f.name, "'"^d.k)
Base.show(io::IO, ::MIME"text/plain", d::SymDerivative) = show(io, d)

struct FunctionRule
    julia::Union{Nothing,Function}
    eval::Union{Nothing,Function}
    diff::Union{Nothing,Function}
end

const FUNCTION_RULES = Dict{Symbol,FunctionRule}()

"""
    SymJ.register_function!(name; julia=nothing, eval=nothing, diff=nothing)

Attach behaviour to the function `name`: a Julia function used for numeric evaluation,
a simplification rule `args -> Union{Nothing,Sym}`, and a partial-derivative rule
`(args, i) -> Union{Nothing,Sym}` (`nothing` for a derivative with no closed form,
which is left unevaluated).
"""
function register_function!(name::Symbol; julia = nothing, eval = nothing, diff = nothing)
    FUNCTION_RULES[name] = FunctionRule(julia, eval, diff)
    name
end

"""
    apply_function(name::Symbol, args::Vector{Sym})

Apply the function `name`, simplifying with its registered rule if any.
"""
function apply_function(name::Symbol, args::Vector{Sym})
    rule = get(FUNCTION_RULES, name, nothing)
    if rule !== nothing
        # Float arguments are values: evaluate first, before symbolic rules such as
        # asin(-u) = -asin(u), which do not hold for Julia's principal values on a
        # branch cut (asin(-17.0 + 0im) is not -asin(17.0 + 0im)).
        if rule.julia !== nothing && _float_args(args)
            # (behind an inference barrier: see _apply_float)
            v = Base.inferencebarrier(_apply_float)(rule, name, args)::Union{Nothing,Sym}
            v === nothing || return v
        end
        # a complex float argument (2.0im, 2.0 + 2.0im): its value too (SymPy #6071)
        if rule.julia !== nothing && length(args) == 1 && (z = _complex_float(args[1])) !== nothing
            v = try Sym(rule.julia(z)) catch err; err isa Union{DomainError,MethodError,InexactError} || rethrow(); nothing end
            v === nothing || return v
        end
        if rule.eval !== nothing
            r = rule.eval(args)
            r === nothing || return r::Sym
        end
    end
    _mkfn_node(name, args)
end

# the value of a registered function at float arguments, or nothing. Called from
# apply_function through Base.inferencebarrier: _evalc and _exact_limit_value reach the
# whole system (limits, signs, solve), and with a direct call every sin(::Sym) inference
# had to infer it all. The REPL's hint inference does that from scratch, with its own
# cache, as one types integrate(sin(x), ... : it ran for minutes, and on Julia 1.13 it
# froze Ctrl-R history search (which waits on the same thread pool).
function _apply_float(rule::FunctionRule, name::Symbol, args::Vector{Sym})
    v = _try_numeric(rule.julia, args)
    if v === nothing
        v = _try_numeric_complex(rule.julia, args)
        # no real value and a complex infinity: a pole, as for exact arguments
        # (gamma(-1.0) gave -Inf - Inf*im; SymPy #16847); at an infinite argument
        # no pole: left unevaluated, as tan(Inf) is (cos(Inf); SymPy #13499)
        if v !== nothing
            if !all(a -> isfinite(_evalc(a)), args)
                v = nothing
            elseif !isfinite(try _evalc(v) catch; 0.0 end)
                throw(DivisionByZero("$name($(join(args, ", "))): a pole"))
            end
        end
    end
    v === nothing && return nothing
    # at ±∞ (the only float arguments): the limit, exact when it is an integer
    # or a rational multiple of π (sinint(∞) = π/2, atan(∞), erf(-∞); SymPy #6602)
    if any(a -> isfloat(a.num) && isinf(_evalc(a)), args) &&
       all(a -> !isfloat(a.num) || isinf(_evalc(a)), args)
        w = _exact_limit_value(v)
        w === nothing || return w
    end
    v::Sym
end

_float_args(args) = all(a -> a.kind == NUMBER, args) && any(a -> isfloat(a.num), args)

# a + b i with a float among a, b (numbers: 2.0im is the product 2.0 (-1)^(1/2)), as a
# ComplexF64; else nothing
function _complex_float(a::Sym)
    isim(t) = t.kind == MUL && length(t.args) == 1 && t.args[1].kind == NUMBER && t.args[1].num == -1 && t.exps[1] == 1 // 2
    if isim(a)
        isfloat(a.num) || return nothing
        return ComplexF64(0, Float64(a.num))
    end
    a.kind == ADD && length(a.args) == 1 && _eq(a.args[1], IM) || return nothing
    b = a.coefs[1]
    (isfloat(a.num) || isfloat(b)) || return nothing
    ComplexF64(Float64(a.num), Float64(b))
end

# a Float64 value that is an integer, a half-integer, or pπ/q (q ≤ 12) to the last bits, as
# that exact number; else nothing (fresnels(∞) = 1/2)
function _exact_limit_value(v::Sym)
    v.kind == NUMBER && isfloat(v.num) || return nothing
    c = _evalc(v)
    (iszero(imag(c)) && isfinite(real(c))) || return nothing
    f = real(c)
    isinteger(f) && abs(f) < 2^52 && return Sym(Int(f))
    isinteger(2f) && abs(f) < 2^51 && return Sym(Int(2f) // 2)
    for q in 1:12
        p = round(f * q / π)
        (p != 0 && abs(p) <= 48 && abs(f - p * π / q) <= 4eps(abs(f))) && return Sym(Rational{Int}(Int(p), q)) * PI
    end
    nothing
end

# outside the real domain (asin(17.0)): the principal complex value
function _try_numeric_complex(f, args)
    try
        return Sym(f(map(a -> complex(Float64(a.num)), args)...))
    catch err
        err isa Union{DomainError,MethodError,InexactError} || rethrow()
        return nothing
    end
end

function _try_numeric(f, args)
    try
        return Sym(f(map(a -> Float64(a.num), args)...))
    catch err
        err isa Union{DomainError,MethodError,InexactError} || rethrow()   # e.g. factorial(2.5)
        return nothing
    end
end

function _fn_operation(name::Symbol)
    f = get(_NAMED_OPS, name, nothing)
    f === nothing || return f
    r = get(FUNCTION_RULES, name, nothing)
    r === nothing || r.julia === nothing ? SymFunction(name) : r.julia
end

## helpers for rules

# Would `-x` print more simply than `x`? Exactly one of `x`, `-x` answers true (or both
# false), which makes it usable to canonicalize odd and even functions. A sum is asked
# about the term it prints first (show.jl's order, which ignores signs), so sin(y - x)
# becomes -sin(x - y), not -sin(-x + y).
function _could_extract_minus(x::Sym)
    k = x.kind
    (k == NUMBER || k == MUL) && return sign(x.num) < 0
    k == ADD && return sign(x.coefs[_first_printed(x)]) < 0
    false
end

# the index of the term of a sum printed first: the highest degree, then the printed form
function _first_printed(x::Sym)
    best, deg = 1, _degree(x.args[1])
    tie = false
    for i in 2:length(x.args)
        d = _degree(x.args[i])
        d > deg ? ((best, deg, tie) = (i, d, false)) : d == deg && (tie = true)
    end
    tie || return best
    if all(a -> a.kind == SYMBOL || _degree(a) != deg, x.args)   # names compare as printed
        best = 0
        for i in eachindex(x.args)
            _degree(x.args[i]) == deg && (best == 0 || isless(x.args[i].name, x.args[best].name)) && (best = i)
        end
        return best
    end
    key = nothing
    for i in eachindex(x.args)
        _degree(x.args[i]) == deg || continue
        s = _sortstring(x.args[i])
        (key === nothing || s < key) && ((best, key) = (i, s))
    end
    best
end

# x == r*π for an exact rational r
function _pi_multiple(x::Sym)
    iszero(x) && isexact(x.num) && return Rational{BigInt}(0)
    x === PI && return Rational{BigInt}(1)
    if x.kind == MUL && length(x.args) == 1 && x.args[1] === PI && _isone(x.exps[1]) && isexact(x.num)
        return _rat(x.num)
    end
    nothing
end

# a/π a constant known to be an integer that is no number (1342^(87236487262873^5) + 4,
# kept as a power for its size: sin and tan of π times it are 0; SymPy #10607)
function _int_times_pi(a::Sym)
    (_hassymbol(a) || !has(a, PI)) && return false
    q = a / PI                                # (no expand: the power is too big for it)
    # (by form only: these rules run when patterns are built at load time, and the
    # assumptions machinery behind _isinteger made that inference explode)
    q.kind != NUMBER && !_hassymbol(q) && !has(q, PI) && _const_int_power(q)
end

# sin(rπ) when it has a simple closed form
function _sinpi(r::Rational{BigInt})
    r = mod(r, 2)
    neg = r >= 1
    neg && (r -= 1)
    r > 1 // 2 && (r = 1 - r)
    # (π/12, π/10, π/5 and their multiples too, as SymPy evaluates them; SymPy #6048)
    v = r == 0 ? ZERO :
        r == 1 // 12 ? (sqrt(Sym(6)) - sqrt(Sym(2))) / 4 :
        r == 1 // 10 ? (sqrt(Sym(5)) - 1) / 4 :
        r == 1 // 6 ? HALF :
        r == 1 // 5 ? sqrt(10 - 2sqrt(Sym(5))) / 4 :
        r == 1 // 4 ? sqrt(Sym(2)) / 2 :
        r == 3 // 10 ? (sqrt(Sym(5)) + 1) / 4 :
        r == 1 // 3 ? sqrt(Sym(3)) / 2 :
        r == 2 // 5 ? sqrt(10 + 2sqrt(Sym(5))) / 4 :
        r == 5 // 12 ? (sqrt(Sym(6)) + sqrt(Sym(2))) / 4 :
        r == 1 // 2 ? ONE : nothing
    v === nothing && return nothing
    neg ? -v : v
end

# tan(rπ) where sin/cos would leave a quotient of radicals (tan(π/12) = 2 - √3)
function _tanpi(r::Rational{BigInt})
    r = mod(r, 1)
    neg = r > 1 // 2
    neg && (r = 1 - r)
    v = r == 1 // 12 ? 2 - sqrt(Sym(3)) :
        r == 5 // 12 ? 2 + sqrt(Sym(3)) :
        r == 1 // 10 ? sqrt(25 - 10sqrt(Sym(5))) / 5 :
        r == 3 // 10 ? sqrt(25 + 10sqrt(Sym(5))) / 5 :
        r == 1 // 5 ? sqrt(5 - 2sqrt(Sym(5))) :
        r == 2 // 5 ? sqrt(5 + 2sqrt(Sym(5))) : nothing
    v === nothing && return nothing
    neg ? -v : v
end

_isfn(x::Sym, name::Symbol) = x.kind == FN && x.name === name && length(x.args) == 1

_odd(name, a) = _could_extract_minus(a) && !_on_cut(name, a) ? -apply_function(name, Sym[-a]) : nothing

# a constant on the branch cut of an odd inverse function, where f(-a) = -f(a) fails for
# the values N gives (Julia's, above the cut: asin(-2) is -π/2 + i acosh(2), not -asin(2);
# found in the review of SymPy #4004): the cut is on the real or imaginary axis, |a| > 1
# (:out) or 0 < |a| < 1 (:in)
const _ODD_CUTS = Dict(:asin => (:real, :out), :atanh => (:real, :out), :acsc => (:real, :in), :acoth => (:real, :in),
                       :asinh => (:imag, :out), :atan => (:imag, :out), :acsch => (:imag, :in), :acot => (:imag, :in))
function _on_cut(name, a::Sym)
    c = get(_ODD_CUTS, name, nothing)
    (c === nothing || a.kind == FN) && return false
    isdefined(@__MODULE__, :_evalc) || return false
    Base.inferencebarrier(_hassymbol)(a) && return false
    v = try Base.inferencebarrier(_evalc)(a)::ComplexF64 catch; return false end
    isfinite(v) || return false
    axis, side = c
    tiny = 1e-14 * max(1.0, abs(v))
    t = axis === :real ? (abs(imag(v)) <= tiny ? abs(real(v)) : return false) :
                         (abs(real(v)) <= tiny ? abs(imag(v)) : return false)
    side === :out ? t > 1 + 1e-12 : 1e-14 < t < 1 - 1e-12
end

# u for a = i u with u a real number (exact), or nothing
function _imag_number(a::Sym)
    _eq(a, IM) && return ONE
    a.kind == MUL && length(a.args) == 1 && _eq(a.args[1]^a.exps[1], IM) && isexact(a.num) || return nothing
    Sym(a.num)
end

# sinh, cosh, tanh of asinh, acosh, atanh: algebraic (principal branches; SymPy #4136)
function _hyp_of_inverse(f::Symbol, a::Sym)
    a.kind == FN && length(a.args) == 1 && a.name in (:asinh, :acosh, :atanh) || return nothing
    u, g = a.args[1], a.name
    if g === :asinh
        return f === :cosh ? sqrt(1 + u^2) : f === :tanh ? u / sqrt(1 + u^2) : nothing
    elseif g === :acosh
        return f === :sinh ? sqrt(u - 1) * sqrt(u + 1) : f === :tanh ? sqrt(u - 1) * sqrt(u + 1) / u : nothing
    end
    f === :sinh ? u / sqrt(1 - u^2) : f === :cosh ? inv(sqrt(1 - u^2)) : nothing
end
_even(name, a) = _could_extract_minus(a) ? apply_function(name, Sym[-a]) : nothing

function _unary(f)
    args -> length(args) == 1 ? f(args[1]) : throw(ArgumentError("expected one argument, got $(length(args))"))
end

## elementary functions

# u - kπ for an integer k ≠ 0 among u's terms: (u - kπ, k), or nothing (sin, cos and tan
# of a sum with such a shift: SymPy #26373, cos(2acos(3/5) - 4π))
function _pi_shift(a::Sym)
    a.kind == ADD || return nothing
    i = findfirst(t -> t === PI, a.args)
    i === nothing && return nothing
    c = a.coefs[i]
    (isexact(c) && _isint(c) && !iszero(c)) || return nothing
    a - Sym(c) * PI, _toint(c)
end

# sin, cos, tan of asin, acos, atan: algebraic (principal branches; SymPy #25685)
function _trig_of_inverse(f::Symbol, a::Sym)
    a.kind == FN && length(a.args) == 1 && a.name in (:asin, :acos, :atan) || return nothing
    u, g = a.args[1], a.name
    ((f === :sin && g === :asin) || (f === :cos && g === :acos) || (f === :tan && g === :atan)) && return u
    g === :asin && return f === :cos ? sqrt(1 - u^2) : f === :tan ? u / sqrt(1 - u^2) : nothing
    g === :acos && return f === :sin ? sqrt(1 - u^2) : f === :tan ? sqrt(1 - u^2) / u : nothing
    f === :sin ? u / sqrt(1 + u^2) : f === :cos ? inv(sqrt(1 + u^2)) : nothing
end

# u with a = i u, u free of im (a product with the factor im), or nothing: for simplify's
# sin(iu) = i sinh(u), cos(iu) = cosh(u), tan(iu) = i tanh(u) (sin(t [0 1; -1 0])
# showed im*sin(im*t); found in the review of SymPy #18631). Not an automatic rule: as
# one, compile-time evaluation of the sin and sinh rules did not finish.
function _im_multiple(a::Sym)
    a.kind == MUL || return nothing
    any(i -> _eq(a.args[i]^a.exps[i], IM), eachindex(a.args)) || return nothing      # (im is (-1)^(1/2))
    u = expand(a / IM)
    # (im inside a product is the factor (-1)^(1/2), which has() does not see: with one
    # left in u, sin and sinh rewrote each other without end)
    _has_node(u, t -> _eq(t, IM) || (t.kind == MUL && any(i -> _eq(t.args[i]^t.exps[i], IM), eachindex(t.args)))) && return nothing
    u
end

register_function!(:sin; julia = sin,
    eval = _unary(a -> begin
        r = _pi_multiple(a)
        r !== nothing && (v = _sinpi(r); v === nothing || return v)
        _int_times_pi(a) && return ZERO
        _isfn(a, :asin) && return a.args[1]
        v = _trig_of_inverse(:sin, a)
        v === nothing || return v
        sh = _pi_shift(a)
        sh === nothing || return iseven(sh[2]) ? sin(sh[1]) : -sin(sh[1])
        _odd(:sin, a)
    end),
    diff = (args, i) -> cos(args[1]))

register_function!(:cos; julia = cos,
    eval = _unary(a -> begin
        r = _pi_multiple(a)
        r !== nothing && (v = _sinpi(r + 1 // 2); v === nothing || return v)
        _isfn(a, :acos) && return a.args[1]
        v = _trig_of_inverse(:cos, a)
        v === nothing || return v
        sh = _pi_shift(a)
        sh === nothing || return iseven(sh[2]) ? cos(sh[1]) : -cos(sh[1])
        _even(:cos, a)
    end),
    diff = (args, i) -> -sin(args[1]))

register_function!(:tan; julia = tan,
    eval = _unary(a -> begin
        r = _pi_multiple(a)
        if r !== nothing
            t = _tanpi(r)
            t === nothing || return t
            s, c = _sinpi(r), _sinpi(r + 1 // 2)
            if s !== nothing && c !== nothing
                iszero(c) && throw(DivisionByZero("tan($a): a pole"))
                return s / c
            end
        end
        _int_times_pi(a) && return ZERO
        _isfn(a, :atan) && return a.args[1]
        v = _trig_of_inverse(:tan, a)
        v === nothing || return v
        sh = _pi_shift(a)
        sh === nothing || return tan(sh[1])
        _odd(:tan, a)
    end),
    diff = (args, i) -> 1 + tan(args[1])^2)

# exact values of asin and atan at the standard points (odd: the negatives follow);
# acos(v) = π/2 - asin(v)
const _ASIN_VALUES = Pair{Sym,Sym}[ZERO => ZERO, HALF => PI / 6, sqrt(Sym(2)) / 2 => PI / 4,
                                   sqrt(Sym(3)) / 2 => PI / 3, ONE => PI / 2]
# the sines of π/12, π/10, π/8, π/5 and their complements, added at the end of the module
# (built from sums and nested radicals, which need code defined later; asin((√6 - √2)/4)
# stayed; SymPy #5342)
function _more_asin_values!()
    push!(_ASIN_VALUES, (sqrt(Sym(6)) - sqrt(Sym(2))) / 4 => PI / 12,
          (sqrt(Sym(6)) + sqrt(Sym(2))) / 4 => 5PI / 12,
          (sqrt(Sym(5)) - 1) / 4 => PI / 10, (sqrt(Sym(5)) + 1) / 4 => 3PI / 10,
          sqrt(2 - sqrt(Sym(2))) / 2 => PI / 8, sqrt(2 + sqrt(Sym(2))) / 2 => 3PI / 8,
          sqrt(10 - 2sqrt(Sym(5))) / 4 => PI / 5, sqrt(10 + 2sqrt(Sym(5))) / 4 => 2PI / 5)
end
const _ATAN_VALUES = Pair{Sym,Sym}[ZERO => ZERO, sqrt(Sym(3)) / 3 => PI / 6, ONE => PI / 4,
                                   sqrt(Sym(3)) => PI / 3, 2 - sqrt(Sym(3)) => PI / 12,
                                   2 + sqrt(Sym(3)) => 5PI / 12, sqrt(Sym(2)) - 1 => PI / 8,
                                   sqrt(Sym(2)) + 1 => 3PI / 8]
function _table_value(table, a::Sym)
    for (v, θ) in table
        _eq(a, v) && return θ
        _eq(a, -v) && return -θ
    end
    nothing
end

# atan(tan c), asin(sin c), acos(cos c) for a real constant c: c moved into the principal
# range by a multiple of the period (atan(tan 5) = 5 - 2π; SymPy #17300); nothing when the
# value is at the range's edge (decided only approximately) or not of this form
function _inv_periodic(name::Symbol, a::Sym)
    (a.kind == FN && length(a.args) == 1) || return nothing
    (name, a.name) in ((:atan, :tan), (:asin, :sin), (:acos, :cos)) || return nothing
    u = a.args[1]
    _hassymbol(u) && return nothing
    v = try _evalc(u) catch; return nothing end
    (isfinite(v) && abs(imag(v)) <= 1e-14 * max(1, abs(v))) || return nothing
    c = real(v)
    near(t, s) = abs(t - s) < 1e-9
    if name === :atan
        k = round(Int, c / π)
        r = c - k * π
        (near(r, π / 2) || near(r, -π / 2)) && return nothing
        return k == 0 ? u : expand(u - k * PI)
    end
    k = round(Int, c / (2π))
    r = c - 2k * π                                          # in [-π, π]
    w = k == 0 ? u : expand(u - 2k * PI)
    (near(abs(r), π) || near(r, 0.0) || near(abs(r), π / 2)) && return nothing
    name === :acos && return r > 0 ? w : -w
    abs(r) < π / 2 ? w : r > 0 ? PI - w : -PI - w
end

register_function!(:asin; julia = asin,
    eval = _unary(a -> something(_table_value(_ASIN_VALUES, a), _inv_periodic(:asin, a), _odd(:asin, a), Some(nothing))),
    diff = (args, i) -> inv(sqrt(1 - args[1]^2)))

register_function!(:acos; julia = acos,
    eval = _unary(a -> (θ = _table_value(_ASIN_VALUES, a); θ === nothing ? _inv_periodic(:acos, a) : PI / 2 - θ)),
    diff = (args, i) -> -inv(sqrt(1 - args[1]^2)))

register_function!(:atan; julia = atan,
    eval = _unary(a -> begin
        (_eq(a, IM) || _eq(a, -IM)) && throw(DivisionByZero("atan($a): a pole"))
        something(_table_value(_ATAN_VALUES, a), _inv_periodic(:atan, a), _odd(:atan, a), Some(nothing))
    end),
    diff = (args, i) -> inv(1 + args[1]^2))

register_function!(:sinh; julia = sinh,
    eval = _unary(a -> something(_of_own_inverse(:sinh, a), _hyp_of_inverse(:sinh, a), Some(iszero(a) ? ZERO : _odd(:sinh, a)))),
    diff = (args, i) -> cosh(args[1]))

register_function!(:cosh; julia = cosh,
    eval = _unary(a -> something(_of_own_inverse(:cosh, a), _hyp_of_inverse(:cosh, a), Some(iszero(a) ? ONE : _even(:cosh, a)))),
    diff = (args, i) -> sinh(args[1]))

register_function!(:tanh; julia = tanh,
    eval = _unary(a -> something(_of_own_inverse(:tanh, a), _hyp_of_inverse(:tanh, a), Some(iszero(a) ? ZERO : _odd(:tanh, a)))),
    diff = (args, i) -> 1 - tanh(args[1])^2)

register_function!(:asinh; julia = asinh,
    # asinh(i u) = i asin(u) at asin's exact values, |u| ≤ 1 (asinh(i) = iπ/2; SymPy #3906)
    eval = _unary(a -> iszero(a) ? ZERO : something(_imag_inverse(:asin, a), _odd(:asinh, a), Some(nothing))),
    diff = (args, i) -> inv(sqrt(args[1]^2 + 1)))

register_function!(:acosh; julia = acosh,
    # acosh(u) = i acos(u) for a real u in [-1, 1], as Julia's acosh(0 + 0im) = iπ/2 (SymPy #3663)
    eval = _unary(a -> isone(a) ? ZERO : _acosh_unit(a)),
    diff = (args, i) -> inv(sqrt(args[1] - 1) * sqrt(args[1] + 1)))

register_function!(:atanh; julia = atanh,
    eval = _unary(a -> begin
        a.kind == NUMBER && isexact(a.num) && abs(a.num) == N1 && throw(DivisionByZero("atanh($a): a pole"))
        # atanh(i u) = i atan(u) at atan's exact values (atanh(i) = iπ/4; SymPy #3906)
        iszero(a) ? ZERO : something(_imag_inverse(:atan, a), _odd(:atanh, a), Some(nothing))
    end),
    diff = (args, i) -> inv(1 - args[1]^2))

# i g(u) for a = i u, u a real number with |u| ≤ 1, when g(u) has an exact value
function _imag_inverse(g::Symbol, a::Sym)
    u = _imag_number(a)
    u === nothing && return nothing
    g === :asin && abs(u.num) > N1 && return nothing
    r = apply_function(g, Sym[u])
    _hasfn(r, g) ? nothing : IM * r
end
function _acosh_unit(a::Sym)
    (a.kind == NUMBER && isexact(a.num) && abs(a.num) <= N1) || return nothing
    r = apply_function(:acos, Sym[a])
    _hasfn(r, :acos) ? nothing : IM * r
end

register_function!(:abs; julia = abs,
    eval = _unary(a -> begin
        a.kind == NUMBER && return Sym(abs(a.num))
        _isfn(a, :abs) && return a
        # a real constant of known sign: ±itself (abs(√2) = √2, abs(-π) = π; SymPy #11193 area)
        if !_hassymbol(a) && !has(a, IM) && !_hasfn(a, :abs)
            s = Base.inferencebarrier(_csign)(a)::Union{Nothing,Int}
            s === nothing || return s >= 0 ? a : -a
        end
        _even(:abs, a)
    end),
    diff = (args, i) -> args[1] / abs(args[1]))

# log(x) and log(b, x) (base first, as in Base.log); log(b, x) is kept unevaluated
register_function!(:log; julia = log,
    eval = function (args)
        if length(args) == 1
            a = args[1]
            isone(a) && return ZERO
            a === E && return ONE
            # log(0) = -∞ (as for floats; then 0*log(0) is NaN, not 0; SymPy #22612)
            a.kind == NUMBER && isexact(a.num) && iszero(a.num) && return Sym(-Inf)
            if a.kind == POW && a.args[1] === E && a.args[2].kind == NUMBER
                return a.args[2]                 # log(e^r) == r for real r
            end
            # log(|u|^q) = q log|u| for a real number q: |u| ≥ 0 (log(|x|²)/log|x| is 2;
            # SymPy #15854)
            if a.kind == POW && _isfn1(a.args[1], :abs) && a.args[2].kind == NUMBER && isreal(a.args[2].num)
                return a.args[2] * log(a.args[1])
            end
            # log(e^z), z a constant: z with its imaginary part moved into (-π, π] by a
            # multiple of 2π (log(e^(iπ/3)) = iπ/3; SymPy #17208)
            if a.kind == POW && a.args[1] === E && !_hassymbol(a.args[2])
                z = a.args[2]
                v = try _evalc(z) catch; ComplexF64(NaN) end
                if isfinite(v)
                    k = round(Int, imag(v) / (2π))
                    r = imag(v) - 2π * k
                    abs(abs(r) - π) > 1e-9 && return k == 0 ? z : expand(z - 2k * PI * IM)
                end
            end
            # the principal value at an exact negative number or ±im: log(-c) = log(c) + iπ
            a.kind == NUMBER && isexact(a.num) && a.num < N0 && return log(Sym(-a.num)) + IM * PI
            _eq(a, IM) && return IM * PI / 2
            _eq(a, -IM) && return -IM * PI / 2
            # a complex constant of modulus 1 with an angle a rational multiple of π:
            # i·angle (log(1/2 + i√3/2) = iπ/3; SymPy #10428)
            if !_hassymbol(a) && has(a, IM) && a.kind in (ADD, MUL)
                v = try _evalc(a) catch; nothing end
                if v isa Number && isfinite(v) && abs(abs(v) - 1) < 1e-12 && !iszero(imag(v))
                    q = rationalize(angle(v) / π; tol = 1e-12)
                    # (checked exactly by expand, no simplify: simplify takes logs itself)
                    if denominator(q) <= 24
                        θ = Sym(q) * PI
                        iszero(expand(cos(θ) + IM * sin(θ) - a)) && return IM * θ
                    end
                end
            end
            return nothing
        end
        length(args) == 2 || throw(ArgumentError("log takes one or two arguments"))
        b, x = args
        # base 1: log x/log 1 is undefined, 0/0 at x = 1 (it was 0; SymPy #5100)
        isone(b) && return isone(x) ? Sym(NaN) : nothing
        isone(x) && return ZERO
        _eq(b, x) && return ONE
        b === E && return apply_function(:log, Sym[x])
        # log(b, b^k) = k for positive rationals b ≠ 1 and x (log(2, 8) = 3; SymPy #20953)
        if b.kind == NUMBER && x.kind == NUMBER && isexact(b.num) && isexact(x.num) && b.num > N0 && x.num > N0 && !isone(b.num)
            k = round(Int, log(Float64(_rat(x.num))) / log(Float64(_rat(b.num))))
            abs(k) <= 10_000 && _rat(b.num)^k == _rat(x.num) && return Sym(k)
            # powers of one integer: log(4, 8) = 3/2 (SymPy #7228)
            q = _log_ratio_int(_rat(x.num), _rat(b.num))
            q === nothing || return Sym(q)
        end
        nothing
    end,
    diff = function (args, i)
        length(args) == 1 && return inv(args[1])
        b, x = args
        i == 2 ? inv(x * log(b)) : -log(x) / (b * log(b)^2)
    end)

# log(x)/log(b) for x = p^m, b = p^n integers > 1 (or reciprocals) of one p: m/n, else nothing
function _log_ratio_int(x::Rational, b::Rational)
    function pk(r::Rational)
        s = 1
        r < 1 && (r = inv(r); s = -1)
        isone(denominator(r)) && numerator(r) > 1 || return nothing
        p, k = _int_root_power(big(numerator(r)))
        (p, s * k)
    end
    (px, kx), (pb, kb) = something(pk(x), (0, 0)), something(pk(b), (0, 0))
    (px == 0 || pb == 0 || px != pb) && return nothing
    kx // kb
end

## reciprocal and inverse-reciprocal functions
#
# These stay as their own functions (sec(x), not 1/cos(x)). Their exact values come
# from the rules of the function they are defined through, accepted only when that
# function no longer appears: sec(π/3) = 1/cos(π/3) = 2 and sec(acos(x)) = 1/x, while
# sec(x) stays sec(x). A pole (sec(π/2)) is a DivisionByZero error.

_hasfn(ex::Sym, name::Symbol) =
    (ex.kind == FN && ex.name === name) || any(a -> _hasfn(a, name), ex.args) || any(a -> _hasfn(a, name), ex.exps)

# f(a) = num(a) / den(a), with num and den among the existing functions (num may be ONE)
# f(f⁻¹(u)) = u, for every u where it is defined (coth(acoth(x)) stayed; SymPy #5346)
const _OWN_INVERSE = Dict(:sec => :asec, :csc => :acsc, :cot => :acot, :sech => :asech, :csch => :acsch,
                          :coth => :acoth, :sinh => :asinh, :cosh => :acosh, :tanh => :atanh)
_of_own_inverse(name::Symbol, a::Sym) =
    a.kind == FN && length(a.args) == 1 && a.name === get(_OWN_INVERSE, name, nothing) ? a.args[1] : nothing

function _ratio_rule(name::Symbol, num::Union{Nothing,Symbol}, den::Symbol, parity::Symbol)
    _unary(function (a)
        v = _of_own_inverse(name, a)
        v === nothing || return v
        d = apply_function(den, Sym[a])
        n = num === nothing ? ONE : apply_function(num, Sym[a])
        if !_hasfn(d, den) && (num === nothing || !_hasfn(n, num))
            iszero(d) || return n / d
            # a pole: an error, as 1/0 is (not a finite-looking value; SymPy #27185)
            throw(DivisionByZero("$name($a): a pole"))
        end
        parity === :odd ? _odd(name, a) : parity === :even ? _even(name, a) : nothing
    end)
end

# f(a) = g(1/a), with g an existing inverse function
function _invrecip_rule(name::Symbol, g::Symbol, parity::Symbol, atzero)
    _unary(function (a)
        # a pole at 0 (acsc, asec, acsch, asech), as csc(0) is one (acsc(0) stayed, and N
        # gave NaN; SymPy #17143)
        iszero(a) && (atzero === nothing ? throw(DivisionByZero("$name(0): a pole")) : return atzero)
        r = apply_function(g, Sym[inv(a)])
        # (not a complex value for a real a: 1/(a + 0im) has a -0.0 imaginary part, the other
        # side of g's cut, so asech(-2) is -2iπ/3, not acosh(-1/2) = 2iπ/3)
        _hasfn(r, g) || (has(r, IM) && !has(a, IM)) || return r
        # :odd_nonzero: odd for a nonzero number only (a symbol may be 0)
        parity === :odd || (parity === :odd_nonzero && a.kind == NUMBER) ? _odd(name, a) : nothing
    end)
end

register_function!(:sec; julia = sec, eval = _ratio_rule(:sec, nothing, :cos, :even),
    diff = (args, i) -> sec(args[1]) * tan(args[1]))
register_function!(:csc; julia = csc, eval = _ratio_rule(:csc, nothing, :sin, :odd),
    diff = (args, i) -> -csc(args[1]) * cot(args[1]))
register_function!(:cot; julia = cot, eval = _ratio_rule(:cot, :cos, :sin, :odd),
    diff = (args, i) -> -1 - cot(args[1])^2)
register_function!(:sech; julia = sech, eval = _ratio_rule(:sech, nothing, :cosh, :even),
    diff = (args, i) -> -sech(args[1]) * tanh(args[1]))
register_function!(:csch; julia = csch, eval = _ratio_rule(:csch, nothing, :sinh, :odd),
    diff = (args, i) -> -csch(args[1]) * coth(args[1]))
register_function!(:coth; julia = coth, eval = _ratio_rule(:coth, :cosh, :sinh, :odd),
    diff = (args, i) -> 1 - coth(args[1])^2)

# as in Base: asec(x) = acos(1/x), acot(x) = atan(1/x), asech(x) = acosh(1/x), ...
# acot and acoth are odd except at 0 (acot(0) = π/2), so a minus sign is taken out only
# of a nonzero number: acot(-x) stays, as x may be 0
register_function!(:asec; julia = asec, eval = _invrecip_rule(:asec, :acos, :none, nothing),
    diff = (args, i) -> inv(args[1]^2 * sqrt(1 - inv(args[1]^2))))
register_function!(:acsc; julia = acsc, eval = _invrecip_rule(:acsc, :asin, :odd, nothing),
    diff = (args, i) -> -inv(args[1]^2 * sqrt(1 - inv(args[1]^2))))
register_function!(:acot; julia = acot, eval = _invrecip_rule(:acot, :atan, :odd_nonzero, PI / 2),
    diff = (args, i) -> -inv(1 + args[1]^2))
register_function!(:asech; julia = asech, eval = _invrecip_rule(:asech, :acosh, :none, nothing),
    diff = (args, i) -> -inv(args[1]^2 * sqrt(inv(args[1]) - 1) * sqrt(inv(args[1]) + 1)))
register_function!(:acsch; julia = acsch, eval = _invrecip_rule(:acsch, :asinh, :odd, nothing),
    diff = (args, i) -> -inv(args[1]^2 * sqrt(1 + inv(args[1]^2))))
# (acoth(0) = -iπ/2: Julia's acoth(u + 0im) is atanh(u) - iπ/2 for |u| < 1, on both sides of
# 0; it was a pole; SymPy #3663)
register_function!(:acoth; julia = acoth, eval = _invrecip_rule(:acoth, :atanh, :odd_nonzero, -IM * PI / 2),
    diff = (args, i) -> inv(1 - args[1]^2))

## real cube root: differs from x^(1/3) (the principal root) for negative x
register_function!(:cbrt; julia = cbrt,
    eval = _unary(function (a)
        (a.kind == NUMBER && isexact(a.num)) || return _odd(:cbrt, a)
        sign(a.num) < 0 ? -pow(-a, Sym(SymNumber(1 // 3))) : pow(a, Sym(SymNumber(1 // 3)))
    end),
    diff = (args, i) -> inv(3 * cbrt(args[1])^2))

## piecewise constant functions. Their derivatives are left unevaluated, since they
## are zero only away from the jumps.

# the value of a constant expression, when safely away from `avoid(v)` (a jump)
function _constant_float(a::Sym, avoid)
    (a.kind == NUMBER || !_hassymbol(a)) || return nothing
    v = try
        N(a)
    catch
        return nothing
    end
    v isa Real || return nothing
    abs(v - avoid(v)) > 1e-10 * max(1, abs(v)) && return v
    # too near a jump for Float64: log(2, 2^360 - 1) + 1 is 361 - 1.4e-109, 361.0 in
    # Float64 (SymPy #26368); decided in BigFloat when the precision separates them
    # (a value is taken when two precisions agree on it, and it holds its fractional part:
    # π^1000 at 1024 bits is an integer value, its floor came out 4e191 off, and π^1000 -
    # floor(π^1000) at 1024 bits is rounding noise of size 1e189; SymPy #6167)
    prev, pprev = nothing, 0
    for p in (256, 1024, 4096, 16384)
        w = _bigvalue(a, p)
        (w === nothing || !iszero(imag(w))) && return nothing
        r = real(w)
        isfinite(r) || (prev = nothing; continue)
        agree = prev !== nothing && abs(r - prev) <= big(2.0)^(64 - pprev ÷ 2) * max(1, abs(r))
        prev, pprev = r, p
        (agree && (iszero(r) || exponent(r) + 64 < p)) || continue
        ok = setprecision(BigFloat, p) do
            abs(r - avoid(r)) > big(2.0)^(32 - p) * max(1, abs(r))
        end
        ok && return r
    end
    nothing
end

register_function!(:sign; julia = sign,
    eval = _unary(function (a)
        a.kind == NUMBER && isexact(a.num) && return Sym(sign(a.num))
        _isfn(a, :sign) && return a
        v = _constant_float(a, _ -> 0.0)
        v === nothing || return Sym(Int(sign(v)))
        _exact_jump_value(a, ZERO) && return ZERO
        # a complex constant: a/|a| (sign(3 + 4i) = 3/5 + 4i/5; SymPy #5763)
        if !_hassymbol(a) && _has_imag_unit(a)
            m = abs(a)
            _isfn(m, :abs) || iszero(m) || return expand(a / m)
        end
        # a positive constant factor dropped (sign(2x) = sign(x); SymPy #4708)
        if a.kind == MUL
            pos = [i for i in eachindex(a.args) if !_hassymbol(a.args[i]) &&
                   (w = _constant_float(a.args[i]^a.exps[i], _ -> 0.0); w !== nothing && w > 0)]
            if !isone(abs(a.num)) || !isempty(pos)
                rest = Sym(sign(a.num)) * prod((a.args[i]^a.exps[i] for i in eachindex(a.args) if !(i in pos)); init = ONE)
                _eq(rest, a) || return sign(rest)
            end
        end
        _odd(:sign, a)
    end))

# whether a constant within rounding of the jump n (one no precision separates) is n
# exactly, by simplify (a disguised 0, 5log5/(6log2 + 6log3 - 5log5 - 6log6) + 1: sign 0,
# floor 0; SymPy #20897, #20865)
function _exact_jump_value(a::Sym, n::Sym)
    _hassymbol(a) && return false
    w = try _evalc(a) catch; return false end
    (isfinite(w) && abs(w - ComplexF64(N(n))) <= 1e-10 * max(1, abs(w))) || return false
    d = try simplify(a - n) catch err; err isa InterruptException && rethrow(); return false end
    iszero(d)
end

# the jumps to stay away from: integers, or for round the half-integers
for (f, g, jump) in ((:floor, floor, round), (:ceil, ceil, round), (:trunc, trunc, round),
                     (:round, round, v -> floor(v) + 0.5))
    register_function!(f; julia = g,
        eval = _unary(function (a)
            a.kind == NUMBER && isexact(a.num) && return Sym(SymNumber(g(BigInt, _rat(a.num))))
            _isfn(a, f) && return a
            # an integer-valued argument (componentwise): floor(ceil(x)) = ceil(x), any
            # rounding of floor/ceil/round/trunc of something likewise (SymPy #11207)
            a.kind == FN && a.name in (:floor, :ceil, :round, :trunc) && length(a.args) == 1 && return a
            # a complex number: its parts (round(2.6 + 1.2i) = 3 + i; SymPy #21513)
            if !_hassymbol(a) && has(a, IM)
                v = try _evalc(a) catch; nothing end
                if v isa Complex && isfinite(v) && !iszero(imag(v))
                    re, ip = real(a), imag(a)
                    (has(re, IM) || has(ip, IM)) && return nothing
                    return apply_function(f, Sym[re]) + IM * apply_function(f, Sym[ip])
                end
            end
            # integer-valued terms out: floor(floor(x) + 1) = floor(x) + 1, as floor(G + r)
            # = G + floor(r) for G a floor, ceil or integer (componentwise for complex
            # values too); not trunc or round (SymPy #18689)
            if f in (:floor, :ceil) && a.kind == ADD
                intv(t) = t.kind == FN && t.name in (:floor, :ceil) && length(t.args) == 1
                ix = [i for i in eachindex(a.args) if intv(a.args[i]) && isexact(a.coefs[i]) && _isint(a.coefs[i])]
                n0 = isexact(a.num) && _isint(a.num) ? Sym(a.num) : ZERO
                if !isempty(ix) || !iszero(n0)
                    G = n0 + _add_many(Sym(a.coefs[i]) * a.args[i] for i in ix)
                    rest = expand(a - G)
                    if !_eq(rest, a)
                        return G + apply_function(f, Sym[rest])
                    end
                end
            end
            # floor(floor(u)/m) = floor(u/m) for an integer m > 0, ceil alike (SymPy #27400)
            if f in (:floor, :ceil) && a.kind == MUL && length(a.args) == 1 && _isone(a.exps[1]) &&
               _isfn(a.args[1], f) && isexact(a.num) && isone(_numerator(a.num)) && _denominator(a.num) > 1
                return apply_function(f, Sym[a.args[1].args[1] * Sym(a.num)])
            end
            v = _constant_float(a, jump)
            v === nothing || return Sym(SymNumber(g(BigInt, v)))
            # (at a jump exactly: the integer itself, for floor, ceil, trunc)
            if f !== :round && !_hassymbol(a)
                w = try _evalc(a) catch; nothing end
                if w isa Number && isfinite(w) && abs(imag(w)) < 1e-10
                    n = Sym(SymNumber(round(BigInt, real(w))))
                    _exact_jump_value(a, n) && return n
                end
            end
            # floor(u + r) = floor(u + r - n) + n for the integer part n of a rational r, so
            # that floor(z + 3/4) - floor(z - 1/4) is 1; ceil alike
            if f in (:floor, :ceil) && a.kind == ADD && isexact(a.num) && !_isint(a.num)
                n = Sym(SymNumber(floor(BigInt, _rat(a.num))))
                iszero(n) || return apply_function(f, Sym[a - n]) + n
            end
            # floor(-v) = -ceil(v), ceil(-v) = -floor(v); trunc and round (ties to even) are odd
            _could_extract_minus(a) || return nothing
            -apply_function(f === :floor ? :ceil : f === :ceil ? :floor : f, Sym[-a])
        end))
end

# Signatures for extending a function of `n` arguments to expressions: each argument
# a `Sym` or a `Number`, at least one a `Sym`. Used by the package extensions.
_symbolic_signatures(n::Int) =
    [Tuple(T) for T in Iterators.product(ntuple(_ -> (Sym, Number), n)...) if any(==(Sym), T)]

for f in (:sin, :cos, :tan, :asin, :acos, :atan, :sinh, :cosh, :tanh, :asinh, :acosh, :atanh, :abs, :log,
          :sec, :csc, :cot, :asec, :acsc, :acot, :sech, :csch, :coth, :asech, :acsch, :acoth,
          :cbrt, :sign)
    @eval Base.$f(x::Sym) = apply_function($(QuoteNode(f)), Sym[x])
end
# (floor, ceil, trunc and round, which take a rounding mode and digits, are in
# mathfunctions.jl)

# functions with an existing symbolic form
Base.log2(x::Sym) = log(Sym(2), x)
Base.log10(x::Sym) = log(Sym(10), x)
Base.log1p(x::Sym) = log(1 + x)
Base.exp2(x::Sym) = Sym(2)^x
Base.exp10(x::Sym) = Sym(10)^x
Base.expm1(x::Sym) = exp(x) - 1
Base.sinpi(x::Sym) = sin(PI * x)
Base.cospi(x::Sym) = cos(PI * x)
Base.sincos(x::Sym) = (sin(x), cos(x))
Base.log(b::Sym, x::Sym) = apply_function(:log, Sym[b, x])
Base.log(b::Number, x::Sym) = log(Sym(b), x)
Base.log(b::Sym, x::Number) = log(b, Sym(x))
