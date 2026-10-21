# Unevaluated function applications (kind FN) and the rules that simplify them.
#
# A function is identified by its name. Behaviour is looked up in `FUNCTION_RULES`:
#   julia: a Julia function for numeric evaluation (and TermInterface.operation)
#   eval:  args -> Union{Nothing,Sym}; return a simplified value or `nothing`
#   diff:  (args, i) -> partial derivative with respect to the i-th argument
# Names with no rule (user functions) stay unevaluated.

"""
    SymFunction(name)

A symbolic function. Calling it builds an unevaluated application:

```julia
f = SymFunction(:f)     # or `@syms f()`
f(x, 2)                 # f(x, 2)
```
"""
struct SymFunction <: AbstractSymbolic
    name::Symbol
end
(f::SymFunction)(args...) = apply_function(f.name, Sym[Sym(a) for a in args])
Base.show(io::IO, f::SymFunction) = print(io, f.name)
Base.show(io::IO, ::MIME"text/plain", f::SymFunction) = show(io, f)

struct FunctionRule
    julia::Union{Nothing,Function}
    eval::Union{Nothing,Function}
    diff::Union{Nothing,Function}
end

const FUNCTION_RULES = Dict{Symbol,FunctionRule}()

"""
    register_function!(name; julia=nothing, eval=nothing, diff=nothing)

Attach behaviour to the function `name`: a Julia function used for numeric evaluation,
a simplification rule `args -> Union{Nothing,Sym}`, and a partial-derivative rule
`(args, i) -> Sym`.
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
        if rule.eval !== nothing
            r = rule.eval(args)
            r === nothing || return r::Sym
        end
        if rule.julia !== nothing && _float_args(args)
            v = _try_numeric(rule.julia, args)
            v === nothing || return v
        end
    end
    _mkfn_node(name, args)
end

_float_args(args) = all(a -> a.kind == NUMBER, args) && any(a -> isfloat(a.num), args)

function _try_numeric(f, args)
    try
        return Sym(f(map(a -> Float64(a.num), args)...))
    catch err
        err isa DomainError || rethrow()
        return nothing
    end
end

_fn_operation(name::Symbol) =
    (r = get(FUNCTION_RULES, name, nothing); r === nothing || r.julia === nothing ? SymFunction(name) : r.julia)

## helpers for rules

# Would `-x` print more simply than `x`? Exactly one of `x`, `-x` answers true (or both
# false), which makes it usable to canonicalize odd and even functions.
function _could_extract_minus(x::Sym)
    k = x.kind
    (k == NUMBER || k == MUL) && return sign(x.num) < 0
    k == ADD && return iszero(x.num) ? sign(x.coefs[1]) < 0 : sign(x.num) < 0
    false
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

# sin(rπ) when it has a simple closed form
function _sinpi(r::Rational{BigInt})
    r = mod(r, 2)
    neg = r >= 1
    neg && (r -= 1)
    r > 1 // 2 && (r = 1 - r)
    v = r == 0 ? ZERO :
        r == 1 // 6 ? HALF :
        r == 1 // 4 ? sqrt(Sym(2)) / 2 :
        r == 1 // 3 ? sqrt(Sym(3)) / 2 :
        r == 1 // 2 ? ONE : nothing
    v === nothing && return nothing
    neg ? -v : v
end

_isfn(x::Sym, name::Symbol) = x.kind == FN && x.name === name && length(x.args) == 1

_odd(name, a) = _could_extract_minus(a) ? -apply_function(name, Sym[-a]) : nothing
_even(name, a) = _could_extract_minus(a) ? apply_function(name, Sym[-a]) : nothing

function _unary(f)
    args -> length(args) == 1 ? f(args[1]) : throw(ArgumentError("expected one argument, got $(length(args))"))
end

## elementary functions

register_function!(:sin; julia = sin,
    eval = _unary(a -> begin
        r = _pi_multiple(a)
        r !== nothing && (v = _sinpi(r); v === nothing || return v)
        _isfn(a, :asin) && return a.args[1]
        _odd(:sin, a)
    end),
    diff = (args, i) -> cos(args[1]))

register_function!(:cos; julia = cos,
    eval = _unary(a -> begin
        r = _pi_multiple(a)
        r !== nothing && (v = _sinpi(r + 1 // 2); v === nothing || return v)
        _isfn(a, :acos) && return a.args[1]
        _even(:cos, a)
    end),
    diff = (args, i) -> -sin(args[1]))

register_function!(:tan; julia = tan,
    eval = _unary(a -> begin
        r = _pi_multiple(a)
        if r !== nothing
            s, c = _sinpi(r), _sinpi(r + 1 // 2)
            s !== nothing && c !== nothing && !iszero(c) && return s / c
        end
        _isfn(a, :atan) && return a.args[1]
        _odd(:tan, a)
    end),
    diff = (args, i) -> 1 + tan(args[1])^2)

register_function!(:asin; julia = asin,
    eval = _unary(a -> iszero(a) ? ZERO : isone(a) ? PI / 2 : _odd(:asin, a)),
    diff = (args, i) -> inv(sqrt(1 - args[1]^2)))

register_function!(:acos; julia = acos,
    eval = _unary(a -> isone(a) ? ZERO : iszero(a) ? PI / 2 : nothing),
    diff = (args, i) -> -inv(sqrt(1 - args[1]^2)))

register_function!(:atan; julia = atan,
    eval = _unary(a -> iszero(a) ? ZERO : isone(a) ? PI / 4 : _odd(:atan, a)),
    diff = (args, i) -> inv(1 + args[1]^2))

register_function!(:sinh; julia = sinh,
    eval = _unary(a -> iszero(a) ? ZERO : _odd(:sinh, a)),
    diff = (args, i) -> cosh(args[1]))

register_function!(:cosh; julia = cosh,
    eval = _unary(a -> iszero(a) ? ONE : _even(:cosh, a)),
    diff = (args, i) -> sinh(args[1]))

register_function!(:tanh; julia = tanh,
    eval = _unary(a -> iszero(a) ? ZERO : _odd(:tanh, a)),
    diff = (args, i) -> 1 - tanh(args[1])^2)

register_function!(:asinh; julia = asinh,
    eval = _unary(a -> iszero(a) ? ZERO : _odd(:asinh, a)),
    diff = (args, i) -> inv(sqrt(args[1]^2 + 1)))

register_function!(:acosh; julia = acosh,
    eval = _unary(a -> isone(a) ? ZERO : nothing),
    diff = (args, i) -> inv(sqrt(args[1] - 1) * sqrt(args[1] + 1)))

register_function!(:atanh; julia = atanh,
    eval = _unary(a -> iszero(a) ? ZERO : _odd(:atanh, a)),
    diff = (args, i) -> inv(1 - args[1]^2))

register_function!(:abs; julia = abs,
    eval = _unary(a -> begin
        a.kind == NUMBER && return Sym(abs(a.num))
        _isfn(a, :abs) && return a
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
            if a.kind == POW && a.args[1] === E && a.args[2].kind == NUMBER
                return a.args[2]                 # log(e^r) == r for real r
            end
            return nothing
        end
        length(args) == 2 || throw(ArgumentError("log takes one or two arguments"))
        b, x = args
        isone(x) && return ZERO
        _eq(b, x) && return ONE
        b === E && return apply_function(:log, Sym[x])
        nothing
    end,
    diff = function (args, i)
        length(args) == 1 && return inv(args[1])
        b, x = args
        i == 2 ? inv(x * log(b)) : -log(x) / (b * log(b)^2)
    end)

for f in (:sin, :cos, :tan, :asin, :acos, :atan, :sinh, :cosh, :tanh, :asinh, :acosh, :atanh, :abs, :log)
    @eval Base.$f(x::Sym) = apply_function($(QuoteNode(f)), Sym[x])
end
Base.log(b::Sym, x::Sym) = apply_function(:log, Sym[b, x])
Base.log(b::Number, x::Sym) = log(Sym(b), x)
Base.log(b::Sym, x::Number) = log(b, Sym(x))
