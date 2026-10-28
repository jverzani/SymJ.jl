# Rewriting between forms (SymPy's rewrite): the functions of an expression written
# through a target function, as trigonometric and hyperbolic functions through
# exponentials, inverse ones through logarithms, exponentials of imaginary arguments
# through sine and cosine, or factorials through gamma. Each target is a table from the
# name of a function to its form in the target, applied bottom up.

# u = a + b im, the symbols taken as real (the conjugate is u with im replaced by -im)
function _imsplit(u::Sym)
    c = subs(u, IM => -IM)
    c === u && return (u, ZERO)
    (expand((u + c) / 2), expand((u - c) / (2IM)))
end

_eip(u) = exp(IM * u)
_eim(u) = exp(-IM * u)

# trigonometric and hyperbolic functions through exponentials
const _TO_EXP = Dict{Symbol,Function}(
    :sin => u -> (_eip(u) - _eim(u)) / (2IM),
    :cos => u -> (_eip(u) + _eim(u)) / 2,
    :tan => u -> (_eip(u) - _eim(u)) / (IM * (_eip(u) + _eim(u))),
    :cot => u -> IM * (_eip(u) + _eim(u)) / (_eip(u) - _eim(u)),
    :sec => u -> 2 / (_eip(u) + _eim(u)),
    :csc => u -> 2IM / (_eip(u) - _eim(u)),
    :sinh => u -> (exp(u) - exp(-u)) / 2,
    :cosh => u -> (exp(u) + exp(-u)) / 2,
    :tanh => u -> (exp(u) - exp(-u)) / (exp(u) + exp(-u)),
    :coth => u -> (exp(u) + exp(-u)) / (exp(u) - exp(-u)),
    :sech => u -> 2 / (exp(u) + exp(-u)),
    :csch => u -> 2 / (exp(u) - exp(-u)),
)

# inverse functions through logarithms (principal branches)
_asin_log(u) = -IM * log(IM * u + sqrt(1 - u^2))
_acosh_log(u) = log(u + sqrt(u - 1) * sqrt(u + 1))
_asinh_log(u) = log(u + sqrt(u^2 + 1))
const _TO_LOG = Dict{Symbol,Function}(
    :asin => _asin_log,
    :acos => u -> PI / 2 - _asin_log(u),
    :atan => u -> IM * (log(1 - IM * u) - log(1 + IM * u)) / 2,
    :acot => u -> IM * (log(1 - IM / u) - log(1 + IM / u)) / 2,
    :asec => u -> PI / 2 - _asin_log(1 / u),
    :acsc => u -> _asin_log(1 / u),
    :asinh => _asinh_log,
    :acosh => _acosh_log,
    :atanh => u -> (log(1 + u) - log(1 - u)) / 2,
    :acoth => u -> (log(1 + 1 / u) - log(1 - 1 / u)) / 2,
    :asech => u -> _acosh_log(1 / u),
    :acsch => u -> _asinh_log(1 / u),
)

# exp(a + b im) = exp(a)(cos b + im sin b)
function _exp_sincos(u)
    a, b = _imsplit(u)
    iszero(b) && return nothing
    exp(a) * (cos(b) + IM * sin(b))
end
const _TO_SINCOS = Dict{Symbol,Function}(
    :tan => u -> sin(u) / cos(u), :cot => u -> cos(u) / sin(u),
    :sec => u -> 1 / cos(u), :csc => u -> 1 / sin(u),
    :exp => _exp_sincos,
)

# exp(a + b im) = (cosh a + sinh a) exp(b im)
function _exp_sinhcosh(u)
    a, b = _imsplit(u)
    iszero(a) && return nothing
    (cosh(a) + sinh(a)) * exp(IM * b)
end
const _TO_SINHCOSH = Dict{Symbol,Function}(
    :tanh => u -> sinh(u) / cosh(u), :coth => u -> cosh(u) / sinh(u),
    :sech => u -> 1 / cosh(u), :csch => u -> 1 / sinh(u),
    :exp => _exp_sinhcosh,
)

# through the tangent of the half angle, t = tan(u/2) (tanh(u/2))
_ht(u) = tan(u / 2)
_hth(u) = tanh(u / 2)
const _TO_TAN = Dict{Symbol,Function}(
    :sin => u -> 2_ht(u) / (1 + _ht(u)^2),
    :cos => u -> (1 - _ht(u)^2) / (1 + _ht(u)^2),
    :sec => u -> (1 + _ht(u)^2) / (1 - _ht(u)^2),
    :csc => u -> (1 + _ht(u)^2) / 2_ht(u),
    :cot => u -> 1 / tan(u),
)
const _TO_TANH = Dict{Symbol,Function}(
    :sinh => u -> 2_hth(u) / (1 - _hth(u)^2),
    :cosh => u -> (1 + _hth(u)^2) / (1 - _hth(u)^2),
    :sech => u -> (1 - _hth(u)^2) / (1 + _hth(u)^2),
    :csch => u -> (1 - _hth(u)^2) / 2_hth(u),
    :coth => u -> 1 / tanh(u),
)

# factorials, binomials and beta through gamma, and back
_gammafn(u) = apply_function(:gamma, Sym[u])
const _TO_GAMMA = Dict{Symbol,Function}(
    :factorial => u -> _gammafn(u + 1),
    :binomial => (n, k) -> _gammafn(n + 1) / (_gammafn(k + 1) * _gammafn(n - k + 1)),
    :beta => (a, b) -> _gammafn(a) * _gammafn(b) / _gammafn(a + b),
)
const _TO_FACTORIAL = Dict{Symbol,Function}(
    :gamma => u -> factorial(u - 1),
    :binomial => (n, k) -> factorial(n) / (factorial(k) * factorial(n - k)),
    :beta => (a, b) -> factorial(a - 1) * factorial(b - 1) / factorial(a + b - 1),
)

const _REWRITE_TARGETS = Dict{Symbol,Dict{Symbol,Function}}(
    :exp => _TO_EXP, :log => _TO_LOG,
    :sin => _TO_SINCOS, :cos => _TO_SINCOS,
    :sinh => _TO_SINHCOSH, :cosh => _TO_SINHCOSH,
    :tan => _TO_TAN, :tanh => _TO_TANH,
    :gamma => _TO_GAMMA, :factorial => _TO_FACTORIAL,
)

"""
    SymJ.rewrite(ex, target)
    SymJ.rewrite(ex, fs => target)

`ex` with its functions written through the function `target`, as SymPy's `rewrite`.
With `fs => target` (`fs` a function or a tuple of them), only those functions are
rewritten. The targets:

- `exp`: trigonometric and hyperbolic functions, `sin(x) = (exp(im*x) - exp(-im*x))/(2im)`,
  `cosh(x) = (exp(x) + exp(-x))/2`;
- `log`: inverse trigonometric and hyperbolic functions (principal branches),
  `asinh(x) = log(x + sqrt(x^2 + 1))`, `atan(x) = im*(log(1 - im*x) - log(1 + im*x))/2`;
- `sin` or `cos`: `tan`, `cot`, `sec`, `csc` through `sin` and `cos`, and exponentials
  with imaginary parts, `exp(a + b*im) = exp(a)(cos(b) + im*sin(b))` (the symbols taken as
  real);
- `sinh` or `cosh`: the same for hyperbolic functions, `exp(x) = cosh(x) + sinh(x)`;
- `tan` (`tanh`): through the tangent of the half angle, `sin(x) = 2t/(1 + t^2)` with
  `t = tan(x/2)`;
- `gamma` (from SpecialFunctions) and `factorial`: factorials, binomials and `beta` through
  `gamma`, and gamma, binomials and `beta` through factorials.

The result is not simplified. Not exported.

```julia
SymJ.rewrite(sin(x), exp)                 # -im*(exp(im*x) - exp(-im*x))/2
SymJ.rewrite(exp(im*x), cos)              # cos(x) + im*sin(x)
SymJ.rewrite(tan(x)^2 + sec(x), sin)      # sin(x)^2/cos(x)^2 + 1/cos(x)
SymJ.rewrite(asinh(x), log)               # log(x + sqrt(x^2 + 1))
SymJ.rewrite(binomial(n, k), factorial)   # factorial(n)/(factorial(k)*factorial(n - k))
SymJ.rewrite(sin(x) + cos(x), sin => exp) # only sin rewritten
```
"""
rewrite(ex, target::Function) = _rewrite_to(Sym(ex), _rewrite_table(target), nothing)
function rewrite(ex, p::Pair)
    fs = p.first isa Union{Tuple,AbstractVector} ? p.first : (p.first,)
    all(f -> f isa Function && !(f isa Sym), fs) ||
        throw(ArgumentError("rewrite: in fs => target, fs is a function or a tuple of functions, not $(p.first)"))
    p.second isa Function || throw(ArgumentError("rewrite: the target $(p.second) is not a function"))
    _rewrite_to(Sym(ex), _rewrite_table(p.second), Set{Symbol}(nameof(f) for f in fs))
end

function _rewrite_table(target::Function)
    target isa Sym && throw(ArgumentError("rewrite: the target is a function, such as exp or sin, not the expression $target"))
    t = get(_REWRITE_TARGETS, nameof(target), nothing)
    t === nothing && throw(ArgumentError("rewrite: no rewriting in terms of $target; the targets are " *
                                         join(sort!(collect(keys(_REWRITE_TARGETS))), ", ")))
    t
end

function _rewrite_to(e::Sym, table::Dict{Symbol,Function}, only)
    _replace_walk(e) do t
        if t.kind == FN
            name, args = t.name, t.args
        elseif t.kind == POW && _eq(t.args[1], E)
            name, args = :exp, t.args[2:2]
        else
            return nothing
        end
        (only === nothing || name in only) || return nothing
        f = get(table, name, nothing)
        f === nothing && return nothing
        applicable(f, args...) || return nothing
        f(args...)
    end
end
