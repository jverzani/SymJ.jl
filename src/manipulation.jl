## traversal

"""
    free_symbols(ex)

The variables appearing in `ex`, sorted by name, numbers in names by value: `x2` before
`x10`, and `@syms xs[1:10]` in the order of the indices.
"""
function free_symbols(ex::Sym)
    s = Set{Sym}()
    _free!(s, ex)
    _sort_by_name!(collect(s))
end

# symbols in the natural order of their names: runs of digits (subscripts too, negative
# after ₋) compare as numbers, so x2 < x10 and xs₋₁ < xs₀ < xs₁ < xs₁₀
_sort_by_name!(v::AbstractVector) = sort!(v; lt = (a, b) -> _natural_lt(string(a), string(b)))

_digitval(c::Char) = '0' <= c <= '9' ? c - '0' : '₀' <= c <= '₉' ? c - '₀' : -1

# the name split into runs of text (Strings) and of digits (BigInts, negative after ₋)
function _name_tokens(s::AbstractString)
    toks = Any[]
    cs = collect(s)
    i = 1
    while i <= length(cs)
        neg = cs[i] == '₋' && i < length(cs) && _digitval(cs[i + 1]) >= 0
        if neg || _digitval(cs[i]) >= 0
            j = neg ? i + 1 : i
            v = big(0)
            while j <= length(cs) && _digitval(cs[j]) >= 0
                v = 10v + _digitval(cs[j])
                j += 1
            end
            push!(toks, neg ? -v : v)
            i = j
        else
            j = i
            while j <= length(cs) && _digitval(cs[j]) < 0 && !(cs[j] == '₋' && j < length(cs) && _digitval(cs[j + 1]) >= 0)
                j += 1
            end
            push!(toks, String(cs[i:(j - 1)]))
            i = j
        end
    end
    toks
end

function _natural_lt(a::AbstractString, b::AbstractString)
    ta, tb = _name_tokens(a), _name_tokens(b)
    for (x, y) in zip(ta, tb)
        x == y && continue
        x isa BigInt && y isa BigInt && return x < y
        x isa BigInt && return true                  # a number before text
        y isa BigInt && return false
        return x < y
    end
    length(ta) == length(tb) ? a < b : length(ta) < length(tb)
end
free_symbols(::Number) = Sym[]

# (internal, not exported: parameters are hidden until the name question is settled)
#     free_parameters(ex)
#
# The parameters appearing in `ex`, sorted by name.
function free_parameters(ex::Sym)
    s = Set{Sym}()
    _free!(s, ex, PARAMETER)
    _sort_by_name!(collect(s))
end
free_parameters(::Number) = Sym[]

function _free!(s, ex::Sym, k::Kind = SYMBOL)
    ex.kind == k && return push!(s, ex)
    if _isderivat(ex)                    # Derivative(f, v, a): v is bound
        t = _free!(Set{Sym}(), ex.args[1], k)
        delete!(t, ex.args[2])
        union!(s, t)
        return _free!(s, ex.args[3], k)
    end
    if _isdefint(ex)                     # Integral(f, v, a, b): v is bound
        t = _free!(Set{Sym}(), ex.args[1], k)
        delete!(t, ex.args[2])
        union!(s, t)
        _free!(s, ex.args[3], k)
        return _free!(s, ex.args[4], k)
    end
    for a in ex.args
        _free!(s, a, k)
    end
    for a in ex.exps
        _free!(s, a, k)
    end
    s
end

"""
    has(ex, x)

Does `x` occur as a subexpression of `ex`? (Matches whole nodes only.)
"""
has(ex::Sym, x::Sym) = _eq(ex, x) || any(a -> has(a, x), ex.args) || any(a -> has(a, x), ex.exps)

## substitution and calling

"""
    subs(ex, x => v, ...)
    subs(ex, d::AbstractDict)

Replace subexpressions of `ex` and re-simplify; the parts no substitution reaches are
kept as they are, not rebuilt. Expressions are callable, so `ex(x => 1, y => 2)` is the
same as `subs(ex, x => 1, y => 2)`; see [`Sym`](@ref) calling for the single-value forms.

A symbol or a function call is replaced where it occurs. A power, sum or product is
also matched algebraically, by identities that hold for every value:

- a power `b^P`: `b^q` with `q/P` an integer `n` is `w^n`, and an exponent with terms is
  split, `b^(nP + r) = w^n b^r`;
- a sum: where every one of its terms occurs the same number `m` of times, `m*w` plus
  the rest;
- a product: an integer power of it inside a product.

```julia
subs(x^4, x^2 => w)                       # w^2
subs(x^3, x^2 => w)                       # x^3 (3/2 is not an integer)
subs(exp(x + 1), exp(-x) => w)            # ℯ/w
subs(2a + 2b + x, a + b => w)             # 2w + x
subs(x^2 * log(x)^2, x * log(x) => w)     # w^2
```
"""
subs(ex::Sym, ps::Pair...) = _subs_top(ex, Dict{Sym,Sym}(Sym(p.first) => Sym(p.second) for p in ps))
subs(ex::Sym, d::AbstractDict) = _subs_top(ex, Dict{Sym,Sym}(Sym(k) => Sym(v) for (k, v) in d))
subs(ex::Number, args...) = Sym(ex)

# whole nodes, unless a target is a power, a sum or a product, also matched algebraically
function _subs_top(ex::Sym, d::Dict{Sym,Sym})
    ck = [k => v for (k, v) in d if _subs_compound(k)]
    isempty(ck) ? _subs(ex, d) : _subs_alg(ex, d, ck)
end
_subs_compound(k::Sym) = k.kind == POW || k.kind == ADD || (k.kind == MUL && !(length(k.args) == 1 && _isone(k.exps[1])))

# A node none of whose parts change is returned as it is (no rebuilding, so a
# substitution that misses most of a large expression allocates little).
function _subs(ex::Sym, d::Dict{Sym,Sym})
    v = get(d, ex, nothing)
    v === nothing || return v
    k = ex.kind
    if k == ADD
        new = _subs_args(ex.args, d)
        new === nothing && return ex
        return _add_many(Iterators.flatten(((Sym(ex.num),),
            (_scale(new[i], ex.coefs[i]) for i in eachindex(ex.args)))))
    elseif k == MUL
        powkeys = any(key -> key.kind == POW, keys(d))   # the dictionary is small
        # a factor b^e is stored as a pair, not a node: match it whole first (only
        # building the node when a power is being substituted)
        if powkeys && any(i -> !_isone(ex.exps[i]) && haskey(d, _mkpow(ex.args[i], ex.exps[i])), eachindex(ex.args))
            r = Sym(ex.num)
            for i in eachindex(ex.args)
                v = !_isone(ex.exps[i]) ? get(d, _mkpow(ex.args[i], ex.exps[i]), nothing) : nothing
                r = mul(r, v === nothing ? pow(_subs(ex.args[i], d), _subs(ex.exps[i], d)) : v)
            end
            return r
        end
        bs, es = _subs_args(ex.args, d), _subs_args(ex.exps, d)
        bs === nothing && es === nothing && return ex
        bs === nothing && (bs = ex.args)
        es === nothing && (es = ex.exps)
        r = Sym(ex.num)
        for i in eachindex(bs)
            r = mul(r, pow(bs[i], es[i]))
        end
        return r
    elseif k == POW
        b, e = _subs(ex.args[1], d), _subs(ex.args[2], d)
        b === ex.args[1] && e === ex.args[2] && return ex
        return pow(b, e)
    elseif k == FN
        ex.name === :Derivative && return _subs_derivative(ex, d)
        _isdefint(ex) && return _subs_defint(ex, d)
        new = _subs_args(ex.args, d)
        new === nothing && return ex
        return apply_function(ex.name, new)
    end
    ex
end
# the parts substituted, or nothing when none changes (allocating only from the first
# change on)
function _subs_args(args::Vector{Sym}, d::Dict{Sym,Sym})
    for i in eachindex(args)
        a = _subs(args[i], d)
        a === args[i] && continue
        new = args[1:(i - 1)]
        push!(new, a)
        for j in (i + 1):length(args)
            push!(new, _subs(args[j], d))
        end
        return new
    end
    nothing
end

## algebraic substitution of compound targets
#
# Each rule is an identity: (b^P)^n = b^(nP) for an integer n, so b^q with q/P an integer
# n is w^n; b^(nP + r) = (b^P)^n b^r, for an exponent with terms (e^(x + 1) with e^x => w
# is e*w, e^x with e^(-x) => w is 1/w); a sum K inside a sum, every term of K there m times
# (m a number), is m*w plus the rest (2a + 2b + x with a + b => w is 2w + x); a product K
# inside a product as K^n, n an integer, is w^n times the rest (x^2 log(x)^2 with
# x log(x) => w is w^2). The parts are substituted first, then the node.

function _subs_alg(ex::Sym, d::Dict{Sym,Sym}, ck)
    v = get(d, ex, nothing)
    v === nothing || return v
    k = ex.kind
    r = if k == ADD
        _add_many(Iterators.flatten(((Sym(ex.num),),
            (_scale(_subs_alg(a, d, ck), c) for (a, c) in zip(ex.args, ex.coefs)))))
    elseif k == MUL
        prod((pow(_subs_alg(b, d, ck), _subs_alg(e, d, ck)) for (b, e) in zip(ex.args, ex.exps)); init = Sym(ex.num))
    elseif k == POW
        pow(_subs_alg(ex.args[1], d, ck), _subs_alg(ex.args[2], d, ck))
    elseif k == FN
        ex.name === :Derivative || _isdefint(ex) ? _subs(ex, d) : apply_function(ex.name, Sym[_subs_alg(a, d, ck) for a in ex.args])
    else
        ex
    end
    for (K, w) in ck
        r = _subs_node(r, K, w)
    end
    r
end

# the node r with the target K => w matched algebraically, or r
function _subs_node(r::Sym, K::Sym, w::Sym)
    _eq(r, K) && return w
    if K.kind == POW
        r.kind == POW && return something(_subs_pow(r.args[1], r.args[2], K, w), r)
        if r.kind == MUL
            out, hit = Sym(r.num), false
            for (b, e) in zip(r.args, r.exps)
                t = _subs_pow(b, e, K, w)
                hit |= t !== nothing
                out = mul(out, t === nothing ? pow(b, e) : t)
            end
            return hit ? out : r
        end
    elseif K.kind == ADD && r.kind == ADD
        return something(_subs_sum(r, K, w), r)
    elseif K.kind == MUL && r.kind == MUL
        return something(_subs_prod(r, K, w), r)
    end
    r
end

# b^q with the target B^P => w: w^n b^(q - nP), or nothing
function _subs_pow(b::Sym, q::Sym, K::Sym, w::Sym)
    B, P = K.args[1], K.args[2]
    _eq(b, B) || return nothing
    n = q / P
    if n.kind == NUMBER
        isinteger(n) && !iszero(n) || return nothing
        return pow(w, n)
    end
    # an exponent with terms: n from the first term of P, the rest without P's terms
    q, P = expand(q), expand(P)
    tP, tq = _subs_terms(P), _subs_terms(q)
    isempty(tP) && return nothing
    c = get(tq, first(tP).first, nothing)
    c === nothing && return nothing
    n = c / first(tP).second
    (isinteger(n) && !iszero(n)) || return nothing
    rest = expand(q - n * P)
    any(t -> haskey(_subs_terms(rest), t), keys(tP)) && return nothing
    mul(pow(w, n), pow(b, rest))
end

# a sum's terms (constant left out) as term => numeric coefficient, in order
function _subs_terms(e::Sym)
    out = Dict{Sym,Sym}()
    e.kind == ADD && return Dict{Sym,Sym}(a => Sym(c) for (a, c) in zip(e.args, e.coefs))
    e.kind == NUMBER && return out
    e.kind == MUL && !isone(e.num) && return Dict{Sym,Sym}(e / Sym(e.num) => Sym(e.num))
    Dict{Sym,Sym}(e => ONE)
end

# the target sum K, every term of it m times in r: m*w + (r - m*K), or nothing
function _subs_sum(r::Sym, K::Sym, w::Sym)
    tr = Dict{Sym,Sym}(a => Sym(c) for (a, c) in zip(r.args, r.coefs))
    m = nothing
    for (a, c) in zip(K.args, K.coefs)
        cr = get(tr, a, nothing)
        cr === nothing && return nothing
        mi = cr / Sym(c)
        m === nothing ? (m = mi) : (_eq(m, mi) || return nothing)
    end
    m === nothing && return nothing
    add(mul(m, w), r - m * K)
end

# the target product K inside r as K^n (n an integer of one sign): w^n * r / K^n, or nothing
function _subs_prod(r::Sym, K::Sym, w::Sym)
    fr = Dict{Sym,Sym}(zip(r.args, r.exps))
    n = nothing
    for (b, e) in zip(K.args, K.exps)
        er = get(fr, b, nothing)
        er === nothing && return nothing
        q = er / e
        (q.kind == NUMBER && isexact(q.num)) || return nothing
        pos = q.num > N0
        ni = pos ? floor(q) : -floor(-q)          # toward 0
        iszero(ni) && return nothing
        if n === nothing
            n = ni
        else
            (n.num > N0) == pos || return nothing
            n = pos ? min(n, ni) : max(n, ni)
        end
    end
    n === nothing && return nothing
    mul(pow(w, n), r / K^n)
end

# a derivative's variable is bound: Derivative(f, v) with v => a symbol renames it, with
# v => a value gives the derivative at a point, Derivative(f, v, a); at a point, only f's
# other symbols and the point are substituted
_isderivat(ex::Sym) = ex.kind == FN && ex.name === :Derivative && length(ex.args) == 3

# a definite integral Integral(f, v, a, b), ∫ₐᵇ f dv (printed Integral(f, (v, a, b)), as
# integrate takes it): v is bound, so only f's other symbols and the limits are substituted
_isdefint(ex::Sym) = ex.kind == FN && ex.name === :Integral && length(ex.args) == 4
function _subs_defint(ex::Sym, d::Dict{Sym,Sym})
    f, v, a, b = ex.args
    df = haskey(d, v) ? Dict(k => w for (k, w) in d if !_eq(k, v)) : d
    _mkfn_node(:Integral, Sym[_subs(f, df), v, _subs(a, d), _subs(b, d)])
end
function _subs_derivative(ex::Sym, d::Dict{Sym,Sym})
    f, v = ex.args[1], ex.args[2]
    inner = haskey(d, v) ? delete!(copy(d), v) : d
    fs = _subs(f, inner)
    if length(ex.args) == 3
        return apply_function(:Derivative, Sym[fs, v, _subs(ex.args[3], d)])
    end
    w = get(d, v, nothing)
    w === nothing && return apply_function(:Derivative, Sym[fs, v])
    w.kind == SYMBOL && return apply_function(:Derivative, Sym[_subs(fs, Dict(v => w)), w])
    apply_function(:Derivative, Sym[fs, v, w])
end

## calling

# One method handles every call form, so dispatch has no ambiguities to resolve:
#   ex()            -> the value of a constant expression as a Julia number
#                      (`Sym(2)() === 2`, `Sym(π)() === π`); ex itself if it has a
#                      free symbol (nothing to substitute). Generic code (e.g.
#                      GruntzLimit) reads numeric leaves this way, and treats a
#                      non-number as a variable. This must not throw: a throw costs
#                      ~40 μs (the backtrace), paid for every symbol leaf read.
#   ex(x => a, ...) -> substitution, returns a Sym
#   ex(a::Sym)      -> substitution into the only variable, returns a Sym
#   ex(a)           -> numeric evaluation at the only variable, returns what Julia
#                      arithmetic on `a` returns (Int, Float64, BigFloat, a dual number, ...);
#                      with parameters present, substitution instead (`(a*x^2)(3) == 9a`)
#   ex(; x = a, ...) -> values for variables and parameters by name: numeric evaluation
#                      when every one is given a non-symbolic value, else substitution.
#                      Names not in ex are ignored; not mixed with the other forms.
function (ex::Sym)(args...; kws...)
    isempty(kws) || (isempty(args) ? (return _call_named(ex, values(kws))) : throw(ArgumentError(
        "call an expression with keywords (`ex(; x = 1)`) or with pairs or a value, not both")))
    isempty(args) && return _constant_value(ex)
    all(a -> a isa Pair, args) && return subs(ex, args...)
    length(args) == 1 || throw(ArgumentError(
        "call an expression with one value (it has at most one free symbol) or with pairs `x => value`"))
    a = only(args)
    a isa Sym || _hasparameter(ex) ? _subs1(ex, Sym(a)) : evaluate(ex, a)
end

# ex(; x = 1, a = 2): the expression's variables and parameters matched by name
function _call_named(ex::Sym, kw::NamedTuple)
    leaves = Dict{Symbol,Sym}()
    for s in Iterators.flatten((free_symbols(ex), free_parameters(ex)))
        s.name in keys(kw) || continue
        haskey(leaves, s.name) && throw(ArgumentError(
            "$(s.name) names both a variable and a parameter of $ex; call it with pairs instead"))
        leaves[s.name] = s
    end
    nleaves = length(free_symbols(ex)) + length(free_parameters(ex))
    vals = [kw[n] for n in keys(leaves)]
    if length(leaves) == nleaves && !any(v -> v isa Sym, vals)
        nt = NamedTuple{Tuple(keys(leaves))}(Tuple(vals))
        T = isempty(vals) ? Int : mapreduce(typeof, promote_type, vals)
        return _ev(ex, _Named(nt), _liftfloat(T), Ref{Union{Nothing,Sym}}(nothing))
    end
    subs(ex, (s => kw[n] for (n, s) in leaves)...)
end

# values for named leaves, for _ev
struct _Named{NT<:NamedTuple}
    vals::NT
end
_namedval(a::_Named, ex::Sym) = getfield(a.vals, ex.name)

function _constant_value(ex::Sym)
    _hassymbol(ex) && return ex
    evaluate(ex, 0)      # the argument is unused, and an Int keeps exact values exact
end

"""
    SymJ.unwrap_const(x)

The value of a constant expression as a Julia number, and anything else as it is: `ex()`
as a function, for `map` and other generic code (as `unwrap_const` in SymbolicUtils,
SymEngine and Giac). Exact values stay exact: `SymJ.unwrap_const(Sym(1)/2) === 1//2`,
`SymJ.unwrap_const(PI) === π`, `SymJ.unwrap_const(sqrt(Sym(2))) === sqrt(2)`; an expression
with a free symbol and a non-`Sym` argument are returned unchanged. Not exported.

```julia
map(SymJ.unwrap_const, [Sym(2), x + 1, Sym(3)/4])   # Any[2, x + 1, 3//4]
```
"""
unwrap_const(x) = x
unwrap_const(ex::Sym) = _constant_value(ex)

"""
    float(ex::Sym)

The value of a constant expression as a `Float64` (a `ComplexF64` when it is not real),
as [`N`](@ref); an expression with a free symbol is returned unchanged, so generic code
calling `float` works on it. For an array, broadcast: `float.(A)`. `float(Sym)` is `Sym`.
"""
Base.float(ex::Sym) = _hassymbol(ex) ? ex : N(ex)
Base.float(::Type{Sym}) = Sym

# does ex contain a symbol? (free_symbols without the allocation)
# does ex contain a variable or parameter (_hassymbol), or a parameter? Plain loops: these
# run before every numeric call of an expression.
_hassymbol(ex::Sym) = _hasleaf(ex, true)
_hasparameter(ex::Sym) = _hasleaf(ex, false)
function _hasleaf(ex::Sym, anynamed::Bool)
    (anynamed ? _isnamed(ex) : ex.kind == PARAMETER) && return true
    for a in ex.args
        _hasleaf(a, anynamed) && return true
    end
    for a in ex.exps
        _hasleaf(a, anynamed) && return true
    end
    false
end

function _subs1(ex::Sym, a::Sym)
    vars = free_symbols(ex)
    isempty(vars) && return ex
    length(vars) == 1 || _too_many_symbols(ex)
    subs(ex, only(vars) => a)
end

_too_many_symbols(ex) = throw(ArgumentError(
    "$ex has free symbols $(join(free_symbols(ex), ", ")); call it with pairs, e.g. ex(x => 1, y => 2)"))

## numeric evaluation with Julia arithmetic

"""
    evaluate(ex, a)

Evaluate `ex`, which may have at most one free symbol, at `a` using Julia arithmetic on
`a`. This is what calling an expression on a non-symbolic value does: if `u = x^2` then
`u(3) === 9`, `u(3.0) === 9.0` and `u(1//2) === 1//4`. Because the arithmetic is generic,
`a` may be any number type (big floats, dual numbers, intervals, ...).

- Exact coefficients are used as Julia integers and rationals, so exact inputs give
  exact results where Julia arithmetic would (`(x/3)(3) == 1//1`).
- When `a` is a floating point type (real or complex), constants such as `1/3`,
  `sqrt(2)` and `π` are computed in that type, so `BigFloat` inputs keep their precision.
- Negative integer powers use `inv` (`(1/x)(2) === 0.5`), and a negative base with a
  fractional power gives the principal complex value (`sqrt(x)(-4) == 2im`).
- Functions need a Julia method (registered with `julia = ...`); user functions without
  one throw an `ArgumentError`.
"""
function evaluate(ex::Sym, a)
    var = Ref{Union{Nothing,Sym}}(nothing)
    _ev(ex, a, _liftfloat(typeof(a)), var)
end

_liftfloat(::Type{T}) where {T<:AbstractFloat} = T
_liftfloat(::Type{Complex{T}}) where {T<:AbstractFloat} = T
_liftfloat(::Type) = nothing

function _jlnum(n::SymNumber)
    n.tag == T_FLOAT && return n.f
    n.tag == T_SMALL && return n.d == 1 ? n.n : n.n // n.d
    r = n.big::Rational{BigInt}
    isone(denominator(r)) ? numerator(r) : r
end

_leaf(n::SymNumber, F) = F === nothing ? _jlnum(n) : convert(F, _jlnum(n))
_constant(c, F) = F === nothing ? c : convert(F, c)

function _ev(ex::Sym, a, F, var)
    k = ex.kind
    k == NUMBER && return _leaf(ex.num, F)
    if k == SYMBOL
        a isa _Named && return _namedval(a, ex)
        if var[] === nothing
            var[] = ex
        elseif !_eq(var[]::Sym, ex)
            throw(ArgumentError("expression has more than one free symbol; call it with pairs, e.g. ex(x => 1, y => 2)"))
        end
        return a
    elseif k == PARAMETER
        a isa _Named && return _namedval(a, ex)
        throw(ArgumentError("cannot evaluate the parameter $ex numerically; substitute a value for it"))
    elseif k == CONSTANT
        ex.name === :pi && return _constant(π, F)
        ex.name === :E && return _constant(ℯ, F)
        throw(ArgumentError("no numeric value for constant $ex"))
    elseif k == ADD
        s = _leaf(ex.num, F)
        for i in eachindex(ex.args)
            t = _ev(ex.args[i], a, F, var)
            s += _isN1(ex.coefs[i]) ? t : _leaf(ex.coefs[i], F) * t
        end
        return s
    elseif k == MUL
        p = _leaf(ex.num, F)
        for i in eachindex(ex.args)
            p *= _evpow(ex.args[i], ex.exps[i], a, F, var)
        end
        return p
    elseif k == POW
        return _evpow(ex.args[1], ex.args[2], a, F, var)
    end
    rule = get(FUNCTION_RULES, ex.name, nothing)
    (rule === nothing || rule.julia === nothing) &&
        throw(ArgumentError("$(ex.name) has no numeric method; register one with `register_function!`"))
    rule.julia(map(x -> _intarg(x) ? _toint(x.num) : _ev(x, a, F, var), ex.args)...)
end

function _evpow(b::Sym, e::Sym, a, F, var)
    b === E && return exp(_ev(e, a, F, var))      # exp(x) is stored as ℯ^x
    bv = _ev(b, a, F, var)
    if e.kind == NUMBER && isexact(e.num)
        if _isint(e.num)
            k = _toint(e.num)
            return k >= 0 ? bv^k : inv(bv)^(-k)
        end
        if e.num == SymNumber(1 // 2)
            return bv isa Real && bv < 0 ? sqrt(complex(bv)) : sqrt(bv)
        end
        ev = convert(F === nothing ? Float64 : F, _jlnum(e.num))
    else
        ev = _ev(e, a, F, var)
    end
    bv isa Real && ev isa Real && bv < 0 && !isinteger(ev) && return complex(bv)^ev
    bv^ev
end

## expansion

"""
    expand(ex)

Distribute products and positive integer powers over sums. Expressions that are
polynomials with numeric coefficients in their symbols are expanded through [`Poly`](@ref).
"""
function expand(ex::Sym)
    if ex.kind == MUL || ex.kind == POW || ex.kind == ADD
        r = _expand_poly(ex)
        r === nothing || return r
    end
    _expand_generic(ex)
end

function _expand_poly(ex::Sym)
    gs = Set{Sym}()
    _numpoly_symbols!(gs, ex) || return nothing
    length(gs) > 16 && return nothing             # keep at least 4 bits per exponent
    g = sort!(collect(gs); lt = _lt)
    try
        return Sym(Poly{SymNumber}(ex, g))
    catch err
        err isa OverflowError || rethrow()
        return nothing                            # degree too large to pack
    end
end

# is ex a polynomial with numeric coefficients in its symbols? (collects the symbols)
function _numpoly_symbols!(gs::Set{Sym}, ex::Sym)
    k = ex.kind
    k == NUMBER && return true
    _isnamed(ex) && (push!(gs, ex); return true)
    k == ADD && return all(t -> _numpoly_symbols!(gs, t), ex.args)
    k == MUL && return all(i -> _nonneg_int(ex.exps[i]) && _numpoly_symbols!(gs, ex.args[i]), eachindex(ex.args))
    k == POW && return _nonneg_int(ex.args[2]) && _numpoly_symbols!(gs, ex.args[1])
    false
end

function _expand_generic(ex::Sym)
    k = ex.kind
    if k == ADD
        return _add_many(Iterators.flatten(((Sym(ex.num),),
            (_scale(expand(ex.args[i]), ex.coefs[i]) for i in eachindex(ex.args)))))
    elseif k == MUL
        r = Sym(ex.num)
        for i in eachindex(ex.args)
            r = _expand_mul(r, expand(_mkpow(ex.args[i], ex.exps[i])))
        end
        return r
    elseif k == POW
        b, e = expand(ex.args[1]), expand(ex.args[2])
        if b.kind == ADD && e.kind == NUMBER && _isint(e.num)
            n = _toint(e.num)
            n > 0 && return _expand_pow(b, n)
            return pow(_expand_pow(b, -n), MINUS_ONE)
        end
        return pow(b, e)
    elseif k == FN
        return apply_function(ex.name, Sym[expand(a) for a in ex.args])
    end
    ex
end
expand(x::Number) = Sym(x)

function _expand_mul(a::Sym, b::Sym)
    (a.kind == ADD || b.kind == ADD) || return mul(a, b)
    ta, ca = _termlist(a)
    tb, cb = _termlist(b)
    d = Dict{Sym,SymNumber}()
    sizehint!(d, length(ta) * length(tb))
    c = N0
    for i in eachindex(ta), j in eachindex(tb)
        c = _accum_sym!(d, c, mul(ta[i], tb[j]), ca[i] * cb[j])
    end
    _add_from_dict(c, d)
end

# a sum as parallel vectors of keys and coefficients (constant term keyed by ONE)
function _termlist(x::Sym)
    if x.kind == ADD
        iszero(x.num) && return (x.args, x.coefs)
        return (_insert(x.args, 1, ONE), _insert(x.coefs, 1, x.num))
    end
    x.kind == NUMBER && return (Sym[ONE], SymNumber[x.num])
    c, t = _split_coef(x)
    (Sym[t], SymNumber[c])
end

# (t_1 + ... + t_k)^n by the multinomial theorem. Prefix products are shared down the
# recursion and powers of each term are tabulated, so each monomial costs about one `mul`.
function _expand_pow(b::Sym, n::Int)
    n == 1 && return b
    ts, cs = _termlist(b)
    k = length(ts)
    k == 1 && return _scale(pow(ts[1], Sym(n)), cs[1]^n)
    tpow = [Sym[pow(ts[i], Sym(a)) for a in 1:n] for i in 1:k]
    cpow = [SymNumber[cs[i]^a for a in 1:n] for i in 1:k]
    d = Dict{Sym,SymNumber}()
    sizehint!(d, binomial(n + k - 1, k - 1))
    c = _multinomial!(d, N0, tpow, cpow, 1, n, N1, ONE)
    _add_from_dict(c, d)
end

function _multinomial!(d, c, tpow, cpow, i, rem, coef, mono)
    if i == length(tpow) || rem == 0
        rem > 0 && (coef *= cpow[i][rem]; mono = mul(mono, tpow[i][rem]))
        return _accum_sym!(d, c, mono, coef)
    end
    for a in 0:rem
        ca = a == 0 ? coef : coef * _binom(rem, a) * cpow[i][a]
        ma = a == 0 ? mono : mul(mono, tpow[i][a])
        c = _multinomial!(d, c, tpow, cpow, i + 1, rem - a, ca, ma)
    end
    c
end

_binom(n::Int, k::Int) = n <= 60 ? SymNumber(binomial(n, k)) : SymNumber(binomial(big(n), k))

## numeric evaluation

"""
    N(ex)

Evaluate `ex` numerically as a `Float64` (or `ComplexF64` when the value is not real).
"""
function N(ex::Sym)
    z = _evalc(ex)
    iszero(imag(z)) ? real(z) : z
end
N(x::Number) = x

function _evalc(ex::Sym)::ComplexF64
    k = ex.kind
    k == NUMBER && return Float64(ex.num)
    if k == CONSTANT
        ex.name === :pi && return π
        ex.name === :E && return ℯ
        throw(ArgumentError("no numeric value for constant $ex"))
    end
    _isnamed(ex) && throw(ArgumentError("cannot evaluate free symbol $ex numerically"))
    if k == ADD
        s = complex(Float64(ex.num))
        for i in eachindex(ex.args)
            s += Float64(ex.coefs[i]) * _evalc(ex.args[i])
        end
        return s
    elseif k == MUL
        p = complex(Float64(ex.num))
        for i in eachindex(ex.args)
            p *= _cpow(ex.args[i], ex.exps[i])
        end
        return p
    elseif k == POW
        return _cpow(ex.args[1], ex.args[2])
    end
    ex.name === :piecewise && return _evalc_piecewise(ex)
    rule = get(FUNCTION_RULES, ex.name, nothing)
    (rule === nothing || rule.julia === nothing) && throw(ArgumentError("no numeric method for $(ex.name)"))
    # an exact integer argument stays an Int (an order, a kind: besselh(ν, 1, x))
    any(_intarg, ex.args) && return _fn_with_ints(rule.julia, ex.args)
    vals = map(_evalc, ex.args)
    # real arguments: call the real method (special functions often have no complex
    # one); fall back to complex for a value outside the real domain (log(-1.0))
    if all(v -> iszero(imag(v)), vals)
        try
            return rule.julia(map(real, vals)...)
        catch err
            err isa DomainError || rethrow()
        end
    end
    rule.julia(vals...)
end

# An exact integer given directly to a function stays an Int when evaluated, as in
# Julia: besselj(2, x) passes 2, and some parameters must be integers (polygamma's
# order, besselh's kind).
_intarg(a::Sym) = a.kind == NUMBER && _isint(a.num) && isexact(a.num) && a.num.tag == T_SMALL

function _fn_with_ints(f, args)
    vals = Any[_intarg(a) ? _toint(a.num) : _evalc(a) for a in args]
    if all(v -> v isa Int || iszero(imag(v)), vals)
        try
            return f(map(v -> v isa Int ? v : real(v), vals)...)
        catch err
            err isa DomainError || rethrow()
        end
    end
    try
        return f(vals...)
    catch err
        err isa DomainError || rethrow()
    end
    # outside the real domain with the integers kept (log(-3)): complex values throughout
    f(map(v -> v isa Int ? complex(float(v)) : v, vals)...)
end

function _cpow(b::Sym, e::Sym)
    z = _evalc(b)
    e.kind == NUMBER && _isint(e.num) && return z^_toint(e.num)
    w = _evalc(e)
    iszero(imag(z)) && real(z) >= 0 && iszero(imag(w)) && return complex(real(z)^real(w))
    z^w
end

function Base.Float64(ex::Sym)
    v = N(ex)
    v isa Real || throw(InexactError(:Float64, Float64, ex))
    v
end

## terms and coefficients

# a term c x₁^p₁ ⋯ x_m^p_m with c free of the xs: (p, c), or nothing when the term
# involves an x otherwise (sin(x), (x + 1)^(1/2)); the powers may be any expressions free
# of the xs
function _monomial(t::Sym, xs::Vector{Sym})
    ps = Sym[ZERO for _ in xs]
    any(x -> has(t, x), xs) || return ps, t
    if t.kind == MUL
        c = Sym(t.num)
        for (b, p) in zip(t.args, t.exps)
            i = findfirst(x -> _eq(b, x), xs)
            if i !== nothing && !any(x -> has(p, x), xs)
                ps[i] += p
            elseif any(x -> has(b, x) || has(p, x), xs)
                return nothing
            else
                c *= b^p
            end
        end
        return ps, c
    end
    b, p = t.kind == POW ? (t.args[1], t.args[2]) : (t, ONE)
    i = findfirst(x -> _eq(b, x), xs)
    (i === nothing || any(x -> has(p, x), xs)) && return nothing
    ps[i] = p
    ps, ONE
end
_terms(e::Sym) = e.kind == ADD ? [Sym(c) * a for (c, a) in zip(e.coefs, e.args)] : [e]

# x^n given as the variable: (x, n)
function _coeff_var(x::Sym, n)
    x.kind == POW && (_isnamed(x.args[1]) || x.args[1].kind == FN) && return x.args[1], x.args[2] * n
    x, Sym(n)
end

"""
    SymJ.coeff(ex, x, n = 1)
    SymJ.coeff(ex, x^n)

The coefficient of `x^n` in `ex`, after expanding: the sum of the factors that multiply
`x^n` exactly in its terms, so `coeff(a*x^2 + b*x + c, x)` is `b` and `n = 0` gives the
terms free of `x`. Terms with `x` in some other way (`sin(x)`) are not counted. `x` may be
a symbol or a function call (`y(t)`), and `n` any power, `-1` or `1/2` too. Not exported.

```julia
@syms a b c x y
SymJ.coeff(a*x^2 + b*x + c, x, 2)           # a
SymJ.coeff((x + y)^3, x^2)                  # 3y
SymJ.coeff(a*x + sin(x) + b, x, 0)          # b
```
"""
function coeff(ex, x::Sym, n = 1)
    v, p = _coeff_var(x, n)
    e = expand(Sym(ex))
    out = ZERO
    for t in _terms(e)
        m = _monomial(t, [v])
        m === nothing && continue
        _eq(expand(m[1][1] - p), ZERO) && (out += m[2])
    end
    e.kind == ADD && !iszero(e.num) && iszero(p) && (out += Sym(e.num))
    out
end

"""
    SymJ.gather(ex, x)
    SymJ.gather(ex, [x, y, ...])

`ex` expanded and its terms gathered by powers of `x` (by monomials in `x, y, ...` for
several): `a*x^2 + b*x^2 + a*x - b*x + c` as `(a + b)*x^2 + (a - b)*x + c`. Terms with `x`
in some other way (`sin(x)`) are left as they are. (SymPy's `collect`; Julia's `collect`
is another function.) Not exported.

```julia
@syms a b c x y
SymJ.gather(a*x^2 + b*x^2 + a*x - b*x + c, x)    # (a + b)*x^2 + (a - b)*x + c
SymJ.gather(x*y + a*x*y + x + b*y, [x, y])       # (a + 1)*x*y + x + b*y
```
"""
gather(ex, x::Sym) = gather(ex, [x])
function gather(ex, xs::AbstractVector)
    vs = Sym[Sym(x) for x in xs]
    all(v -> _isnamed(v) || v.kind == FN, vs) ||
        throw(ArgumentError("gather: $xs are not symbols or function calls"))
    e = expand(Sym(ex))
    groups = Dict{Vector{Sym},Sym}()
    order = Vector{Sym}[]
    rest = e.kind == ADD ? Sym(e.num) : ZERO
    for t in _terms(e)
        m = _monomial(t, vs)
        if m === nothing
            rest += t
            continue
        end
        ps, c = m
        haskey(groups, ps) || push!(order, ps)
        groups[ps] = get(groups, ps, ZERO) + c
    end
    out = rest
    for ps in order
        out += groups[ps] * prod((v^p for (v, p) in zip(vs, ps)); init = ONE)
    end
    out
end

## Horner form

"""
    SymJ.horner(p, x)
    SymJ.horner(p, [x, y, ...])
    SymJ.horner(p)

The polynomial `p` in nested (Horner) form in `x`: `2x^3 - x + 5` as `(2x^2 - 1)x + 5`
(`((2x)x - 1)x + 5`, a run of zero coefficients as one power), so it evaluates with
the fewest multiplications. With several variables (all of `p`'s, in order, when none are
given), each coefficient is put in Horner form in the next ones. The coefficients may be
any expressions free of the variables. Julia's `evalpoly` evaluates numbers this way.
Not exported.

```julia
@syms x y
SymJ.horner(2x^3 - x + 5, x)                  # (2x^2 - 1)*x + 5
SymJ.horner(x^2*y + x*y^2 + x + 1)            # x*(x*y + y^2 + 1) + 1
```
"""
horner(p, x::Sym) = horner(p, [x])
horner(p) = horner(p, free_symbols(Sym(p)))
function horner(p, xs::AbstractVector)
    e = expand(Sym(p))
    isempty(xs) && return e
    x = Sym(first(xs))
    _isnamed(x) || throw(ArgumentError("horner: $x is not a symbol"))
    has(e, x) || return horner(e, xs[2:end])
    P = try Poly{Sym}(e, [x]) catch
        throw(ArgumentError("horner: $p is not a polynomial in $x"))
    end
    cs = Dict{Int,Sym}(ev[1] => c for (ev, c) in zip(exponent_vectors(P), coefficients(P)))
    n = maximum(keys(cs))
    h = horner(cs[n], xs[2:end])
    for i in (n - 1):-1:0
        c = get(cs, i, ZERO)
        h = h * x + (iszero(c) ? ZERO : horner(c, xs[2:end]))
    end
    h
end

## floats as rationals

"""
    rationalize(ex::Sym; tol = nothing)
    rationalize(T, ex::Sym; tol = nothing)

`ex` with each floating-point number in it replaced by a nearby rational
(`rationalize(T, f; tol)` from Julia, each float's own `eps` when `tol` is not given):
`rationalize(0.5x + 1.25)` is `x/2 + 5/4`. Exact numbers are left alone. The result has
only rationals, not recognized constants such as `√2` or `π`.
"""
Base.rationalize(ex::Sym; tol = nothing) = rationalize(Int, ex; tol)
Base.rationalize(::Type{T}, ex::Sym; tol = nothing) where {T<:Integer} = _rationalize(ex, T, tol)

_rationalize(e::Sym, T, tol) = _mapfloats(e) do f
    SymNumber(tol === nothing ? rationalize(T, f) : rationalize(T, f; tol))
end

# e with each finite float f (coefficients and exponents too) replaced by r(f)
function _mapfloats(r, e::Sym)
    n(x::SymNumber) = isfloat(x) && isfinite(Float64(x)) ? Sym(r(Float64(x))) : Sym(x)
    R(a::Sym) = _mapfloats(r, a)
    k = e.kind
    k == NUMBER && return n(e.num)
    isempty(e.args) && return e
    k == ADD && return n(e.num) + _add_many(n(c) * R(a) for (c, a) in zip(e.coefs, e.args))
    k == MUL && return n(e.num) * prod((R(b)^R(p) for (b, p) in zip(e.args, e.exps)); init = ONE)
    k == POW && return R(e.args[1])^R(e.args[2])
    _rebuild_node(e, Sym[R(a) for a in e.args])
end

## convenience constructors

"""
    @syms x y f()
    @syms xs[1:10] A[0:2, 0:2] fs[1:3]()

Create symbols (and, with `f()`, a [`SymFunction`](@ref)) and bind them to the names.
`xs[r₁, r₂, …]`, with ranges (evaluated where the macro is called, so `xs[1:n]` works),
binds an array of symbols named by their indices: `@syms xs[0:3]` gives the vector
`[xs₀, xs₁, xs₂, xs₃]` (indexed from 1 as usual; the offset is in the names), and
`A[1:2, 1:2]` a matrix of `A₁_₁, A₁_₂, …`. `fs[1:3]()` gives functions. The names may
also be separated by commas, `@syms x, y`.
"""
macro syms(xs...)
    # commas: a single tuple of declarations
    length(xs) == 1 && xs[1] isa Expr && xs[1].head === :tuple && (xs = Tuple(xs[1].args))
    ex = Expr(:block)
    names = Symbol[]
    for x in xs
        if x isa Symbol
            push!(ex.args, :($(esc(x)) = Sym($(QuoteNode(x)))))
            push!(names, x)
        elseif x isa Expr && x.head === :call && length(x.args) == 1 && x.args[1] isa Symbol
            f = x.args[1]
            push!(ex.args, :($(esc(f)) = SymFunction($(QuoteNode(f)))))
            push!(names, f)
        elseif x isa Expr && x.head === :ref && x.args[1] isa Symbol && length(x.args) >= 2
            push!(ex.args, _syms_array(x, false))
            push!(names, x.args[1])
        elseif x isa Expr && x.head === :call && length(x.args) == 1 && x.args[1] isa Expr &&
               x.args[1].head === :ref && x.args[1].args[1] isa Symbol && length(x.args[1].args) >= 2
            push!(ex.args, _syms_array(x.args[1], true))
            push!(names, x.args[1].args[1])
        else
            throw(ArgumentError("@syms: unsupported argument $x (use x, f(), xs[1:n] or fs[1:n]())"))
        end
    end
    push!(ex.args, Expr(:tuple, map(esc, names)...))
    ex
end

# xs[r...] (fs[r...]()): the array built at run time from the ranges
_syms_array(x::Expr, fn::Bool) =
    :($(esc(x.args[1])) = $(_symarray)($(QuoteNode(x.args[1])), ($(map(esc, x.args[2:end])...),), $fn))

const _SUBSCRIPTS = Dict{Char,Char}('-' => '₋', '0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄',
                                    '5' => '₅', '6' => '₆', '7' => '₇', '8' => '₈', '9' => '₉')
_subscript(i::Integer) = join(_SUBSCRIPTS[c] for c in string(i))

function _symarray(base::Symbol, ranges::Tuple, fn::Bool)
    all(r -> r isa AbstractRange{<:Integer}, ranges) ||
        throw(ArgumentError("@syms: $base[...] takes integer ranges, such as $base[1:5]"))
    map(Iterators.product(ranges...)) do ind
        name = Symbol(base, join(map(_subscript, ind), "_"))
        fn ? SymFunction(name) : Sym(name)
    end
end

# (internal, not exported: parameters are hidden until the name question is settled)
#     @parameters a b
#
# Create symbolic parameters (see [`parameter`](@ref)) and bind them to the names.
macro parameters(xs...)
    all(x -> x isa Symbol, xs) || throw(ArgumentError("@parameters: names only"))
    ex = Expr(:block)
    for x in xs
        push!(ex.args, :($(esc(x)) = parameter($(QuoteNode(x)))))
    end
    push!(ex.args, Expr(:tuple, map(esc, xs)...))
    ex
end

"""
    symbols("x y z")

A tuple of symbols named by a space- or comma-separated string.
"""
symbols(s::AbstractString) = Tuple(Sym(Symbol(n)) for n in split(s, r"[\s,]+"; keepempty = false))
