# TermInterface v2 contract. Sums, products, powers and function applications are call
# expressions; numbers, symbols and constants are leaves.

TermInterface.isexpr(x::Sym) = x.kind == ADD || x.kind == MUL || x.kind == POW || x.kind == FN
TermInterface.iscall(x::Sym) = isexpr(x)

function TermInterface.operation(x::Sym)
    k = x.kind
    k == ADD && return (+)
    k == MUL && return (*)
    k == POW && return (^)
    k == FN && return _fn_operation(x.name)
    throw(ArgumentError("$x is not a call expression"))
end
TermInterface.head(x::Sym) = operation(x)

# The operation of a function whose numeric method has another name (_num_factorial for
# factorial): a function of the function's own name, so that a consumer matching by name
# (GruntzLimit's factorial, lambertw; SymPy #19379, #18992, #18969) finds it. On Syms it
# builds the node, on numbers it evaluates.
module _NamedOps
import ..SymJ
for f in (:factorial, :binomial, :lambertw, :harmonic, :li, :heaviside, :dirac)
    @eval function $f(args...)
        any(a -> a isa SymJ.Sym, args) && return SymJ.apply_function($(QuoteNode(f)), SymJ.Sym[SymJ.Sym(a) for a in args])
        SymJ.FUNCTION_RULES[$(QuoteNode(f))].julia(args...)
    end
end
end
const _NAMED_OPS = Dict{Symbol,Function}(f => getfield(_NamedOps, f) for f in
    (:factorial, :binomial, :lambertw, :harmonic, :li, :heaviside, :dirac))

# For POW and FN the node's own vector is returned; treat it as read-only.
function TermInterface.arguments(x::Sym)
    k = x.kind
    if k == ADD
        out = Sym[_with_coef(x.coefs[i], x.args[i]) for i in eachindex(x.args)]
        iszero(x.num) || push!(out, Sym(x.num))
        return out
    elseif k == MUL
        out = Sym[_mkpow(x.args[i], x.exps[i]) for i in eachindex(x.args)]
        _isN1(x.num) || pushfirst!(out, Sym(x.num))
        return out
    end
    (k == POW || k == FN) && return x.args
    throw(ArgumentError("$x is not a call expression"))
end
TermInterface.children(x::Sym) = arguments(x)

# Arguments in display order: deterministic and independent of the hash-based storage.
function TermInterface.sorted_arguments(x::Sym)
    args = arguments(x)
    x.kind == ADD || x.kind == MUL ? _display_sort(args) : args
end
TermInterface.sorted_children(x::Sym) = sorted_arguments(x)

function TermInterface.maketerm(::Type{<:Sym}, head, args, metadata)
    head === (+) && return _add_many(args)
    head === (*) && return foldl(mul, (Sym(a) for a in args); init = ONE)
    Sym(head(map(Sym, args)...))
end
