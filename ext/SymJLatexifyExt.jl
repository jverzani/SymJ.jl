# LaTeX output through Latexify: `latexify(ex)` for expressions, relations and arrays of
# expressions, and a `text/latex` display for them, so notebooks and the documentation
# show results as typeset mathematics. As SymPyCore does, the expression's printed form
# (valid Julia) is handed to Latexify, which typesets it; relations, printed with ≪ and
# friends, are rebuilt with Julia's operators first.
module SymJLatexifyExt

using SymJ
using SymJ: Relation
using Latexify

# the printed form as a Julia expression, with ℯ and im as the letters e and i, and big
# integer literals (parsed as @int128_str or @big_str calls) as integers
_expr(x::Sym) = _letters(Meta.parse(string(x)))
_letters(s::Symbol) = s === :ℯ ? :e : s === :im ? :i : s
function _letters(e::Expr)
    e.head === :macrocall && e.args[end] isa String && return parse(BigInt, e.args[end])
    Expr(e.head, map(_letters, e.args)...)
end
_letters(v) = v

const _OPS = Dict(:eq => :(=), :lt => :<, :le => :≤, :gt => :>, :ge => :≥, :ne => :≠)

Latexify.@latexrecipe function _(x::Sym)
    return _expr(x)
end

Latexify.@latexrecipe function _(r::Relation)
    op = _OPS[r.op]
    return op === :(=) ? Expr(:(=), _expr(r.lhs), _expr(r.rhs)) : Expr(:call, op, _expr(r.lhs), _expr(r.rhs))
end

# Each display falls back to the printed form (in \mathtt) if Latexify cannot typeset
# it, so showing a result never fails.
function _latex(f, x)
    try
        f()
    catch
        "\$\\mathtt{" * replace(string(x), r"([\\{}_^#%&$~])" => s"\\\1") * "}\$"
    end
end

Base.show(io::IO, ::MIME"text/latex", x::Union{Sym,Relation}) = print(io, _latex(() -> latexify(x), x))
function Base.show(io::IO, ::MIME"text/latex", x::AbstractArray{Sym})
    # inline (\$…\$) rather than an equation environment, which not every renderer takes
    isempty(x) && return print(io, "\$\\left[ \\; \\right]\$")
    print(io, _latex(() -> latexify(x; env = :inline), x))
end
# a tuple of expressions (several results at once)
function Base.show(io::IO, ::MIME"text/latex", t::Tuple{Union{Sym,Relation},Vararg{Union{Sym,Relation}}})
    print(io, _latex(() -> "\$\\left( " * join((latexify(u; env = :raw) for u in t), ", \\; ") * " \\right)\$", t))
end

end
