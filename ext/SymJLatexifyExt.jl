# LaTeX output through Latexify: `latexify(ex)` for expressions and relations, and a
# `text/latex` display for SymJ's results, so notebooks and the documentation show them
# as typeset mathematics: expressions, relations, exact numbers, polynomials, series,
# arrays of these, solutions of systems (Dicts and vectors of them), sets of intervals,
# regions of inequality systems, cases, and tuples of any of them. Julia's own numbers
# are left alone.
#
# As SymPyCore does, an expression's printed form (valid Julia) is handed to Latexify,
# which typesets it; relations, printed with ≪ and friends, are rebuilt with Julia's
# operators first. A display that cannot be typeset falls back to the printed form.
module SymJLatexifyExt

using SymJ
using SymJ: Relation, SymSet, SymInterval, CylindricalSet, Cases, SymSeries, Membership,
            _NotFact, NumberSet
using Latexify

# the printed form as a Julia expression, with ℯ and im as the letters e and i, and big
# integer literals (parsed as @int128_str or @big_str calls) as integers; a piecewise
# function as an if chain (typeset as cases), heaviside and dirac as θ and δ
_expr(x) = _letters(Meta.parse(string(x)))
_letters(s::Symbol) = s === :ℯ ? :e : s === :im ? :i : s
function _letters(e::Expr)
    e.head === :macrocall && e.args[end] isa String && return parse(BigInt, e.args[end])
    if e.head === :call && e.args[1] === :piecewise
        return _ifchain(map(_letters, e.args[2:end]))
    elseif e.head === :call && e.args[1] === :DifferenceRoot
        # the recurrence and its initial values as equations
        eq(a) = a isa Expr && a.head === :call && a.args[1] in (:~, :(=>)) ?
            Expr(:call, :(==), map(_letters, a.args[2:end])...) : _letters(a)
        # (built here: Latexify reads a leading a == b argument as a keyword)
        parts = [Latexify.latexraw(eq(a)) for a in e.args[2:end]]
        return Expr(:call, :DifferenceRoot, Latexify.LaTeXString(join(parts, raw",\; ")))
    elseif e.head === :call && e.args[1] in (:heaviside, :dirac)
        return Expr(:call, e.args[1] === :heaviside ? :θ : :δ, map(_letters, e.args[2:end])...)
    elseif e.head === :call && haskey(_RELOPS, e.args[1])
        return Expr(:call, _RELOPS[e.args[1]], map(_letters, e.args[2:end])...)
    end
    Expr(e.head, map(_letters, e.args)...)
end
const _RELOPS = Dict(:≪ => :<, :≦ => :≤, :≫ => :>, :≧ => :≥, :~ => :(==), :≶ => :≠)

# pieces c => v (Exprs), the last with condition true
function _ifchain(pairs, head = :if)
    c, v = pairs[1].args[2], pairs[1].args[3]
    rest = pairs[2:end]
    c === true && return v
    Expr(head, _texcondexpr(c), v, isempty(rest) ? Latexify.LaTeXString("\\text{undefined}") : _ifchain(rest, :elseif))
end
# a conjunction: a lower and an upper bound chained, a ≤ x < b, else the conditions
# with ∧
function _texcondexpr(c)
    c isa Expr && c.head === :tuple || return c
    if length(c.args) == 2 && all(a -> a isa Expr && a.head === :call, c.args)
        (l, u) = c.args
        flip = Dict(:> => :<, :≥ => :≤)
        if l.args[1] in keys(flip) && u.args[1] in (:<, :≤) && l.args[2] == u.args[2]
            return Expr(:comparison, l.args[3], flip[l.args[1]], l.args[2], u.args[1], u.args[3])
        end
    end
    Expr(:call, :∧, c.args...)
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

## the LaTeX of each kind of result, without delimiters

_raw(x) = String(latexify(x; env = :raw))
_tex(x::Union{Sym,Relation}) = _raw(x)
_tex(n::SymNumber) = _raw(Sym(n))
_tex(p::Poly) = _raw(Sym(p))
_tex(s::SymSeries) = _raw(_expr(s))
# a vector as a bracketed list, [-2, 2] (a matrix as an array)
_tex(x::AbstractVector{Sym}) = isempty(x) ? "\\left[ \\; \\right]" : "\\left[ " * join(map(_raw, x), ", \\; ") * " \\right]"
_tex(x::AbstractArray{Sym}) = isempty(x) ? "\\left[ \\; \\right]" : _raw(x)
_tex(x::AbstractArray{SymNumber}) = _tex(map(Sym, x))

# a solution of a system: {x = 2, y = 1}, in the order of the unknowns' names
function _tex(d::AbstractDict{Sym,Sym})
    ks = sort!(collect(keys(d)); by = string)
    "\\left\\{ " * join((_raw(k) * " = " * _raw(d[k]) for k in ks), ", \\; ") * " \\right\\}"
end
_tex(v::AbstractVector{<:AbstractDict{Sym,Sym}}) = isempty(v) ? "\\left[ \\; \\right]" :
    "\\left[ \\begin{array}{l} " * join(map(_tex, v), " \\\\ ") * " \\end{array} \\right]"

# (-∞, -1) ∪ [2, 3] ∪ {5}
_end(e, lower) = e === nothing ? (lower ? "-\\infty" : "\\infty") : _raw(e)
function _tex(I::SymInterval)
    I.lo !== nothing && I.hi !== nothing && I.lo == I.hi && return "\\left\\{ " * _raw(I.lo) * " \\right\\}"
    (I.lclosed ? "\\left[ " : "\\left( ") * _end(I.lo, true) * ", " * _end(I.hi, false) *
        (I.rclosed ? " \\right]" : " \\right)")
end
_tex(S::SymSet) = isempty(S) ? "\\emptyset" : join(map(_tex, S.intervals), " \\cup ")

# the region of a system of inequalities, as it prints: each unknown between bounds
function _tex(C::CylindricalSet)
    C.feasible || return "\\emptyset"
    parts = String[]
    for (x, lo, hi) in reverse(C.levels)
        uniform(bs) = isempty(bs) || all(b -> b[2] == bs[1][2], bs)
        if uniform(lo) && uniform(hi)
            side(bs, lower) = isempty(bs) ? (lower ? "-\\infty" : "\\infty") : length(bs) == 1 ? _raw(bs[1][1]) :
                (lower ? "\\max" : "\\min") * "\\left( " * join((_raw(b[1]) for b in bs), ", ") * " \\right)"
            push!(parts, _raw(x) * " \\in " * (isempty(lo) || lo[1][2] ? "\\left( " : "\\left[ ") * side(lo, true) *
                         ", " * side(hi, false) * (isempty(hi) || hi[1][2] ? " \\right)" : " \\right]"))
        else
            append!(parts, [_raw(b[1]) * (b[2] ? " < " : " \\leq ") * _raw(x) for b in lo])
            append!(parts, [_raw(x) * (b[2] ? " < " : " \\leq ") * _raw(b[1]) for b in hi])
        end
    end
    append!(parts, [_raw(v) * " = " * _raw(e) for (v, e) in C.equalities])
    join(parts, ", \\; ")
end

# cases: the answers, each with its conditions
const _SETTEX = ("\\mathbb{Z}", "\\mathbb{Q}", "\\mathbb{R}", "\\mathbb{C}")
_texcond(c::Relation) = _raw(c)
_texcond(c::_NotFact) = "\\lnot \\left( " * _texcond(c.f) * " \\right)"
_texcond(m::Membership) = _raw(m.x) * (m.negated ? " \\notin " : " \\in ") * _SETTEX[m.set.rank]
_texcond(c) = "\\text{" * string(c) * "}"
function _texconds(conds)
    isempty(conds) && return "\\text{always}"
    parts, used = String[], falses(length(conds))
    for (i, c) in enumerate(conds)              # a lower and an upper bound chained
        used[i] && continue
        used[i] = true
        if c isa Relation && c.op in (:gt, :ge)
            j = findfirst(k -> !used[k] && conds[k] isa Relation && conds[k].op in (:lt, :le) &&
                               isequal(conds[k].lhs, c.lhs), eachindex(conds))
            if j !== nothing
                used[j] = true
                rel(op) = op in (:lt, :gt) ? " < " : " \\leq "
                push!(parts, _raw(c.rhs) * rel(c.op) * _raw(c.lhs) * rel(conds[j].op) * _raw(conds[j].rhs))
                continue
            end
        end
        push!(parts, _texcond(c))
    end
    join(parts, " \\land ")
end
# an error case shows its message (the series diverges), not the exception's type
_errmsg(v) = hasproperty(v, :msg) && v.msg isa AbstractString ? v.msg : sprint(showerror, v)
_texanswer(v) = v isa Exception ? "\\text{" * replace(_errmsg(v), r"[\\{}_^#%&$~]" => "") * "}" :
    applicable(_tex, v) ? _tex(v) : "\\text{" * string(v) * "}"
_tex(cs::Cases) = "\\begin{cases} " *
    join((_texanswer(v) * " & " * _texconds(conds) for (conds, v) in cs.cases), " \\\\ ") * " \\end{cases}"

const _TeX = Union{Sym,Relation,SymNumber,Poly,SymSeries,AbstractArray{Sym},AbstractArray{SymNumber},
                   AbstractDict{Sym,Sym},AbstractVector{<:AbstractDict{Sym,Sym}},SymSet,CylindricalSet,Cases}
_tex(t::Tuple) = "\\left( " * join(map(_tex, t), ", \\; ") * " \\right)"

## display

# falls back to the printed form (in \mathtt), so showing a result never fails
function _latex(x)
    try
        "\$" * _tex(x) * "\$"
    catch
        "\$\\mathtt{" * replace(string(x), r"([\\{}_^#%&$~])" => s"\\\1") * "}\$"
    end
end

Base.show(io::IO, ::MIME"text/latex", x::_TeX) = print(io, _latex(x))
Base.show(io::IO, ::MIME"text/latex", t::Tuple{_TeX,Vararg{_TeX}}) = print(io, _latex(t))

end
