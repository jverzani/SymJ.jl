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
            _NotFact, NumberSet, SymNumber, SolutionSet, ImageSet, ConditionSet
using Latexify

# the printed form as a Julia expression, with ℯ and im as the letters e and i, and big
# integer literals (parsed as @int128_str or @big_str calls) as integers; a piecewise
# function as an if chain (typeset as cases), heaviside and dirac as θ and δ
_expr(x) = _letters(Meta.parse(string(x)))
_letters(s::Symbol) = s === :ℯ ? :e : s === :im ? :i : _greek_name(s)
# a name spelled with a Greek letter's name, or with several subscripts: beta_13_2 as
# \beta_{13,2} (it was beta_{13\_2}; SymPy #4576)
const _GREEK = Set(["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota",
    "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi",
    "chi", "psi", "omega", "Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi", "Sigma", "Upsilon",
    "Phi", "Psi", "Omega"])
function _greek_name(s::Symbol)
    str = String(s)
    # an array element from @syms A[1:2, 1:2]: A₁_₂ as A_{1,2}, x₋₁ as x_{-1} (Latexify made
    # each subscript character _1, A_1__{2}, which is not LaTeX and stopped the display)
    m = match(r"^([^₀-₉₋]+?)([₀-₉₋]+(?:_[₀-₉₋]+)*)$", str)
    if m !== nothing
        idx = [map(c -> c == '₋' ? '-' : Char('0' + (c - '₀')), p) for p in split(m.captures[2], '_')]
        b = String(m.captures[1])
        base = b in _GREEK ? (b == "omicron" ? "o" : "\\" * b) : Latexify.latexraw(Symbol(b))
        return Latexify.LaTeXString(base * "_{" * join(idx, ",") * "}")
    end
    parts = split(str, '_')
    (any(isempty, parts) || !(parts[1] in _GREEK || length(parts) > 2)) && return s
    base = parts[1] in _GREEK ? (parts[1] == "omicron" ? "o" : "\\" * parts[1]) : parts[1]
    Latexify.LaTeXString(length(parts) == 1 ? base : base * "_{" * join(parts[2:end], ",") * "}")
end
# a float with an exponent in standard form, 1.0e20 as 1.0 \cdot 10^{20} (SymPy #4429)
function _letters(v::Float64)
    t = string(v)
    occursin('e', t) || return v
    m, x = split(t, 'e')
    Latexify.LaTeXString(m * " \\cdot 10^{" * string(parse(Int, x)) * "}")
end
_rawx(a) = Latexify.latexraw(_letters(a))
# Sum, Integral and Derivative nodes in their notation (they were typeset as calls;
# SymPy #4767, #4438)
function _calculus(e::Expr)
    f, as = e.args[1], e.args[2:end]
    if f === :Sum && length(as) == 4
        return Latexify.LaTeXString("\\sum_{" * _rawx(as[2]) * "=" * _rawx(as[3]) * "}^{" * _rawx(as[4]) * "} " * _rawx(as[1]))
    elseif f === :Integral && length(as) == 2
        # (printed Integral(f, x) or Integral(f, (x, a, b)))
        v, lim = as[2], ""
        if v isa Expr && v.head === :tuple && length(v.args) == 3
            v, lim = v.args[1], "_{" * _rawx(v.args[2]) * "}^{" * _rawx(v.args[3]) * "}"
        end
        return Latexify.LaTeXString("\\int" * lim * " " * _rawx(as[1]) * " \\, d" * _rawx(v))
    elseif f === :Derivative && length(as) == 2
        vars = Any[as[2]]
        F = as[1]
        while F isa Expr && F.head === :call && F.args[1] === :Derivative && length(F.args) == 3
            pushfirst!(vars, F.args[3])
            F = F.args[2]
        end
        partial = !(F isa Expr && F.head === :call && length(F.args) == 2)
        d = partial ? "\\partial" : "d"
        n = length(vars)
        num = n == 1 ? d : d * "^{" * string(n) * "}"
        den = join((d * " " * _rawx(v) * (k == 1 ? "" : "^{" * string(k) * "}") for (v, k) in
                    [(v, count(==(v), vars)) for v in unique(vars)]), " ")
        return Latexify.LaTeXString("\\frac{" * num * "}{" * den * "} " * _rawx(F))
    end
    nothing
end
function _letters(e::Expr)
    e.head === :macrocall && e.args[end] isa String && return parse(BigInt, e.args[end])
    if e.head === :call && e.args[1] in (:Sum, :Integral, :Derivative)
        r = _calculus(e)
        r === nothing || return r
    end
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
    elseif e.head === :call && e.args[1] === :/ && length(e.args) == 3 && (n = _negated(e.args[2])) !== nothing
        # the sign of a quotient in front: -\frac{x}{2}, not \frac{ - x}{2} (SymPy #3973)
        return Expr(:call, :-, _letters(Expr(:call, :/, n, e.args[3])))
    elseif e.head === :call && e.args[1] === :conj && length(e.args) == 2
        # a conjugate as an overline, \overline{\Psi} (SymPy #2779's review)
        return Latexify.LaTeXString("\\overline{" * Latexify.latexraw(_letters(e.args[2])) * "}")
    elseif e.head === :call && haskey(_FNTEX, e.args[1])
        return Expr(:call, _FNTEX[e.args[1]], map(_letters, e.args[2:end])...)
    elseif e.head === :call && e.args[1] === :^ && length(e.args) == 3 && e.args[2] isa Expr &&
           e.args[2].head === :call && e.args[2].args[1] in (:Sum, :Integral, :Derivative)
        # a power of a derivative (sum, integral) in parentheses: (n'(r))², not d/dr of
        # n(r)² (SymPy #3179)
        b = _letters(e.args[2])
        b isa Latexify.LaTeXString && (b = Latexify.LaTeXString("\\left(" * b * "\\right)"))
        return Expr(:call, :^, b, _letters(e.args[3]))
    end
    Expr(e.head, map(_letters, e.args)...)
end
const _RELOPS = Dict(:≪ => :<, :≦ => :≤, :≫ => :>, :≧ => :≥, :~ => :(==), :≶ => :≠)
# -a for a numerator that is a negation (-x), a negative number (-1) or a product with a
# negative leading number (-3x): a; else nothing
function _negated(n)
    n isa Real && n < 0 && return -n
    n isa Expr && n.head === :call || return nothing
    n.args[1] === :- && length(n.args) == 2 && return n.args[2]
    if n.args[1] === :* && length(n.args) >= 3 && n.args[2] isa Real && n.args[2] < 0
        c = -n.args[2]
        return isone(c) ? (length(n.args) == 3 ? n.args[3] : Expr(:call, :*, n.args[3:end]...)) : Expr(:call, :*, c, n.args[3:end]...)
    end
    nothing
end
# functions named like a Greek letter whose symbol is the capital: Γ(x) (not γ, Euler's
# constant), the beta function B(x, y) (SymPy #3911)
const _FNTEX = Dict(:gamma => :Γ, :beta => Latexify.LaTeXString("\\mathrm{B}"),
                    :fresnels => Latexify.LaTeXString("\\mathrm{S}"), :fresnelc => Latexify.LaTeXString("\\mathrm{C}"))

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
    e = _expr(x)
    # (a LaTeXString alone, a Sum or a name, in a block: it gets its $ delimiters)
    return e isa Latexify.LaTeXString ? Expr(:block, e) : e
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
# equations (a system, or dsolve's answer for one), one per line
_tex(v::AbstractVector{<:Relation}) = isempty(v) ? "\\left[ \\; \\right]" :
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

# solveset's sets: points, families over ℤ, unsolved parts, joined by ∪
_tex(f::ImageSet) = "\\left\\{ " * _raw(f.expr) * " \\;\\middle|\\; " * join(map(_raw, f.vars), ", ") *
    " \\in \\mathbb{Z}" * join((", \\; " * _texcond(c) for c in f.conds)) * " \\right\\}"
_tex(c::ConditionSet) = "\\left\\{ " * _raw(c.x) * " \\in " * _SETTEX[c.domain.rank] * " \\;\\middle|\\; " * _raw(c.eq) * " \\right\\}"
function _tex(S::SolutionSet)
    parts = String[]
    isempty(S.points) || push!(parts, "\\left\\{ " * join(map(_raw, S.points), ", ") * " \\right\\}")
    append!(parts, map(_tex, S.families))
    append!(parts, map(_tex, S.unsolved))
    isempty(parts) ? "\\emptyset" : join(parts, " \\cup ")
end

const _TeX = Union{Sym,Relation,SymNumber,Poly,SymSeries,AbstractArray{Sym},AbstractArray{SymNumber},
                   AbstractDict{Sym,Sym},AbstractVector{<:AbstractDict{Sym,Sym}},AbstractVector{<:Relation},
                   SymSet,CylindricalSet,Cases,SolutionSet}
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
