# Symbolic equations and inequalities.
#
# `==`, `<`, `≤` keep returning Bool (structural equality is what Dicts, `in` and generic
# code need). Relations are written with `~` (equations; the convention of Symbolics.jl,
# ModelingToolkit and SimpleExpressions) and with the operators and function names of
# CommonEq.jl, defined here without depending on it:
#
#   =  : lhs ~ rhs,  lhs ⩵ rhs (\Equal),   Eq(lhs, rhs)
#   <  : a ≪ b (\ll),                      Lt(a, b)
#   ≤  : a ≦ b (\leqq),                    Le(a, b)
#   >  : a ≫ b (\gg),                      Gt(a, b)
#   ≥  : a ≧ b (\geqq),                    Ge(a, b)
#   ≠  : a ≶ b (\lessgtr), a ≷ b,          Ne(a, b)
#
# A relation is not an algebraic quantity (two inequalities are not added), so it is its
# own small type rather than a Sym node.

"""
    Relation

A symbolic equation or inequality between two expressions: `x^2 + 1 ~ 3`, `x ≪ 1`
(`Lt(x, 1)`), `x ≦ 1`, `x ≫ 1`, `x ≧ 1`, `x ≶ 1` (`Ne`). See [`lhs`](@ref),
[`rhs`](@ref).

Calling a relation substitutes into both sides, as calling an expression does
(`(x^2 ~ 4)(x => 2)`, `(x^2 ~ 4)(2)`); when no symbols remain the result is a `Bool`.
`r()` evaluates a relation between constants.
"""
struct Relation <: AbstractSymbolic
    op::Symbol            # :eq, :lt, :le, :gt, :ge, :ne
    lhs::Sym
    rhs::Sym
end

"""
    lhs(r::Relation), rhs(r::Relation)

The two sides of a relation.
"""
lhs(r::Relation) = r.lhs
rhs(r::Relation) = r.rhs

const _REL_INFO = Dict(:eq => ("~", ==), :lt => ("≪", <), :le => ("≦", <=),
                       :gt => ("≫", >), :ge => ("≧", >=), :ne => ("≶", !=))

"""
    Eq(a, b), Lt(a, b), Le(a, b), Gt(a, b), Ge(a, b), Ne(a, b)

The relations `a = b`, `a < b`, `a ≤ b`, `a > b`, `a ≥ b`, `a ≠ b` between expressions
(the names of CommonEq.jl). Infix forms: `~` or `⩵`, `≪`, `≦`, `≫`, `≧`, `≶` or `≷`.
"""
function Eq end
function Lt end
function Le end
function Gt end
function Ge end
function Ne end

for (f, op) in ((:Eq, :eq), (:Lt, :lt), (:Le, :le), (:Gt, :gt), (:Ge, :ge), (:Ne, :ne))
    for (S, T) in ((Sym, Sym), (Sym, Number), (Number, Sym))
        @eval $f(a::$S, b::$T) = Relation($(QuoteNode(op)), Sym(a), Sym(b))
    end
end
for (S, T) in ((Sym, Sym), (Sym, Number), (Number, Sym))
    @eval Base.:~(a::$S, b::$T) = Eq(a, b)
end

# `<`, `≤`, `>`, `≥` compare expressions as numbers (a Bool), so they work for
# constants; for anything with a symbol, an error pointing to the symbolic operator.
const _SYMBOLIC_OP = Dict(:< => ("≪", "\\ll", "Lt"), :<= => ("≦", "\\leqq", "Le"),
                          :> => ("≫", "\\gg", "Gt"), :>= => ("≧", "\\geqq", "Ge"))
function _compare(op::Symbol, a::Sym, b::Sym)
    if _hassymbol(a) || _hassymbol(b)
        infix, tab, f = _SYMBOLIC_OP[op]
        throw(ArgumentError("`$a $op $b`: `$op` compares numbers and returns a Bool; for a " *
                            "symbolic inequality use `$infix` (typed $tab) or `$f(a, b)`"))
    end
    va, vb = N(a), N(b)
    (va isa Real && vb isa Real) || throw(ArgumentError("`$a $op $b`: not real numbers"))
    op === :< ? va < vb : op === :<= ? va <= vb : op === :> ? va > vb : va >= vb
end
for (S, T) in ((Sym, Sym), (Sym, Real), (Real, Sym))
    @eval Base.:<(a::$S, b::$T) = _compare(:<, Sym(a), Sym(b))
    @eval Base.:<=(a::$S, b::$T) = _compare(:<=, Sym(a), Sym(b))
    @eval Base.:>(a::$S, b::$T) = _compare(:>, Sym(a), Sym(b))
    @eval Base.:>=(a::$S, b::$T) = _compare(:>=, Sym(a), Sym(b))
    @eval Base.isless(a::$S, b::$T) = _compare(:<, Sym(a), Sym(b))
end

const ⩵ = Eq
const ≪ = Lt
const ≦ = Le
const ≫ = Gt
const ≧ = Ge
const ≶ = Ne
const ≷ = Ne

Base.:(==)(a::Relation, b::Relation) = a.op === b.op && a.lhs == b.lhs && a.rhs == b.rhs
Base.hash(r::Relation, h::UInt) = hash(r.rhs, hash(r.lhs, hash(r.op, h)))

function Base.show(io::IO, r::Relation)
    _show(io, r.lhs, 0)
    print(io, " ", _REL_INFO[r.op][1], " ")
    _show(io, r.rhs, 0)
end
Base.show(io::IO, ::MIME"text/plain", r::Relation) = show(io, r)

## TermInterface: a call of the relation's function on its two sides

const _REL_FUN = Dict(:eq => (~), :lt => Lt, :le => Le, :gt => Gt, :ge => Ge, :ne => Ne)
TermInterface.isexpr(::Relation) = true
TermInterface.iscall(::Relation) = true
TermInterface.operation(r::Relation) = _REL_FUN[r.op]
TermInterface.head(r::Relation) = operation(r)
TermInterface.arguments(r::Relation) = Sym[r.lhs, r.rhs]
TermInterface.children(r::Relation) = arguments(r)
TermInterface.maketerm(::Type{Relation}, head, args, metadata) = head(args...)

## calling: substitute into both sides; a relation between constants is decided

# (a side called with numbers evaluates to a plain number: wrapped back as a Sym)
(r::Relation)(args...) = _decide(Relation(r.op, Sym(r.lhs(args...)), Sym(r.rhs(args...))))
(r::Relation)() = _decide(r)

# Bool when both sides are constants, else the relation. Exact numbers compare exactly;
# otherwise numerically (equality to a relative tolerance).
function _decide(r::Relation)
    r.lhs == r.rhs && return _REL_INFO[r.op][2](0, 0)    # identical sides: whatever they are
    (_hassymbol(r.lhs) || _hassymbol(r.rhs)) && return r
    d = r.lhs - r.rhs
    _, cmp = _REL_INFO[r.op]
    if d.kind == NUMBER && isexact(d.num)
        return cmp(sign(d.num), 0)
    end
    v = N(d)
    v isa Real || return r.op === :eq ? false : r.op === :ne ? true :
        throw(ArgumentError("$r: the sides do not compare as real numbers"))
    scale = max(1.0, abs(N(r.lhs)), abs(N(r.rhs)))
    abs(v) <= 1e-12 * scale ? cmp(0, 0) : cmp(v, 0)
end

# the sides of a relation as one expression compared with zero. Only for relations: a
# bare expression is never taken to mean `ex ~ 0` (solving requires an equation).
_as_zero(r::Relation) = r.lhs - r.rhs
