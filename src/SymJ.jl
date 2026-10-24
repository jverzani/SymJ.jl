"""
    SymJ

Symbolic manipulation with a single concrete, hash-carrying node type, eager
canonicalization of sums, products and powers, exact numbers, and support for the
TermInterface v2 contract. See `docs/src/design.md` for the design.
"""
module SymJ

using TermInterface
using RuntimeGeneratedFunctions
using ScopedValues: ScopedValue, with
import GruntzLimit
using GruntzLimit: GruntzError
import CommonEq: Eq, Lt, Le, Gt, Ge, Ne, ⩵, ≪, ≦, ≫, ≧, ≶, ≷
RuntimeGeneratedFunctions.init(@__MODULE__)
import TermInterface: isexpr, iscall, head, children, operation, arguments, maketerm

export AbstractSymbolic, Sym, SymNumber, SymFunction
export @syms, symbols, subs, evaluate, lambdify, lambdify_expr, expand, free_symbols, has, N
export register_function!
export numer_denom, cancel, factor_squarefree, squarefree_part, factor, apart, resultant, discriminant, groebner, reduced
export solve, lambertw, series, limit, GruntzError, simplify, refine
export Assumptions, assuming, @assume, @assuming, ℤ, ℚ, ℝ, ℂ, Integers, Rationals, Reals, Complexes
export Relation, lhs, rhs, Eq, Lt, Le, Gt, Ge, Ne, ⩵, ≪, ≦, ≫, ≧, ≶, ≷
export Poly, gens, nterms, coefficients, exponent_vectors, degree, divides, divexact

include("numbers.jl")
include("types.jl")
include("arithmetic.jl")
include("functions.jl")
include("mathfunctions.jl")
include("lambertw.jl")
include("calculus.jl")
include("manipulation.jl")
include("poly.jl")
include("polydiv.jl")
include("polygcd.jl")
include("modgcd.jl")
include("rational.jl")
include("squarefree.jl")
include("factor.jl")
include("vanhoeij.jl")
include("eez.jl")
include("ratfun.jl")
include("apart.jl")
include("iteration.jl")
include("terminterface.jl")
include("relations.jl")
include("solve.jl")
include("intervals.jl")
include("linsolve.jl")
include("ineqsystems.jl")
include("assumptions.jl")
include("cases.jl")
include("limit.jl")
include("simplify.jl")
include("series.jl")
include("resultant.jl")
include("groebner.jl")
include("nlsolve.jl")
include("lambdify.jl")
include("show.jl")
include("precompile.jl")

end
