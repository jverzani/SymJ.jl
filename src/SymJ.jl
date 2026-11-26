"""
    SymJ

Symbolic manipulation with a single concrete, hash-carrying node type, eager
canonicalization of sums, products and powers, exact numbers, and support for the
TermInterface v2 contract. See `docs/src/design/` for the design.
"""
module SymJ

using TermInterface
using RuntimeGeneratedFunctions
using ScopedValues: ScopedValue, with
import GruntzLimit
import LinearAlgebra
using GruntzLimit: GruntzError
import CommonEq: Eq, Lt, Le, Gt, Ge, Ne, ⩵, ≪, ≦, ≫, ≧, ≶, ≷
RuntimeGeneratedFunctions.init(@__MODULE__)
import TermInterface: isexpr, iscall, head, children, operation, arguments, maketerm

export Sym
export @syms, subs, lambdify, expand, free_symbols, has, N
export cancel, factor, apart, resultant, discriminant, groebner, reduced
export solve, lambertw, series, limit, simplify, refine
export assuming, @assume, @assuming, ℤ, ℚ, ℝ, ℂ, Integers, Rationals, Reals, Complexes
export Relation, lhs, rhs, Eq, Lt, Le, Gt, Ge, Ne, ⩵, ≪, ≦, ≫, ≧, ≶, ≷
export Poly, gens, coefficients, exponent_vectors, degree, divides

include("numbers.jl")
include("types.jl")
include("arithmetic.jl")
include("functions.jl")
include("mathfunctions.jl")
include("lambertw.jl")
include("fresnel.jl")
include("calculus.jl")
include("manipulation.jl")
include("quadrature.jl")
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
include("decompose.jl")
include("iteration.jl")
include("terminterface.jl")
include("relations.jl")
include("solve.jl")
include("quintic.jl")
include("intervals.jl")
include("linsolve.jl")
include("ineqsystems.jl")
include("assumptions.jl")
include("cases.jl")
include("solveset.jl")
include("limit.jl")
include("simplify.jl")
include("complexparts.jl")
include("linalg.jl")
include("radsimp.jl")
include("normalforms.jl")
include("piecewise.jl")
include("enclose.jl")
include("summation.jl")
include("recurrence.jl")
include("gammasimp.jl")
include("ode.jl")
include("laplace.jl")
include("patterns.jl")
include("simprules.jl")
include("rewrite.jl")
include("identify.jl")
include("definite.jl")
include("series.jl")
include("resultant.jl")
include("minpoly.jl")
include("groebner.jl")
include("nlsolve.jl")
include("realroots.jl")
include("chebzeros.jl")
include("lambdify.jl")
include("fsolve.jl")
# SpecialFunctions, HypergeometricFunctions and PolyLog on expressions
include("specialfunctions.jl")
include("hypergeometric.jl")
include("polylog.jl")
include("show.jl")
_more_asin_values!()
include("precompile.jl")
# (after the precompile workload: on Julia 1.12 any submodule defined before it, even an
# empty one, brought back #230 — with Test or JSON loaded first, the first show of an
# expression re-inferred SymJ's printing code for minutes, `@syms x` at the REPL too)
include("commands.jl")

end
