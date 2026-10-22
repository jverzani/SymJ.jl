"""
    SymJ

Symbolic manipulation with a single concrete, hash-carrying node type, eager
canonicalization of sums, products and powers, exact numbers, and support for the
TermInterface v2 contract. See `docs/src/design.md` for the design.
"""
module SymJ

using TermInterface
import TermInterface: isexpr, iscall, head, children, operation, arguments, maketerm

export AbstractSymbolic, Sym, SymNumber, SymFunction
export @syms, symbols, subs, evaluate, expand, free_symbols, has, N
export register_function!
export numer_denom, cancel
export Poly, gens, nterms, coefficients, exponent_vectors, degree, divides, divexact

include("numbers.jl")
include("types.jl")
include("arithmetic.jl")
include("functions.jl")
include("mathfunctions.jl")
include("calculus.jl")
include("manipulation.jl")
include("poly.jl")
include("polydiv.jl")
include("polygcd.jl")
include("modgcd.jl")
include("rational.jl")
include("iteration.jl")
include("terminterface.jl")
include("show.jl")

end
