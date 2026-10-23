# Symbolics (SymbolicUtils terms), with GruntzLimit for limits and RischIntegrate for
# integrals, as for SymJ (Symbolics has neither); Nemo and Groebner load its factoring
# and polynomial solving
const T_START = time()
using Symbolics, Nemo, Groebner
import GruntzLimit, RischIntegrate, TermInterface
import SpecialFunctions, PolyLog, HypergeometricFunctions
using Symbolics: unwrap
using TermInterface: iscall, operation, arguments
const NAMES = (:x, :y, :z, :w, :t, :a, :b, :c, :d, :n,
               :v1, :v2, :v3, :v4, :v5, :v6, :v7, :v8, :v9, :v10, :v11, :v12)
const VAR = Dict{Symbol,Any}(v => Symbolics.variable(v) for v in NAMES)
# exact π and ℯ as in GruntzLimit's tests: identity(π) and exp of a constant 1
const ONE = SymbolicUtils.Const{SymbolicUtils.SymReal}(1)
VAR[:π] = Num(SymbolicUtils.maketerm(typeof(unwrap(VAR[:x])), identity, Any[π], nothing))
VAR[:ℯ] = Num(exp(ONE))
EXACT(k) = Num(k)
num(e) = e isa Num ? e : Num(e)

const UNAVAILABLE = (:gcd,)                       # no polynomial gcd for Symbolics expressions
EXPAND(f) = Symbolics.expand(f)
# factor_use_nemo gives a repeated factor as a power
function FACTOR(f)
    u, fs = Symbolics.factor_use_nemo(f)
    out = map(fs) do g
        e = unwrap(g)
        k = iscall(e) && operation(e) === (^) ? SymbolicUtils.unwrap_const(arguments(e)[2]) : nothing
        k isa Integer ?
            (Num(arguments(e)[1]), Int(k)) : (g, 1)
    end
    (u, out)
end
CANCEL(f) = simplify_fractions(f)
APART(f, x) = partial_frac_decomposition(f, x)
DIFF(f, x, k) = expand_derivatives((Differential(x)^k)(f))
SOLVE(l, r, x) = symbolic_solve(num(l) ~ num(r), x)
LINSOLVE(eqs, vars) = symbolic_linear_solve([num(l) ~ num(r) for (l, r) in eqs], vars)
LIMIT(f, x, x0, dir) = GruntzLimit.gruntz_limit(unwrap(f), unwrap(x), x0 isa Num ? unwrap(x0) : x0; dir)
INTEGRATE(f, x) = RischIntegrate.integrate(unwrap(f), unwrap(x))
# no answer: nothing (solve, partial fractions), or roots_of placeholders
UNSOLVED(op, out) = out === nothing || occursin("roots_of", string(out))

function EVAL(ex, vals)
    ex isa RischIntegrate.IntegralWithRootSums && return ex(vals[Symbol(string(ex.var))])
    ex isa Num || ex isa SymbolicUtils.BasicSymbolic || return ex isa Number ? ex : nothing
    # real values first; if that fails (log of a negative, say), complex ones, as the
    # other systems' evaluation does
    for cplx in (false, true)
        u = try
            d = Dict{Any,Any}(unwrap(VAR[v]) => (cplx ? complex(p) : p) for (v, p) in vals)
            SymbolicUtils.unwrap_const(SymbolicUtils.substitute(unwrap(ex), d; fold = Val(true)))
        catch
            nothing
        end
        u isa Number && return u
    end
    nothing
end
function NTERMS(ex)
    e = unwrap(ex)
    iscall(e) && operation(e) === (+) ? length(arguments(e)) : 1
end
DENDEG(ex, x) = Symbolics.degree(Symbolics.denominator(num(ex)), x)
CLEARCACHE() = nothing
VERSIONS() = [("Symbolics", string(pkgversion(Symbolics))), ("SymbolicUtils", string(pkgversion(SymbolicUtils))),
              ("GruntzLimit", string(pkgversion(GruntzLimit))), ("RischIntegrate", string(pkgversion(RischIntegrate)))]
include("core.jl")
