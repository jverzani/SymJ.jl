# SymJ, with GruntzLimit for limits and RischIntegrate for integrals
const T_START = time()
using SymJ
import GruntzLimit, RischIntegrate, TermInterface
import SpecialFunctions, PolyLog, HypergeometricFunctions   # load SymJ's extensions for special functions
const NAMES = (:x, :y, :z, :w, :t, :a, :b, :c, :d, :n,
               :v1, :v2, :v3, :v4, :v5, :v6, :v7, :v8, :v9, :v10, :v11, :v12)
const VAR = Dict{Symbol,Any}(v => Sym(v) for v in NAMES)
VAR[:π], VAR[:ℯ] = Sym(π), Sym(ℯ)
EXACT(n) = Sym(n)
sym(e) = e isa Sym ? e : Sym(e)
gens(f...) = sort!(collect(union(free_symbols.(sym.(f))...)); by = string)

EXPAND(f) = expand(f)
function FACTOR(f)
    unit, fs = factor(Poly(f, gens(f)...))
    (Sym(unit), [(Sym(p), m) for (p, m) in fs])
end
GCD(f, g) = gcd(f, g)
CANCEL(f) = cancel(f)
APART(f, x) = apart(f, x)
DIFF(f, x, n) = foldl((e, _) -> diff(e, x), 1:n; init = f)
SOLVE(l, r, x) = solve(sym(l) ~ sym(r), x)
LINSOLVE(eqs, vars) = (d = solve([sym(l) ~ sym(r) for (l, r) in eqs], vars); [d[v] for v in vars])
LIMIT(f, x, x0, dir) = GruntzLimit.gruntz_limit(f, x, x0; dir)
INTEGRATE(f, x) = RischIntegrate.integrate(f, x)
UNSOLVED(op, out) = false                          # failures throw

function EVAL(ex, vals)
    ex isa Number && return ex
    # an answer with RootSums (or ₂F₁ kept numeric) is called at a value of its variable;
    # the other values must be those of INTPARAMS (it has no parameters to substitute)
    ex isa RischIntegrate.IntegralWithRootSums &&
        return ex(vals[Symbol(string(ex.var))])
    u = N(subs(ex, (VAR[v] => p for (v, p) in vals)...))
    u isa Number ? u : nothing
end
NTERMS(ex) = ex isa Sym && TermInterface.iscall(ex) && TermInterface.operation(ex) === (+) ?
    length(TermInterface.arguments(ex)) : 1
DENDEG(ex, x) = (d = numer_denom(sym(ex))[2]; degree(Poly(d, x), x))
CLEARCACHE() = nothing
VERSIONS() = [("SymJ", string(pkgversion(SymJ))), ("GruntzLimit", string(pkgversion(GruntzLimit))),
              ("RischIntegrate", string(pkgversion(RischIntegrate)))]
include("core.jl")
