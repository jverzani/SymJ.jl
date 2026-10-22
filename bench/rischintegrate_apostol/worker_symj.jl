const T_START = time()
using SymJ
import RischIntegrate
import SpecialFunctions, PolyLog          # as in env_risch; loads SymJ's extensions
const VAR = Dict{Symbol,Any}(v => Sym(v) for v in (:x, :t, :a, :b, :c, :d, :z, :n))
INTEGRATE(f, x) = RischIntegrate.integrate(f, x)
SOLVED(out) = true  # failures throw
# substitute the values, then evaluate numerically
function EVAL(ex, vals)
    isa(ex, Number) && return complex(ex)
    u = N(subs(ex, (VAR[v] => p for (v, p) in vals)...))
    isa(u, Number) ? complex(u) : nothing
end
include("core.jl")
