# SymPy (through SymPyPythonCall), for reference
const T_START = time()
using SymPyPythonCall
const S = SymPyPythonCall.sympy
@syms x y z w t a b c d n v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12
const VAR = Dict{Symbol,Any}(v => getfield(Main, v) for v in
    (:x, :y, :z, :w, :t, :a, :b, :c, :d, :n, :v1, :v2, :v3, :v4, :v5, :v6, :v7, :v8, :v9, :v10, :v11, :v12))
VAR[:π], VAR[:ℯ] = SymPyPythonCall.PI, SymPyPythonCall.E
EXACT(k) = Sym(k)
toint(m) = Int(N(m))

EXPAND(f) = expand(f)
FACTOR(f) = ((u, fs) = S.factor_list(f); (u, [(g, toint(m)) for (g, m) in fs]))
GCD(f, g) = S.gcd(f, g)
CANCEL(f) = cancel(f)
APART(f, x) = apart(f, x)
DIFF(f, x, k) = diff(f, x, k)
SOLVE(l, r, x) = solve(Eq(l, r), x)
LINSOLVE(eqs, vars) = (d = solve([Eq(l, r) for (l, r) in eqs], vars); [d[v] for v in vars])
limpoint(x0) = x0 isa Real && isinf(x0) ? (x0 > 0 ? oo : -oo) : x0
LIMIT(f, x, x0, dir) = limit(f, x => limpoint(x0); dir = dir === :+ ? "+" : dir === :- ? "-" : "+-")
INTEGRATE(f, x) = integrate(f, x)
# what SymPy leaves unevaluated, or answers with no limit (zoo: complex infinity)
function UNSOLVED(op, out)
    s = string(out)
    op === :integrate && return occursin("Integral(", s)
    op === :limit && return occursin(r"Limit\(|zoo|nan|AccumBounds", s)
    false
end

# (a Sym is a Number: a Julia number is anything else)
function EVAL(ex, vals)
    ex isa Sym || return ex isa Number ? ex : nothing
    u = N(subs(ex, (VAR[v] => p for (v, p) in vals)...))
    u isa Number && !(u isa Sym) ? u : nothing
end
NTERMS(ex) = ex isa Sym && Bool(ex.is_Add) ? length(ex.args) : 1
DENDEG(ex, x) = toint(S.degree(S.denom(ex), x))
CLEARCACHE() = S.core.cache.clear_cache()
VERSIONS() = [("SymPy", string(S.__version__)), ("SymPyPythonCall", string(pkgversion(SymPyPythonCall)))]
include("core.jl")
