# Giac (through Giac.jl), for reference
const T_START = time()
using Giac
@giac_var x y z w t a b c d n v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12
const VAR = Dict{Symbol,Any}(v => getfield(Main, v) for v in
    (:x, :y, :z, :w, :t, :a, :b, :c, :d, :n, :v1, :v2, :v3, :v4, :v5, :v6, :v7, :v8, :v9, :v10, :v11, :v12))
VAR[:π], VAR[:ℯ] = giac_eval("pi"), giac_eval("exp(1)")
EXACT(k) = giac_eval(string(k))
gx(e) = e isa GiacExpr ? e : giac_eval(string(e))
cmd(c, args...) = invoke_cmd(c, args...)
toint(m) = Int(to_julia(m))

EXPAND(f) = cmd(:expand, f)
function FACTOR(f)                                 # [f1, m1, f2, m2, ...], constants included
    L = to_julia(cmd(:factors, f))
    (1, [(L[2k-1], toint(L[2k])) for k in 1:length(L)÷2])
end
GCD(f, g) = cmd(:gcd, f, g)
CANCEL(f) = cmd(:normal, f)
APART(f, x) = cmd(:partfrac, f, x)
DIFF(f, x, k) = cmd(:diff, f, x, k)
SOLVE(l, r, x) = to_julia(cmd(:solve, cmd(:equal, gx(l), gx(r)), x))
LINSOLVE(eqs, vars) = to_julia(cmd(:linsolve, cmd(:makevector, (cmd(:equal, gx(l), gx(r)) for (l, r) in eqs)...),
                                   cmd(:makevector, vars...)))
# limits as in GruntzLimit's benchmark: direction 1 (from the right) or -1; ±infinity
const SIDE = Dict(:+ => 1, :- => -1)
atinf(x0) = x0 isa Real && isinf(x0)
function LIMIT(f, x, x0, dir)
    p = atinf(x0) ? giac_eval(x0 > 0 ? "+infinity" : "-infinity") : gx(x0)
    out = atinf(x0) || dir === :both ? cmd(:limit, f, x, p) : cmd(:limit, f, x, p, SIDE[dir])
    s = string(out)
    s == "+infinity" ? Inf : s == "-infinity" ? -Inf : out
end
INTEGRATE(f, x) = cmd(:integrate, f, x)
function UNSOLVED(op, out)
    s = string(out)
    op === :integrate && return occursin("integrate", s)
    op === :limit && return !(out isa Real) && occursin(r"limit|undef|bounded_function|^infinity$", s)
    false
end

# values rounded to 12 digits: Giac misreads longer floats (see RischIntegrate's benchmark)
r12(p::Real) = round(p, sigdigits = 12)
r12(p::Complex) = complex(r12(real(p)), r12(imag(p)))
function EVAL(ex, vals)
    ex isa Number && return ex
    u = to_julia(cmd(:evalf, Giac.substitute(gx(ex), Dict(VAR[v] => r12(p) for (v, p) in vals))))
    u isa Number ? u : nothing
end
NTERMS(ex) = ex isa GiacExpr && string(cmd(:sommet, ex)) == "'+'" ? toint(cmd(:nops, ex)) : 1
DENDEG(ex, x) = toint(cmd(:degree, cmd(:denom, gx(ex)), x))
CLEARCACHE() = nothing
VERSIONS() = [("Giac.jl", string(pkgversion(Giac)))]
include("core.jl")
