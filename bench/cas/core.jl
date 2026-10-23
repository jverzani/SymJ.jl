# Shared by the workers. Before including, a worker defines
#   VAR::Dict{Symbol,Any}    its variables (see ALLVARS) and the constants :π, :ℯ
#   the operations           EXPAND(f), FACTOR(f) -> (unit, [(g, m), ...]), GCD(f, g),
#                            CANCEL(f), APART(f, x), DIFF(f, x, n), SOLVE(lhs, rhs, x) -> roots,
#                            LINSOLVE(eqs, vars) -> values in the order of vars,
#                            LIMIT(f, x, x0, dir), INTEGRATE(f, x)
#   UNSOLVED(op, out)::Bool  did it return the problem unevaluated (or no answer)?
#   EVAL(ex, vals)           the value of `ex` with each variable v replaced by vals[v],
#                            as a number, or `nothing`
#   NTERMS(ex)               the number of terms of a sum (1 if not a sum)
#   DENDEG(ex, x)            the degree in x of the denominator
#   CLEARCACHE()             forget cached results (before the repeat timing)
# Usage: julia --project=ENV worker_X.jl START OUTFILE [STOP] (run by supervise.jl)
const T_LOADED = time()

const PROBLEMS = include(joinpath(@__DIR__, "problems.jl"))
const ALLVARS = (:x, :y, :z, :w, :t, :a, :b, :c, :d, :n,
                 :v1, :v2, :v3, :v4, :v5, :v6, :v7, :v8, :v9, :v10, :v11, :v12)

# the term from a problem's expression, with this worker's variables. Numeric constants
# such as sqrt(2), log(6) or 1/2 would be computed in floating point by Julia before the
# system sees them, so their arguments go through the worker's EXACT (exact numbers).
const BUILDMOD = Module()
for v in (ALLVARS..., :π, :ℯ)
    Core.eval(BUILDMOD, :(const $v = $(VAR[v])))
end
Core.eval(BUILDMOD, :(const EXACT = $EXACT))
isnum(a) = a isa Number && !(a isa Bool)
function exactify(ex)
    ex isa Expr || return ex
    args = map(exactify, ex.args)
    if ex.head === :call
        f = args[1]
        if f === :/ && all(isnum, args[2:end])
            args[2] = :(EXACT($(args[2])))
        elseif !(f in (:+, :-, :*, :^, ://, :/)) && length(args) >= 2 && all(isnum, args[2:end])
            args[2] = :(EXACT($(args[2])))
        end
    end
    Expr(ex.head, args...)
end
build(ex) = ex isa Union{Expr,Symbol} ? Core.eval(BUILDMOD, exactify(ex)) : ex

## references: the problem's expressions evaluated in plain Julia, in BigFloat

# vals: Symbol => number. Complex values are allowed (for an integrand at a point where
# it is complex).
function refval(ex, vals; prec = 256)
    setprecision(BigFloat, prec) do
        m = Module()
        Core.eval(m, :(const π = $(big(π)); const ℯ = $(big(ℯ))))
        for (v, p) in vals
            Core.eval(m, :(const $v = $(p isa Complex ? Complex{BigFloat}(p) : BigFloat(p))))
        end
        Core.eval(m, ex)
    end
end
f64(u) = u isa Real ? Float64(u) : ComplexF64(u)

# values for every variable: two points, distinct values, and a different shift for each
# variable between them (so that x - y, say, changes too)
point(k) = Dict(v => 0.31 + 0.071 * i + (k - 1) * (0.137 + 0.029 * i) for (i, v) in enumerate(ALLVARS))
const PTS = (point(1), point(2))
# parameters of equations and limits (positive, as the limit answers assume), and the
# integrals' values as in RischIntegrate's benchmark
const EQPARAMS = Dict(:a => 1.7, :b => 2.3, :c => 0.4)
const INTPARAMS = Dict(v => 1.0 + 0.137 * i for (i, v) in enumerate((:x, :t, :a, :b, :c, :d, :z, :n)))

close_to(u, r; tol = 1e-6) = isfinite(u) && isfinite(r) && abs(u - r) <= tol * max(1, abs(r))
value(out, vals) = (u = EVAL(out, vals); u === nothing ? nothing : complex(u))

# does out equal the expression ex (the input, or a reference answer) at both points?
function same_value(out, ex)
    for p in PTS
        u = value(out, p)
        u === nothing && return nothing
        close_to(u, complex(f64(refval(ex, p)))) || return false
    end
    true
end

## checks: :verified, :partial (the value is right, the form is not what was asked),
## :wrong, or :unverified (the answer could not be evaluated)

verdict(ok) = ok === nothing ? :unverified : ok ? :verified : :wrong

function check(op, args, expect, out)
    if op === :expand
        v = same_value(out, args[1])
        v === true || return verdict(v)
        NTERMS(out) == expect ? :verified : :partial
    elseif op === :factor
        unit, fs = out
        nonconst = 0
        for p in PTS
            u = value(unit, p)
            u === nothing && return :unverified
            for (g, m) in fs
                gv = value(g, p)
                gv === nothing && return :unverified
                u *= gv^m
            end
            close_to(u, complex(f64(refval(args[1], p)))) || return :wrong
        end
        # irreducible factors: count those that are not constants
        nonconst = sum((m for (g, m) in fs if value(g, PTS[1]) != value(g, PTS[2])); init = 0)
        nonconst == expect ? :verified : :partial
    elseif op === :gcd
        r = map(PTS) do p
            u = value(out, p)
            u === nothing ? nothing : u / complex(f64(refval(expect, p)))
        end
        any(isnothing, r) && return :unverified
        abs(r[1]) > 1e-12 && close_to(r[1], r[2]) ? :verified : :wrong
    elseif op === :cancel
        v = same_value(out, args[1])
        v === true || return verdict(v)
        DENDEG(out, VAR[:x]) == expect ? :verified : :partial
    elseif op === :apart
        v = same_value(out, args[1])
        v === true || return verdict(v)
        NTERMS(out) == expect ? :verified : :partial
    elseif op === :diff
        f, xs, n = args
        for p in PTS
            u = value(out, p)
            u === nothing && return :unverified
            close_to(u, complex(f64(nth_derivative(f, xs, n, p)))) || return :wrong
        end
        :verified
    elseif op === :solve
        lhs, rhs, xs = args
        reals = ComplexF64[]
        for r in out
            u = value(r, EQPARAMS)
            u === nothing && return :unverified
            vals = merge(EQPARAMS, Dict(xs => u))
            res = try
                complex(f64(refval(:($lhs - $rhs), vals))) / max(1, abs(complex(f64(refval(lhs, vals)))))
            catch
                return :unverified
            end
            abs(res) <= 1e-8 || return :wrong
            abs(imag(u)) <= 1e-8 * max(1, abs(u)) && !any(s -> abs(s - u) <= 1e-8, reals) && push!(reals, u)
        end
        length(reals) == expect ? :verified : :partial
    elseif op === :linsolve
        eqs, vars = args
        vals = copy(EQPARAMS)
        for (v, s) in zip(vars, out)
            u = value(s, EQPARAMS)
            u === nothing && return :unverified
            vals[v] = u
        end
        all(((l, r),) -> abs(complex(f64(refval(:($l - $r), vals)))) <= 1e-8, eqs) ? :verified : :wrong
    elseif op === :limit
        occursin("NaN", string(out)) && return :wrong
        u = value(out, EQPARAMS)
        u === nothing && return :unverified
        abs(imag(u)) <= 1e-8 * max(1, abs(u)) || return :wrong
        v, r = real(u), Float64(refval(expect, EQPARAMS))
        isinf(r) ? (v == r ? :verified : :wrong) : close_to(v, r; tol = 1e-8) ? :verified : :wrong
    elseif op === :integrate
        # the central difference of the answer against the integrand, at two points where
        # the integrand is real if there are such points (on a branch cut a complex value
        # can take either side), else at 0.31 and 0.57 with complex values
        f, xs = args
        h = 1e-5
        at(x, cplx) = Dict(v => (v === xs ? x : INTPARAMS[v]) * (cplx ? 1 + 0im : 1) for v in keys(INTPARAMS))
        realval(x) = try
            u = f64(refval(f, at(x, false)))
            u isa Real && isfinite(u) ? u : nothing
        catch
            nothing
        end
        xr = filter(x -> realval(x) !== nothing, [0.31, 0.57, 0.13, 0.77, -0.21, 0.93, -0.47, 1.37, 2.3])
        pts = length(xr) >= 2 ? [(x, false) for x in xr[1:2]] : [(0.31, true), (0.57, true)]
        for (x0, cplx) in pts
            vals = try
                (complex(f64(refval(f, at(x0, cplx)))), value(out, at(x0 + h, cplx)), value(out, at(x0 - h, cplx)))
            catch
                nothing
            end
            vals === nothing && return :unverified
            any(isnothing, vals) && return :unverified
            fv, Fp, Fm = vals
            (isfinite(fv) && isfinite(Fp) && isfinite(Fm)) || return :unverified
            abs((Fp - Fm) / (2h) - fv) <= 1e-4 * max(1, abs(fv)) || return :wrong
        end
        :verified
    else
        error("unknown op $op")
    end
end

# d^n f/dx^n at the point, by a central difference in high precision (exact to well
# beyond Float64 for these smooth functions)
function nth_derivative(f, xs, n, p)
    setprecision(BigFloat, 4096) do
        h = big(10.0)^-60
        s = sum((-1)^k * binomial(n, k) * refval(f, merge(p, Dict(xs => big(p[xs]) + (n / 2 - k) * h)); prec = 4096)
                for k in 0:n)
        s / h^n
    end
end

## running a problem

function operate(op, args)
    op === :expand && return EXPAND(args[1])
    op === :factor && return FACTOR(args[1])
    op === :gcd && return GCD(args[1], args[2])
    op === :cancel && return CANCEL(args[1])
    op === :apart && return APART(args[1], args[2])
    op === :diff && return DIFF(args[1], args[2], args[3])
    op === :solve && return SOLVE(args[1], args[2], args[3])
    op === :linsolve && return LINSOLVE(args[1], args[2])
    op === :limit && return LIMIT(args[1], VAR[:x], args[2], args[3])
    op === :integrate && return INTEGRATE(args[1], args[2])
    error("unknown op $op")
end

# the arguments built with this worker's variables; factor and gcd inputs expanded
# (not timed: factoring starts from the expanded polynomial)
function inputs(op, args)
    op === :factor && return (EXPAND(build(args[1])),)
    op === :gcd && return (EXPAND(build(args[1])), EXPAND(build(args[2])))
    op === :apart && return (build(args[1]), VAR[args[2]])
    op === :diff && return (build(args[1]), VAR[args[2]], args[3])
    op === :solve && return (build(args[1]), build(args[2]), VAR[args[3]])
    op === :linsolve && return ([(build(l), build(r)) for (l, r) in args[1]], [VAR[v] for v in args[2]])
    op === :limit && return (build(args[1]), args[2] isa Real && isinf(args[2]) ? args[2] : build(args[2]), args[3])
    op === :integrate && return (build(args[1]), VAR[args[2]])
    (build(args[1]),)
end

clean(s) = replace(s, r"\s+" => " ")[1:min(end, 400)]
show_out(out) = out isa Tuple && length(out) == 2 && out[2] isa AbstractVector ?   # factor
    string(out[1], " · ", join(("($g)" * (m == 1 ? "" : "^$m") for (g, m) in out[2]), " · ")) :
    out isa AbstractVector ? "[" * join(string.(out), ", ") * "]" : string(out)

# Run once (the first call, compilation included), then, if that took under 15 s, again after
# clearing any cache (the repeat call)
function run_one(io, i)
    area, op, args, expect, _ = PROBLEMS[i]
    println(io, "S\t$i\t", time()); flush(io)
    t1, t2, out = NaN, NaN, nothing
    t0 = time()
    status, ver, s = try
        a = inputs(op, args)
        t0 = time(); out = operate(op, a); t1 = time() - t0
        if t1 < 15
            CLEARCACHE()
            t0 = time(); operate(op, a); t2 = time() - t0
        end
        if UNSOLVED(op, out)
            (:unsolved, :none, show_out(out))
        else
            v = try check(op, args, expect, out) catch; :unverified end
            (:solved, v, show_out(out))
        end
    catch err
        isnan(t1) && (t1 = time() - t0)
        (:error, :none, sprint(showerror, err))
    end
    println(io, "R\t$i\t$status\t$ver\t$t1\t$t2\t", clean(s)); flush(io)
end

let start = parse(Int, ARGS[1]), outfile = ARGS[2], stop = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : length(PROBLEMS)
    open(outfile, "a") do io
        println(io, "L\t", T_LOADED - T_START, "\t", start); flush(io)
        start == 1 && foreach(VERSIONS()) do (k, v)
            println(io, "V\t", k, "\t", v)
        end
        foreach(i -> run_one(io, i), start:stop)
        println(io, "D"); flush(io)
    end
end
