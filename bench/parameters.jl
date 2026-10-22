# Regression benchmark for the core operations, used to check that adding symbolic
# parameters costs nothing measurable. Minimum times (steadier than medians for
# spotting regressions).
#
#   julia --project=bench bench/parameters.jl save NAME       # results to bench/results_NAME.tsv
#   julia --project=bench bench/parameters.jl compare NAME    # this checkout against NAME
#
# Single workloads vary by about ±10% between identical runs, and a saved baseline
# drifts with machine load. For a verdict, put the old commit in a git worktree with
# its own environment and run old and new alternately (`save A1`, `save B1`, `save A2`,
# ...), comparing each workload's best time; see docs/src/design.md, "Parameters".
using BenchmarkTools
using SymJ
using SymJ: Sym

const x, y, z, w = Sym(:x), Sym(:y), Sym(:z), Sym(:w)

build_sum(n) = (acc = 0 * x; for i in 1:n; acc += i * x^i * y^(n - i); end; acc)
build_prod(n) = (acc = 1 + 0 * x; for i in 1:n; acc *= (x + i)^i; end; acc)

const BIG = expand((x + y + z + 1)^8)                             # 165 terms
const MIXED = sin(x)^2 + 3x^5 - 2x + exp(-x / 2) + sqrt(Sym(2)) * x / 3 + log(1 + x^2) / (x + 1)
const RATIONAL = (x^3 - y^3) / (x - y) + 1 / (x + y) - x / (x^2 - y^2)
const FACT = expand((x^2 + 1)^3 * (x - 3)^2 * (2x + 5) * (x^4 - x + 1))
const PG1 = Poly((x + y + z + w + 1)^4 * ((x - y + 2z * w + 1)^2 + y), x, y, z, w)
const PG2 = Poly((x + y + z + w + 1)^4 * ((y * z - x + 3)^2 + w), x, y, z, w)
const LAMB = lambdify(MIXED)
const TERMS = SymJ.arguments(BIG)
const HASPARAMS = isdefined(SymJ, :parameter)              # absent before parameters were added
const PA, PB = HASPARAMS ? (SymJ.parameter(:a), SymJ.parameter(:b)) : (Sym(:a), Sym(:b))
const PEX = PA * x^2 + PB * x + sin(PA * x)
const PBIG = expand((x + PA + PB + 1)^6)
const PLAMB = lambdify(PEX, x, PA, PB)
const R15 = Ref(1.5)                                              # not a literal: no constant folding

const SUITE = [
    "x + y" => () -> x + y,
    "x * y" => () -> x * y,
    "x^2" => () -> x^2,
    "build_sum(200)" => () -> build_sum(200),
    "build_prod(30)" => () -> build_prod(30),
    "expand (x+y+z+1)^15" => () -> expand((x + y + z + 1)^15),
    "expand (a+x)^100" => () -> expand((y + x)^100),
    "diff MIXED" => () -> diff(MIXED, x),
    "diff BIG" => () -> diff(BIG, x),
    "subs BIG" => () -> subs(BIG, x => y + 1),
    "free_symbols BIG" => () -> free_symbols(BIG),
    "hash+== BIG terms" => () -> (d = Dict{Sym,Int}(); for t in TERMS; d[t] = 1; end; BIG == BIG + 0),
    "call MIXED(1.5)" => () -> MIXED(R15[]),
    "call BIG(x => 2)" => () -> BIG(x => 2),
    "N(sqrt(2) + π)" => () -> N(sqrt(Sym(2)) + Sym(π)),
    "lambdified call" => () -> LAMB(R15[]),
    "Poly(BIG)" => () -> Poly(BIG, x, y, z),
    "Sym(Poly)" => () -> Sym(PG1),
    "gcd Poly 4 vars" => () -> gcd(PG1, PG2),
    "factor univariate" => () -> factor(FACT),
    "cancel RATIONAL" => () -> cancel(RATIONAL),
    "string BIG" => () -> string(BIG),
    # with parameters (no baseline: these are new)
    "PARAM a*x + b" => () -> PA * x + PB,
    "PARAM expand (a+b+x+1)^8" => () -> expand((PA + PB + x + 1)^8),
    "PARAM call PEX(x => 3)" => () -> PEX(x => 3),
    "PARAM free_parameters" => () -> HASPARAMS ? SymJ.free_parameters(PBIG) : nothing,
    "PARAM lambdified call" => () -> PLAMB(R15[], R15[], R15[]),
]

function run_suite()
    out = Pair{String,Float64}[]
    for (name, f) in SUITE
        f()
        b = @benchmark $f() samples = 500 seconds = 3
        push!(out, name => minimum(b).time)
        println(rpad(name, 24), lpad(BenchmarkTools.prettytime(minimum(b).time), 12))
    end
    out
end

resultfile(name) = joinpath(@__DIR__, "results_$name.tsv")

function compare(name)
    base = Dict(split(l, '\t')[1] => parse(Float64, split(l, '\t')[2]) for l in eachline(resultfile(name)))
    now = run_suite()
    ratios = Float64[]
    println("\n", rpad("benchmark", 24), lpad(name, 12), lpad("now", 12), lpad("ratio", 8))
    for (k, v) in now
        haskey(base, k) || (println(rpad(k, 24), lpad("", 12), lpad(BenchmarkTools.prettytime(v), 12)); continue)
        r = v / base[k]
        push!(ratios, r)
        println(rpad(k, 24), lpad(BenchmarkTools.prettytime(base[k]), 12), lpad(BenchmarkTools.prettytime(v), 12),
                lpad(string(round(r; digits = 2)), 8), r > 1.15 ? "   <-- slower" : "")
    end
    # single workloads vary by about ±10% between identical runs on this machine; the
    # geometric mean over the suite is the steadier summary
    println("\ngeometric mean ratio: ", round(exp(sum(log, ratios) / length(ratios)); digits = 3),
            "   (flag: a workload > 1.15)")
end

mode, name = length(ARGS) >= 2 ? (ARGS[1], ARGS[2]) : ("run", "")
if mode == "save"
    r = run_suite()
    open(resultfile(name), "w") do io
        for (k, v) in r
            println(io, k, '\t', v)
        end
    end
    println("saved ", resultfile(name))
elseif mode == "compare"
    compare(name)
else
    run_suite()
end
