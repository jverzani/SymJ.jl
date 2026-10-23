# julia build.jl
# Builds results/cas_comparison.html from results/*.tsv and template.html
const HERE = @__DIR__
const RESULTS = joinpath(HERE, "results")
const SYSTEMS = [
    ("symj", "SymJ", "with GruntzLimit (limits) and RischIntegrate (integrals)"),
    ("symbolics", "Symbolics", "SymbolicUtils terms; GruntzLimit and RischIntegrate as for SymJ"),
    ("sympy", "SymPy", "through SymPyPythonCall"),
    ("giac", "Giac", "through Giac.jl"),
]

# a small JSON writer, enough for the page's data
json(x::Nothing) = "null"
json(x::Bool) = string(x)
json(x::Real) = isfinite(x) ? string(x) : "null"
json(x::Symbol) = json(string(x))
function json(s::AbstractString)
    io = IOBuffer(); print(io, '"')
    for c in s
        c == '"' ? print(io, "\\\"") : c == '\\' ? print(io, "\\\\") :
        c == '<' ? print(io, "\\u003c") : c < ' ' ? print(io, "\\u", string(UInt16(c), base = 16, pad = 4)) : print(io, c)
    end
    print(io, '"'); String(take!(io))
end
json(v::Union{AbstractVector,Tuple}) = "[" * join(json.(v), ",") * "]"
json(d::AbstractDict) = "{" * join((json(string(k)) * ":" * json(v) for (k, v) in d), ",") * "}"

# how a problem reads on the page
showex(e) = replace(string(e), "Inf" => "∞")
function describe(op, args)
    op === :expand && return "expand " * showex(args[1])
    op === :factor && return "factor " * showex(args[1])
    op === :gcd && return "gcd(" * showex(args[1]) * ", " * showex(args[2]) * ")"
    op === :cancel && return "cancel " * showex(args[1])
    op === :apart && return "partial fractions of " * showex(args[1]) * " in " * string(args[2])
    op === :diff && return (args[3] == 1 ? "d/d$(args[2]) " : "d^$(args[3])/d$(args[2])^$(args[3]) ") * showex(args[1])
    op === :solve && return "solve " * showex(args[1]) * " = " * showex(args[2]) * " for " * string(args[3])
    op === :linsolve && return "solve " * join((showex(l) * " = " * showex(r) for (l, r) in args[1]), ", ") *
                               " for " * join(args[2], ", ")
    op === :limit && return "lim x → " * showex(args[2]) * (args[3] === :both ? "" : string(args[3])) * " of " * showex(args[1])
    op === :integrate && return "∫ " * showex(args[1]) * " d" * string(args[2])
end

problems = [Dict("i" => i, "area" => area, "op" => op, "src" => src, "p" => describe(op, args))
            for (i, (area, op, args, _, src)) in enumerate(include(joinpath(HERE, "problems.jl")))]
N = length(problems)
areas = unique(p["area"] for p in problems)

versions = Dict{String,String}("Julia" => string(VERSION))
systems = map(SYSTEMS) do (key, name, how)
    cells = Vector{Any}(nothing, N)
    loads, wall, limit = Float64[], nothing, nothing
    for line in eachline(joinpath(RESULTS, "results_$key.tsv"))
        p = split(line, '\t')
        if p[1] == "R"
            t1, t2 = parse(Float64, p[5]), parse(Float64, p[6])
            isnan(t1) && (t1 = 0.0)                    # not attempted (not available)
            cells[parse(Int, p[2])] = Dict("s" => p[3], "v" => p[4], "t" => round(t1, digits = 4),
                                           "r" => isnan(t2) ? nothing : round(t2, digits = 4),
                                           "o" => length(p) > 6 ? p[7] : "")
        elseif p[1] == "L"
            push!(loads, parse(Float64, p[2]))
        elseif p[1] == "T"
            wall, limit = parse(Float64, p[2]), parse(Float64, p[3])
        elseif p[1] == "V"
            versions[p[2]] = p[3]
        end
    end
    missing = findall(isnothing, cells)
    isempty(missing) || error("$key: no result for problems $missing")
    Dict("key" => key, "name" => name, "how" => how, "cells" => cells,
         "load" => isempty(loads) ? nothing : round(first(loads), digits = 2),
         "restarts" => max(length(loads) - 1, 0),
         "wall" => wall === nothing ? nothing : round(wall, digits = 1), "limit" => limit)
end

commits = Dict{String,String}()
f = joinpath(RESULTS, "commits.tsv")
isfile(f) && for line in eachline(f)
    k, v = split(line, '\t'); commits[k] = v
end
date = Libc.strftime("%Y-%m-%d", isfile(f) ? mtime(f) : time())

data = Dict("problems" => problems, "areas" => areas, "systems" => systems, "date" => date,
            "limit" => maximum(something(s["limit"], 0) for s in systems),
            "versions" => versions, "commits" => commits)
html = replace(read(joinpath(HERE, "template.html"), String), "/*DATA*/null" => json(data))
write(joinpath(RESULTS, "cas_comparison.html"), html)

# a text summary
kind(c) = c["s"] == "solved" ? c["v"] : c["s"] in ("timeout", "crash") ? "timeout" : c["s"]   # or "unavailable"
println(rpad("area", 20), join((rpad(s["name"], 30) for s in systems)))
for a in areas
    idx = [p["i"] for p in problems if p["area"] == a]
    print(rpad("$a ($(length(idx)))", 20))
    for s in systems
        c = s["cells"][idx]
        v = count(x -> kind(x) == "verified", c)
        t = sum(x -> x["t"], c)
        print(rpad("$v verified, $(round(t, digits = 1)) s", 30))
    end
    println()
end
