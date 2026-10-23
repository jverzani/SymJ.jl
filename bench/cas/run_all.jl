# julia run_all.jl [LIMIT]
# Runs the systems one after another (a fresh process each) with a limit of LIMIT
# seconds per problem (default 60), then builds the page.
const HERE = @__DIR__
const LIMIT = get(ARGS, 1, "60")
mkpath(joinpath(HERE, "results"))

# the git commit of each package used from a local checkout, and whether it has
# uncommitted changes (other than this benchmark's results). Found before anything is
# written.
function commit(dir)
    try
        c = readchomp(`git -C $dir rev-parse --short HEAD`)
        changes = readchomp(`git -C $dir status --porcelain --untracked-files=no -- . ':(exclude)bench/cas/results'`)
        c * (isempty(changes) ? "" : " (with uncommitted changes)")
    catch
        "unknown"
    end
end
checkouts = [("SymJ", joinpath(HERE, "..", "..")), ("GruntzLimit", joinpath(HERE, "..", "..", "..", "GruntzLimit")),
             ("RischIntegrate", joinpath(HERE, "..", "..", "..", "RischIntegrate"))]
open(joinpath(HERE, "results", "commits.tsv"), "w") do io
    foreach(((name, dir),) -> println(io, name, "\t", commit(dir)), checkouts)
end

for (name, env) in (("symj", "env_symj"), ("symbolics", "env_symbolics"), ("sympy", "env_sympy"), ("giac", "env_giac"))
    @info "running $name"
    run(`julia $(joinpath(HERE, "supervise.jl")) $name $env $LIMIT`)
end
run(`julia $(joinpath(HERE, "build.jl"))`)
