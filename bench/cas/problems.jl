# All the problems, as (area, op, args, expect, source); see problems_basic.jl for the
# operations and what `expect` means
[
    include(joinpath(@__DIR__, "problems_basic.jl"));
    [("Limits", :limit, (f, x0, dir), ans, src) for (f, x0, dir, ans, src) in include(joinpath(@__DIR__, "problems_limits.jl"))];
    [("Integrals", :integrate, (f, x), nothing, "Apostol $i") for (i, (f, x)) in enumerate(include(joinpath(@__DIR__, "problems_integrals.jl")))]
]
