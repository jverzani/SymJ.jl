using SymJ
using Documenter, MaterialDocs
using Latexify      # results in the examples are shown typeset (SymJLatexifyExt)

DocMeta.setdocmeta!(SymJ, :DocTestSetup, :(using SymJ); recursive=true)

# The reference, one page per topic in the tutorial's order. A page lists the source
# files whose docstrings it shows: they appear in the order of the files, and within a
# file by kind (module, type, function, macro), then by name. `first` puts some names at
# the top, in the order given (each must be documented in one file only; a name with
# docstrings in several files, such as `solve` or `diff`, stays with its files). `see`
# links to names documented on other pages (a docstring is shown on one page only).
#
# A new docstring in a listed file needs nothing here. A new source file with docstrings
# is added to a page; until it is, `makedocs` fails (`checkdocs = :all`) naming the
# docstrings not shown anywhere. The pages are written to src/reference/ at each build
# (not committed).
const REFERENCE = [
    (title = "Expressions and evaluation", file = "expressions",
     files = ["SymJ.jl", "types.jl", "numbers.jl", "functions.jl", "arithmetic.jl", "relations.jl",
              "manipulation.jl", "identify.jl", "lambdify.jl"],
     first = ["SymJ", "@syms", "symbols", "Sym", "SymFunction", "subs", "N"]),
    (title = "Polynomials", file = "polynomials",
     files = ["poly.jl", "polydiv.jl", "polygcd.jl", "squarefree.jl", "factor.jl", "decompose.jl", "resultant.jl", "groebner.jl",
              "realroots.jl"]),
    (title = "Rational functions", file = "rational",
     files = ["rational.jl", "apart.jl"]),
    (title = "Simplification", file = "simplification",
     files = ["simplify.jl", "gammasimp.jl", "rewrite.jl"],
     first = ["simplify", "refine"]),
    (title = "Assumptions", file = "assumptions",
     files = ["assumptions.jl", "cases.jl"]),
    (title = "Solving equations", file = "solving",
     files = ["solve.jl", "linsolve.jl", "intervals.jl", "ineqsystems.jl"],
     see = ["rsolve", "dsolve"]),
    (title = "Sums and recurrences", file = "sums",
     files = ["summation.jl", "recurrence.jl"]),
    (title = "Matrices", file = "matrices",
     files = ["linalg.jl"]),
    (title = "Piecewise functions", file = "piecewise",
     files = ["piecewise.jl"]),
    (title = "Calculus", file = "calculus",
     files = ["limit.jl", "calculus.jl", "series.jl", "laplace.jl"]),
    (title = "Differential equations", file = "ode",
     files = ["ode.jl"]),
    (title = "Pattern matching", file = "patterns",
     files = ["patterns.jl"],
     first = ["pattern", "@pattern", "rule", "@rule"]),
]

function write_reference(dir)
    rm(dir; force = true, recursive = true)
    mkpath(dir)
    meta = "```@meta\nCurrentModule = SymJ\nEditURL = nothing\n```\n"
    for p in REFERENCE
        isempty(p.files) && error("reference page $(p.title): no files (an empty list would mean all of them)")
        first_ = get(p, :first, String[])
        io = IOBuffer()
        println(io, "# ", p.title, "\n\n", meta)
        see = get(p, :see, String[])
        isempty(see) || println(io, "See also ", join(("[`$n`](@ref)" for n in see), ", "), ".\n")
        isempty(first_) || println(io, "```@docs\n", join(first_, "\n"), "\n```\n")
        objs = join(("getfield(SymJ, Symbol(\"$n\"))" for n in first_), ", ")
        println(io, "```@autodocs\nModules = [SymJ]\nPages = [", join(repr.(p.files), ", "), "]")
        isempty(first_) || println(io, "Filter = t -> !any(o -> t === o, ($objs,))")
        println(io, "```")
        write(joinpath(dir, p.file * ".md"), take!(io))
    end
    write(joinpath(dir, "index.md"), """
        # Reference

        $meta
        Every documented name, by page. The pages follow the tutorial's order.

        ```@index
        Pages = [$(join((repr(p.file * ".md") for p in REFERENCE), ", "))]
        ```
        """)
end
write_reference(joinpath(@__DIR__, "src", "reference"))

makedocs(;
    modules=[SymJ],
    authors="jverzani <jverzani@gmail.com> and contributors",
    sitename="SymJ.jl",
    format=Material3(;
        theme=:default,
        dark_mode=:toggle,
        prettyurls=get(ENV, "CI", nothing) == "true",
        edit_link="main",
    ),
    checkdocs=:all,
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Implementation" => "implementation.md",
        "Design notes" => ["design/index.md", "design/algebra.md", "design/calculus.md"],
        "Reference" => ["reference/index.md"; ["reference/$(p.file).md" for p in REFERENCE]],
    ],
)

deploydocs(;
    repo="github.com/jverzani/SymJ.jl",
    devbranch="main",
)
