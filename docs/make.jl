using SymJ
using Documenter, MaterialDocs
using Latexify      # results in the examples are shown typeset (SymJLatexifyExt)

DocMeta.setdocmeta!(SymJ, :DocTestSetup, :(using SymJ); recursive=true)

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
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Implementation" => "implementation.md",
        "Design notes" => "design.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/jverzani/SymJ.jl",
    devbranch="main",
)
