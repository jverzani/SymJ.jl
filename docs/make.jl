using SymJ
using Documenter

DocMeta.setdocmeta!(SymJ, :DocTestSetup, :(using SymJ); recursive=true)

makedocs(;
    modules=[SymJ],
    authors="jverzani <jverzani@gmail.com> and contributors",
    sitename="SymJ.jl",
    format=Documenter.HTML(;
        canonical="https://jverzani.github.io/SymJ.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jverzani/SymJ.jl",
    devbranch="main",
)
