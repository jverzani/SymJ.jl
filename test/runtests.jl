using SymJ
using SymJ: PI, E, IM, ADD, MUL, POW, FN, NUMBER, kind
using SymJ: @parameters, parameter, isparameter, free_parameters   # internal for now
using TermInterface
using Test
import ForwardDiff
using Latexify: latexify

@syms x y z f() g()
_add_many_test(f, k, n) = sum(subs(f, k => i) for i in 1:n)

@testset "SymJ.jl" begin
    # one file per area, each compiled on its own (one huge testset function took LLVM
    # minutes, or hung, on Julia 1.13)
    for file in ("numbers.jl", "functions.jl", "expressions.jl", "polynomials.jl", "simplify.jl", "solve.jl",
                 "calculus.jl", "sums.jl", "odes.jl", "linalg.jl", "compile.jl", "review.jl")
        include(file)
    end

    @testset "ambiguities" begin
        @test isempty(Test.detect_ambiguities(SymJ; recursive = true))
    end

    include("extensions.jl")
end
