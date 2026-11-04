# Exact numbers and canonical arithmetic (included by runtests.jl)

@testset "SymNumber" begin
    a = SymNumber(typemax(Int))
    @test (a + a) == SymNumber(2 * big(typemax(Int)))
    @test (a + a) - a == a
    @test (a + a - a).tag == SymJ.T_SMALL        # demotes back to Int
    @test SymNumber(1 // 3) + SymNumber(2 // 3) == SymNumber(1)
    @test SymNumber(2 // 3)^-2 == SymNumber(9 // 4)
    @test SymNumber(3)^100 == SymNumber(big(3)^100)
    @test !isequal(SymNumber(1), SymNumber(1.0))
    @test SymNumber(1) == SymNumber(1.0)
    @test SymNumber(1 // 2) < SymNumber(2 // 3)
end

@testset "canonical sums" begin
    @test x + x == 2x
    @test x - x == 0
    @test x + y == y + x
    @test (x + y) + z == x + (y + z)
    @test 2x + 3y - 2x == 3y
    @test kind(x + y) == ADD
    @test 2 * (x + y) == 2x + 2y          # numbers distribute over a sum
    @test (x + y) - (x + y) == 0
    @test x + 1 - 1 == x
    @test sum([x, y, x, 1, -1]) == 2x + y
end

@testset "canonical products and powers" begin
    @test x * x == x^2
    @test x * y == y * x
    @test x^2 * x^3 == x^5
    @test x * y / y == x
    @test x^2 / x^2 == 1
    @test (x * y)^2 == x^2 * y^2
    @test (x^2)^3 == x^6
    @test kind((x^2)^(1 // 2)) == POW      # not simplified to x
    @test (x^(1 // 2))^2 == x
    @test kind((x + y)^2) == POW          # no automatic expansion
    @test x^0 == 1
    @test x^1 == x
    @test 1^x == 1
    @test exp(x) * exp(y) == exp(x + y)
    @test exp(log(x)) == x
end

@testset "exact numbers" begin
    # an exact base to an infinite power that is 0 is an exact 0 (SymPy #27605)
    @test Sym(0)^Sym(Inf) === SymJ.ZERO
    @test (Sym(1) / 2)^Sym(Inf) === SymJ.ZERO
    @test Sym(2)^Sym(-Inf) === SymJ.ZERO
    @test Sym(2)^Sym(Inf) == Inf
    @test sqrt(Sym(4)) == 2
    @test sqrt(Sym(8)) == 2 * sqrt(Sym(2))
    @test sqrt(Sym(2)) * sqrt(Sym(2)) == 2
    @test sqrt(Sym(2))^3 == 2 * sqrt(Sym(2))
    @test Sym(4)^(1 // 4) == sqrt(Sym(2))
    @test Sym(27)^(2 // 3) == 9
    @test Sym(2)^(-1 // 2) == sqrt(Sym(2)) / 2
    @test sqrt(Sym(1 // 4)) == Sym(1 // 2)
    @test sqrt(Sym(-1)) == IM
    @test IM * IM == -1
    @test IM^4 == 1
    @test sqrt(Sym(-4)) == 2IM
    @test Sym(2)^100 == Sym(big(2)^100)
    @test Sym(1) / 3 == Sym(1 // 3)
    # a power too large to compute stays a power (SymPy #28112)
    @test kind(Sym(5)^(Sym(5)^(Sym(5)^5))) == POW
    @test Sym(0)^0 == 1
    @test Sym(-1)^(Sym(10)^30 + 1) == -1
    @test Sym(1 // 2)^(-3) == 8
    # numer_denom of a long sum stays shallow (deep nesting overflowed the stack, #27957)
    let depth(e) = isempty(e.args) ? 0 : 1 + maximum(depth, e.args)
        @test depth(numer_denom(sum(1 / (x + k) for k in 1:100))[1]) < 10
        @test cancel(sum(1 / (x + k) for k in 1:20) - (numer_denom(sum(1 / (x + k) for k in 1:20))[1] /
                                                      numer_denom(sum(1 / (x + k) for k in 1:20))[2])) == 0
    end
    @test N(sqrt(Sym(2))) ≈ sqrt(2)
    @test N(IM) == im
    @test Sym(π) == PI
    @test N(Sym(π) / 4) ≈ π / 4
    # checked in higher precision where Float64 loses the value (SymPy #28602, #28579, #28280)
    let h = Sym(13) / Sym(10)^46, e = -177682458096 * sqrt(Sym(3)) - 84161668320 * sqrt(7 - 4sqrt(Sym(3))) +
            48590761860 * sqrt(Sym(3)) * sqrt(7 - 4sqrt(Sym(3))) + 307755045036
        @test N(log1p(h)) ≈ 1.3e-45
        @test N(exp(h) - 1) ≈ 1.3e-45
        @test N(log(1 / (Sym(10)^-22 + 1))) ≈ -1e-22
        @test N(e) ≈ 1.2263867743179892e-6
        @test 20 + log(1 + Sym(10)^-20) > 20
        @test 1 + Sym(10)^-20 > 1
        @test !(exp(Sym(10)^-20) <= 1)
    end
    @test Sym(0.5) + Sym(0.25) == Sym(0.75)
end

@testset "one-argument + and *" begin
    @test +(x) == x                      # as for numbers (RischIntegrate builds *(a))
    @test *(x) == x
end

@testset "unwrap_const, float" begin
    @test SymJ.unwrap_const(Sym(2)) === 2
    @test SymJ.unwrap_const(Sym(1) / 2) === 1 // 2
    @test SymJ.unwrap_const(PI) === π
    @test SymJ.unwrap_const(sqrt(Sym(2))) === sqrt(2)
    @test SymJ.unwrap_const(Sym(2.5)) === 2.5
    e = x + 1
    @test SymJ.unwrap_const(e) === e
    @test SymJ.unwrap_const(3) === 3
    @test SymJ.unwrap_const("a") == "a"
    @test map(SymJ.unwrap_const, [Sym(2), x + 1]) == Any[2, x + 1]
    @test float(Sym(1) / 4) === 0.25
    @test float(PI) === Float64(π)
    @test float(sqrt(Sym(2)) + 1) ≈ sqrt(2) + 1
    @test float(2 + 3IM) === 2.0 + 3.0im
    @test float(e) === e
    @test float(Sym) === Sym
    @test float.([Sym(1) / 2, PI]) == [0.5, Float64(π)]
    @test isequal(float.([x, Sym(1) / 2]), Any[x, 0.5])
end

@testset "iteration (like Number)" begin
    for ex in (x, x + y, sin(x), Sym(2), Poly(x + 1, x))
        @test collect(ex) == fill(ex)
        @test length(ex) == 1 && size(ex) == () && ndims(ex) == 0 && !isempty(ex)
        @test eltype(ex) == typeof(ex)
        @test [e for e in ex] == fill(ex)
        @test first(ex) == ex && last(ex) == ex && only(ex) == ex
        @test ex[] == ex && ex[1] == ex && ex[1, 1] == ex
        @test_throws BoundsError ex[2]
        @test ex in ex
        (a,) = ex
        @test a == ex
    end
    @test sum(x + y) == x + y        # sum over a one-element iterable
    @test [x; y] == [x, y]
    @test x .+ [1, 2] == [x + 1, x + 2]
end
