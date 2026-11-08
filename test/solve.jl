# Solving equations (included by runtests.jl)

@testset "solve" begin
    # a candidate from squaring is dropped when it fails for every sampled parameter
    # value (SymPy #27381); one valid for some parameter values is kept
    let a = Sym(:a), b = Sym(:b), x = Sym(:x)
        @test solve(a * (2sqrt(1 - b^2) + 1) - b * (2sqrt(1 - a^2) + 1) ~ 0, b) == [a]
        @test solve(sqrt(x) ~ a, x) == [a^2]
        @test length(solve(sqrt(x + 1) ~ x - a, x)) == 2
    end
    # a system's components where a divisor of its solutions vanishes (SymPy #27008, #27458)
    let (x0, x1, x2, x3, x4, x5) = Sym.([:x0, :x1, :x2, :x3, :x4, :x5]), (x, y, z) = Sym.([:x, :y, :z])
        r = solve([x0 * x1 - x2 ~ 0, x3 + x4 - x5 ~ 0, x0 - x5 ~ 1], [x0, x1, x2, x3, x4, x5])
        @test length(r) == 2
        @test any(s -> s[x5] == -1 && s[x0] == 0 && s[x2] == 0 && s[x1] == x1, r)
        @test length(solve([x * y ~ 0, x + z ~ 1], [x, y, z])) == 2
    end
    # unknowns that are function calls or derivatives (implicit differentiation)
    let (a, b, c, d) = symbols("a b c d"), u = SymFunction(:u), v = SymFunction(:v)
        ex2 = diff(x * u(x) - (a * x^3 + b * x^2 + c * x + d), x)
        du = diff(u(x), x)
        @test solve(ex2 ~ 0, du) == [(3a * x^2 + 2b * x + c - u(x)) / x]
        @test solve(diff(x^2 + u(x)^2 - 1, x) ~ 0, du) == [-x / u(x)]
        @test solve(x * u(x) + sin(x) ~ 1, u(x)) == [(1 - sin(x)) / x]
        @test solve(du + u(x) ~ x, u(x)) == [x - du]                 # the derivative is held
        @test solve([u(x) + v(x) ~ 3, u(x) - v(x) ~ 1], [u(x), v(x)]) == [Dict(u(x) => 2, v(x) => 1)]
        @test solve(u(x)^2 ~ 4, u(x); assume = @assume(u(x) ≫ 0)) == [2]
        @test (@assuming u(x) ≫ 0 begin solve(u(x)^2 ~ 4, u(x)) end) == [2]
        @test_throws ArgumentError solve(x ~ 1, x + 1)
    end
    # several logarithms of x: combined (rational coefficients), candidates checked
    @test solve(log(x) + log(x - 1) ~ log(Sym(6)), x) == [3]       # not -2 (complex logs)
    @test solve(2log(x) - log(x + 2) ~ 0, x) == [2]
    @test solve(log(x) + log(x + 3) ~ 2log(Sym(2)), x) == [1]
    @test solve(log(x^2) - log(x) ~ log(Sym(5)), x) == [5]
    @test solve(log(x) - log(x - 1) ~ 1, x) == [Sym(ℯ) / (Sym(ℯ) - 1)]
    @test_throws ArgumentError solve(log(x) + sin(x) ~ 1, x)
    @test_throws ArgumentError solve(x * log(x) + log(x - 1) ~ 1, x)
    @test solve(exp(x) ~ -1, x) == [IM * PI]
    # N outside the real domain, with an exact integer argument
    @test N(log(Sym(-3))) ≈ log(complex(-3.0)) && N(log(Sym(-2)) + log(Sym(-3)) - log(Sym(6))) ≈ 2π * im
    @parameters a b c
    satisfies(eq, rs) = all(r -> abs(N(subs(lhs(eq) - rhs(eq), x => r))) < 1e-8 * max(1, abs(N(r))), rs)
    cases = [(2x + 3 ~ 0, 1), (x^2 - 2 ~ 0, 2), (x^2 ~ -1, 2), (x^2 + x + 1 ~ 0, 2), (x^4 ~ 1, 4),
             (x^5 ~ 2, 5), (x^3 - 2x - 5 ~ 0, 3), (x^3 - 3x + 1 ~ 0, 3), (x^4 + x + 1 ~ 0, 4),
             (x^4 - 4x^3 + x + 7 ~ 0, 4), (x^4 + 1 ~ 0, 4), (x^4 - 10x^2 + 1 ~ 0, 4), (x^5 - x - 1 ~ 0, 5),
             ((x^2 - 1) / (x - 1) ~ 0, 1), (1 / x ~ 0, 0), (x^6 - 1 ~ 0, 6), (0.5x^2 - 1 ~ 0, 2),
             ((x - 1)^3 * (x + 2) ~ 0, 2), (x^3 ~ x, 3), (2x^2 ~ x + 1, 2)]
    for (eq, nroots) in cases
        r = solve(eq, x)
        @test length(r) == nroots
        @test satisfies(eq, r)
    end
    @test solve(x^2 - 2 ~ 0, x) == [-sqrt(Sym(2)), sqrt(Sym(2))]      # exact, sorted
    @test solve(x^2 ~ -1, x) == [-IM, IM]
    @test solve(2x^2 ~ x + 1, x) == [Sym(-1 // 2), Sym(1)]
    @test solve(x^3 - 3x + 1 ~ 0, x) == [2cos(10PI / 9), 2cos(4PI / 9), 2cos(2PI / 9)]  # real form
    @test all(r -> imag(N(r)) == 0, solve(x^3 - 3x + 1 ~ 0, x))
    # a triple root, and Δ0 = 0, without Cardano's 0/0 (coefficients in im and y)
    @test solve((x - 2IM * y)^3 ~ 0, x) == [2IM * y]
    @test solve((x - 2IM * y)^3 * (x + y) ~ 0, x) == [-y, 2IM * y]
    @test solve(x^3 + 3x^2 + 3x + 1 - y ~ 0, x)[1] == -(-y)^(1 // 3) - 1
    @test all(r -> !SymJ._hassymbol(r) && N(r) isa Number, solve(x^5 - x - 1 ~ 0, x))  # numeric (a + b*im)
    # symbolic coefficients
    # (cases where a coefficient vanishes)
    cs = solve(a * x + b ~ 0, x)
    @test cs isa SymJ.Cases
    @test cs(a => 2, b => 4) == [-2]
    @test cs(a => 0, b => 1) == Sym[]
    @test last(cs.cases).second == [-b / a]
    @test solve(a * x + b ~ 0, x; assume = (a ≫ 0,)) == [-b / a]
    @test solve(x^2 ~ a, x) == [-sqrt(a), sqrt(a)]
    cs = solve(exp(a * x) ~ 2, x)
    @test cs(a => 0) == Sym[]
    @test cs(a => 1) == [log(Sym(2))]
    @test solve(x / a ~ 1, x)(a => 0) == Sym[]
    @test solve(x / a ~ 1, x)(a => 3) == [3]
    cs = solve(a * x^2 + b * x + c ~ 0, x)
    @test cs(a => 0, b => 2, c => -4) == [2]
    @test cs(a => 0, b => 0, c => 1) == Sym[]
    q = last(cs.cases).second
    @test length(q) == 2 && all(r -> cancel(subs(a * x^2 + b * x + c, x => r)) == 0 ||
                                 abs(N(subs(a * x^2 + b * x + c, x => r, a => 2, b => 3, c => -5))) < 1e-12, q)
    @test solve(x^2 - a^2 ~ 0, x) == [-a, a] || Set(solve(x^2 - a^2 ~ 0, x)) == Set([-a, a])
    # errors
    @test_throws ArgumentError solve(x^2 + 1, x)          # an equation is required
    @test string(solve(x ≪ 1, x)) == "(-∞, 1)"            # inequalities: sets
    @test_throws ArgumentError solve(x ~ x, x)
    # through one function of x: substitute, solve, invert; extraneous roots dropped
    @test solve(sin(x) ~ 1 // 2, x) == [PI / 6, 5PI / 6]
    @test solve(cos(x) ~ 0, x) == [PI / 2, 3PI / 2]
    @test solve(tan(x) ~ 1, x) == [PI / 4]
    # the solutions in [0, 2π], each once (0 brings 2π)
    @test solve(cos(x) ~ 1, x) == [0, 2PI]
    @test solve(cos(x) ~ -1, x) == [PI]
    @test solve(cos(x) ~ -1 // 2, x) == [2PI / 3, 4PI / 3]
    @test solve(sin(x) ~ -1, x) == [3PI / 2]
    @test solve(sin(x) ~ 0, x) == [0, PI, 2PI]
    @test solve(sin(x) ~ -1 // 2, x) == [7PI / 6, 11PI / 6]
    @test solve(diff(sin(x) - x, x) ~ 0, x) == [0, 2PI]
    @test solve(sin(2x) ~ 1, x) == [PI / 4]
    @test solve(sin(x)^2 ~ 1 // 4, x) == [PI / 6, 5PI / 6, 7PI / 6, 11PI / 6]
    @test solve(exp(x) ~ 2, x) == [log(Sym(2))] && solve(log(x) ~ 2, x) == [exp(Sym(2))]
    @test solve(exp(x) ~ a, x) == [log(a)]
    @test solve(sqrt(x) ~ 3, x) == [Sym(9)] && solve(sqrt(x) ~ -1, x) == Sym[]
    @test solve(2^x ~ 8, x) == [Sym(3)]                     # an exact integer log
    # powers of one exponential, also with integer bases that are powers of one integer
    @test solve(exp(2x) - 3exp(x) + 2 ~ 0, x) == [Sym(0), log(Sym(2))]
    @test solve(exp(x) + exp(-x) ~ 2, x) == [Sym(0)]
    @test solve(4^x - 3 * 2^x + 2 ~ 0, x) == [Sym(0), Sym(1)] && solve(8^x ~ 2 * 4^x, x) == [Sym(1)]
    # a square root mixed with x: squared, extraneous roots dropped
    @test solve(sqrt(x + 1) ~ x - 1, x) == [Sym(3)]                # 0 is extraneous
    @test solve(sqrt(x) + x ~ 6, x) == [Sym(4)] && solve(x - sqrt(x) ~ 0, x) == [Sym(0), Sym(1)]
    # sin, cos, tan, … of one argument: through t = tan(x/2), solutions in [0, 2π) (SymPy #29745)
    @test solve(sin(x) + cos(x) ~ 1, x) == [0, PI / 2]
    # x = π (t = ∞) only where it is a solution (SymPy #28213)
    @test solve(sin(x) + cos(x) ~ -1, x) == [PI, 3PI / 2]
    @test !any(r -> isequal(r, PI), solve(sin(x) + a * cos(x) ~ a, x))
    @test !any(r -> isequal(r, PI), solve(sin(x) + 2cos(x) * abs(a) ~ 1, x))
    # a float or symbolic exponent: the principal root (SymPy #28772)
    @test solve(x^a ~ b, x) == [b^(1 / a)]
    @test N(solve(x^1.11 ~ 2, x)[1])^1.11 ≈ 2
    # abs for a real unknown: ±r, or the cases of the argument's sign (SymPy #28742)
    @test solve(abs(x) ~ 1, x; assume = (x ∈ ℝ,)) == [-1, 1]
    @test solve(abs(x - 2) ~ 3, x; assume = (x ∈ ℝ,)) == [-1, 5]
    @test solve(abs(2x + 1) ~ x + 4, x; assume = (x ∈ ℝ,)) == [-5 // 3, 3]
    @test solve(abs(x - 1) + abs(x + 1) ~ 4, x; assume = (x ∈ ℝ,)) == [-2, 2]
    @test solve(abs(x) ~ a, x; assume = (x ∈ ℝ, a ∈ ℝ))(a => -1) == Sym[]
    @test_throws ArgumentError solve(abs(x) ~ 1, x)                   # a circle in ℂ
    # a root where a denominator vanishes for some parameter values: dropped in that case
    cs = solve(sin(x) / (cos(x) - a) ~ 0, x)
    @test cs(a => 1) == [PI]
    @test cs(a => -1) == [0, 2PI]
    @test cs(a => 0) == [0, PI, 2PI]
    @test solve(x / (x - a) ~ 0, x)(a => 0) == Sym[]
    # a vanishing leading coefficient in tan(x/2): x = π is a root then (b + c = 0)
    let (b, c) = symbols("b c"), cs = solve(a * sin(x) + b * cos(x) ~ c, x)
        @test cs(a => 1, b => 2, c => -2)[1] == PI
        @test cs(a => 0, b => 2, c => -2) == [PI]
        @test length(cs(a => 1, b => 1, c => 1)) == 2
    end
    # leading coefficients in one function of x: c cos(x) = -c, a eˣ = 1
    let c = Sym(:c)
        @test solve(c * cos(x) ~ -c, x)(c => 2) == [PI]
        @test_throws ArgumentError solve(c * cos(x) ~ -c, x)(c => 0)      # every x
    end
    @test solve(a * exp(x) ~ 1, x)(a => 0) == Sym[]
    @test solve(2sin(x) - 2sqrt(Sym(3)) * cos(x) + sqrt(Sym(3)) * tan(x) - 3 ~ 0, x) == [PI / 3, 5PI / 6, 7PI / 6, 4PI / 3]
    @test solve(tan(x) + cot(x) ~ 2, x) == [PI / 4, 5PI / 4]
    @test solve(sin(2x) + cos(2x) ~ 1, x) == [0, PI / 4]
    @test_throws ArgumentError solve(sin(x) + cos(2x) ~ 1, x)    # two arguments: not yet
    # Lambert W: (αx + β) e^(γx + δ) = K, both real branches when there are two
    @test solve(x * exp(x) ~ 1, x) == [lambertw(Sym(1))] && solve(x ~ exp(-x), x) == [lambertw(Sym(1))]
    @test Set(solve(x * exp(x) ~ -1 // 5, x)) == Set([lambertw(Sym(-1 // 5)), lambertw(Sym(-1 // 5), -1)])
    @test solve(x * exp(x) ~ 0, x) == [Sym(0)]
    for eq in (2^x ~ 3x, exp(x) ~ x + 2, (2x + 1) * exp(3x) ~ 4)
        rs = solve(eq, x)
        @test !isempty(rs) && satisfies(eq, rs)
    end
    @test_throws ArgumentError solve(x * exp(x^2) ~ 1, x)       # not a Lambert form

    # linear systems: one Dict in a vector; free variables omitted; [] when inconsistent
    holds(eqs, d) = all(eq -> cancel(subs(lhs(eq) - rhs(eq), d)) == 0 ||
                              abs(N(subs(lhs(eq) - rhs(eq), d))) < 1e-12, eqs)
    @test solve([x + y ~ 3, x - y ~ 1], [x, y]) == [Dict(x => 2, y => 1)]
    # not polynomial but triangular: one unknown at a time (SymPy #28064)
    @test solve([exp(x) ~ 2, x + y ~ 3], [x, y]) == [Dict(x => log(Sym(2)), y => 3 - log(Sym(2)))]
    @test solve([sin(x) ~ 0, x * y ~ 1], [x, y]) == [Dict(x => PI, y => 1 / PI), Dict(x => 2PI, y => 1 / (2PI))]
    let (T, p, h, c) = symbols("T p h c")
        @test solve([T * c ~ h, log(p) ~ 2T], [T, p])(c => 2) == [Dict(T => h / 2, p => exp(h))]
    end
    @test solve([(x^2 - 1) / (x - 1) ~ 2, y ~ x], [x, y]) == Dict{Sym,Sym}[]   # x = 1 is not in the domain
    @test solve([(x^2 - 1) / (x - 1) ~ 3, x + y ~ 3], [x, y]) == [Dict(x => 2, y => 1)]
    sys3 = [2x + 3y - z ~ 1, x - y + 2z ~ 4, 3x + y + z ~ 6]
    @test solve(sys3, [x, y, z]) == [Dict(x => Sym(11 // 5), y => Sym(-1), z => Sym(2 // 5))]
    @test solve([x + y + z ~ 1, x - z ~ 0], [x, y, z]) == [Dict(x => z, y => 1 - 2z, z => z)]
    @test solve([x + y ~ 1, x + y ~ 2], [x, y]) == Dict{Sym,Sym}[]
    @test solve([x + y ~ 1, 2x + 2y ~ 2], [x, y]) == [Dict(x => 1 - y, y => y)]
    @test_throws ArgumentError solve([x + y ~ 1, x - y ~ 2], [x, y]; assume = (x ≫ 0,))   # not yet for systems
    @test solve((x + y ~ 3, x - y ~ 1), (x, y)) == [Dict(x => 2, y => 1)]       # tuples
    for (eqs, vs) in (([a * x + y ~ 1, x - y ~ b], [x, y]), ([a * x + b * y ~ 1, b * x - a * y ~ 0], [x, y]),
                      ([x / 2 + y / 3 ~ 1, x ~ y], [x, y]), ([x + y ~ sin(a), x - y ~ cos(a)], [x, y]),
                      ([x + 2y + z ~ a, y - z ~ b, x + z ~ 0], [x, y, z]))
        ds = solve(eqs, vs)
        ds isa SymJ.Cases && (ds = last(ds.cases).second)       # the generic case
        @test length(ds) == 1 && length(ds[1]) == length(vs) && holds(eqs, ds[1])
    end
    @test last(solve([a * x + b * y ~ 1, b * x - a * y ~ 0], [x, y]).cases).second == [Dict(x => a / (a^2 + b^2), y => b / (a^2 + b^2))]
    # pivots that can vanish give cases (SymPy #30224: a = 0 kept)
    cs = solve([a * x + y ~ 1, x - y ~ 2], [x, y])
    @test cs(a => 0) == [Dict(x => 3, y => 1)]
    @test cs(a => -1) == Dict{Sym,Sym}[]
    @test cs(a => 1) == [Dict(x => Sym(3 // 2), y => Sym(-1 // 2))]
    cs = solve([x + y ~ 1, x + y ~ a], [x, y])
    @test cs(a => 1) == [Dict(x => 1 - y, y => y)]
    @test cs(a => 2) == Dict{Sym,Sym}[]
    cs = solve([a * x + b * y ~ 0, x^2 + y^2 ~ 1], [x, y])
    @test Set(cs(a => 0, b => 1)) == Set([Dict{Sym,Sym}(x => -1, y => 0), Dict{Sym,Sym}(x => 1, y => 0)])
    d = only(solve([0.5x + y ~ 1, x - y ~ 0.25], [x, y]))
    @test N(d[x]) ≈ 5 / 6 && N(d[y]) ≈ 7 / 12

    # nonlinear systems: through a lex Gröbner basis
    sameset(A, B) = length(A) == length(B) && all(a -> any(b -> keys(a) == keys(b) &&
                                                     all(k -> iszero(expand(a[k] - b[k])), keys(a)), B), A)
    @test sameset(solve([x^2 + y^2 ~ 1, x ~ y], [x, y]),
                  [Dict(x => sqrt(Sym(2)) / 2, y => sqrt(Sym(2)) / 2), Dict(x => -sqrt(Sym(2)) / 2, y => -sqrt(Sym(2)) / 2)])
    @test sameset(solve([x * y ~ 1, x + y ~ 3], [x, y]),
                  [Dict(x => (3 + sqrt(Sym(5))) / 2, y => (3 - sqrt(Sym(5))) / 2),
                   Dict(x => (3 - sqrt(Sym(5))) / 2, y => (3 + sqrt(Sym(5))) / 2)])
    @test length(solve([x^2 ~ 1, y^2 ~ 1], [x, y])) == 4                     # not in shape position
    # an algebraic number in the coefficients is not a parameter (SymPy #29750)
    let (s, c) = (Sym(:s), Sym(:c)), r3 = sqrt(Sym(3))
        eqs = [2s * c - 2r3 * c^2 + r3 * s - 3c ~ 0, s^2 + c^2 - 1 ~ 0]
        sols = solve(eqs, [s, c])
        @test length(sols) == 4
        @test all(d -> all(e -> iszero(simplify(subs(e.lhs - e.rhs, d...))), eqs), sols)
        @test_throws ArgumentError groebner([e.lhs for e in eqs], s, c)
    end
    # (sqrt(2)*sqrt(3) unmerged: compared after simplify)
    @test sameset(map(d -> Dict(k => simplify(v) for (k, v) in d), solve([x^2 + y^2 ~ 1, x ~ sqrt(Sym(2)) * y], [x, y])),
                  [Dict(x => sqrt(Sym(6)) / 3, y => sqrt(Sym(3)) / 3), Dict(x => -sqrt(Sym(6)) / 3, y => -sqrt(Sym(3)) / 3)])
    # float solutions refined by Newton's method: small components keep their relative
    # accuracy (SymPy #29140)
    let
        @syms H c1 c2 c3 NaOH Na OH
        K1, K2, Kw = Sym(10)^(-189 // 20), Sym(10)^(-64 // 5), Sym(10)^(-14)
        eqs = [K1 * c1 ~ H * c2, K2 * c2 ~ H * c3, OH * Na ~ NaOH, H * OH ~ Kw,
               Na + H ~ OH + c2 + 2c3, c1 + c2 + c3 ~ 9 // 2000, NaOH + Na ~ 1 // 8]
        sols = solve(eqs, [H, c1, c2, c3, NaOH, Na, OH])
        phys = only(filter(s -> all(v -> N(s[v]) > 0, (H, c1, c2, c3, NaOH, Na, OH)), sols))
        @test isapprox(N(phys[H]), 9.458052236307555e-14; rtol = 1e-10)
        @test isapprox(N(phys[H] * phys[OH]), 1e-14; rtol = 1e-10)
    end
    @test sameset(solve([x^2 + y^2 ~ 4, (x - 1)^2 + y^2 ~ 4], [x, y]),
                  [Dict(x => 1 // 2, y => sqrt(Sym(15)) / 2), Dict(x => 1 // 2, y => -sqrt(Sym(15)) / 2)])
    @test isempty(solve([x^2 + y^2 ~ 1, x^2 + y^2 ~ 4], [x, y]))
    @test sameset(solve([x^2 + y^2 ~ 1], [x, y]), [Dict(x => sqrt(1 - y^2), y => y), Dict(x => -sqrt(1 - y^2), y => y)])   # y free
    @test sameset(solve([1 / x + 1 / y ~ 1, x ~ y], [x, y]), [Dict(x => 2, y => 2)])         # rational
    @test sameset(solve([y ~ x^2 / (x - 1), x * (x - 1) ~ 0], [x, y]), [Dict(x => 0, y => 0)])   # pole dropped
    @test length(solve([x + y + z ~ 0, x * y + y * z + z * x ~ 0, x * y * z ~ 1], [x, y, z])) == 6
    @test length(solve([x^2 + y^2 ~ a^2, x ~ y], [x, y])) == 2                # a parameter
    for (eqs, vs, n) in (([x^2 + x * y + y^2 ~ 3, x^2 - y^2 ~ x - 1], [x, y], 4),
                         ([x^2 + y^2 + z^2 ~ 6, x * y ~ 2, y + z ~ 3], [x, y, z], 4),
                         ([0.5x^2 + y^2 ~ 1.5, x ~ 2y], [x, y], 2),                     # floats
                         ([y ~ x^5 - x - 1, y ~ 0], [x, y], 5),                        # numeric roots
                         ([x + 2y + 2z ~ 1, x^2 + 2y^2 + 2z^2 ~ x, 2x * y + 2y * z ~ y], [x, y, z], 4))
        ds = solve(eqs, vs)
        @test length(ds) == n && all(d -> holds(eqs, d), ds)
    end
    @test_throws ArgumentError solve([sin(x) ~ y, x ~ y^2], [x, y])
    @test solve([1 / x ~ 1], [x]) == [Dict(x => 1)]
    # the call's form fixes the result's: one equation and one unknown, or vectors of both
    @test solve(x ~ 2, x) == [Sym(2)]
    @test solve([x ~ 2], [x]) == [Dict(x => 2)]
    @test_throws ArgumentError solve(x ~ 2, [x])
    @test_throws ArgumentError solve([x ~ 2], x)
    @test_throws ArgumentError solve((x ~ 2,), x)
    @test_throws ArgumentError solve(x + y ~ 2, x, y)
    @test string(solve([x ≪ 1], [x])) == "(-∞, 1)"                # a one-unknown system

    # < and friends compare constants as numbers; with symbols, a pointer to ≪
    @test sqrt(Sym(2)) < 2 && Sym(π) > 3 && !(2 ≥ sqrt(Sym(5))) && Sym(1 // 3) <= 0.5
    @test sort([Sym(3), sqrt(Sym(2)), Sym(π)]) == [sqrt(Sym(2)), Sym(3), Sym(π)]
    @test_throws ArgumentError x + y < 2
    @test_throws ArgumentError x >= 1
    @test occursin("≪", (try x < 1 catch e; sprint(showerror, e) end))
    @test occursin("≧", (try x >= 1 catch e; sprint(showerror, e) end))

    # systems of inequalities: one unknown (sets), several (cylindrical, linear)
    @test string(solve([x ≫ 1, x ≪ 3], [x])) == "(1, 3)"
    @test string(solve((x^2 - 1) / (x - 1) ≫ 0, x)) == "(-1, 1) ∪ (1, ∞)"   # 1 is not in the domain
    @test string(solve((x^3 - x) / (x^2 - x) ≧ 0, x)) == "[-1, 0) ∪ (0, 1) ∪ (1, ∞)"
    @test string(solve((y + 1)^2 - y * (y + 2) ≫ 0, x)) == "(-∞, ∞)"   # constant once cancelled
    @test string(solve([x^2 ≪ 4, x ≧ 0], [x])) == "[0, 2)"
    @test string(solve([x^2 ~ 4, x ≫ 0], [x])) == "{2}"
    @test string(solve([x ≶ 2, x^2 ≦ 9], [x])) == "[-3, 2) ∪ (2, 3]"
    C2 = solve([x + y ≪ 2, y ≪ 3], [x, y])
    @test C2 isa SymJ.CylindricalSet && string(C2) == "y ∈ (-∞, 3), x ∈ (-∞, -y + 2)"
    T2 = solve([x ≧ 0, y ≧ 0, x + y ≦ 1], [x, y])
    @test string(T2) == "y ∈ [0, 1], x ∈ [0, -y + 1]"
    @test string(solve([x ≧ 0, y ≧ 0, x + y ≦ 1], [y, x])) == "x ∈ [0, 1], y ∈ [0, -x + 1]"
    @test Dict(x => 0.2, y => 0.3) in T2 && !(Dict(x => 0.8, y => 0.3) in T2) && Dict(x => 0, y => 1) in T2
    @test isempty(solve([x + y ≪ 1, x + y ≫ 2], [x, y]))
    S3 = solve([x + y + z ~ 1, x ≧ 0, y ≧ 0, z ≧ 0], [x, y, z])
    @test string(S3) == "z ∈ [0, 1], y ∈ [0, -z + 1], x = -y - z + 1"
    @test Dict(x => 0.5, y => 0.25, z => 0.25) in S3 && !(Dict(x => 0.5, y => 0.5, z => 0.25) in S3)
    @test string(solve([x ≧ 0, y ≧ 0, z ≧ 0, x + 2y + 3z ≦ 6], [x, y, z])) ==
          "z ∈ [0, 2], y ∈ [0, -3*z/2 + 3], x ∈ [0, -2*y - 3*z + 6]"
    # calling a set: a Bool with every unknown given, else the slice in the rest
    @test C2(x => 1, y => 0.5) === true && C2(x => 1, y => 4) === false
    @test string(C2(y => 1)) == "(-∞, 1)" && C2(y => 1) isa SymJ.SymSet && C2(y => 1)(0)
    T3 = solve([x ≧ 0, y ≧ 0, z ≧ 0, x + 2y + 3z ≦ 6], [x, y, z])
    @test string(T3(z => 1)) == "y ∈ [0, 3/2], x ∈ [0, -2*y + 3]"
    @test string(T3(z => 1, y => 0)) == "[0, 3]" && isempty(T3(x => 7))
    @test T3(x => 1, y => 1, z => 1) === true && T3(x => 2, y => 1, z => 1) === false
    @test isempty(solve([x + y ≪ 1, x + y ≫ 2], [x, y])(x => 0))
    @test_throws ArgumentError C2(z => 1)                                # not an unknown
    S1 = solve(x^2 ≫ 1, x)
    @test S1(2) && !S1(0) && S1(x => 2) && S1(Sym(3) // 2) && S1.var == x
    @test_throws ArgumentError S1(y => 2)
    @test_throws ArgumentError solve([x * y ≪ 1, x ≫ 0], [x, y])        # nonlinear: not yet
    @test_throws ArgumentError solve([x ≶ y, x ≫ 0], [x, y])           # not convex

    # inequalities: unions of intervals with exact endpoints
    S(r) = string(solve(r, x))
    @test S(x^2 ≫ 1) == "(-∞, -1) ∪ (1, ∞)" && S(x^2 ≦ 1) == "[-1, 1]"
    @test S(x^2 ≦ 0) == "{0}" && S(x^2 ≪ 0) == "∅" && S(x^2 + 1 ≫ 0) == "(-∞, ∞)"
    @test S((x - 1) / (x + 2) ≧ 0) == "(-∞, -2) ∪ [1, ∞)" && S(1 / x ≪ 0) == "(-∞, 0)"
    @test S(x^3 - x ≧ 0) == "[-1, 0] ∪ [1, ∞)" && S(x^2 - 2 ≪ 0) == "(-sqrt(2), sqrt(2))"
    @test S(x ≶ 3) == "(-∞, 3) ∪ (3, ∞)" && S((x - 1)^2 ≫ 0) == "(-∞, 1) ∪ (1, ∞)"
    @test S(x^2 ≧ x) == "(-∞, 0] ∪ [1, ∞)"
    P2s, T = solve(x^2 ≫ 1, x), solve(x ≦ 3, x)
    @test P2s isa SymJ.SymSet && !(0 in P2s) && 2 in P2s && -1.5 in P2s && !(-1 in P2s)
    @test string(P2s ∪ T) == "(-∞, ∞)" && string(P2s ∩ T) == "(-∞, -1) ∪ (1, 3]"
    @test string(setdiff(T, P2s)) == "[-1, 1]" && isempty(setdiff(T, T))
    # ⊆, with a set or one of its intervals
    @test solve(x^2 ≪ 1, x) ⊆ T
    @test !(T ⊆ P2s)
    @test P2s.intervals[1] ⊆ P2s
    @test !(solve(x ≦ 1, x) ⊆ P2s.intervals[1])
    @test sqrt(Sym(2)) in solve(x^2 ≦ 2, x) && !(sqrt(Sym(2)) in solve(x^2 ≪ 2, x))
    @test_throws ArgumentError solve(a * x ≫ 1, x)               # parameters: not yet
    @test_throws ArgumentError solve(sin(x) ≫ 0, x)
    @test_throws ArgumentError solve(x^5 + a * x + 1 ~ 0, x)  # degree 5, symbolic

    # decomposable polynomials of degree above 4: exact, through the components
    for p in ((x^2 + x)^3 - 2, (x^3 + x)^2 - 3, (x^2 + 1)^4 + (x^2 + 1) + 1, (x^3 - 3x)^3 - 3(x^3 - 3x) - 1, x^6 - 5x^3 + 1)
        p = expand(p)
        rs = solve(p ~ 0, x)
        @test length(rs) == length(SymJ._xcoeffs_sym(p, x)) - 1
        @test all(r -> abs(N(subs(p, x => r))) < 1e-8, rs)
        @test !any(r -> SymJ.isnumber(r) && SymJ.isfloat(r.num), rs)
    end
    @test (-1 + sqrt(1 + 4 * Sym(2)^(1 // 3))) / 2 in solve(expand((x^2 + x)^3 - 2) ~ 0, x)

    # products: solved factor by factor
    cs = solve(diff((a / x)^x, x) ~ 0, x)                           # (a/x)^x is never 0
    @test last(cs.cases).second == [a / SymJ.E]
    @test first(cs.cases).first == Any[a ~ 0]                       # (0/x)^x is constant
    @test solve(diff(x^x, x) ~ 0, x) == [1 / SymJ.E]
    @test solve(x * exp(x) ~ 0, x) == [0]
    @test solve((x - 1)^2 * exp(2x) ~ 0, x) == [1]
    @test solve(2^x * (x + 1) ~ 0, x) == [-1]
    @test solve(exp(x) * (exp(x) - 3) ~ 0, x) == [log(Sym(3))]
    @test solve(x^2 * cos(x) ~ 0, x) == [0, PI / 2, 3PI / 2]
    @test solve((x - 1) * log(x + 2) ~ 0, x) == [-1, 1]
    @test solve(diff(x * exp(-x), x) ~ 0, x) == [1]                 # a common factor
    @test solve(x * log(x) ~ 0, x) == [1]                           # log(0) is not defined
    @test solve(sin(x) * log(x) ~ 0, x) == [1, PI, 2PI]
    @test solve(sqrt(x) * (x - 1) ~ 0, x) == [0, 1]
    # removable singularities: cancel removes x - 1, but x = 1 is not in the domain
    @test solve((x^2 - 1) / (x - 1) ~ 2, x) == Sym[]
    @test solve((x^2 - 1) / (x - 1) ~ 3, x) == [2]
    @test solve(tan(x) * cos(x) ~ 0, x) == [0]                      # π/2: tan is not defined
    let (a,) = symbols("a")
        @test solve((x^2 - a^2) / (x - a) ~ 2a, x) == Sym[]
    end
    # rational powers of x as powers of one root of x (t = x^(1/10): t³ + t² = 3); the
    # root kept unexpanded, whose numeric check otherwise loses it to cancellation
    let r = solve(x^(3 // 10) + x^(1 // 5) ~ 3, x)
        @test length(r) == 1
        @test N(r[1]) ≈ 4.997466445630026
    end
    # a high degree (u = x^(3/100): degree 181) has numeric roots, not NaN
    let r = solve(x^(543 // 100) + x^(222 // 100) ~ 5033 // 100, x)
        @test N(r[1]) ≈ 2.0204868226226256
        @test all(v -> abs(N(v^(543 // 100) + v^(222 // 100)) - 50.33) < 1e-8, r)
    end
    # several radicals: squared away one at a time, extraneous roots dropped
    @test solve(sqrt(x) + sqrt(x + 5) ~ 5, x) == [4]
    @test solve(sqrt(x + 1) - sqrt(x - 1) ~ 1, x) == [5 // 4]
    @test solve(sqrt(x) - sqrt(x + 1) ~ 2, x) == Sym[]
    @test solve(sqrt(x^2 + 1) ~ sqrt(2x^2 - 3), x) == [-2, 2]
    let (y,) = symbols("y")             # implicit differentiation: where dy/dx = 0
        cs = solve(3sqrt(x^2 + y^2) - x * sqrt(x^2 + y^2) - x * sqrt((3 - x)^2 + (3 - y)^2) ~ 0, x)
        @test last(cs.cases).second == [y, 3y / (2y - 3)]
        @test cs(y => 3 // 2) == [3 // 2]                   # 3y/(2y - 3) has no value there
    end
end

@testset "solve with assumptions" begin
    @syms a b n
    # a float double root split by rounding into a complex pair is taken as real, with a
    # warning (SymPy #30346)
    let eq = -1.30616509926855e-8 * x^3 + 0.00273103448275862 * x ~ 0.480662068965517
        rs = @test_logs (:warn, r"taken as real") solve(eq, x; assume = (x ≧ 0, x ≦ 264))
        @test length(rs) == 1
        @test isapprox(N(only(rs)), 264; rtol = 1e-10)
    end
    @test solve(x^2 ~ 4, x; assume = (x ≫ 0,)) == [2] == solve(x^2 ~ 4, x; assume = @assume(x > 0))
    @test solve(x^2 ~ 2, x; assume = @assume(x in Rationals)) == Sym[]
    @test (@assuming x > 0 begin
        solve(x^2 ~ 4, x)
    end) == [2]
    @test @assuming(x > 0, x != 2, begin solve(x^2 ~ 4, x) end) == Sym[]
    @test (@assuming x >= 1 begin
        @assuming x < 3 begin solve(x^2 ~ 9, x) end              # nested: both apply
    end) == Sym[]
    @test solve(x^2 ~ 4, x; assume = Assumptions(x ≶ 2)) == [-2]
    @test solve(x^2 ~ -1, x; assume = (x ∈ ℝ,)) == Sym[]
    @test solve(x^3 ~ 8, x; assume = (x ∈ ℝ,)) == [2]
    @test solve(x^2 - 5x + 6 ~ 0, x; assume = (x ≫ 2,)) == [3]
    @test solve(x^2 ~ a, x; assume = (x ≫ 0, a ≫ 0)) == [sqrt(a)]
    @test solve(x^2 ~ a, x; assume = (x ≫ 0, a ~ 9)) == [3]        # a substituted
    @test solve(x^2 ~ -1, x) == [-IM, IM]                          # none: unchanged
    assuming(x ≫ 0) do
        @test solve(x^2 ~ 4, x) == [2]
        @test solve(x^2 ~ 4, x; assume = (x ≪ 1,)) == Sym[]       # the keyword adds
        @test string(solve(x^2 ≫ 1, x)) == "(1, ∞)"
    end
    # cases: a sign split when a is real, the fact and its negation when not
    cs = solve(x^2 ~ a, x; assume = (x ∈ ℝ, a ∈ ℝ))
    @test cs isa SymJ.Cases && length(cs) == 3
    @test cs[1] == (Any[a ≪ 0] => Sym[]) && cs[2] == (Any[a ~ 0] => [Sym(0)])
    @test cs(a => 4) == [-2, 2] && cs(a => 0) == [0] && cs(a => -1) == Sym[]
    cs = solve(x^2 ~ a, x; assume = (x ≫ 0,))
    @test length(cs) == 2 && cs[1] == (Any[a ≫ 0] => [sqrt(a)])
    @test cs(a => 4) == [2] && cs(a => -1) == Sym[] && cs[a => 4] == cs[1]
    cs = solve(a * x ~ 1, x; assume = (x ≫ 0,))
    @test cs(a => 2) == [1 // 2]
    @test cs(a => 0) == Sym[]
    @test cs(a => -1) == Sym[]
    # drop_empty: the cases with no solutions dropped, one left given as its solutions
    @test solve(x^2 ~ a, x; assume = (x ≫ 0,), drop_empty = true) == [sqrt(a)]
    @test solve(x^2 ~ 4, x; assume = (x ≫ 0,), drop_empty = true) == [2]
    @test solve(x^2 ~ -4, x; assume = (x ∈ ℝ,), drop_empty = true) == Sym[]
    cs = solve(x^2 ~ a, x; assume = (x ∈ ℝ, a ∈ ℝ), drop_empty = true)
    @test cs isa SymJ.Cases
    @test length(cs) == 2
    let (x0, v0, α, t) = symbols("x0 v0 α t")                   # the projectile's t(x)
        @test solve(x ~ x0 + v0 * cos(α) * t, t; assume = (x0 ∈ ℝ, v0 ∈ ℝ, α ∈ ℝ, t ≫ 0), drop_empty = true) ==
              [(x - x0) / (v0 * cos(α))]
    end
    cs = solve(x^2 ~ 9, x; assume = (x ≫ a,))
    @test [c.second for c in cs] == [[-3, 3], [3], [3], Sym[]]
    @test cs(a => -5) == [-3, 3] && cs(a => 0) == [3] && cs(a => 4) == Sym[]
    cs = solve(x^2 - 3x + 2 ~ 0, x; assume = (n ∈ ℤ, x ≪ n))
    @test cs(n => 1) == Sym[] && cs(n => 2) == [1] && cs(n => 5) == [1, 2]
    # rationality of constants: radicals, π and ℯ
    @test solve(x^2 ~ 2, x; assume = (x ∈ ℚ,)) == Sym[] == solve(x^2 ~ 2, x; assume = (x ∈ ℤ,))
    @test solve(4x^2 ~ 1, x; assume = (x ∈ ℚ,)) == [-1 // 2, 1 // 2]
    @test solve(4x^2 ~ 1, x; assume = (x ∈ ℤ,)) == Sym[]
    @test solve(x^2 - x - 1 ~ 0, x; assume = (x ∈ ℚ,)) == Sym[]
    @test solve(x^3 ~ 2, x; assume = (x ∈ ℚ,)) == Sym[]          # complex roots too
    s2, s3, s6 = sqrt(Sym(2)), sqrt(Sym(3)), sqrt(Sym(6))
    @test SymJ._const_rational(s2 * s3 * s6) === true && SymJ._const_rational(s2 * s3) === false
    @test SymJ._const_rational(2 + 3 * PI) === false && SymJ._const_rational(PI * E) === nothing
    @test SymJ._const_rational(s2 + s3) === nothing                # several radicals: unknown
    @test !(s2 ∈ ℚ) && !(PI ∈ ℤ) && (s2 * s3 * s6 ∈ ℤ)
    cs = solve(2x ~ n, x; assume = (x ∈ ℤ,))
    @test length(cs) == 2 && cs[1].second == [n / 2]
    @test occursin("a > 0  =>  [sqrt(a)]", sprint(show, MIME"text/plain"(), solve(x^2 ~ a, x; assume = (x ≫ 0,))))
    # inequalities: within the unknown's set
    @test string(solve(x^2 ≫ 1, x; assume = (x ≪ 0,))) == "(-∞, -1)"
    # a fact about a product: no recursion (a's sign through the bounds of a*b)
    P = @assume(a * b > 0, a in Reals, b in Reals)
    @test SymJ._signs(a, P) == 0x07 && SymJ._signs(a * b, P) == 0x04
    @test solve(x^2 ~ a * b, x; assume = P) == [-sqrt(a * b), sqrt(a * b)]
    # a product's zero split one factor at a time; linear conditions solved for a symbol
    @test SymJ._solved_for(Relation(:lt, a - b, Sym(0))) in (a ≪ b, b ≫ a)
    @test SymJ._zero_factor(a * b, Assumptions()) in (a, b)
    cs = solve(x^2 ~ 9, x; assume = @assume(x > a - b))
    @test length(cs) == 4 && cs(a => 0, b => 5) == [-3, 3] && cs(a => 5, b => 0) == Sym[]
    @test occursin("a - 3 < b < a + 3", sprint(show, MIME"text/plain"(), cs))
    @test_throws ArgumentError solve(x^2 ~ 4, x; assume = (x ≫ 1, x ≪ 0))   # contradictory
    # a predicate on a symbol is false (not known), and the message says so
    @test !isodd(x)
    @test occursin("isodd(n)", (try Assumptions(isodd(x)) catch e; sprint(showerror, e) end))
    # membership of ±x + k (ℤ) or c x + d (ℚ, ℝ, ℂ) is membership of x (SymPy #29082)
    @test iszero(refine(sin(PI * x); assume = (-x ∈ ℤ,)))
    @test iszero(refine(sin(PI * x); assume = (x + 3 ∈ ℤ,)))
    @test iszero(refine(imag(x); assume = (2x - 1 ∈ ℝ,)))
    @test_throws ArgumentError Assumptions(2x ∈ ℤ)          # says less than x ∈ ℤ
    @test_throws ArgumentError Assumptions(x^2 ∈ ℝ)
    # b^q, b < 0 and q real, not an integer, is not real (SymPy #29335)
    @test SymJ._isreal(x^y, Assumptions(x ≪ 0, y ∈ ℚ, y ∉ ℤ)) === false
    @test SymJ._isreal(x^y, Assumptions(x ≪ 0, y ∈ ℚ)) === nothing
end
