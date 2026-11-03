# Polynomials and rational functions (included by runtests.jl)

@testset "Poly" begin
    # an expression as a generator: its symbols must not be in a coefficient (SymPy #29080)
    @test_throws ArgumentError Poly(x + 1 / x, 1 / x)
    # division by a constant; by a polynomial in the generators: divrem
    @test Poly(x^2 + 2x + 1, x) / 3 == Poly(x^2 / 3 + 2x / 3 + 1 // 3, x)
    @test Poly(x^2 + 2x + 1, x) // 3 == Poly(x^2 / 3 + 2x / 3 + 1 // 3, x)
    @test_throws ArgumentError Poly(x^2 + 1, x) / x
    # gcd with symbolic (algebraic) coefficients, and with a Julia number
    @test gcd(Poly(sqrt(Sym(2)) * x^2 + sqrt(Sym(2)) * x, x), Poly(x^2 - 1, x)) == Poly(x + 1, x)
    @test gcd(Poly(x^2 - 2, x), Poly(x - sqrt(Sym(2)), x)) == Poly(x - sqrt(Sym(2)), x)
    @test lcm(Poly(sqrt(Sym(2)) * x, x), Poly(x + 1, x)) == Poly(x^2 + x, x)
    @test gcd(2x, 0) == 2x
    @test lcm(2x, 3) == 6x
    @test gcd(6, 4x) == 2
    @test_throws ArgumentError Poly(x^2 + sin(x), sin(x))
    @test degree(Poly(sin(x)^2 + 2sin(x), sin(x)), sin(x)) == 2
    @syms a b
    p = Poly((x + y + 1)^3, x, y)
    @test p isa Poly{SymNumber}
    @test nterms(p) == 10
    @test Sym(p) == SymJ._expand_generic((x + y + 1)^3)
    @test gens(p) == [x, y]
    @test degree(p) == 3 && degree(p, y) == 3
    @test exponent_vectors(Poly(x^2 * y + 3, x, y)) == [[2, 1], [0, 0]]
    @test coefficients(Poly(2x^2 - 1 // 2, x)) == [SymNumber(2), SymNumber(-1 // 2)]

    # products of integer polynomials: dense slices and the hash table, with Int128
    # or 192-bit accumulation, agree with expression arithmetic
    let z = Sym(:z), f = 1 + x + y + z, c = 2^62 - 1
        P = Poly(f, x, y, z)
        @test Sym(P^6 * (P^6 + 1)) == SymJ._expand_generic(f^6 * (f^6 + 1))     # dense
        @test Sym(Poly(x^9 - y^7 + 3z, x, y, z) * Poly(x^8 + 5y^3 * z - 1, x, y, z)) ==
              SymJ._expand_generic((x^9 - y^7 + 3z) * (x^8 + 5y^3 * z - 1))   # sparse
        # Int64 coefficients whose sums pass 2^127: 192-bit accumulation, both signs
        h = sum(c * x^i * y^j for i in 0:4, j in 0:4)
        H = Poly(h, x, y)
        @test Sym(H * H) == SymJ._expand_generic(h * h)                       # dense
        @test Sym(H * -H) == SymJ._expand_generic(-h * h)
        s = sum(c * x^(7i) * y^(3i) for i in 0:9) - c * z^20
        S = Poly(s, x, y, z)
        @test Sym(S * S) == SymJ._expand_generic(s * s)                       # sparse
    end

    # arithmetic agrees with expression arithmetic
    q = Poly(x - 2y + 3, x, y)
    for r in (p + q, p - q, p * q, q^4, -q, 2 * p, p + 1, p * x)
        @test r isa Poly
    end
    @test Sym(p * q) == expand(Sym(p) * Sym(q))
    @test Sym(q^4) == SymJ._expand_generic((x - 2y + 3)^4)
    @test p - p == zero(p) && iszero(p - p)
    @test isone(q^0)
    @test p * q == q * p
    @test (p + q) * q == p * q + q * q

    # coefficient paths: Int128 fast path, BigInt fallback, rationals, floats
    for ex in ((x + 2^40 * y + 1)^5, (x + typemax(Int) * y)^3, (x / 3 + y / 2 - 1)^6, (0.5x + y)^3)
        @test Sym(Poly(ex, x, y)) == SymJ._expand_generic(ex)
    end
    @test coefficients(Poly((x + 2^40)^4, x))[end] == SymNumber(big(2)^160)

    # powers: multinomial theorem (sparse bases) and repeated multiplication (dense)
    rep(p, k) = foldl(*, fill(p, k))
    for (b, k) in ((Poly(x + y, x, y), 100), (Poly(x - 3y, x, y), 7), (Poly(x / 2 + y - 1, x, y), 9),
                   (Poly(big(2)^70 * x + 1, x), 5), (Poly(sum(x^i for i in 0:9), x), 20),
                   (Poly(x * y + z + 1, x, y, z), 12))
        @test b^k == rep(b, k)
    end
    @test coefficients(Poly(x + y, x, y)^100)[51] == SymNumber(binomial(big(100), 50))
    @test Poly(Sym(:a) * x + 1, x)^3 == rep(Poly(Sym(:a) * x + 1, x), 3)   # Sym coefficients

    # differing generators are unified
    @test Poly(x + 1, x) * Poly(y - 1, y) == Poly(x * y - x + y - 1, x, y)

    # symbolic coefficients
    s = Poly(a * x^2 + b * x + 3, x)
    @test s isa Poly{Sym}
    @test coefficients(s) == [a, b, Sym(3)]
    @test expand(Sym(s * s)) == expand((a * x^2 + b * x + 3)^2)
    @test_throws ArgumentError Poly{SymNumber}(a * x, x)
    @test_throws ArgumentError Poly(sin(x), x)
    @test_throws ArgumentError Poly(sqrt(x), x)

    @test diff(p, x) == Poly(diff(Sym(p), x), x, y)
    @test diff(p, z) == zero(p)
    @test p(x => 1, y => 1) == 27
    @test string(Poly(x^2 - 2x * y + 1, x, y)) == "Poly(x^2 - 2*x*y + 1, x, y)"

    # Exponents beyond the uniform packing (16 generators: 4 bits each, max 14). A
    # field too small gets a new layout in the same word (fields sized per generator);
    # past the word, construction and gcd, factor, cancel, ... widen to UInt128 or
    # BigInt, while direct arithmetic (type-stable) raises.
    g = [Sym(Symbol(:v, i)) for i in 1:16]
    h = Poly(g[1] + 1, g...)
    @test h.layout == SymJ._uniform(16, 4)
    @test nterms(h^14) == 15 && (h^14).layout == h.layout
    h15 = h^15                                       # g[1] needs 5 bits: new layout
    @test SymJ._monotype(h15) === UInt64 && h15.layout != h.layout
    @test nterms(h15) == 16 && Sym(h15) == expand((g[1] + 1)^15) && h15 == h^14 * h
    @test hash(h15) == hash(Poly((g[1] + 1)^15, g...))
    hs = Poly(sum(g) + 1, g...)
    @test_throws SymJ._PackOverflow hs^15            # 16 fields of 5 bits: not in a UInt64
    @test occursin("widen", sprint(showerror, SymJ._PackOverflow(15)))
    @test SymJ._monotype(Poly(prod(g)^15, g...)) === UInt128
    @test SymJ._monotype(Poly(prod(g)^(2^40), g...)) === BigInt
    x40 = Poly(g[1]^(2^40) + g[16], g...)           # one high degree among 16: one UInt64
    @test SymJ._monotype(x40) === UInt64 && degree(x40) == 2^40 && Sym(x40) == g[1]^(2^40) + g[16]
    w = convert(SymJ.Poly{SymNumber,UInt128}, h^14)
    @test w == h^14 && hash(w) == hash(h^14) && w * convert(SymJ.Poly{SymNumber,UInt128}, h) == h15
    @test Poly(g[1]^100 * g[2] * g[3], g...) == Poly(g[1]^100 * g[2], g...) * Poly(g[3], g...)
    @test diff(Poly(g[1]^30, g...), g[1]) == Poly(30 * g[1]^29, g...)
    @test gcd(Poly(g[1]^20 - 1, g...), Poly(g[1]^10 - 1, g...)) == Poly(g[1]^10 - 1, g...)
    # division whose remainder outgrows the inputs' degrees (x^k by x + y^5)
    qd, rd = divrem(Poly(g[1]^4, g...), Poly(g[1] + g[2]^5, g...))
    @test rd == Poly(g[2]^20, g...) && qd * Poly(g[1] + g[2]^5, g...) + rd == Poly(g[1]^4, g...)
    big14 = prod(g[3:16])^(2^40)                     # BigInt words
    num, den = Poly(g[1]^4 + big14, g...), Poly(g[1] + g[2]^5, g...)
    qd, rd = divrem(num, den)
    @test SymJ._monotype(num) === BigInt
    @test qd * den + rd == num && rd == Poly(g[2]^20 + big14, g...)
    # more than 64 generators
    g70 = [Sym(Symbol(:u, i)) for i in 1:70]
    P70 = Poly(sum(g70)^2, g70...)
    @test nterms(P70) == 70 * 71 ÷ 2 && Sym(P70) == expand(sum(g70)^2)
    @test gcd(P70, Poly(sum(g70) * g70[1], g70...)) == Poly(sum(g70), g70...)
    @test expand((g[1] + 1)^20) == SymJ._expand_generic((g[1] + 1)^20)  # unaffected (1 generator)
end

@testset "Poly division and gcd" begin
    P(e) = Poly(e, x, y, z)
    # small deterministic pseudo-random polynomials (LCG; no extra test deps)
    state = Ref(UInt64(12345))
    rnd(k) = (state[] = state[] * 6364136223846793005 + 1442695040888963407; Int((state[] >> 33) % k))
    randpoly(nt, d) = P(sum((rnd(7) - 3) * x^rnd(d) * y^rnd(d) * z^rnd(2) for _ in 1:nt) + rnd(5) + 1)

    # long division
    q, r = divrem(Poly(x^3 - 1, x), Poly(x - 1, x))
    @test q == Poly(x^2 + x + 1, x) && iszero(r)
    q, r = divrem(Poly(x^2 + 1, x), Poly(2x, x))
    @test q == Poly(x / 2, x) && r == Poly(Sym(1), x)
    @test_throws DivideError divrem(P(x), zero(P(x)))
    @test divexact(P((x + y)^3), P(x + y)) == P((x + y)^2)
    @test_throws ArgumentError divexact(P(x + 1), P(y))
    @test divides(P(x * y + x), P(y + 1)) == (true, P(x))

    for _ in 1:30
        a, b = randpoly(6, 4), randpoly(3, 3)
        q, r = divrem(a, b)
        @test q * b + r == a
        @test !any(m -> SymJ._mono_divides(b.monos[1], m, b.layout, 3), r.monos)
    end

    # gcd
    @test gcd(P((x + 1) * (x - 2)^2 * (x + y)), P((x + 1) * (x + y) * (y - 3))) == P((x + 1) * (x + y))
    @test gcd(P(2x + 2), P(4x + 4)) == P(2x + 2)              # integer content kept
    @test gcd(P(x / 2 + 1 // 2), P(x^2 - 1)) == P(x + 1)      # rational input: monic
    @test gcd(P(-x - 1), P(x + 1)) == P(x + 1)                # positive leading coefficient
    @test isone(gcd(P(x^2 + y^2), P(x + y)))
    @test gcd(zero(P(x)), P(3x + 6)) == P(3x + 6)
    @test gcd(P(Sym(6)), P(Sym(4))) == P(Sym(2))
    @test gcd(Poly(x^2 - 1, x), Poly(y * x + y, y, x)) == Poly(x + 1, x, y)  # gens unified
    @test_throws ArgumentError gcd(P(0.5x), P(x))
    @test lcm(P(x^2 - 1), P(x^2 + 2x + 1)) == P((x - 1) * (x + 1)^2)

    for _ in 1:20
        a, b, g = randpoly(3, 3), randpoly(3, 3), randpoly(2, 2)
        A, B = a * g, b * g
        G = gcd(A, B)
        @test G == gcd(B, A)
        @test divides(A, G)[1] && divides(B, G)[1]
        @test divides(G, g)[1]
        @test isone(gcd(divexact(A, G), divexact(B, G)))
    end

    # modular gcd agrees with the PRS reference implementation
    for _ in 1:40
        a, b, g = randpoly(4, 4), randpoly(4, 3), randpoly(3, 3)
        @test gcd(a * g, b * g) == SymJ._gcd(a * g, b * g)
        @test gcd(a, b) == SymJ._gcd(a, b)
    end
    # large coefficients: several primes and Chinese remaindering
    g = P(big(3)^40 * x^2 * y - big(7)^25 * z + 2^40)
    A, B = g * P(x^3 + 5y - big(11)^20), g * P(y^2 * z - 3x + 1)
    @test gcd(A, B) == g
    @test gcd(P(big(2)^100 * (x + 1)), P(big(6)^80 * (x + 1) * (y + 2))) == P(big(2)^80 * (x + 1))
    # the first prime used, 2^31 - 1, is unlucky here: the cofactors agree mod p
    @test gcd(P((x + y) * x), P((x + y) * (x + (2^31 - 1)))) == P(x + y)
    @test gcd(P((x + y) * (z + 1) * x), P((x + y) * (z + 1) * (x + (2^31 - 1) * y))) == P((x + y) * (z + 1))
    # Zippel (opt-in): the same results, sparse or not, falling back when its
    # interpolation cannot apply (a content in the other variables, a dense gcd)
    w = Sym(:w)
    P4(e) = Poly(e, x, y, z, w)
    for (g, c1, c2) in ((x^10 * y^3 * z + 2y^7 * z^5 * w - x^2 * w^8 + 3y^9 + x * y * z * w,
                         x^9 * z^4 + y^8 * w^2 + 5, y^10 * w^3 - x^6 * w^5 + 7z),
                        ((y + z) * x^7 + y^9 * w^4 + z^6 * w^5 + 1, x^4 * w + y^5 + 2, (w^3 - y) * x^5 + z^8),
                        ((y^2 + 1) * (x^3 + z), x + z * w, x^2 - w),
                        ((x + y + z + w + 1)^4, (x - y + 2z * w + 1)^2 + y, (y * z - x + 3)^2 + w),
                        (big(3)^20 * x * y + big(5)^15 * z^4 * w - 7, x - y * z, y^3 + w))
        G = P4(expand(g))
        A, B = G * P4(c1), G * P4(c2)
        @test gcd(A, B; method = :zippel) == gcd(A, B) == G
    end
    saved = state[]                  # leave the later random tests unchanged
    for _ in 1:20
        a, b, g = randpoly(4, 4), randpoly(4, 3), randpoly(3, 3)
        @test gcd(a * g, b * g; method = :zippel) == gcd(a * g, b * g)
    end
    state[] = saved
    @test gcd(P(x // 2 + y // 3), P((3x + 2y) * z); method = :zippel) == P(x + 2y // 3)
    @test_throws ArgumentError gcd(P(x), P(y); method = :sparse)
    # :auto (the default): dense gcds stay with Brown, sparse or many-variable ones
    # take Zippel; the results agree
    Gd = P4(expand((x + y + z + w + 1)^6))
    Ad, Bd = Gd * P4(x - y + 2z + 1), Gd * P4(y * z - x + 3)
    @test !SymJ._use_zippel(Ad, Bd)
    vs = [Sym(Symbol(:v, i)) for i in 1:12]
    gv = expand(sum(vs)^2 + vs[1] * vs[2])
    Av, Bv = Poly(expand(gv * (vs[1] - vs[3])), vs...), Poly(expand(gv * (vs[2] + vs[4] - 1)), vs...)
    @test SymJ._use_zippel(Av, Bv)
    @test gcd(Av, Bv) == gcd(Av, Bv; method = :zippel) == Poly(gv, vs...)
    @test !SymJ._use_zippel(P(x^2 * y * z - 1), P(x * y * z^2 + 1))      # tiny: dense
    @test occursin(":auto", sprint(showerror, try gcd(P(x), P(y); method = :sparse) catch e; e end))
    # the keyword through the functions built on gcd (the same results)
    zg = x^10 * y^3 * z + 2y^7 * z^5 * w - x^2 * w^8 + 3y^9 + x * y * z * w
    za, zb = expand(zg * (x^9 * z^4 + y^8 * w^2 + 5)), expand(zg * (y^10 * w^3 - x^6 * w^5 + 7z))
    # (the sign follows the order of the generators, which is by hash and so differs
    # between Julia versions)
    @test gcd(za, zb; method = :zippel) == gcd(za, zb) &&
          (iszero(expand(gcd(za, zb) - zg)) || iszero(expand(gcd(za, zb) + zg)))
    @test lcm(za, zb; method = :zippel) == lcm(za, zb)
    @test lcm(P4(za), P4(zb); method = :zippel) == lcm(P4(za), P4(zb))
    @test cancel(za / zb; gcd_method = :zippel) == cancel(za / zb)
    zs = expand(zg^2 * (x * y - z))
    @test factor_squarefree(zs; gcd_method = :zippel) == factor_squarefree(zs)
    @test factor_squarefree(P4(zs); gcd_method = :zippel) == factor_squarefree(P4(zs))
    @test squarefree_part(zs; gcd_method = :zippel) == squarefree_part(zs)
    zf = expand((x * y + z^2 + 1)^2 * (x - y * w))
    @test factor(zf; gcd_method = :zippel) == factor(zf)
    @test factor(P4(zf); gcd_method = :zippel) == factor(P4(zf))
    @test apart(1 / ((x - y) * (x + y * z)), x; gcd_method = :zippel) == apart(1 / ((x - y) * (x + y * z)), x)
    @test apart(1 / (x^2 - 1); gcd_method = :zippel) == apart(1 / (x^2 - 1))
    @test cancel(Sym(3); gcd_method = :zippel) == 3 && factor(6; gcd_method = :zippel) == 6
    @test_throws ArgumentError cancel(za / zb; gcd_method = :sparse)
    @test_throws ArgumentError factor(zf; gcd_method = :sparse)

    # univariate, high degree
    u, v = Poly((x + 1)^30 * (x - 2)^20, x), Poly((x + 1)^25 * (x + 3)^20, x)
    @test gcd(u, v) == Poly((x + 1)^25, x)

    # squarefree factorization
    prodf(u, fs) = foldl((r, fm) -> r * fm.first^fm.second, fs; init = P(Sym(u)))
    u, fs = factor_squarefree(P((x - 1)^3 * (x + 2)^2 * (x + 3)))
    @test u == 1 && fs == [P(x + 3) => 1, P(x + 2) => 2, P(x - 1) => 3]
    @test factor_squarefree(P(6 * (x + 1)^2)) == (SymNumber(6), [P(x + 1) => 2])
    @test factor_squarefree(P(-(x - 1)^2)) == (SymNumber(-1), [P(x - 1) => 2])
    @test factor_squarefree(P((x / 2 + 1 // 2)^2)) == (SymNumber(1 // 4), [P(x + 1) => 2])
    @test factor_squarefree(P(x * (y + 1)^2)) == (SymNumber(1), [P(x) => 1, P(y + 1) => 2])
    @test factor_squarefree(P(Sym(5))) == (SymNumber(5), Pair{Poly{SymNumber},Int}[])
    @test_throws ArgumentError factor_squarefree(zero(P(x)))
    @test squarefree_part(P((x + y)^2 * (x - y) * (y + 1)^3)) == P((x + y) * (x - y) * (y + 1))
    for _ in 1:15
        fs0 = [randpoly(2, 2) for _ in 1:3]
        p = P(Sym(rnd(5) + 1)) * fs0[1] * fs0[2]^2 * fs0[3]^3
        u, fs = factor_squarefree(p)
        @test prodf(u, fs) == p
        @test issorted(last.(fs)) && allunique(last.(fs))
        for (f, _) in fs                     # squarefree: gcd(f, ∂f/∂x, ∂f/∂y, ∂f/∂z) = 1
            @test isone(foldl(gcd, (diff(f, v) for v in (x, y, z)); init = f))
            @test sign(f.coefs[1]) > 0
        end
        for i in eachindex(fs), j in (i + 1):length(fs)   # pairwise coprime
            @test isone(gcd(fs[i].first, fs[j].first))
        end
    end
    # expressions
    @test factor_squarefree(x^3 - x^2 - x + 1) == (x + 1) * (x - 1)^2
    @test factor_squarefree(2x^2 + 4x + 2) == 2 * (x + 1)^2
    @test factor_squarefree((x^2 - 1) / (x^2 + 2x + 1)^2) == (x - 1) / (x + 1)^3
    @test factor_squarefree(sin(x)^2 + 2sin(x) + 1) == (sin(x) + 1)^2
    @test factor_squarefree(x / (x + 0.5)) == x / (x + 0.5)   # inexact: unchanged
    @test squarefree_part(x^3 - x^2 - x + 1) == (x + 1) * (x - 1)

    # on expressions
    # beyond the 64-bit packing: the slower path, not a refusal or a silent no-op
    @syms a1 a2 a3 a4 a5 a6 a7 a8
    @test gcd(x^130 * a1 - a2 + a3 + a4 + a5 + a6 + a7 + a8, x^65 + a1) == 1
    wg = x^70 + a1 * a2 + a3^2
    @test gcd(expand(wg * (x - a4 + a5 * a6)), expand(wg * (x + a7 - a8))) == wg
    @test cancel((x^40 - a1^40) / (x - a1) + a2 + a3 + a4 + a5 + a6 + a7 + a8 + y + z) ==
          expand(sum(x^k * a1^(39 - k) for k in 0:39)) + a2 + a3 + a4 + a5 + a6 + a7 + a8 + y + z
    @test factor(expand((x^40 + a1) * (x - a2) * (a3 + a4 + a5 + a6 + a7 + a8 + y + z))) ==
          (x - a2) * (x^40 + a1) * (a3 + a4 + a5 + a6 + a7 + a8 + y + z)

    @test gcd(x^2 - 1, x^2 + 2x + 1) == x + 1
    @test gcd(x^2 * y - y, x * y + y) == x * y + y
    @test lcm(x^2 - 1, x + 1) == x^2 - 1
    @test_throws ArgumentError gcd(sin(x), x)
end

@testset "numer_denom and cancel" begin
    # square roots of rationals: cancelled through their conjugates, kept when the
    # denominator would not get simpler
    s2, s3 = sqrt(Sym(2)), sqrt(Sym(3))
    @test cancel((x^2 + 2s2 * x + 2) / (x + s2)) == x + s2
    @test cancel((x^2 - 2) / (x - s2)) == x + s2
    @test cancel((x^2 - 5 - 2s2 * s3) / (x - s2 - s3)) == x + s2 + s3
    @test cancel((x^2 - 3) / ((x + s3) * (x + 2))) == (x - s3) / (x + 2)
    @test iszero(cancel(cancel((x + s2 * y) / (x^2 - 2y^2)) - 1 / (x - s2 * y)))
    @test cancel(1 / (x - s2)) == 1 / (x - s2)                     # not (x + √2)/(x² - 2)
    @test cancel((x + s2) / (x^2 + 1)) == (x + s2) / (x^2 + 1)
    @test cancel((sqrt(x) + 1) / (x - 1)) == (sqrt(x) + 1) / (x - 1)  # sqrt(x) is no number
    @test numer_denom(x / y + 1) == (x + y, y)
    @test numer_denom(3x / 4) == (3x, Sym(4))
    @test numer_denom(1 / sqrt(x)) == (Sym(1), sqrt(x))
    @test numer_denom(x) == (x, Sym(1))
    @test numer_denom(Sym(2) / 3) == (Sym(2), Sym(3))

    @test cancel((x^2 - 1) / (x + 1)) == x - 1
    @test cancel((x^2 - 1) / (x^2 + 2x + 1)) == (x - 1) / (x + 1)
    @test cancel(1 / x + 1 / x^2) == (x + 1) / x^2
    @test cancel((2x + 2) / (4x)) == (x + 1) / (2x)
    @test cancel((x * y + x) / (y^2 - 1)) == x / (y - 1)
    @test cancel((x - 1) / (1 - x)) == -1
    @test cancel(1 / (x + 1) - 1 / (x - 1)) == -2 / (x^2 - 1)
    @test cancel((x / 2 + 1) / (x + 2)) == Sym(1 // 2)
    @test cancel(x * (x + 1)) == x^2 + x                   # polynomials come back expanded
    @test cancel((sin(x)^2 - 1) / (sin(x) + 1)) == sin(x) - 1
    @test cancel((x - 1) / (sqrt(x) - 1)) == (x - 1) / (sqrt(x) - 1)  # radicals are independent
    @test cancel(x / (x + 0.5)) == x / (x + 0.5)           # floats: unchanged
    @test cancel(Sym(3) / 6) == Sym(1 // 2)
    @test cancel(2) == Sym(2)

    # cancel never changes the value
    for e in ((x^3 - y^3) / (x - y) + 1 / (x + y), (x^2 * y - y) / (x * y + y) - x / (y + 1),
              (1 / x - 1 / y) / (x - y), ((x + 1)^3 - 1) / x)
        c = cancel(e)
        for (a, b) in ((2, 3), (-1 // 2, 5), (7, -4))
            @test e(x => a, y => b) == c(x => a, y => b)
        end
    end
    @test cancel((x^3 - y^3) / (x - y)) == x^2 + x * y + y^2
    @test cancel((1 / x - 1 / y) / (x - y)) == -1 / (x * y)
end

@testset "factor and apart" begin
    P1(e) = Poly(e, x)
    P3(e) = Poly(e, x, y, z)
    prodf(p, (u, fs)) = foldl((r, fm) -> r * fm.first^fm.second, fs; init = Poly(Sym(u), gens(p)...))
    @test factor(P1(x^4 - 1)) == (SymNumber(1), [P1(x + 1) => 1, P1(x - 1) => 1, P1(x^2 + 1) => 1])
    @test factor(P1(6x^2 + 5x + 1)) == (SymNumber(1), [P1(2x + 1) => 1, P1(3x + 1) => 1])
    @test factor(P1(2x^3 - 6x)) == (SymNumber(2), [P1(x) => 1, P1(x^2 - 3) => 1])
    @test factor(P1(-(x^2 + 1)^2 * (x - 3) / 4)) == (SymNumber(-1 // 4), [P1(x - 3) => 1, P1(x^2 + 1) => 2])
    @test length(factor(P1(x^5 - x - 1))[2]) == 1                 # irreducible
    @test length(factor(P1(x^4 + 1))[2]) == 1                     # irreducible over Q, splits mod every p
    @test Set(factor(P1(x^4 + 4))[2]) == Set([P1(x^2 - 2x + 2) => 1, P1(x^2 + 2x + 2) => 1])
    @test [degree(g) for (g, _) in factor(P1(x^60 - 1))[2]] == [1, 1, 2, 2, 2, 4, 4, 4, 8, 8, 8, 16]
    @test length(factor(P1(prod(x - i for i in -15:15)))[2]) == 31
    big2 = P1((big(10)^20 * x^3 + 7x + big(3)^40) * (x^4 - big(2)^60 * x + 1))
    @test prodf(big2, factor(big2)) == big2 && length(factor(big2)[2]) == 2
    # Swinnerton-Dyer: irreducible of degree 16, many factors mod every prime
    s2, s3, s5, s7 = sqrt(Sym(2)), sqrt(Sym(3)), sqrt(Sym(5)), sqrt(Sym(7))
    sd = expand(prod(x + a + b + c + d for a in (s2, -s2), b in (s3, -s3), c in (s5, -s5), d in (s7, -s7)))
    @test length(factor(P1(sd))[2]) == 1
    # multivariate
    @test Set(factor(P3(x^2 * y^2 - 1))[2]) == Set([P3(x * y - 1) => 1, P3(x * y + 1) => 1])
    @test factor(P3(x^3 + y^3))[2] == [P3(x + y) => 1, P3(x^2 - x * y + y^2) => 1]
    @test factor(P3(x * y + x + y + 1))[2] == [P3(x + 1) => 1, P3(y + 1) => 1]
    @test length(factor(P3(x^6 - y^6))[2]) == 4
    @test factor(P3((x + y + z + 1)^2 * (x - y) * (x * y * z - 1)))[2] ==
          [P3(x - y) => 1, P3(x * y * z - 1) => 1, P3(x + y + z + 1) => 2]
    @test_throws ArgumentError factor(zero(P1(x)))

    # EEZ: against Kronecker substitution (the reference) on random products
    w = Sym(:w)
    P4(e) = Poly(e, x, y, z, w)
    state4 = Ref(UInt64(77))
    rnd4(k) = (state4[] = state4[] * 6364136223846793005 + 1442695040888963407; Int((state4[] >> 33) % k))
    rp4() = P4(sum((rnd4(5) - 2) * x^rnd4(3) * y^rnd4(2) * z^rnd4(2) for _ in 1:3) + rnd4(3) + 1)
    for _ in 1:12
        p = rp4() * rp4() * rp4()
        _, sq = factor_squarefree(p)
        for (q, _) in sq
            act = [i for i in 1:4 if SymJ._deg_in_all(q, i) > 0]
            length(act) > 1 || continue
            @test Set(SymJ._factor_eez(q, act)) == Set(SymJ._factor_kronecker(q, act))
        end
    end
    # a point that splits an irreducible factor: recombination recovers the factors
    P2(e) = Poly(e, x, y)
    for (e, us, want) in ((x^2 - y, [P2(x - 2), P2(x + 2)], [x^2 - y]),
                          ((x^2 - y) * (x + y), [P2(x - 2), P2(x + 2), P2(x + 4)], [x + y, x^2 - y]))
        qs = SymJ._shift_var(P2(expand(e)), 2, 4)
        gs = [SymJ._shift_var(g, 2, -4) for g in SymJ._eez_split(qs, us, 1, [2])]
        @test Set(Sym.(gs)) == Set(want)
    end
    # van Hoeij recombination: the same factors as subset search, and the
    # Swinnerton-Dyer polynomials (irreducible, many factors mod every prime)
    densez(S) = (v = zeros(BigInt, degree(S) + 1);
                 for (e, c) in zip(exponent_vectors(S), coefficients(S)); v[e[1]+1] = big(SymJ._numerator(c)); end; v)
    function vh_vs_subsets(e)
        f = densez(P1(expand(e)))
        p, mf = SymJ._choose_prime(f)
        vh = SymJ._vanhoeij(f, mf, p)
        B = big(2)^(length(f) - 1) * (isqrt(sum(c -> c^2, f)) + 1)
        k, pk = 1, big(p)
        while pk <= 2 * abs(f[end]) * B; k += 1; pk *= p; end
        zs = SymJ._recombine(f, SymJ._hensel_lift(SymJ._zmod(f, pk), mf, p, k), pk)
        length(mf), vh, zs
    end
    for e in (prod(x - i for i in 1:12), prod(x^2 + i for i in 1:6), x^30 - 1,
              (x^4 - 10x^2 + 1) * (x^2 - 2) * (x - 1) * (x^3 - 2))
        nmod, vh, zs = vh_vs_subsets(e)
        @test vh !== nothing && Set(vh) == Set(zs)
    end
    function swinnerton_dyer(k)
        S = x
        for q in (2, 3, 5, 7, 11)[1:k]
            r = sqrt(Sym(q))
            S = expand(expand(S(x => x + r)) * expand(S(x => x - r)))
        end
        P1(S)
    end
    for k in (4, 5)
        f = densez(swinnerton_dyer(k))
        p, mf = SymJ._choose_prime(f)
        k == 5 && @test length(mf) >= SymJ._VANHOEIJ_MIN_FACTORS  # the lattice path
        @test SymJ._vanhoeij(f, mf, p) == [f]                     # irreducible
    end
    @test length(factor(swinnerton_dyer(5))[2]) == 1

    # Wang's leading-coefficient distribution, at the point y = 0 where its test holds
    @test SymJ._wang_nondivisors(BigInt[2, 5], big(1), big(1)) == BigInt[2, 5]
    @test SymJ._wang_nondivisors(BigInt[3, 6], big(2), big(1)) === nothing   # 6 shares 2 and 3
    qw = P2(expand(((y + 2) * x + 1) * ((y + 5) * x^2 + 3)))
    gw = SymJ._eez_wang(qw, [P2(2x + 1), P2(5x^2 + 3)], 1, [2])
    @test gw !== nothing && Set(Sym.(gw)) == Set(expand.([(y + 2) * x + 1, (y + 5) * x^2 + 3]))
    # beyond Kronecker substitution's reach (image degree ~ 10*9*6*6)
    big4 = P4(expand((x^5 * y^3 + z^4 * w^2 - 7) * (x^4 - y^5 * z + w^3 + 2) * (x * y * z * w + 1)))
    @test_throws ArgumentError SymJ._factor_kronecker(big4, [1, 2, 3, 4])
    u4, f4 = factor(big4)
    @test prodf(big4, (u4, f4)) == big4 && length(f4) == 3
    # content in the main variable is factored separately
    @test Set(factor(P4((y^2 - 1) * (x * y + z)))[2]) == Set([P4(y - 1) => 1, P4(y + 1) => 1, P4(x * y + z) => 1])

    # random products: exact product, irreducible factors, right multiplicities
    state = Ref(UInt64(2024))
    rnd(k) = (state[] = state[] * 6364136223846793005 + 1442695040888963407; Int((state[] >> 33) % k))
    rp() = P3(sum((rnd(5) - 2) * x^rnd(3) * y^rnd(2) for _ in 1:3) + rnd(4) + 1)
    for _ in 1:12
        p = rp() * rp()^2 * rp()
        iszero(p) && continue
        u, fs = factor(p)
        @test prodf(p, (u, fs)) == p
        @test all(fm -> length(factor(fm.first)[2]) == 1 && sign(fm.first.coefs[1]) > 0, fs)
    end
    for _ in 1:12
        p = P1(sum((rnd(9) - 4) * x^i for i in 0:rnd(6)) + 1) * P1(sum((rnd(9) - 4) * x^i for i in 0:rnd(6)) + 1)
        iszero(p) && continue
        u, fs = factor(p)
        @test prodf(p, (u, fs)) == p
        @test all(fm -> length(factor(fm.first)[2]) == 1, fs)
    end

    # expressions
    @test factor(x^4 - 1) == (x + 1) * (x - 1) * (x^2 + 1)
    @test factor(x^2 * y^2 - 1) == (x * y + 1) * (x * y - 1)
    @test factor((x^2 - 1) / (x^2 + 2x + 1)) == (x - 1) / (x + 1)
    @test factor(sin(x)^2 - 1) == (sin(x) - 1) * (sin(x) + 1)
    @test factor(x / (x + 0.5)) == x / (x + 0.5)
    @test factor(12) == 12

    # partial fractions
    @test apart(1 / (x^2 - 1)) == (1 // 2) / (x - 1) - (1 // 2) / (x + 1)
    @test apart((x^3 + 1) / (x^2 - 1)) == x + 1 / (x - 1)
    @test apart(1 / (x^3 + x)) == 1 / x - x / (x^2 + 1)
    @test apart((3x + 5) / (x + 1)^2) == 3 / (x + 1) + 2 / (x + 1)^2
    @test apart((2x^2 + 3) / (x^3 - 2)) == (2x^2 + 3) / (x^3 - 2)   # irreducible denominator
    @test apart(x^2 + 1) == x^2 + 1
    for e in (1 / ((x - 1)^2 * (x + 1)), (x^4 + 1) / (x^2 * (x^2 + 1)^2), 1 / (x^4 - 1),
              (x / 2 + 1) / (3x^2 - 3), (x^5 + 2) / ((x^2 + x + 1)^2 * (x - 2)^3 * x))
        @test cancel(apart(e) - e) == 0
    end
    # symbolic coefficients: over ℚ(params), with variables or parameters
    @parameters a b c
    @test apart(1 / (x^2 - a^2), x) == 1 / (2a * (x - a)) - 1 / (2a * (x + a))
    @test apart(1 / (x^2 + a), x) == 1 / (x^2 + a)                 # irreducible over ℚ(a)
    for e in (1 / ((x - a) * (x - b)), (x + 1) / ((x - a)^2 * (x + b)), x^3 / ((x - a) * (x - 1)),
              1 / (a * x * (x + b)), (x^2 + c) / ((x - a)^2 * (x^2 + b)), 1 / ((x - y) * (x + y)),
              (a * x + 1) / ((x - a) * (x^2 - b)), 1 / ((x - a) * (x - 2a) * (x - 3a)))
        r = apart(e, x)
        @test cancel(r - e) == 0
    end
    @test length(SymJ.arguments(apart(1 / ((x - a) * (x - b) * (x - c)), x))) == 3
    @test_throws ArgumentError apart(x / y)                         # two symbols: give x
end

@testset "RatFun" begin
    R = SymJ.RatFun
    P(e) = Poly(e, x, y)
    r1 = SymJ._ratfun(P(x^2 - y^2), P(2x - 2y))            # reduces to (x + y)/2
    @test SymJ._ispoly(r1) && r1.num == P(x / 2 + y / 2)
    r2 = SymJ._ratfun(P(x), P(-2y))                         # integer num/den, positive den
    @test r2.den == P(2y) && r2.num == P(-x)
    @test r1 + r2 == R(P(x / 2 + y / 2)) + SymJ._ratfun(P(-x), P(2y))
    @test iszero(r1 - r1) && isone(r1 / r1) && r1 * inv(r1) == one(R)
    @test (r2^2) == r2 * r2 && r2^-1 == inv(r2)
    @test R(SymNumber(3)) + R(P(x)) == R(P(x + 3))           # constants mix with polynomials
    @test Sym(r2) == -x / (2y) && zero(R) == R(SymNumber(0))
    @test_throws DivideError inv(zero(R))
end

@testset "Gröbner bases" begin
    @syms gx gy gz gr ga
    same(A, B) = length(A) == length(B) && all(iszero(expand(a - b)) for (a, b) in zip(A, B))
    @test same(groebner([gx^2 + gy^2 - 1, gx - gy], gx, gy), [gx - gy, gy^2 - 1 // 2])
    @test same(groebner([gx^2 + gy^2 - 1, gx - gy], gx, gy; order = :grevlex), [gy^2 - 1 // 2, gx - gy])
    @test same(groebner([gx * gy - 1, gx^2 + gy^2 - 4], gx, gy), [gx + gy^3 - 4gy, gy^4 - 4gy^2 + 1])
    @test same(groebner([gx * gy - 1, gx^2 + gy^2 - 4], gx, gy; order = :grevlex),
               [gy^3 + gx - 4gy, gx^2 + gy^2 - 4, gx * gy - 1])
    @test same(groebner([gx + gy + gz, gx * gy + gy * gz + gz * gx, gx * gy * gz - 1], gx, gy, gz),
               [gx + gy + gz, gy^2 + gy * gz + gz^2, gz^3 - 1])                        # cyclic-3
    @test same(groebner([gx + 2gy + 2gz - 1, gx^2 + 2gy^2 + 2gz^2 - gx, 2gx * gy + 2gy * gz - gy], gx, gy, gz),
               [gx - 60gz^3 + 158gz^2 / 7 + 8gz / 7 - 1, gy + 30gz^3 - 79gz^2 / 7 + 3gz / 7,
                gz^4 - 10gz^3 / 21 + gz^2 / 84 + gz / 84])                              # katsura-2
    @test groebner([gx - 1, gx - 2], gx) == [1]                                          # no common zero
    @test same(groebner([gy - gx^2, gz - gx^3], gz, gy, gx), [gz - gx^3, gy - gx^2])     # infinitely many
    @test same(groebner([gy - gx^2, gz - gx^3], gx, gy, gz; order = :grevlex), [gx^2 - gy, gx * gy - gz, gy^2 - gx * gz])
    @test same(groebner([gx^2 + gy^2 - gr^2, gx - gy], gx, gy), [gx - gy, gy^2 - gr^2 / 2])   # a parameter
    @test same(groebner([ga * gx - 1, gx * gy - ga], gx, gy), [gx - 1 / ga, gy - ga^2])
    @test same(groebner([gx^2 - 2, gy - gx]), [gx - gy, gy^2 - 2])                      # all the symbols
    @test same(groebner([gx^2 + gy^2 ~ 1, gx ~ gy], gx, gy), [gx - gy, gy^2 - 1 // 2])   # equations
    @test same(groebner([gx^300 - 1, gx^200 - 1], gx), [gx^100 - 1])
    vs = [Sym(Symbol("gv$i")) for i in 1:7]                                            # wider words
    @test same(groebner([vs[i]^200 - vs[i+1] for i in 1:6], vs...; order = :grevlex)[end:end], [vs[6]^200 - vs[7]])
    # FGLM agrees with lex computed directly, and the inputs are in the ideal
    for (F, xs) in (([gx^2 + 2gy * gz - 3, gy^2 - gx * gz + 1, gz^2 + gx * gy - 2], [gx, gy, gz]),
                    ([gx + gy + gz, gx * gy + gy * gz + gz * gx, gx * gy * gz - 1], [gx, gy, gz]))
        @test same(groebner(F, xs...), SymJ._groebner(F, xs, Sym[], :lex; fglm = false))
        for ord in (:lex, :grlex, :grevlex)
            G = groebner(F, xs...; order = ord)
            @test all(f -> iszero(reduced(f, G, xs...; order = ord)[2]), F)
        end
    end
    # with parameters, FGLM over Q(params) agrees with lex computed directly
    let F = [gx^2 + gy^2 - ga, gx * gy - gr]
        @test all(iszero, simplify.(groebner(F, gx, gy) .- SymJ._groebner(F, [gx, gy], [ga, gr], :lex; fglm = false)))
    end
    # lex with parameters through grevlex and FGLM, not a lex Buchberger that hangs
    # (SymPy #30609); at numeric values of the parameters, the basis there
    let gb = Sym(:gb), gc = Sym(:gc)
        F1 = (1 + gx * gy)^3 * gz + gy^2 * (1 + gx * gy) * (4 + 3gx * gy)
        F2 = gy + 3gx * (1 + gx * gy)^2 * gz + 3gx * gy^2 * (4 + 3gx * gy)
        F3 = 2gx - 3gx^2 * gy - gx^3 * gz
        t = @elapsed G = groebner([F1, F2 - gb, F3 - gc], gx, gy, gz)
        @test t < 30
        @test all(iszero, simplify.(subs.(G, gb => 2, gc => 3) .- groebner([F1, F2 - 2, F3 - 3], gx, gy, gz)))
    end
    # division
    G = groebner([gx^2 + gy^2 - 1, gx - gy], gx, gy)
    q, r = reduced(gx^3 + gy, G, gx, gy)
    @test r == 3gy / 2 && iszero(expand(q[1] * G[1] + q[2] * G[2] + r - (gx^3 + gy)))
    q, r = reduced(gx^2 * gy - gr * gx, [gx * gy - gr, gy^2 - 1], gx, gy; order = :grevlex)
    @test iszero(expand(q[1] * (gx * gy - gr) + q[2] * (gy^2 - 1) + r - (gx^2 * gy - gr * gx)))
    @test_throws ArgumentError groebner([gx], gx; order = :revlex)
    @test_throws ArgumentError groebner([gx^2 - 0.5], gx)
    @test_throws ArgumentError groebner([sin(gx) - gy], gx, gy)
    @test_throws ArgumentError groebner([gx < gy], gx, gy)
end

@testset "resultants" begin
    @syms ra rb rc rd re rf rg rh ry rz
    for m in (:auto, :modular, :subresultant)
        R(p, q) = resultant(p, q, x; method = m)
        @test R(x^2 - 2, x^2 - 3) == 1
        @test R(x^2 + ry^2 - 1, x - ry) == 2ry^2 - 1
        @test R(ra * x^2 + rb * x + rc, 2ra * x + rb) == 4ra^2 * rc - ra * rb^2
        # res(Π(x - aᵢ), Π(x - cⱼ)) = Π(aᵢ - cⱼ)
        @test iszero(expand(R((x - ra) * (x - rb), (x - rc) * (x - rd)) - (ra - rc) * (ra - rd) * (rb - rc) * (rb - rd)))
        @test R((x - 1) * (x + 2), (x - 1) * (x^2 + 5)) == 0                  # a common root
        @test R(x^2 / 2 - 1 // 3, x / 5 + 1) == Sym(73 // 150)
        @test R(x^3 + ry, x^2 * ry - rz) == ry^5 - rz^3
        @test R(Sym(3), x^2 + 1) == 9 && R(x^2 + 1, Sym(3)) == 9 && R(Sym(2), Sym(3)) == 1
        @test R(ry * x^2 + x + 1, x - ry) == ry^3 + ry + 1                     # lc vanishes at ry = 0
        @test R(x - ry^100, x^3 - 2) == ry^300 - 2
        @test R(Sym(0), x + 1) == 0
        # res(B, A) = (-1)^(mn) res(A, B)
        @test R(x^2 + x + 1, x^3 - ry) == R(x^3 - ry, x^2 + x + 1)          # mn even
        @test R(x + ry, x^3 - 2) == -R(x^3 - 2, x + ry)
    end
    # the methods agree on random dense and sparse inputs
    st = Ref(UInt64(7))                       # a small LCG: no Random dependency
    rnd(r) = (st[] = st[] * 6364136223846793005 + 1442695040888963407; first(r) + Int((st[] >> 33) % length(r)))
    rp(vs, d, k) = sum(rnd(-9:9) * prod(v^rnd(0:d) for v in vs) for _ in 1:k)
    for (vs, d, k) in (([x], 12, 10), ([x, ry], 5, 10), ([x, ry, rz], 3, 8))
        for _ in 1:3
            p, q = rp(vs, d, k), rp(vs, d, k)
            has(p, x) && has(q, x) || continue
            r1 = resultant(p, q, x; method = :modular)
            @test iszero(expand(r1 - resultant(p, q, x; method = :subresultant)))
        end
    end
    @test nterms(Poly(resultant(ra*x^3 + rb*x^2 + rc*x + rd, re*x^3 + rf*x^2 + rg*x + rh, x),
                      ra, rb, rc, rd, re, rf, rg, rh)) == 34                  # the known 34 terms
    # polynomials, coefficients in functions, floats
    P = resultant(Poly(x^2 + ry^2 - 1, x, ry), Poly(x - ry, x, ry), x)
    @test P isa Poly && Sym(P) == 2ry^2 - 1
    @test Sym(resultant(Poly(x^2 - 2, x), Poly(x - 3, x))) == 7
    @test resultant(x^2 - sin(ry), x - cos(ry), x) == cos(ry)^2 - sin(ry)
    @test N(resultant(x^2 - 2.0, x - 1.5, x)) ≈ 0.25
    @test_throws ArgumentError resultant(sin(x), x - 1, x)
    @test_throws ArgumentError resultant(x, x, x; method = :fast)
    # discriminants
    @test discriminant(ra * x^2 + rb * x + rc, x) == rb^2 - 4ra * rc
    @test discriminant(x^3 + ra * x + rb, x) == -4ra^3 - 27rb^2
    @test discriminant((x - 1)^2 * (x + 1), x) == 0
    @test discriminant(3x + 1, x) == 1
    @test Sym(discriminant(Poly(x^2 - 5, x))) == 20
    @test_throws ArgumentError discriminant(Sym(4), x)
end

@testset "minimal_polynomial" begin
    MP = SymJ.minimal_polynomial
    s2, s3 = sqrt(Sym(2)), sqrt(Sym(3))
    @test MP(Sym(3) / 4, x) == 4x - 3
    @test MP(s2, x) == x^2 - 2
    @test MP(s2 + s3, x) == x^4 - 10x^2 + 1
    @test MP((1 + sqrt(Sym(5))) / 2, x) == x^2 - x - 1
    @test MP(IM, x) == x^2 + 1
    @test MP(1 / (1 + s2), x) == x^2 + 2x - 1
    @test MP(sqrt(2 + s3), x) == x^4 - 4x^2 + 1
    @test MP(cos(PI / 7), x) == 8x^3 - 4x^2 - 4x + 1
    @test MP(sin(PI / 5), x) == 16x^4 - 20x^2 + 5
    @test MP(tan(PI / 12), x) == x^2 - 4x + 1
    @test MP((-1)^(Sym(2) // 5), x) == x^4 + x^3 + x^2 + x + 1     # a root of unity
    @test MP(Sym(-8)^(Sym(1) // 3), x) == x^2 - 2x + 4                # the principal cube root
    @test MP(cbrt(Sym(-8)), x) == x + 2                               # the real one
    for e in (Sym(2)^(1 // 3) + IM, sqrt(Sym(2)) + sqrt(Sym(3)) + sqrt(Sym(5)), Sym(2)^(1 // 4) + Sym(3)^(1 // 3))
        p = MP(e, x)
        @test abs(N(subs(p, x => e))) < 1e-8
    end
    @test_throws ArgumentError MP(PI, x)
    @test_throws ArgumentError MP(log(Sym(2)), x)
    @test_throws ArgumentError MP(Sym(0.5), x)
    @test_throws ArgumentError MP(x + 1, x)
end

@testset "decompose" begin
    @syms a
    D = SymJ.decompose
    @test D(x^4 + 2x^3 - x - 1, x) == [x^2 - x - 1, x^2 + x]
    @test D(expand((x^2 + x)^3 - 2), x) == [x^3 - 2, x^2 + x]
    @test D(expand((x^3 - 3x)^3 - 3(x^3 - 3x) - 1), x) == [x^3 - 3x - 1, x^3 - 3x]
    p = expand(2((x^2 + 1)^2 + 3)^2 - 5)
    gs = D(p, x)
    @test length(gs) == 3
    @test expand(subs(gs[1], x => subs(gs[2], x => gs[3]))) == p
    @test D(expand((x^2 + a * x)^2 + a), x) == [x^2 + a, x^2 + a * x]   # symbolic coefficients
    @test D(x^5 + x + 1, x) == [x^5 + x + 1]                             # prime degree
    @test D(x^6 + x + 1, x) == [x^6 + x + 1]                             # indecomposable
    @test D(x^2 + 1, x) == [x^2 + 1]
    @test D(Sym(3), x) == [3]
    @test_throws ArgumentError D(sin(x), x)
end

@testset "real root isolation" begin
    RR, RI = SymJ.real_roots, SymJ.root_intervals
    @test RR(x^3 - 2x - 5, x) == [2.0945514815423265] && RR(x^2 + 1, x) == Float64[] && RR(Sym(5), x) == Float64[]
    @test RR((x - 1)^2 * (x^2 - 2), x) == [-sqrt(2), 1.0, 1.0, sqrt(2)]
    @test RR(prod(x - k for k in 1:20), x) == Float64.(1:20)                   # Wilkinson
    @test RR(x * (x - 1 // 2) * (x + 3 // 7), x) == [-3 / 7, 0.0, 0.5]
    @test RR(x^5 - x - 1, x; T = BigFloat)[1] ≈ big"1.167303978261418684256045899854842180720560371525489039140082449275651903429519"
    I = RI((x - 1)^2 * (x^2 - 2), x)
    @test length(I) == 3 && I[2] == (1, 1, 2) && all(t -> t[1] < t[2] || t[1] == t[2], I)
    @test all(t -> t[2] - t[1] < 1 // 100, RI(x^2 - 2, x; width = 1 // 100))
    # an exact root next to an interval of the same factor (its endpoint)
    rs = [-11 // 7, 0, 5 // 7, -27 // 7, -2, -23 // 7, 18 // 7]
    @test RR(prod(x - r for r in rs), x) == sort(Float64.(rs))
    @test_throws ArgumentError RR(x^2 - sqrt(Sym(2)), x)
    @test_throws ArgumentError RR(Sym(0), x)
    @test_throws ArgumentError RR(sin(x), x)
end
