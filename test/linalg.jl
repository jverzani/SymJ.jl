# Linear algebra (included by runtests.jl)

@testset "linear algebra" begin
    @syms a b c d λ
    LA = SymJ.LinearAlgebra
    # det: Bareiss for polynomials, elimination over rational functions otherwise
    @test LA.det([a b; c d]) == a * d - b * c && LA.det([x 1; 1 x]) == x^2 - 1 && LA.det(Sym[1 2; 3 4]) == -2
    @test LA.det(Sym[1 2; 2 4]) == 0 && LA.det(Matrix{Sym}(undef, 0, 0)) == 1
    @test N(LA.det(Sym[1.5 2; 2 3])) ≈ 0.5                             # floating point entries
    @test N.(LA.inv(Sym[1.5 2; 2 3])) ≈ [6 -4; -4 3]
    @test LA.det([1/x 1; 1 x]) == 0 && LA.det([sqrt(Sym(2)) 1; 1 sqrt(Sym(2))]) == 1
    G = [Sym(Symbol("g", i, j)) for i in 1:5, j in 1:5]
    @test length(LA.det(G).args) == 120                                # all 5! terms
    @test LA.det([Sym(1) / (i + j - 1) for i in 1:5, j in 1:5]) == 1 // 266716800000   # Hilbert
    @test LA.det([x 1; 1 x] - λ * LA.I(2)) == expand((x - λ)^2 - 1)  # characteristic polynomial
    # inv: the adjugate for polynomial entries, Gauss-Jordan over rational functions else
    @test inv([a b; c d]) == [d -b; -c a] / (a * d - b * c) && inv(Sym[1 2; 3 4]) == [-2 1; 3//2 -1//2]
    @test inv([x 1; 1 x]) == [x -1; -1 x] / (x^2 - 1) && inv(Sym[2;;]) == Sym[1//2;;]
    @test_throws LA.SingularException inv(Sym[1 2; 2 4])
    @test_throws LA.SingularException inv([1/x 1; 1 x])
    # dot and norm, entrywise (symbols are complex)
    @test LA.dot(x, a) == conj(x) * a
    @test LA.dot([x, a], [1, 2]) == conj(x) + 2conj(a)
    @test LA.dot([1, 2], [x, a]) == x + 2a
    @test LA.dot([x, a, Sym(1)], [Sym(2), Sym(3), b])(x => 2, a => 3, b => 5) == 18
    @test LA.norm(Sym[3, 4]) == 5
    @test LA.norm([x, a]) == sqrt(abs(x)^2 + abs(a)^2)
    @test LA.norm([x, a], 1) == abs(x) + abs(a)
    @test LA.norm(Sym[3, -4], Inf) == 4
    @test LA.norm(Sym[1, 1], 3) == Sym(2)^(1 // 3)
    @test LA.norm(x) == abs(x)
    @test LA.normalize(Sym[3, 4]) == Sym[3 // 5, 4 // 5]
    @test refine(LA.dot([x, a], [x, a]); assume = @assume(x in Reals, a in Reals)) == x^2 + a^2
    @test refine(LA.norm([x, a]); assume = @assume(x in Reals, a in Reals)) == sqrt(x^2 + a^2)
    # nine facts, as many as Assumptions has fields
    vs = [Sym(Symbol(:v, i)) for i in 1:9]
    @test assuming(() -> refine(sqrt(vs[9]^2)), (v ∈ ℝ for v in vs)...) == abs(vs[9])
    G3 = G[1:3, 1:3]
    @test all(iszero, map(simplify, inv(G3) * G3 - LA.I))
    @test map(simplify, inv([sin(x) cos(x); -cos(x) sin(x)])) == [sin(x) -cos(x); cos(x) sin(x)]
    # lu: a LinearAlgebra.LU, first-nonzero pivoting
    F = LA.lu([a b; c d])
    @test F isa LA.LU && F.L == [1 0; c/a 1] && F.p == [1, 2] && LA.det(F) == a * d - b * c
    @test map(simplify, F.L * F.U - [a b; c d][F.p, :]) == zeros(Sym, 2, 2)
    P = LA.lu(Sym[0 1; 1 0])
    @test P.p == [2, 1] && P.U == Sym[1 0; 0 1]
    @test_throws LA.SingularException LA.lu(Sym[1 2; 2 4])
    @test LA.lu(Sym[1 2; 2 4]; check = false).info == 2
    # \ exactly
    @test [a b; c d] \ [x, y] == [(d * x - b * y) / (a * d - b * c), (a * y - c * x) / (a * d - b * c)]
    @test Sym[1 2; 3 4] \ Sym[5, 6] == [-4, 9 // 2]
    # matrices with I and with expressions as scalars
    @test [x 1; 1 x] - LA.I == [x-1 1; 1 x-1] && 2x * [x 1; 1 x] == [2x^2 2x; 2x 2x^2]
    @test [x 1; 1 x] / x == [1 1/x; 1/x 1]
    @test_throws ArgumentError λ * LA.I
    # rref, rank, nullspace (generic: symbolic pivots nonzero)
    @test SymJ.rref(Sym[1 2 3; 4 5 6; 7 8 9]) == Sym[1 0 -1; 0 1 2; 0 0 0] && LA.rank(Sym[1 2 3; 4 5 6; 7 8 9]) == 2
    @test LA.nullspace(Sym[1 2 3; 4 5 6; 7 8 9]) == Sym[1; -2; 1;;]
    @test LA.rank([a b; 2a 2b]) == 1 && LA.nullspace([a b; 2a 2b]) == [-b/a; 1;;] && SymJ.rref([a 1; 1 1]) == Sym[1 0; 0 1]
    # algebraic constants are numbers, not independent atoms: sqrt(2)^2 - 2 is zero
    s2 = sqrt(Sym(2))
    @test LA.rank([s2 1; 2 s2]) == 1 && LA.det([s2 1; 2 s2]) == 0 && LA.nullspace([s2 1; 2 s2]) == [-s2/2; 1;;]
    @test_throws LA.SingularException inv([s2 1; 2 s2])
    @test inv([s2 1; 1 s2]) == [s2 -1; -1 s2]
    # eigenvalues with multiplicity; eigenvectors; eigen for diagonalizable matrices
    @test LA.eigvals(Sym[2 1; 1 2]) == [1, 3] && LA.eigvals(Sym[1 1; 0 1]) == [1, 1]
    @test LA.eigvals(Sym[1 2 0; 2 1 0; 0 0 3]) == [-1, 3, 3] && Set(LA.eigvals([a b; b a])) == Set([a + b, a - b])
    @test LA.eigvals(Sym[0 -1; 1 0]) == [-IM, IM] && LA.eigvecs(Sym[3 0; 0 3]) == Sym[1 0; 0 1]
    @test size(LA.eigvecs(Sym[1 1; 0 1])) == (2, 1)                     # defective
    @test_throws ArgumentError LA.eigen(Sym[1 1; 0 1])
    for A in (Sym[2 1; 1 2], [a b; b a], Sym[2 0 0; 0 3 4; 0 4 9], [x 1; 1 x], Sym[1 2; 3 4], Sym[0 -1; 1 0])
        F = LA.eigen(A)
        @test F isa LA.Eigen && all(iszero, map(simplify, A * F.vectors - F.vectors * LA.Diagonal(F.values)))
    end
    @test SymJ.charpoly([a b; c a], x) == expand((x - a)^2 - b * c)
    # without division, so no spurious roots from an uncancelled denominator (SymPy #28212, #28793)
    let M = Sym[sqrt(Sym(2)) PI; sqrt(Sym(3)) E], C = Sym[IM -IM 0; -IM IM -IM; 0 -IM IM]
        @test length(LA.eigvals(M)) == 2
        @test sort(real.(N.(LA.eigvals(M)))) ≈ sort(real.(LA.eigvals(Float64.(N.(M)))))
        @test length(LA.eigvals(C)) == 3
        @test sort(imag.(N.(LA.eigvals(C)))) ≈ [1 - sqrt(2), 1, 1 + sqrt(2)]
        @test SymJ.charpoly(Sym[1 2 3; 4 5 6; 7 8 10], x) == x^3 - 16x^2 - 12x + 3
    end
    # matrices of square roots: eigvals hung in elimination (SymPy #28606, #28507)
    let M = Sym[sqrt(Sym(3)) 4 sqrt(Sym(19)); sqrt(Sym(62)) sqrt(Sym(14)) 5; sqrt(Sym(7)) sqrt(Sym(59)) 5],
        M5 = Sym[45 15 20 18 sqrt(Sym(5))/2; 15 85 35 22 Sym(11)//37; 20 35 15 40 3sqrt(Sym(3))/7;
                 Sym(9)//2 2 4 12 sqrt(Sym(7))/5; Sym(25)//13 sqrt(Sym(7))/4 2 1 1]
        # (rounded: the real parts of a conjugate pair differ in the last digits by platform)
        by(z) = (round(real(z); digits = 8), round(imag(z); digits = 8))
        num(A) = sort(complex.(LA.eigvals(ComplexF64.(N.(A)))); by)
        @test sort(complex.(N.(LA.eigvals(M))); by) ≈ num(M)
        @test sort(complex.(N.(LA.eigvals(M5))); by) ≈ num(M5)
        @test size(LA.eigvecs(Sym[sqrt(Sym(2)) PI; sqrt(Sym(3)) E])) == (2, 2)
    end
    # powers of im with the sign outside: im^3 = 1/im = -im
    @test IM^3 == -IM && 1 / IM == -IM && IM^-1 + IM == 0 && (-1)^(-1 // 3 * one(Sym)) == -(-1)^(2 // 3 * one(Sym))
    @test simplify(1 / (1 + IM)) == (1 - IM) / 2
    @test LA.eigvals([cos(x) -sin(x); sin(x) cos(x)]) == [cos(x) - IM * sin(x), cos(x) + IM * sin(x)]
    # QR: Gram-Schmidt, exact
    for A in (Sym[1 1; 1 0], Sym[3 1; 4 2], Sym[1 2; 3 4; 5 6])
        Q, R = LA.qr(A)
        @test all(iszero, map(simplify, Q * R - A)) && map(simplify, transpose(Q) * Q) == Matrix{Sym}(LA.I(size(A, 2)))
        @test all(i -> iszero(R[i, 1]), 2:size(R, 1))
    end
    @test LA.qr(Sym[3 1; 4 2]).Q == Sym[3//5 -4//5; 4//5 3//5] && LA.qr(Sym[3 1; 4 2]).R == Sym[5 11//5; 0 2//5]
    @test_throws ArgumentError LA.qr(Sym[1 2; 2 4])
    # SVD: exact, checked numerically where simplify can't (nested radicals)
    F = LA.svd(Sym[3 0; 4 5])
    @test F.S == [3sqrt(Sym(5)), sqrt(Sym(5))] && map(simplify, F.U * LA.Diagonal(F.S) * F.Vt) == Sym[3 0; 4 5]
    @test LA.svdvals(Sym[1 1; 0 0]) == [sqrt(Sym(2)), 0] && LA.svdvals(Sym[2 0; 0 -3]) == [3, 2]
    numf(M) = map(e -> ComplexF64(N(e)), M)
    for A in (Sym[1 2; 3 4], Sym[1 im; 0 1], Sym[2 1 0; 1 2 1], Sym[1 1; 1 1; 0 0], Sym[1 0 1; 0 1 1])
        F = LA.svd(A)
        U, S, Vt = numf(F.U), numf(F.S), numf(F.Vt)
        @test isapprox(U * LA.Diagonal(S) * Vt, numf(A); atol = 1e-12) && isapprox(U' * U, LA.I; atol = 1e-12) &&
              isapprox(Vt * Vt', LA.I; atol = 1e-12) && isapprox(real.(S), LA.svdvals(numf(A)); atol = 1e-12)
    end
    F = LA.svd(Sym[1 1; 1 1; 0 0]; full = true)
    @test size(F.U) == (3, 3) && map(simplify, transpose(F.U) * F.U) == Sym[1 0 0; 0 1 0; 0 0 1]
    @test_throws ArgumentError LA.svd(Sym[1 2 3; 4 5 6; 7 8 10])          # a cubic: refused, not slow
    # symbolic singular values refined under the assumptions
    @syms a b
    @test LA.svdvals(Sym[a 0; 0 b]) == [sqrt(b^2), sqrt(a^2)]
    @test LA.svdvals(Sym[a 0; 0 b]; assume = @assume(a > 0, b > 0)) == [b, a]
    @test assuming(() -> LA.svdvals(Sym[a 0; 0 b]), @assume(a ∈ ℝ, b ∈ ℝ)) == [abs(b), abs(a)]
    @test (@assuming a > 0 begin LA.svd(Sym[a a; a a]).S end) == [2a, 0]
    # Jordan form
    for A in (Sym[1 1; 0 1], Sym[2 1 0; 0 2 0; 0 0 3], Sym[5 4 2 1; 0 1 -1 -1; -1 -1 3 0; 1 1 -1 2], Sym[2 1; 1 2])
        P, J = SymJ.jordan_form(A)
        @test all(iszero, map(simplify, P * J * inv(P) - A))
    end
    @test SymJ.jordan_form(Sym[5 4 2 1; 0 1 -1 -1; -1 -1 3 0; 1 1 -1 2])[2] == Sym[1 0 0 0; 0 2 0 0; 0 0 4 1; 0 0 0 4]
    @test LA.eigvals(Sym[2 0 0; 0 3 4; 0 4 9]) == [1, 2, 11]                 # constants sorted
    # matrix functions through the Jordan form
    @syms t n
    @test exp(t * Sym[0 1; -1 0]) == [cos(t) sin(t); -sin(t) cos(t)]
    @test exp(Sym[1 1; 0 1]) == [E E; 0 E] && exp(t * Sym[2 1; 0 2]) == [exp(2t) t*exp(2t); 0 exp(2t)]
    @test sin(Sym[0 1; 0 0]) == Sym[0 1; 0 0] && cos(Sym[0 1; 0 0]) == Sym[1 0; 0 1]
    Fib = Sym[1 1; 1 0]^n
    @test [round(Int, N(subs(Fib[1, 2], n => k))) for k in 1:10] == [1, 1, 2, 3, 5, 8, 13, 21, 34, 55]
    @test isapprox(N.(subs.(exp(t * Sym[1 2; 3 4]), t => 0.3)), exp(0.3 * [1 2; 3 4]); rtol = 1e-10)
    # sqrt, log and non-integer powers (stdlib's Schur route needs floats)
    @test sqrt(Sym[4 1; 0 4]) == Sym[2 1//4; 0 2]
    @test Sym[4 1; 0 4]^(3//2) == Sym[8 3; 0 8]
    @test Sym[4 1; 0 4]^(2//1) == Sym[16 8; 0 16]
    @test log(Sym[2 1; 0 2]) == [log(Sym(2)) 1//2; 0 log(Sym(2))]
    @test sqrt(Sym[0 0; 0 4]) == Sym[0 0; 0 2]
    let x = Sym(:x), S = Sym[1 x; x 2]^(1//2)
        @test isapprox(N.(subs.(S * S, x => 3)), [1.0 3; 3 2]; atol = 1e-10)
    end
    @test_throws ArgumentError sqrt(Sym[0 1; 0 0])
    @test_throws ArgumentError log(Sym[0 1; 0 0])
    # pseudoinverse from the rank factorization
    @test LA.pinv(Sym[1 2; 2 4]) == Sym[1//25 2//25; 2//25 4//25]
    @test LA.pinv(Sym[1 2 3; 4 5 6]) == Sym[-17//18 4//9; -1//9 1//9; 13//18 -2//9]
    @test LA.pinv(zeros(Sym, 2, 3)) == zeros(Sym, 3, 2)
    let a = Sym(:a), M = Sym[a 1; 2a 2]
        @test map(cancel, M * LA.pinv(M) * M - M) == zeros(Sym, 2, 2)
    end
    let x = Sym(:x)
        @test N.(subs.(LA.pinv(Sym[1 x; 2 2x; 0 1]), x => 2)) ≈ LA.pinv([1 2.0; 2 4; 0 1])
    end
end

@testset "Smith and Hermite normal forms" begin
    diagonal(D) = [D[i, i] for i in 1:minimum(size(D))]
    # integers: the invariant factors, U A V = D with U, V unimodular
    A = Sym[2 4 4; -6 6 12; 10 -4 -16]
    D, U, V = SymJ.smith_normal_decomp(A)
    @test diagonal(D) == Sym[2, 6, 12]
    @test U * A * V == D
    @test abs(SymJ.LinearAlgebra.det(U)) == 1
    @test abs(SymJ.LinearAlgebra.det(V)) == 1
    @test diagonal(SymJ.smith_normal_form(Sym[12 6 4; 3 9 6; 2 16 14])) == Sym[1, 10, 30]
    @test diagonal(SymJ.smith_normal_form(Sym[1 2 3; 4 5 6; 7 8 9])) == Sym[1, 3, 0]
    @test SymJ.smith_normal_form(Sym[6 4; 0 0; 2 8]) == Sym[2 0; 0 20; 0 0]
    # polynomials in one variable: monic invariant factors
    P = Sym[x^2 - 1 x + 1; x - 1 x^2]
    D2, U2, V2 = SymJ.smith_normal_decomp(P)
    @test diagonal(D2) == [Sym(1), expand((x^2 - 1)^2)]
    @test iszero(expand.(U2 * P * V2 - D2))
    @test diagonal(SymJ.smith_normal_form(Sym[x 1; 0 x])) == [Sym(1), x^2]
    # Hermite, row form: U A = H, pivots positive, entries above them reduced
    B = Sym[2 3 6 2; 5 6 1 6; 8 3 1 1]
    H, W = SymJ.hermite_normal_decomp(B)
    @test H == Sym[1 0 50 -11; 0 3 28 -2; 0 0 61 -13]
    @test W * B == H
    @test abs(SymJ.LinearAlgebra.det(W)) == 1
    @test SymJ.hermite_normal_form(Sym[4 6; 2 3; 6 9]) == Sym[2 3; 0 0; 0 0]
    @test_throws ArgumentError SymJ.smith_normal_form(Sym[x 1; 1 y])
    @test_throws ArgumentError SymJ.smith_normal_form(Sym[sin(x) 1; 0 1])
end


@testset "Cholesky, LDLᵀ, permanent, rotation matrices" begin
    @syms a b c d θ
    LA = SymJ.LinearAlgebra
    A = Sym[a 1; 1 a]
    # L as the transpose of U, symbols real as for U (L conjugated them: conj(sqrt(a)))
    C = LA.cholesky(A)
    @test C.L == Sym[sqrt(a) 0; 1/sqrt(a) sqrt(a - 1/a)]
    @test simplify.(C.L * C.U) == A
    L, U = C
    @test U == C.U
    # check = false, NoPivot, Hermitian and Symmetric take the same path (they went to
    # LinearAlgebra's generic code: wrong, or a MethodError)
    @test LA.cholesky(A; check = false).U == C.U
    @test LA.cholesky(Sym[4 2; 2 3], LA.NoPivot()).U == Sym[2 1; 0 sqrt(Sym(2))]
    @test LA.cholesky(LA.Hermitian(Sym[4 2; 2 3])).L == Sym[2 0; 1 sqrt(Sym(2))]
    @test LA.cholesky(LA.Symmetric(A)).L == C.L
    @test !LA.issuccess(LA.cholesky(Sym[1 2; 2 1]; check = false))
    @test LA.cholesky(Sym[4 2; 2 3]) \ Sym[2, 1] == Sym[1 // 2, 0]
    @test simplify.(C \ Sym[1, 0]) == Sym[a / (a^2 - 1), -1 / (a^2 - 1)]
    Ci = LA.cholesky(Sym[2 IM; -IM 2])
    @test simplify.(Ci.L * Ci.U) == Sym[2 IM; -IM 2]
    @test_throws ArgumentError LA.cholesky(Sym[1 2; 3 4])
    # LDLᵀ without square roots (SymPy #28477)
    F = LA.ldlt(Sym[4 2; 2 3])
    @test Matrix(F.L) == Sym[1 0; 1//2 1]
    @test F.D == LA.Diagonal(Sym[4, 2])
    @test F \ Sym[2, 1] == Sym[1 // 2, 0]
    G = LA.ldlt(A)
    @test G.D == LA.Diagonal(Sym[a, a - 1/a])
    @test simplify.(G.L * G.D * transpose(Matrix(G.L))) == A
    @test LA.ldlt(Sym[1 2; 2 1]).D == LA.Diagonal(Sym[1, -3])          # indefinite is fine
    @test_throws LA.ZeroPivotException LA.ldlt(Sym[0 1; 1 0])
    # Symmetric(A) prints (symmetric(::Sym, ::Symbol) was missing)
    @test occursin("a  1", sprint(show, MIME"text/plain"(), LA.Symmetric(A)))
    # the permanent, rectangular and empty too (SymPy #28539)
    per = SymJ.permanent
    @test per(Sym[a b; c d]) == a * d + b * c
    @test per(Sym[1 2 3; 4 5 6; 7 8 9]) == 450
    @test per(Sym[1 2 3; 4 5 6]) == 58
    @test per(Sym[1 4; 2 5; 3 6]) == 58
    @test per(zeros(Sym, 0, 3)) == 1
    @test per(zeros(Sym, 3, 0)) == 1
    # rotation matrices, right-handed (SymPy #24936, #24528)
    R = SymJ.rotation_matrix
    @test R(θ) == Sym[cos(θ) -sin(θ); sin(θ) cos(θ)]
    @test R(θ, :z) == Sym[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
    @test R(θ, 1) == Sym[1 0 0; 0 cos(θ) -sin(θ); 0 sin(θ) cos(θ)]
    @test R(PI / 2, :y) * Sym[0, 0, 1] == Sym[1, 0, 0]                  # z to x about y
    @test R(PI / 2, :x) * Sym[0, 1, 0] == Sym[0, 0, 1]
    @test R(2PI / 3, [1, 1, 1]) * Sym[1, 0, 0] == Sym[0, 1, 0]          # about a vector
    @test R(θ, [0, 0, 5]) == R(θ, 3)
    @test simplify.(R(θ, [1, 2, 2]) * transpose(R(θ, [1, 2, 2]))) == Matrix{Sym}(LA.I(3))
    @test_throws ArgumentError R(θ, 4)
    @test_throws ArgumentError R(θ, [0, 0, 0])
end
