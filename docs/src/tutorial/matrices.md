# Matrices

```@setup tut
using SymJ
using RischIntegrate
@syms x y a b c
```

Matrices of expressions are Julia arrays of `Sym`. With
`LinearAlgebra`, `det`, `inv`, `lu` and `\` work exactly, pivoting (for exact entries) on
the first nonzero entry rather than comparing magnitudes:

```@example tut
using LinearAlgebra
@syms d
A = [a b; c d]
det(A), inv(A)
```

```@example tut
F = lu(A)
F.L, F.U
```

```@example tut
A \ [x, y]
```

Row reduction (`SymJ.rref`, not exported, as it is not a LinearAlgebra
function), rank, null spaces, and eigenvalues (with multiplicity) and
eigenvectors:

```@example tut
SymJ.rref(Sym[1 2 3; 4 5 6; 7 8 9]), nullspace(Sym[1 2 3; 4 5 6; 7 8 9])
```

```@example tut
eigvals(Sym[1 2; 3 4]), eigvecs([a b; b a])
```

`dot` and `norm` work entry by entry. Symbols are complex, so `dot`
conjugates its first argument; for real symbols, `refine` under the
assumptions removes the conjugates:

```@example tut
u = [x, y, 1]
(
dot(u, [1, 2, 3]),
norm(u),
refine(dot(u, u); assume = @assume(x in Reals, y in Reals))
)
```

QR factorization by exact Gram–Schmidt, and functions of matrices
through the Jordan form (`SymJ.jordan_form`), such as the exponential
that solves `x' = A x`:

```@example tut
Q, R = qr(Sym[3 1; 4 2])
Q, R
```

The singular value decomposition, exactly, for small matrices:

```@example tut
F = svd(Sym[3 0; 4 5])
F.S, F.U, F.Vt
```

```@example tut
@syms t n
exp(t * Sym[0 1; -1 0]), exp(t * Sym[2 1; 0 2])
```

```@example tut
Sym[1 1; 1 0]^n
```

`A - I` works as for numbers; a symbolic multiple of the identity is `λ*I(n)`:

```@example tut
@syms λ
det([x 1; 1 x] - λ * I(2))
```
