# Feature audit: SymJ against SymPy and Giac

What SymPy and Giac provide, area by area, and where SymJ stands. Written 2026-10-23
against SymPy 1.14 and Giac (through Giac.jl 0.14), with SymJ 0.1.0 plus GruntzLimit
(limits) and RischIntegrate (integrals). Every SymPy name below was checked to exist in
SymPy 1.14, and every Giac command in Giac.jl's command table; SymJ's column was checked
against its source and by running it.

**Legend:** ✅ has · ◐ partial · ❌ missing · ★ planned (series, assumptions, a basic
`simplify`, summation, differential equations)

## Numbers and the core

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Exact integers and rationals of any size | `Integer`, `Rational` | built in | ✅ `SymNumber` (Int fast path, BigInt fallback) |
| Exact radicals, π, ℯ, exact function values | `sqrt`, `pi`, `E` | built in | ✅ `sqrt(Sym(8)) == 2sqrt(2)`, `sin(π/3) == sqrt(3)/2` |
| Combining and denesting radicals | `radsimp`, `sqrtdenest` | `simplify` | ❌ `sqrt(2)*sqrt(3)` stays; `cancel` handles square roots of rationals |
| Complex numbers: `re`, `im`, `conj`, `arg` | `re`, `im`, `conjugate`, `arg` | `re`, `im`, `conj`, `arg` | ◐ `im` is exact; no real/imaginary parts (they need assumptions) |
| Symbolic infinity and undefined | `oo`, `zoo`, `nan` | `infinity`, `undef` | ❌ limits return Julia's `Inf` |
| Numeric evaluation | `evalf(n)` | `evalf` | ✅ `N` (Float64), calling with `BigFloat`, dual numbers, … |
| Assumptions on symbols | `symbols(positive=True)`, `ask`, `refine`, `Q` | `assume` | ❌ ★ (to be passed per call, not stored on symbols) |
| Parameters (named constants) | — | — | ◐ implemented, hidden until the shared-name question is settled |

## Expressions

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Canonical form as built, fast construction | ✅ | ✅ | ✅ (3–30x faster than SymEngine on the benchmarks) |
| Substitution | `subs` | `subst` | ✅ `subs`, calling with pairs or keywords |
| Compiled numeric functions | `lambdify` | — | ✅ `lambdify` (as fast as hand-written code) |
| Generic term interface | — | — | ✅ TermInterface |
| Pattern matching and rewrite rules | `Wild`, `match`, `replace` | rules | ❌ (TermInterface-based rewriters can work on SymJ terms) |
| Collect by powers, coefficients | `collect`, `coeff` | `collect`, `coeff` | ◐ through `Poly` (`coefficients`, `exponent_vectors`) |
| Horner form | `horner` | `horner` | ❌ |
| Common subexpressions | `cse` | — | ◐ `lambdify` shares them; no user-level `cse` |
| Rewriting between forms (exp ↔ trig, …) | `ex.rewrite(...)` | `exp2pow`, `pow2exp`, `halftan`, `sincos`, `tan2sincos` | ❌ |
| Floats to exact numbers | `nsimplify` | — | ❌ |

## Simplification ★

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| General simplification | `simplify` | `simplify`, `simplifier` | ❌ ★ |
| Trigonometric | `trigsimp`, `exptrigsimp` | `tlin`, `tcollect`, `texpand` | ❌ |
| Powers and logarithms | `powsimp`, `powdenest`, `logcombine` | `lncollect`, `lnexpand` | ❌ (`exp(x)*exp(y) == exp(x + y)` is canonical) |
| Gamma and combinatorial | `gammasimp`, `combsimp` | — | ❌ |
| Rational functions | `cancel`, `together`, `apart`, `ratsimp` | `normal`, `ratnormal`, `partfrac`, `propfrac` | ✅ `cancel`, `apart`, `numer_denom` (◐ no separate `together`) |
| Expansion | `expand`, `expand_trig`, `expand_log` | `expand`, `texpand`, `lnexpand` | ✅ `expand` (fast, multinomial); ❌ trig and log expansion |

## Polynomials

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Polynomial type and arithmetic | `Poly` | built in | ✅ `Poly`, packed monomials |
| Division, gcd, lcm | `div`, `rem`, `gcd`, `lcm` | `quo`, `rem`, `gcd`, `lcm` | ✅ modular gcd, Zippel chosen automatically (faster than FLINT on many benchmarks) |
| Factor over ℚ, any number of variables | `factor`, `factor_list` | `factor` | ✅ Zassenhaus, van Hoeij, Wang's EEZ |
| Squarefree decomposition | `sqf_list` | `sqrfree` | ✅ `factor_squarefree`, `squarefree_part` |
| Factor over extensions, ℤ[i], ℤ/p | `factor(…; extension, gaussian, modulus)` | `cfactor`, factor mod p | ❌ |
| Resultant, discriminant | `resultant`, `discriminant` | `resultant`, `discriminant` | ✅ modular or subresultant, chosen per input |
| Gröbner bases | `groebner` | `gbasis` | ✅ Buchberger with Gebauer–Möller and sugar; lex through FGLM; `reduced` |
| Roots: exact to degree 4, numeric beyond | `roots`, `nroots` | `solve`, `proot` | ✅ through `solve` (cubics with three real roots in trigonometric form) |
| Real root isolation and counting | `real_roots`, `count_roots`, `sturm` | `realroot`, `sturm` | ◐ numeric roots only |
| Interpolation | `interpolate` | `lagrange` | ❌ |
| Orthogonal and cyclotomic polynomials | `legendre`, `chebyshevt`, `hermite`, `cyclotomic_poly` | `legendre`, `tchebyshev1`, `hermite`, `cyclotomic` | ❌ |
| Minimal polynomial; composition | `minimal_polynomial`, `compose`, `decompose` | — | ❌ |

## Equations

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Polynomial equations | `solve`, `solveset` | `solve`, `csolve` | ✅ |
| One function of the unknown, Lambert W, several logarithms | `solve` | `solve` | ✅ |
| Trigonometric equations: all solutions | `solveset` (`ImageSet`) | `solve` | ◐ principal solutions in one period |
| Linear systems (symbolic coefficients) | `linsolve`, `solve` | `linsolve` | ✅ |
| Nonlinear systems | `nonlinsolve`, `solve` | `solve` | ❌ (needs resultants or Gröbner bases) |
| Inequalities: one unknown; linear systems | `solve`, `solveset` | `solve` | ✅ unions of intervals; cylindrical description for linear systems |
| Nonlinear inequality systems | partial | partial | ❌ |
| Numeric root finding, any equation | `nsolve` | `fsolve` | ◐ polynomials only |
| Diophantine equations | `diophantine` | `isolve` | ❌ |
| Recurrences | `rsolve` | `rsolve`, `seqsolve` | ❌ |
| Ordinary differential equations | `dsolve`, `classify_ode`, `checkodesol` | `desolve` | ❌ ★ |
| Partial differential equations | `pdsolve` | — | ❌ |

## Calculus

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Derivatives, higher and mixed | `diff` | `diff` | ✅ `diff(f, x, n)`, `diff(f, x, y, …)` |
| Limits, one-sided, at ∞ | `limit` (Gruntz) | `limit` | ✅ GruntzLimit (66/66 of its benchmark; answers by case of a parameter, which neither offers) |
| Series, Taylor, asymptotic, `O` | `series`, `fps`, `O` | `series`, `taylor` | ✅ Taylor, Laurent, Puiseux, at ∞, log terms (no `fps`, no `O` in expressions) |
| Padé approximants | — | `pade` | ❌ |
| Indefinite integrals | `integrate` | `integrate` | ✅ RischIntegrate (175/175 of the Apostol set, more than either) |
| Definite and improper integrals | `integrate(f, (x, a, b))` | `integrate(f, x, a, b)` | ◐ by evaluating an antiderivative; no improper integrals or check for singularities |
| Numeric integration | `Integral(...).evalf()` | `romberg`, `gaussquad` | ❌ (QuadGK on a `lambdify`'d function) |
| Sums of series | `summation`, `Sum` | `sum` | ❌ ★ |
| Products | `product`, `Product` | `product` | ❌ |
| Laplace, Fourier, Mellin, Z transforms | `laplace_transform`, `inverse_laplace_transform`, `fourier_transform`, `mellin_transform` | `laplace`, `ilaplace`, `fourier`, `ifourier`, `ztrans`, `invztrans` | ❌ |
| Residues | `residue` | `residue` | ❌ |
| Vector calculus | `sympy.vector` | gradient and related | ◐ `diff` with broadcasting; no named operators |

## Linear algebra

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Matrices of expressions: products | `Matrix` | built in | ✅ Julia arrays of `Sym` |
| Determinant, inverse, solving | `Matrix` methods `det`, `inv`, `LUsolve` | `det`, `inverse`, `linsolve` | ❌ Julia's generic LU pivots by comparing values, which SymJ refuses for expressions |
| Row reduction, null space, rank | `Matrix` methods `rref`, `nullspace`, `rank` | `rref`, `ker`, `rank` | ❌ |
| Eigenvalues, eigenvectors, Jordan form | `Matrix` methods `eigenvals`, `eigenvects`, `jordan_form` | `eigenvals`, `eigenvects`, `jordan` | ❌ |
| LU, QR, SVD | `Matrix` methods `LUdecomposition`, `QRdecomposition` | `lu`, `qr`, `svd` | ❌ |

## Functions

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Elementary functions with exact values and rules | ✅ | ✅ | ✅ (Base's mathematical functions) |
| Special functions | `gamma`, `zeta`, `besselj`, `erf`, `polylog`, `hyper`, `meijerg`, … | many | ✅ through extensions (SpecialFunctions, HypergeometricFunctions, PolyLog); ❌ Meijer G |
| Piecewise, Heaviside, Dirac delta | `Piecewise`, `Heaviside`, `DiracDelta` | `piecewise`, `Heaviside`, `Dirac` | ❌ |
| Combinatorial functions | `factorial`, `binomial`, `fibonacci`, `harmonic`, `bell`, … | `factorial`, `comb`, … | ◐ `factorial`, `binomial` (a polynomial for an integer `k`) |

## Number theory

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Integer factorization, primality, next prime | `factorint`, `isprime`, `nextprime` | `ifactor`, `isprime`, `nextprime` | ❌ (Primes.jl does this in Julia) |
| Totient, modular inverse, CRT, continued fractions | `totient`, `mod_inverse`, `continued_fraction` | `euler`, `iegcd`, `ichinrem`, `dfc` | ❌ |

## Sets, logic, output

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Solution sets: intervals, unions | `Interval`, `FiniteSet`, `ImageSet` | sets | ◐ `SymSet`, `CylindricalSet` from `solve` (not exported) |
| Boolean logic, satisfiability | `sympy.logic`, `satisfiable` | — | ❌ |
| Text output | ✅ | ✅ | ✅ |
| LaTeX output | `latex` | `latex` | ❌ |
| Code generation (C, Fortran, …) | `sympy.codegen` | — | ❌ (`lambdify` covers Julia) |
| Plotting | `plot` | `plot` | ✅ expressions are callable, so Plots.jl and others plot them directly |

## Probably out of scope for SymJ

Areas both systems cover that belong in separate packages rather than a symbolic kernel:
statistics and probability (`sympy.stats`; Giac's `mean`, `stddev`, `normald`), physics
(`sympy.physics`), geometry (`sympy.geometry`, Giac's geometry), combinatorics and
permutation groups (`sympy.combinatorics`), tensors and differential geometry
(`sympy.tensor`, `sympy.diffgeom`), graph theory (Giac), holonomic functions
(`sympy.holonomic`).

## Gaps worth considering, beyond the planned five

The planned work (★) is series, assumptions, a basic `simplify`, summation and
differential equations. Of the rest, in rough order of value for calculus-level use:

1. **Symbolic linear algebra**: determinant, inverse, solving, row reduction and
   eigenvalues for matrices of expressions. Matrices are common in coursework, and
   SymJ's exact linear-system solver (Gauss–Jordan over `RatFun`) is most of the work.
2. **LaTeX output**, for notebooks and documents. Small, and much asked for.
3. **Piecewise functions, Heaviside and Dirac**, which also matter for Laplace
   transforms and differential equations.
4. **Definite integrals done properly**: improper integrals through GruntzLimit, and
   detecting singularities inside the interval, building on RischIntegrate's
   antiderivatives.
5. **Laplace transforms**, natural companions to the differential equation work.
6. **Resultants**, then **Gröbner bases**, which open nonlinear systems of equations.
7. **Real and imaginary parts and `conj`**, once assumptions exist.
8. **Numeric root finding for any equation**, for when `solve` has no closed form.
9. **Utilities**: `collect`, `coeff`, Horner form, a user-level `cse`, `nsimplify`.
10. **Rewriting between forms** (exp ↔ trig, half-angle), which a `simplify` will want
    anyway.
