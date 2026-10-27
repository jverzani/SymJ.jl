# Feature audit: SymJ against SymPy and Giac

What SymPy and Giac provide, area by area, and where SymJ stands. Written 2026-10-23
against SymPy 1.14 and Giac (through Giac.jl 0.14), with SymJ 0.1.0 plus GruntzLimit
(limits) and RischIntegrate (integrals). Every SymPy name below was checked to exist in
SymPy 1.14, and every Giac command in Giac.jl's command table; SymJ's column was checked
against its source and by running it. Updated 2026-10-25: series, assumptions, `limit`
with assumptions, `refine` and `simplify`, complex parts and LaTeX output are done.

**Legend:** ✅ has · ◐ partial · ❌ missing · ★ planned (all done: series, assumptions, a basic
`simplify`, summation, differential equations)

## Numbers and the core

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Exact integers and rationals of any size | `Integer`, `Rational` | built in | ✅ `SymNumber` (Int fast path, BigInt fallback) |
| Exact radicals, π, ℯ, exact function values | `sqrt`, `pi`, `E` | built in | ✅ `sqrt(Sym(8)) == 2sqrt(2)`, `sin(π/3) == sqrt(3)/2` |
| Combining and denesting radicals | `radsimp`, `sqrtdenest` | `simplify` | ✅ through `simplify`: combining (`sqrt(2)*sqrt(3) = sqrt(6)`), denesting square roots of `a + b sqrt(c)`, rationalizing denominators; ◐ no deeper denesting |
| Complex numbers: `re`, `im`, `conj`, `arg` | `re`, `im`, `conjugate`, `arg` | `re`, `im`, `conj`, `arg` | ✅ `real`, `imag`, `conj`, `angle`, `reim` (exact for constants; identities for every complex value as built; parts under assumptions through `refine`); `adjoint` for matrices of expressions |
| Symbolic infinity and undefined | `oo`, `zoo`, `nan` | `infinity`, `undef` | ◐ limits return `Sym(Inf)`; no symbolic infinity in expressions |
| Numeric evaluation | `evalf(n)` | `evalf` | ✅ `N` (Float64), calling with `BigFloat`, dual numbers, … |
| Assumptions on symbols | `symbols(positive=True)`, `ask`, `refine`, `Q` | `assume` | ✅ `Assumptions`, `@assume`, `assuming`/`@assuming`, per call or per block, not stored on symbols; used by `solve`, `limit`, `refine` and `simplify` (`sqrt(x^2) = x` for `x > 0`); ◐ not yet by `integrate`, systems of equations, or integer congruences (`mod(n, 3) == 1`) |
| Parameters (named constants) | — | — | ◐ implemented, hidden until the shared-name question is settled |

## Expressions

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Canonical form as built, fast construction | ✅ | ✅ | ✅ (3–30x faster than SymEngine on the benchmarks) |
| Substitution | `subs` | `subst` | ✅ `subs`, calling with pairs or keywords |
| Compiled numeric functions | `lambdify` | — | ✅ `lambdify` (as fast as hand-written code) |
| Generic term interface | — | — | ✅ TermInterface |
| Pattern matching and rewrite rules | `Wild`, `match`, `replace` | rules | ✅ `SymJ.pattern`/`SymJ.@pattern` with `match`, `eachmatch`, `occursin`: associative-commutative, slots `~x`, `~!x`, `~~x`, `~~~x`, head slots, predicates; `SymJ.@rule`, `SymJ.rule`, `SymJ.rewrite`, `replace` (stages 1–4 of `notes/pattern_matching.md`: rule sets indexed by head; trig rules among `simplify`'s candidates) |
| Collect by powers, coefficients | `collect`, `coeff` | `collect`, `coeff` | ◐ through `Poly` (`coefficients`, `exponent_vectors`) |
| Horner form | `horner` | `horner` | ❌ |
| Common subexpressions | `cse` | — | ◐ `lambdify` shares them; no user-level `cse` |
| Rewriting between forms (exp ↔ trig, …) | `ex.rewrite(...)` | `exp2pow`, `pow2exp`, `halftan`, `sincos`, `tan2sincos` | ◐ inside `simplify` (exponentials ↔ trig and hyperbolic, whichever is smaller); no user-level `rewrite` |
| Floats to exact numbers | `nsimplify` | — | ❌ |

## Simplification ★

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| General simplification | `simplify` | `simplify`, `simplifier` | ✅ ★ `simplify` (the smallest of several transformations: rational, powers, trigonometric, logarithms, radicals, expansion, factoring) and `refine` (identities under assumptions) |
| Trigonometric | `trigsimp`, `exptrigsimp` | `tlin`, `tcollect`, `texpand` | ✅ through `simplify` (`SymJ.simplify_trig`, `SymJ.simplify_hyperbolic`): normal form modulo `sin² + cos² = 1` (`cosh² - sinh² = 1`), multiple and sum angles, exponentials ↔ trig and hyperbolic |
| Powers and logarithms | `powsimp`, `powdenest`, `logcombine` | `lncollect`, `lnexpand` | ✅ through `simplify` and `refine`, under assumptions: powers combined and denested, logarithms combined or expanded (`exp(x)*exp(y) == exp(x + y)` is canonical) |
| Gamma and combinatorial | `gammasimp`, `combsimp` | — | ✅ `SymJ.simplify_gamma` (also inside `simplify`): integer shifts combined, ratios as binomials, reflection and duplication formulas |
| Rational functions | `cancel`, `together`, `apart`, `ratsimp` | `normal`, `ratnormal`, `partfrac`, `propfrac` | ✅ `cancel`, `apart`, `numer_denom` (◐ no separate `together`) |
| Expansion | `expand`, `expand_trig`, `expand_log` | `expand`, `texpand`, `lnexpand` | ✅ `expand` (fast, multinomial); `SymJ.expand_trig` (sums and integer multiples, trigonometric and hyperbolic), `SymJ.expand_log` (under assumptions, or `force = true`) |

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
| Nonlinear systems | `nonlinsolve`, `solve` | `solve` | ◐ polynomial and rational systems, through a lex Gröbner basis; not transcendental ones |
| Inequalities: one unknown; linear systems | `solve`, `solveset` | `solve` | ✅ unions of intervals; cylindrical description for linear systems |
| Solutions under assumptions (domain, range, parameter cases) | `solveset(…, domain)`, assumptions on symbols | `assume` | ✅ `solve(…; assume)`: roots filtered by the unknown's facts; a `Cases` result when they depend on open facts. Not yet for systems |
| Equations in several unrelated functions of the unknown | `solve` | `solve` | ❌ e.g. `sin(x) + cos(x) ~ 1` |
| Nonlinear inequality systems | partial | partial | ❌ |
| Numeric root finding, any equation | `nsolve` | `fsolve` | ◐ polynomials: `SymJ.real_roots` (correctly rounded) and `SymJ.root_intervals` (exact isolation); other equations through Roots.jl or NonlinearSolve.jl, which take SymJ expressions as they are |
| Diophantine equations | `diophantine` | `isolve` | ❌ |
| Recurrences | `rsolve` | `rsolve`, `seqsolve` | ✅ `SymJ.rsolve`: linear, polynomial coefficients, by Hyper (hypergeometric solutions, combined, A and B over algebraic numbers too), right sides that are sums of hypergeometric terms; d'Alembertian solutions by reduction of order (harmonic numbers, iterated sums), any right side; Liouvillian solutions by sections, piecewise in `mod(n, m)` (the double factorial); otherwise the recurrence itself, `DifferenceRoot`, evaluated at integers (Franel); ◐ no arithmetic on those |
| Ordinary differential equations | `dsolve`, `classify_ode`, `checkodesol` | `desolve` | ◐ ★ `SymJ.dsolve`: first order (linear, separable, Bernoulli, exact with integrating factors, homogeneous, Riccati), linear with constant coefficients of any order (undetermined coefficients, variation of parameters), Cauchy–Euler, `y` missing, autonomous second order (`y'' = f(y, y')`); initial conditions (`y'(0) => 1`); linear systems with constant coefficients, of any order; series solutions (`series = n`: Taylor at an ordinary point, nonlinear too; Frobenius at a regular singular point); Laplace transforms for constant coefficients (chosen for `heaviside`/`dirac` right sides, or `method = :laplace`), systems too; every solution checked; several branches of `y'` give a solution each |
| Partial differential equations | `pdsolve` | — | ❌ |

## Calculus

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Derivatives, higher and mixed | `diff` | `diff` | ✅ `diff(f, x, n)`, `diff(f, x, y, …)` |
| Limits, one-sided, at ∞ | `limit` (Gruntz) | `limit` | ✅ `limit`, through GruntzLimit (66/66 of its benchmark; answers by case of a parameter, which neither offers; assumptions settle cases) |
| Series, Taylor, asymptotic, `O` | `series`, `fps`, `O` | `series`, `taylor` | ✅ Taylor, Laurent, Puiseux, at ∞, log terms (no `fps`, no `O` in expressions) |
| Padé approximants | — | `pade` | ❌ |
| Indefinite integrals | `integrate` | `integrate` | ✅ RischIntegrate (175/175 of the Apostol set, more than either) |
| Definite and improper integrals | `integrate(f, (x, a, b))` | `integrate(f, x, a, b)` | ◐ by evaluating an antiderivative; no improper integrals or check for singularities |
| Numeric integration | `Integral(...).evalf()` | `romberg`, `gaussquad` | ❌ (QuadGK on a `lambdify`'d function) |
| Sums of series | `summation`, `Sum` | `sum` | ✅ ★ `SymJ.summation`: Faulhaber, Gosper, partial fractions (harmonic, polygamma, `ζ(2n)`), binomial theorem, hypergeometric series, Zeilberger (definite hypergeometric sums: `Σ binomial(n,k)^2`, Vandermonde; recurrences of order up to 3 by Hyper, else d'Alembertian solutions, `Σ (-1)^(k+1) binomial(n,k)/k = H_n`, else interlaced ones, `Σ (-1)^k binomial(n,k)^3` (Dixon); bounds other than the summand's own); harmonic numbers by parts (`Σ k H_k`), generalized ones `H(n, m)` (`Σ 1/k^m`, `Σ H_k/k = (H_n² + H_n^(2))/2`; not Euler sums such as `Σ H_k/k²`); infinite sums by `limit`, convergence by cases; ◐ unevaluated `Sum` when the recurrence's solutions are not d'Alembertian (`Σ binomial(n,k)^3`) |
| Products | `product`, `Product` | `product` | ❌ |
| Laplace, Fourier, Mellin, Z transforms | `laplace_transform`, `inverse_laplace_transform`, `fourier_transform`, `mellin_transform` | `laplace`, `ilaplace`, `fourier`, `ifourier`, `ztrans`, `invztrans` | ◐ `SymJ.laplace_transform`, `SymJ.inverse_laplace_transform`: polynomials, exponentials, sines and cosines, `heaviside`/`dirac` shifts; inverse of rational functions times `e^(-cs)`. Not yet: Fourier, Mellin, Z |
| Residues | `residue` | `residue` | ❌ |
| Vector calculus | `sympy.vector` | gradient and related | ◐ `diff` with broadcasting; no named operators |

## Linear algebra

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Matrices of expressions: products | `Matrix` | built in | ✅ Julia arrays of `Sym` |
| Matrix functions | `Matrix.exp`, `analytic_func` | `exp` (matrix), `matpow` | ✅ `exp`, `sin`, `cos`, `A^n` (symbolic `n`) through the Jordan form |
| Determinant, inverse, solving | `Matrix` methods `det`, `inv`, `LUsolve` | `det`, `inverse`, `linsolve` | ✅ `det` (Bareiss for polynomials), `inv` (adjugate or Gauss–Jordan over rational functions), `A \ b`, exactly with first-nonzero pivoting |
| Row reduction, null space, rank | `Matrix` methods `rref`, `nullspace`, `rank` | `rref`, `ker`, `rank` | ✅ `SymJ.rref`, `nullspace`, `rank` (generic: symbolic pivots nonzero) |
| Eigenvalues, eigenvectors, Jordan form | `Matrix` methods `eigenvals`, `eigenvects`, `jordan_form` | `eigenvals`, `eigenvects`, `jordan` | ✅ `eigvals` (with multiplicity), `eigvecs`, `eigen` (diagonalizable), `SymJ.charpoly`, `SymJ.jordan_form` |
| LU, QR, SVD | `Matrix` methods `LUdecomposition`, `QRdecomposition` | `lu`, `qr`, `svd` | ✅ `lu` (a `LinearAlgebra.LU`), `qr` (exact Gram–Schmidt), `svd` and `svdvals` (a `LinearAlgebra.SVD`; squared singular values of degree at most 2 over the rationals) |

## Functions

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Elementary functions with exact values and rules | ✅ | ✅ | ✅ (Base's mathematical functions) |
| Special functions | `gamma`, `zeta`, `besselj`, `erf`, `polylog`, `hyper`, `meijerg`, … | many | ✅ through extensions (SpecialFunctions, HypergeometricFunctions, PolyLog); ❌ Meijer G |
| Piecewise, Heaviside, Dirac delta | `Piecewise`, `Heaviside`, `DiracDelta` | `piecewise`, `Heaviside`, `Dirac` | ✅ `SymJ.piecewise`, `SymJ.heaviside`, `SymJ.dirac`: substitution, `diff`, `limit` by side, `refine`, `lambdify`, integration interval by interval (continuous) and of `dirac`; rewriting `abs`, `sign`, `max`, `min`; from `Cases`; ◐ no simplification across pieces (sums of piecewise functions are not combined), breakpoints with several parameters not ordered |
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
| LaTeX output | `latex` | `latex` | ✅ with Latexify loaded: `latexify`, and typeset display of results (notebooks, the documentation) |
| Code generation (C, Fortran, …) | `sympy.codegen` | — | ❌ (`lambdify` covers Julia) |
| Plotting | `plot` | `plot` | ✅ expressions are callable, so Plots.jl and others plot them directly |

## Probably out of scope for SymJ

Areas both systems cover that belong in separate packages rather than a symbolic kernel:
statistics and probability (`sympy.stats`; Giac's `mean`, `stddev`, `normald`), physics
(`sympy.physics`), geometry (`sympy.geometry`, Giac's geometry), combinatorics and
permutation groups (`sympy.combinatorics`), tensors and differential geometry
(`sympy.tensor`, `sympy.diffgeom`), graph theory (Giac), holonomic functions
(`sympy.holonomic`).

## Gaps worth considering

The planned work (★) is done: series, assumptions, a basic `simplify`, summation, and
differential equations (`SymJ.dsolve`: first-order and linear methods, linear systems,
autonomous reduction with its branches, Riccati, series solutions and Laplace
transforms). Of the rest, in rough order of value for calculus-level use:

1. **Symbolic linear algebra**: feature complete for now. `det`, `inv`, `lu`, `\`,
   `SymJ.rref`, `rank`, `nullspace`, eigenvalues and eigenvectors, `qr`, `svd`,
   `SymJ.jordan_form` and matrix functions (`exp`, `sin`, `cos`, `A^n`) are done. The SVD
   refuses when `AᴴA` has an irreducible cubic (or higher) factor, whose roots have no
   usable closed form.
2. **Definite hypergeometric sums and recurrences**: Zeilberger and Hyper are done
   (`SymJ.rsolve`; second- and third-order Zeilberger recurrences solved by their
   hypergeometric solutions, Binet's formula from `Σ binomial(n-k, k)`), with right
   sides that are sums of hypergeometric terms, the inhomogeneous recurrences of
   bounds other than the summand's own (`Σ_{k=0}^{n-1} binomial(n, k)`), and A, B
   factors over algebraic numbers; d'Alembertian solutions by reduction of order
   (`harmonic(n)`, iterated sums, any right side), and Liouvillian ones by sections
   (piecewise in `mod(n, m)`), for `rsolve` and Zeilberger; otherwise `rsolve` returns
   the recurrence itself (`DifferenceRoot`, evaluated at integers). Open: arithmetic on
   those (holonomic closure properties), and returning one from `summation` (which keeps
   `Sum` for `Σ binomial(n, k)^3`).
3. **Piecewise functions, Heaviside and Dirac**: done, with `SymJ.piecewise_fold`
   combining sums and products of them into one (and `simplify` doing so when smaller).
   Breakpoints with several parameters are ordered by the assumptions, and split into
   cases (a `Cases`) when the assumptions do not order them.
4. **Definite integrals done properly**: improper integrals through `limit`, and
   detecting singularities inside the interval, building on RischIntegrate's
   antiderivatives.
5. **Laplace transforms**: done for the elementary transforms and `heaviside`/`dirac`
   shifts (and used by `dsolve`). Open: a table of special-function transforms, and
   Fourier, Mellin and Z transforms.
6. **Numeric root finding**: real roots of polynomials are isolated exactly
   (`SymJ.root_intervals`, `SymJ.real_roots`); for other equations Roots.jl and
   NonlinearSolve.jl work on SymJ expressions directly (callable, or through
   `lambdify`), so no solver of its own is planned.
7. **Utilities**: `collect`, `coeff`, Horner form, a user-level `cse`, `nsimplify`, and a
   user-level `rewrite` (exp ↔ trig, half-angle; `simplify` already rewrites inside).
8. **Assumptions, the remaining pieces**: integer congruences, `assume` for systems of
   equations, and assumptions in `integrate`.

Done since the first audit: trigonometric and logarithmic expansion, gamma and combinatorial simplification, LaTeX output (a Latexify extension), resultants and Gröbner
bases (polynomial systems are solved), real and imaginary parts and `conj`, and
simplification under assumptions (`refine`, `simplify`).
