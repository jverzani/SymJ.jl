# Feature audit: SymJ against SymPy and Giac

What SymPy and Giac provide, area by area, and where SymJ stands. Written 2026-10-23
against SymPy 1.14 and Giac (through Giac.jl 0.14), with SymJ 0.1.0 plus GruntzLimit
(limits) and RischIntegrate (integrals). Every SymPy name below was checked to exist in
SymPy 1.14, and every Giac command in Giac.jl's command table; SymJ's column was checked
against its source and by running it. Updated 2026-10-25: series, assumptions, `limit`
with assumptions, `refine` and `simplify`, complex parts and LaTeX output are done.
Updated 2026-11-01: the feature gaps found in a review of SymPy issues #28909–#30902 are
added (and on 2026-11-02 those of #27909–#28908), each with its SymPy issue numbers. Updated 2026-11-04: those of #26909–#27908 and #25409–#26908 are added; on 2026-11-05 those of #23909–#25408; the review's correctness fixes are tracked
separately. Names checked against SymPy 1.14 and Giac.jl's command table as before.

**Legend:** ✅ has · ◐ partial · ❌ missing · ★ planned (all done: series, assumptions, a basic
`simplify`, summation, differential equations)

## Numbers and the core

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Exact integers and rationals of any size | `Integer`, `Rational` | built in | ✅ `SymNumber` (Int fast path, BigInt fallback) |
| Exact radicals, π, ℯ, exact function values | `sqrt`, `pi`, `E` | built in | ✅ `sqrt(Sym(8)) == 2sqrt(2)`, `sin(π/3) == sqrt(3)/2`; no `τ = 2π` constant (requested, SymPy #26454) |
| Combining and denesting radicals | `radsimp`, `sqrtdenest` | `simplify` | ✅ through `simplify`: combining (`sqrt(2)*sqrt(3) = sqrt(6)`), denesting square roots of `a + b sqrt(c)`, rationalizing denominators; ◐ no deeper denesting |
| Complex numbers: `re`, `im`, `conj`, `arg` | `re`, `im`, `conjugate`, `arg` | `re`, `im`, `conj`, `arg` | ✅ `real`, `imag`, `conj`, `angle`, `reim` (exact for constants; identities for every complex value as built; parts under assumptions through `refine`); `adjoint` for matrices of expressions |
| Symbolic infinity and undefined | `oo`, `zoo`, `nan` | `infinity`, `undef` | ◐ limits return `Sym(Inf)`; no symbolic infinity in expressions; no accumulation bounds (`AccumBounds`; SymPy #28576) |
| Numeric evaluation | `evalf(n)` | `evalf` | ✅ `N` (Float64), calling with `BigFloat`, dual numbers, …; no precision tracking: Float64 cancellation is not detected (SymPy #28537, #28514) |
| Assumptions on symbols | `symbols(positive=True)`, `ask`, `refine`, `Q` | `assume` | ✅ `Assumptions`, `@assume`, `assuming`/`@assuming`, per call or per block, not stored on symbols; used by `solve`, `limit`, `refine` and `simplify` (`sqrt(x^2) = x` for `x > 0`); ◐ not yet by `integrate`, systems of equations, or integer congruences (`mod(n, 3) == 1`); no prime or composite predicate (SymPy #27446, #27441), no reasoning about unknown functions (congruence, SymPy #28154); round 5 (SymPy #23909–#25408): transcendental and algebraic predicates (#25237) |
| Parameters (named constants) | — | — | ◐ implemented, hidden until the shared-name question is settled |

## Expressions

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Canonical form as built, fast construction | ✅ | ✅ | ✅ (3–30x faster than SymEngine on the benchmarks) |
| Substitution | `subs` | `subst` | ✅ `subs`, calling with pairs or keywords |
| Compiled numeric functions | `lambdify` | — | ✅ `lambdify` (as fast as hand-written code); no `constants` keyword (substitute first; SymPy #25491) |
| Generic term interface | — | — | ✅ TermInterface |
| Pattern matching and rewrite rules | `Wild`, `match`, `replace` | rules | ✅ `SymJ.pattern`/`SymJ.@pattern` with `match`, `eachmatch`, `occursin`: associative-commutative, slots `~x`, `~!x`, `~~x`, `~~~x`, head slots, predicates; `SymJ.@rule`, `SymJ.rule`, `SymJ.rewrite`, `replace` (rules, or `old => new` pairs: patterns, predicates, subexpressions; SymPy's `replace`) (stages 1–4 of `notes/pattern_matching.md`: rule sets indexed by head; trig rules among `simplify`'s candidates); round 5 (SymPy #23909–#25408): no `epath` path selection (#25139) |
| Collect by powers, coefficients | `collect`, `coeff` | `collect`, `coeff` | ✅ `SymJ.gather` (one variable or several, any powers, function calls), `SymJ.coeff`; `Poly` for polynomials |
| Horner form | `horner` | `horner` | ✅ `SymJ.horner`, in one variable or recursively in several |
| Common subexpressions | `cse` | — | ✅ `SymJ.cse` (repeated nodes, as `lambdify` shares them; not parts of larger sums) |
| Rewriting between forms (exp ↔ trig, …) | `ex.rewrite(...)` | `exp2pow`, `pow2exp`, `halftan`, `sincos`, `tan2sincos` | ✅ `SymJ.rewrite(ex, target)` (and `fs => target`): through `exp`, `log` (inverse functions), `sin`/`cos`, `sinh`/`cosh`, half-angle `tan`/`tanh`, `gamma`, `factorial`; `simplify` also moves between exponentials and trig, whichever is smaller; ❌ special-function targets: `li` ↔ `Shi`/`Chi`, `erfi` ↔ `expint`, Airy ↔ Bessel (SymPy #30352, #30349); targets `abs` (`sign(x)`; #27323) and `Integral` (Si, Ci, Fresnel, Ei, `expint`; #26144, #26143, #26134); `acos(sin(x))` and the like as periodic piecewise functions (#25720, #25719); round 5 (SymPy #23909–#25408): trig at rational multiples of π into radicals (`rewrite(cos(π/17), sqrt)`, #24560), a `cot` target and hyperbolic ↔ trig (`coth` via `cot(ix)`, #24199), trig to hypergeometric, Bessel, Mathieu (#24176); a form preferring negative exponents in output (#24834) |
| Floats to exact numbers | `nsimplify` | — | ✅ `rationalize(ex)`: the floats in an expression as rationals; `SymJ.nsimplify(ex; constants)`: floats recognized as rationals, algebraic numbers (degree ≤ 4, in radicals) and rational combinations of constants (PSLQ) |
| EML (exp-minus-log) form | — (requested, SymPy #29688) | — | ❌ |
| Noncommutative symbols, indexed objects, sequences | `symbols(commutative=False)`, `IndexedBase`, `SeqFormula` | — | ❌ (SymPy #28001, #28037, #28226); noncommutative: #27891, #27091, #26864, #26194, #26114 (free algebras), #26041, #25572; indexed objects, `KroneckerDelta`, `LeviCivita`: #27544, #27359, #26942, #26109, #25861; round 5 (SymPy #23909–#25408): `IndexedBase` (#24032), sequence objects `RecursiveSeq` (#24191) |

## Simplification ★

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| General simplification | `simplify` | `simplify`, `simplifier` | ✅ ★ `simplify` (the smallest of several transformations: rational, powers, trigonometric, logarithms, radicals, expansion, factoring) and `refine` (identities under assumptions) |
| Trigonometric | `trigsimp`, `exptrigsimp` | `tlin`, `tcollect`, `texpand` | ✅ through `simplify` (`SymJ.simplify_trig`, `SymJ.simplify_hyperbolic`): normal form modulo `sin² + cos² = 1` (`cosh² - sinh² = 1`), multiple and sum angles, exponentials ↔ trig and hyperbolic; ❌ trigonometric factorization (SymPy #28214) |
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
| Factor over extensions, ℤ[i], ℤ/p | `factor(…; extension, gaussian, modulus)` | `cfactor`, factor mod p | ❌ also no algebraic number fields, irreducibility tests or `gcdex`/`invert` over them (SymPy #30845, #30613, #30196, #30158), and no polynomials over ℤ/p (#29715); `AlgebraicNumber` arithmetic (#26743), `factor` with `extension` (#25592); round 5 (SymPy #23909–#25408): exact square roots in ground domains (`is_square`, #25094), polynomial n-th roots (#25087) |
| Resultant, discriminant | `resultant`, `discriminant` | `resultant`, `discriminant` | ✅ modular or subresultant, chosen per input |
| Gröbner bases | `groebner` | `gbasis` | ✅ Buchberger with Gebauer–Möller and sugar; lex through FGLM; `reduced` |
| Roots: exact to degree 4, numeric beyond | `roots`, `nroots` | `solve`, `proot` | ✅ through `solve` (cubics with three real roots in trigonometric form); ❌ no `RootOf`/`CRootOf` objects (SymPy #26740, #25965) |
| Real root isolation and counting | `real_roots`, `count_roots`, `sturm` | `realroot`, `sturm` | ◐ numeric roots only; ❌ counting complex roots in a region (SymPy #26202) |
| Interpolation | `interpolate` | `lagrange` | ❌ |
| Orthogonal and cyclotomic polynomials | `legendre`, `chebyshevt`, `hermite`, `cyclotomic_poly` | `legendre`, `tchebyshev1`, `hermite`, `cyclotomic` | ❌ (Legendre of symbolic degree: SymPy #27298) |
| Minimal polynomial; composition | `minimal_polynomial`, `compose`, `decompose` | — | ✅ `SymJ.minimal_polynomial` (radicals, `im`, trig of rational multiples of π; by resultants), `SymJ.decompose` (complete functional decomposition, characteristic 0; `solve` uses it for degree > 4); composition is `subs` |
| Vieta's formulas | `viete` | — | ◐ `expand` of the product of linear factors; no dedicated helper (SymPy #26237) |

## Equations

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Polynomial equations | `solve`, `solveset` | `solve`, `csolve` | ✅ |
| One function of the unknown, Lambert W, several logarithms | `solve` | `solve` | ✅; ◐ not Lambert W of `x*exp(x^2/2)` (SymPy #26105); ❌ inverting `factorial`, `gamma` (#25493) |
| Trigonometric equations: all solutions | `solveset` (`ImageSet`, `imageset`) | `solve` | ◐ principal solutions in one period; no image sets (SymPy #29675, #29435 and its duplicates #29419, #29418, #29416) |
| Linear systems (symbolic coefficients) | `linsolve`, `solve` | `linsolve` | ✅ |
| Nonlinear systems | `nonlinsolve`, `solve` | `solve` | ◐ polynomial and rational systems, through a lex Gröbner basis; not transcendental ones; transcendental systems such as `2^x + 2^y ~ 10` are an error (SymPy #30849, #29315); trigonometric polynomial systems (SymPy #28804), systems with `log` (#28659, #28064), equations in `real(x)`, `imag(x)` (#28851); with `log` or trig (#27233, #26070); in `z` and `conj(z)` (#27382, #27191, #26301); since round 5, `abs` (sign cases) and square roots (new unknowns) are taken; round 5 (SymPy #23909–#25408): trigonometric systems (#24787), exponential systems (#24060, #24824) |
| Inequalities: one unknown; linear systems | `solve`, `solveset` | `solve` | ✅ unions of intervals; cylindrical description for linear systems; ❌ with parameters (`reduce_inequalities`; SymPy #30688, #30598, #27779, #27589); ❌ transcendental and `abs` inequalities: exponential (#27589), trigonometric (#27573, #25942), logarithmic (#25697), `abs` (#25983, #25738); round 5 (SymPy #23909–#25408): integer unknowns (`x ∈ ℤ`), and assumptions for inequality systems (#24753, #24626) |
| Solutions under assumptions (domain, range, parameter cases) | `solveset(…, domain)`, assumptions on symbols | `assume` | ✅ `solve(…; assume)`: roots filtered by the unknown's facts; a `Cases` result when they depend on open facts. Not yet for systems; round 5 (SymPy #23909–#25408): integer solutions of polynomial systems (#24791) |
| Equations in several unrelated functions of the unknown | `solve` | `solve` | ◐ sin, cos, tan, cot, sec, csc of one argument, through `tan(x/2)` (`sin(x) + cos(x) ~ 1`; SymPy #29745); ❌ different arguments, `sin(4x) + sin(3x) ~ 0` (SymPy #30401), `sin(x)*cos(x/2) ~ 1` (#25954, #25941); exponentials with trig or hyperbolic functions (#26015); several radicals, `cbrt(x^2-1)` and `cbrt(x)` (#25711, #25709) |
| Equations with `floor`, `sign` (interval solutions) | `solve`, `nonlinsolve` | not checked | ◐ `c*floor(u) + d ~ 0` (and `ceil`) gives its solution set, as an inequality does (SymPy #29026, #29007, #28970, #25001); ❌ `max`, `min` and piecewise left sides cannot be inverted (SymPy #28331) |
| Congruences: equations in `mod` | `solveset` (with `Mod`), `solve_congruence` | — | ❌ cannot invert `mod` (SymPy #27562, #26632, #25883); round 5 (SymPy #23909–#25408): #24115 |
| Nonlinear inequality systems | partial | partial | ❌ (sample points of sign conditions, CAD-like: SymPy #26177) |
| Solving for an expression; rearranging | `solve(…, exclude)`, `solve` for `Vo/Vin` | — | ❌ solving for a ratio of unknowns (#24148); fixed-point forms `x = g(x)` (#24081) |
| Numeric root finding, any equation | `nsolve` | `fsolve` | ◐ polynomials: `SymJ.real_roots` (correctly rounded) and `SymJ.root_intervals` (exact isolation); other equations through Roots.jl or NonlinearSolve.jl, which take SymJ expressions as they are (no `nsolve`: SymPy #27799, #25897, #25895, #25447); round 5 (SymPy #23909–#25408): #24727 |
| Diophantine equations | `diophantine` | `isolve` | ❌ (SymPy #29269, #26788, #26783, #26045, #25775); round 5 (SymPy #23909–#25408): #24128 |
| Recurrences | `rsolve` | `rsolve`, `seqsolve` | ✅ `SymJ.rsolve`: linear, polynomial coefficients, by Hyper (hypergeometric solutions, combined, A and B over algebraic numbers too), right sides that are sums of hypergeometric terms; d'Alembertian solutions by reduction of order (harmonic numbers, iterated sums), any right side; Liouvillian solutions by sections, piecewise in `mod(n, m)` (the double factorial); otherwise the recurrence itself, `DifferenceRoot`, evaluated at integers (Franel); ◐ no arithmetic on those |
| Ordinary differential equations | `dsolve`, `classify_ode`, `checkodesol` | `desolve` | ◐ ★ `SymJ.dsolve`: first order (linear, separable, Bernoulli, exact with integrating factors, homogeneous, Riccati), linear with constant coefficients of any order (undetermined coefficients, variation of parameters), Cauchy–Euler, `y` missing, autonomous second order (`y'' = f(y, y')`); initial conditions (`y'(0) => 1`); linear systems with constant coefficients, of any order; series solutions (`series = n`: Taylor at an ordinary point, nonlinear too; Frobenius at a regular singular point); Laplace transforms for constant coefficients (chosen for `heaviside`/`dirac` right sides, or `method = :laplace`), systems too; every solution checked; several branches of `y'` give a solution each; ❌ Lie symmetries (`infinitesimals`; SymPy #30644, #30643) and equations nonlinear in y′, other nonlinear first-order classes (#25773); ❌ power-series (Frobenius, Bessel-type) methods for variable coefficients, `x^2*u + u'' ~ 0` (SymPy #28438), reductions to Bessel's equation (#26275, #25819); ❌ differential-algebraic systems (#28390, #26717); round 5 (SymPy #23909–#25408): Picard iteration (#24119) |
| Partial differential equations | `pdsolve` | — | ❌ (canonical forms of second-order PDEs: SymPy #27980; separation of variables `pde_separate`: #27516, #26697; #25869) |

## Calculus

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Derivatives, higher and mixed | `diff` | `diff` | ✅ `diff(f, x, n)`, `diff(f, x, y, …)`; ❌ of symbolic order `n` (SymPy #29389, #29102); round 5 (SymPy #23909–#25408): with respect to an expression or function value, `diff(exp(x), exp(x))`, `diff(L, Derivative(Y(x), x))` (#25355, #25252); functions of a matrix argument (#24229) |
| Limits, one-sided, at ∞ | `limit` (Gruntz) | `limit` | ✅ `limit`, through GruntzLimit (66/66 of its benchmark; answers by case of a parameter, which neither offers; assumptions settle cases); ❌ multivariable (SymPy #30871); elliptic integrals, modified Bessel functions (#26250, #26040), `floor` of an argument tending to ∞ (#27773); round 5 (SymPy #23909–#25408): multivariable limits (#24225) |
| Series, Taylor, asymptotic, `O` | `series`, `fps`, `O` | `series`, `taylor` | ✅ Taylor, Laurent, Puiseux, at ∞, log terms (no `fps`, SymPy #29027, #26824, #26626, #26039; no coefficient of `x^n` for symbolic `n`, #27723; no `O` in expressions); one variable only: no multivariate (total-degree) Taylor series (SymPy #28583); no series of unknown functions (#26858), real branches of `abs` (#25682); asymptotic series at ±∞ with `e^g`, `cos h`, `sin h` factors (erf, Ei, Si, Ci, Bessel, Stirling; #27109, #26257, #26207, #26210, #26096); round 5 (SymPy #23909–#25408): series reversion (#24251, #24249); integrating expressions with `O` terms (#24477) |
| Padé approximants | — | `pade` | ❌ (requested, SymPy #27665) |
| Indefinite integrals | `integrate` | `integrate` | ✅ RischIntegrate (175/175 of the Apostol set, more than either); ❌ antiderivatives in special functions beyond `erf`/`erfi` and polylogarithms: hypergeometric, Fresnel, elliptic, Owen's T (SymPy #30374, #29910, #29900, #30070, #30014, #30110; RischIntegrate's); round 5 (SymPy #23909–#25408): integrals of unevaluated `Sum`s, termwise (#25405); Fresnel antiderivatives, `cos(x^2)` (#25093); `atan2` integrands (#24004) |
| Definite and improper integrals | `integrate(f, (x, a, b))` | `integrate(f, x, a, b)` | ✅ split at the integrand's singular points, one-sided limits there and at ±∞, `±∞` or an error when divergent, `Cases` on parameters (`∫₀^∞ e^(-ax)`); assumptions pick antiderivative forms through RischIntegrate's sign oracle (needs a RischIntegrate with it); with no antiderivative, Beta integrals over [0, π/2] (∫ sinᵃx cosᵇx; SymPy #30664); ❌ other Beta forms such as ∫₀¹ xᵃ(1-x)ᵇ (#29751), and differentiation under the integral sign (#30869, #30677, #26249); integral representations of special functions, digamma and ζ (#26523, #26512); round 5 (SymPy #23909–#25408): integrals giving Bessel functions, `∫₀^π e^(z cos θ) dθ` (#25194); oscillatory improper integrals, `sinc` Fourier pairs (#25137); symbolic powers `n` over [0, π/2] (#24139) |
| Numeric integration | `Integral(...).evalf()` | `romberg`, `gaussquad` | ❌ (QuadGK on a `lambdify`'d function); no unevaluated `Integral` to evaluate numerically when there is no antiderivative (SymPy #29833) |
| Sums of series | `summation`, `Sum` | `sum` | ✅ ★ `SymJ.summation`: Faulhaber, Gosper, partial fractions (harmonic, polygamma, `ζ(2n)`), binomial theorem, hypergeometric series, Zeilberger (definite hypergeometric sums: `Σ binomial(n,k)^2`, Vandermonde; recurrences of order up to 3 by Hyper, else d'Alembertian solutions, `Σ (-1)^(k+1) binomial(n,k)/k = H_n`, else interlaced ones, `Σ (-1)^k binomial(n,k)^3` (Dixon); bounds other than the summand's own); harmonic numbers by parts (`Σ k H_k`), generalized ones `H(n, m)` (`Σ 1/k^m`, `Σ H_k/k = (H_n² + H_n^(2))/2`; not Euler sums such as `Σ H_k/k²`); infinite sums by `limit`, convergence by cases; ◐ unevaluated `Sum` when the recurrence's solutions are not d'Alembertian (`Σ binomial(n,k)^3`); ❌ bilateral sums, from -∞ (SymPy #27786); no numeric evaluation of an unevaluated `Sum` (#27765); round 5 (SymPy #23909–#25408): products of symbolic `Sum`s not merged into a multiple sum by `expand` (#25078) |
| Products | `product`, `Product` | `product` | ❌ |
| Laplace, Fourier, Mellin, Z transforms | `laplace_transform`, `inverse_laplace_transform`, `fourier_transform`, `mellin_transform` | `laplace`, `ilaplace`, `fourier`, `ifourier`, `ztrans`, `invztrans` | ◐ `SymJ.laplace_transform`, `SymJ.inverse_laplace_transform`: polynomials, exponentials, sines and cosines, `heaviside`/`dirac` shifts; inverse of rational functions times `e^(-cs)`. Not yet: Fourier (SymPy #26053), Mellin, Z (Z transforms: SymPy #29695, not in 1.14); `floor` (staircase) and held transforms (#26999), a convolution-theorem helper (#26023); round 5 (SymPy #23909–#25408): Fourier transforms (#24101, #25137) |
| Fourier series | `fourier_series` | `fourier_an`, `fourier_bn`, `fourier_cn` | ❌ (SymPy #30183, #29336); also #28755, #28724 |
| Residues | `residue` | `residue` | ❌ (SymPy #27156) |
| Vector calculus | `sympy.vector`, `Matrix.jacobian`, `hessian` | `grad`, `divergence`, `curl`, `laplacian`, `hessian` | ◐ `diff` with broadcasting; no named operators, no `jacobian`/`hessian` helpers (SymPy #30186), no coordinate systems (#29107, #27427), no functions of vector arguments (#25644); round 5 (SymPy #23909–#25408): `vector_integrate` over regions, line and flux integrals (#24190, #24123), dyadics and `Del` (#24129) |
| Range and monotonicity of a function | `function_range` (`sympy.calculus.util`), `is_increasing`, `is_strictly_increasing` | `maximize`, `minimize`, `fMax`, `fMin` | ❌ (SymPy #30505, #30499, #30386, #30293, #30582); also #28578, #28577, #28571, #28486, #26871, #26742, #26518, #26076; stationary points (`stationary_points`, #26077), convexity (`is_convex`, #26427); round 5 (SymPy #23909–#25408): `maximum`/`minimum` of an expression over a domain (#25187, #24134) |
| Calculus of variations | `euler_equations` | `euler_lagrange` | ❌ (SymPy #29123) |
| Memoized derivatives and integrals | — | — | ❌ no user-facing cache for `diff`, `integrate` (SymPy #27640) |

## Linear algebra

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Matrices of expressions: products | `Matrix` | built in | ✅ Julia arrays of `Sym` |
| Matrix functions | `Matrix.exp`, `analytic_func` | `exp` (matrix), `matpow` | ✅ `exp`, `sin`, `cos`, `sqrt`, `log`, `A^n` (symbolic `n`) and `A^p` (non-integer `p`) through the Jordan form |
| Determinant, inverse, solving | `Matrix` methods `det`, `inv`, `LUsolve` | `det`, `inverse`, `linsolve` | ✅ `det` (Bareiss for polynomials), `inv` (adjugate or Gauss–Jordan over rational functions), `A \ b`, exactly with first-nonzero pivoting |
| Row reduction, null space, rank | `Matrix` methods `rref`, `nullspace`, `rank` | `rref`, `ker`, `rank` | ✅ `SymJ.rref`, `nullspace`, `rank` (generic: symbolic pivots nonzero); no right-hand-side option reducing only the coefficient columns (SymPy #25661) |
| Linear equations to matrix form | `linear_eq_to_matrix` | — | ❌ (SymPy #25423) |
| Eigenvalues, eigenvectors, Jordan form | `Matrix` methods `eigenvals`, `eigenvects`, `jordan_form` | `eigenvals`, `eigenvects`, `jordan` | ✅ `eigvals` (with multiplicity), `eigvecs`, `eigen` (diagonalizable), `SymJ.charpoly`, `SymJ.jordan_form` |
| LU, QR, SVD | `Matrix` methods `LUdecomposition`, `QRdecomposition` | `lu`, `qr`, `svd` | ✅ `lu` (a `LinearAlgebra.LU`), `qr` (exact Gram–Schmidt), `svd` and `svdvals` (a `LinearAlgebra.SVD`; squared singular values of degree at most 2 over the rationals) |
| Smith normal form | `smith_normal_form`, `smith_normal_decomp` (`sympy.matrices.normalforms`) | `smith`, `ismith` | ❌ (SymPy #29139); Hermite normal form (`hermite_normal_form`) too (SymPy #28601); unimodular transforms (#25759) |
| Pseudoinverse | `Matrix.pinv` | — | ✅ `pinv` from the rank factorization (generic rank, as for `rref`) (SymPy #30759) |
| Schur, Cholesky decompositions (symbolic) | `Matrix.cholesky`, `LDLdecomposition` | `cholesky` | ❌ `schur` of a `Sym` matrix is a MethodError (SymPy #28255); no symbolic Cholesky (#28477) |
| Permanent | `Matrix.per` | `permanent` | ❌ (SymPy #28539) |
| Rotation matrices | `rot_axis1`, `rot_axis2`, `rot_axis3` | — | ❌ (SymPy #24936, #24528) |
| Rank conditions with parameters | — | — | ✅ `rank(A; generic = false)`: `Cases` where it drops (SymPy #28232, #24613) |
| Matrix expressions: `MatrixSymbol`, block matrices, matrix calculus, array expressions | `MatrixSymbol`, `BlockMatrix`, `sympy.tensor.array.expressions` | — | ❌ explicit matrices of `Sym` only (SymPy #30900, #30902, #30860, #30527, #30387, #30278, #30272, #30271, #30270, #30269, #30268, #29837, #29758, #29355, #29292, #28916); round 2: SymPy #28838, #28834, #28823, #28799, #28767, #28766, #28721, #28708, #28641, #28640, #28634, #28633, #28632, #28622, #28620, #28618, #28617, #28615, #28523, #28277, #28264, #28225 (block diagonal); rounds 3–4: SymPy #27854, #27406, #26950, #26946, #26585, #26561, #26527, #26491, #26479, #26288, #26229, #25602, #25554; round 5 (SymPy #23909–#25408): #25118, #24980, #24859, #24824, #24735 (`block_collapse`), #24725, #24717, #24697, #24696, #24689 (`cse` of matrix expressions), #24637, #24507 (matrix equations), #24149, #24131, #24021, #23931 (matrix calculus) |

## Functions

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Elementary functions with exact values and rules | ✅ | ✅ | ✅ (Base's mathematical functions) |
| Special functions | `gamma`, `zeta`, `besselj`, `erf`, `polylog`, `hyper`, `meijerg`, … | many | ✅ through SpecialFunctions, HypergeometricFunctions, PolyLog (dependencies); ❌ Meijer G; ❌ Fresnel integrals (`fresnels`, `fresnelc`), associated Legendre (`assoc_legendre`), `hyperexpand`, `besselsimp`, Owen's T (not in SymPy 1.14) (SymPy #30064, #29645, #30708, #29905, #30345; Bessel naming #29706); also Lerch Φ `lerchphi` and `expand_func` (#27276), spherical harmonics `Ynm`, `Znm` (#26906), `hyperexpand` of Meijer G (#26525), Fresnel (#26306), `Shi`, `Chi` (#26208), generalized `expint(ν, z)` (#25845); round 5 (SymPy #23909–#25408): Meijer G and `hyperexpand` (#25123, #24374), Jacobi elliptic functions `sn`, `cn`, `dn` (#24494), regularized incomplete beta `betainc_regularized` (#24186), Euler numbers and polynomials `euler` (#24156), spherical Bessel functions (#23970) |
| Piecewise, Heaviside, Dirac delta | `Piecewise`, `Heaviside`, `DiracDelta` | `piecewise`, `Heaviside`, `Dirac` | ✅ `SymJ.piecewise`, `SymJ.heaviside`, `SymJ.dirac`: substitution, `diff`, `limit` by side, `refine`, `lambdify`, integration interval by interval (continuous) and of `dirac`; rewriting `abs`, `sign`, `max`, `min`; from `Cases`; ◐ no simplification across pieces (sums of piecewise functions are not combined), breakpoints with several parameters not ordered; ❌ `SingularityFunction` (SymPy #26700) |
| Combinatorial functions | `factorial`, `binomial`, `fibonacci`, `harmonic`, `bell`, … | `factorial`, `comb`, … | ◐ `factorial`, `binomial` (a polynomial for an integer `k`) |

## Number theory

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Integer factorization, primality, next prime, prime counting | `factorint`, `isprime`, `nextprime`, `primepi` | `ifactor`, `isprime`, `nextprime` | ❌ (Primes.jl does this in Julia) (prime counting: SymPy #29125; `factorint`: #27509, #27424; `perfect_power`, squarefree tests: #27238, #27271); round 5 (SymPy #23909–#25408): `factorint`/`ecm` (#25239), `perfect_power` (#25138), `primeomega` (#24165) |
| Totient, modular inverse, CRT, continued fractions | `totient`, `mod_inverse`, `continued_fraction` | `euler`, `iegcd`, `ichinrem`, `dfc` | ❌ (discrete logarithm `discrete_log`: SymPy #27612, #26490); round 5 (SymPy #23909–#25408): continued fractions and convergents of quadratic irrationals (#25060) |

## Sets, logic, output

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Solution sets: intervals, unions | `Interval`, `FiniteSet`, `ImageSet` | sets | ◐ `SymSet`, `CylindricalSet` from `solve` (not exported); no image sets (SymPy #29675, #29435); no interval arithmetic (`SetExpr`, `Interval**oo`; #27743); round 5 (SymPy #23909–#25408): symbolic membership `Contains` (#25207, #25115), Cartesian powers `Reals^n` (#25128) |
| Boolean logic, satisfiability | `sympy.logic`, `satisfiable` | — | ❌ (SymPy #30861, #30729, #30719, #30714); normal forms `to_dnf`/`to_cnf`, `Xor`, `simplify_logic` with don't-cares: #28393, #28369; #27820, #27310, #27309, #27273, #26985 (`to_anf`), #26418, #25873, #25827 (`Xor`), #25451; round 5 (SymPy #23909–#25408): #25218 (`to_anf`), #25115; conjunctions of relations, `simplify(Le(x, y) & Ne(x, y))` (#24195) |
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
(`sympy.holonomic`), linear programming (SymPy's `linprog`, Giac's `lpsolve` and
`simplex_reduce`; SymPy #29382, #29368, #28747, #28104, #28089), quaternions (`Quaternion`; #30560, #28556, #27393, #24943, #24925; octonions #25637), symbolic curve fitting (#24255),
finite fields (polynomials and matrices over GF(p): SymPy #28156, #28121, #27663, #25445), discrete
transforms (`fft`; #28169).

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
4. **Definite integrals**: done (singular points, limits, divergence, cases on
   parameters, and Cauchy principal values with `principal_value = true`). Open: limits
   of special functions the limit code does not know beyond `erf`, `erfc`, `erfi`,
   `sinint`, `cosint` at ±∞.
5. **Laplace transforms**: done for the elementary transforms and `heaviside`/`dirac`
   shifts (and used by `dsolve`). Open: a table of special-function transforms, and
   Fourier, Mellin and Z transforms.
6. **Numeric root finding**: real roots of polynomials are isolated exactly
   (`SymJ.root_intervals`, `SymJ.real_roots`); for other equations Roots.jl and
   NonlinearSolve.jl work on SymJ expressions directly (callable, or through
   `lambdify`), so no solver of its own is planned.
7. **Utilities**: `SymJ.gather` (SymPy's `collect`), `SymJ.coeff` and `rationalize(ex)`
   are done, and `SymJ.horner`, `SymJ.cse`, `SymJ.rewrite(ex, target)`,
   `SymJ.nsimplify` (constant recognition).
8. **Assumptions, the remaining pieces**: integer congruences, `assume` for systems of
   equations, and assumptions in `integrate`. (`@assuming refine = true` refines the
   block's value, since construction does not read the facts.) Integer assumptions could
   also narrow conditions: for an integer `x < 1000`, `x ≥ 999` is `x = 999`, after which
   a piecewise's piece at that point merges with the next one (SymPy #27627; the merge
   itself is done, the narrowing is not).
9. **Asymptotic series with exponential factors** — done (2026-11-04, branch
   asymptotic-series): a `SymSeries` at ±∞ carries pieces `e^g·(series)`,
   `e^g cos(h)·(series)` and `e^g sin(h)·(series)`; `exp`, `sin`, `cos` of growing
   arguments, `erf`, `erfc`, `erfi`, `erfcx`, `expinti`, `expint`, `sinint`, `cosint`,
   the four Bessel functions and Stirling's series for `gamma` (`loggamma` as an ordinary
   series). Not yet: a reciprocal of a sum with such factors (`1/(1 + e^x)`), limits
   through these series (#26040), and finite points (essential singularities stay errors).
   The original note: (decided 2026-11-04, for the next
   round of feature work). At ±∞, `series(erf(x), x, Inf)` is an error for now: past
   its 1, erf's expansion is `e^(-x²)` times a series in 1/x, exponentially small, so as
   a power series in 1/x it would be `1 + O(x^-n)` for every n (SymPy #27109). The plan
   is SymPy's form, `1 - e^(-x²)(1/(√π x) - 1/(2√π x³) + …)`: a `SymSeries` whose
   terms may carry a factor `e^(g(x))`, through addition, multiplication and printing.
   Answering `1 + O(…)` (dropping exponentially small terms) was rejected: correct but
   empty, and a change of convention everywhere. With it come the other asymptotic
   expansions at ∞ the reviews found missing: `erfc`, `Ei`/`expinti`, `Si` and `Ci`
   (#26257, #26207), the modified Bessel functions `besseli`/`besselk` (#26210; also
   their limits at ∞, #26040), and Stirling's series for `gamma` (#26096).
10. **Gaps from the seventh SymPy review (#20909–#22408)**, none a wrong answer:
   - series: `log(log(x))` at 0, a logarithmic leading term inside log (#21227);
     irrational exponents, `exp(w^(1/π))` (#21358); a small exponential inside an
     oscillating argument at ∞, `sin(1/x + e^(-x)) - sin(1/x)` (#21938).
   - `dsolve`: Riccati equations not in the recognized forms, `u' = (x + u)²` by the
     substitution w = x + u (#21503).
   - definite integrals without an elementary antiderivative by contour or series
     methods: `∫₀^π e^(cos x) cos(sin x) dx = π` (#21030).
   - algebraic numbers: `atan` of a radical that is `tan(kπ/n)` not recognized, so
     `cos(25 atan(…))` stays (#21324); denesting square roots of sums of several surds,
     `√(19 + 3√2 + 6√3 - √10)` (#21599).
11. **Gaps from the eighth SymPy review (#19409–#20908)**, none a wrong answer:
   - solve: `angle(u) = r` (a ray, not points; #19506); trigonometric functions of
     different multiples of one angle, `tan(x/6) = cos(x/2)` (#19488).
   - simplify: rational functions over an algebraic extension, a derivative in
     `√(x² + x + 1)` that is `(√(x² + x + 1) - 1)/x` modulo the minimal polynomial (#20236);
     trigonometric quotients collapsed to one multiple angle, `64/sin(2x)⁶` (#20091);
     √(a - √(a²+b²))·√(a + √(a²+b²)) products under assumptions (#20251).
   - definite integrals: contour methods for `∫ log(1 + x²)/(9 + (7 - x)²)` over ℝ (#20171);
     a numeric fallback on finite intervals for non-elementary integrands (round-8 extra);
     incomplete gamma for `∫_y^∞ e^(-x³) dx` (#20685, in the order given); |sec x| over a
     period (its breakpoints are poles).
   - `subs` results are not always in the canonical form of the same expression built
     directly, `-(-√2(x - 1) + x + 1)/(x - 1)` (#20616).
   - Also from the first pass: Diophantine equations, function_range/maximum over a
     domain, periodicity, symbolic derivative orders, factoring over GF(p), Legendre
     polynomials, incomplete/regularized beta and gamma, factorial2, special matrices
     (Vandermonde, Toeplitz), divide-and-conquer recurrences, nonlinear ODE systems,
     matrix expressions, `apart` with an extension or RootSum, noncommutative symbols,
     ImageSet/ProductSet, `collect`, `diff` of equations, integer factorization of
     numbers (`factor(12)`).

Done since the first audit: trigonometric and logarithmic expansion, gamma and combinatorial simplification, LaTeX output (a Latexify extension), resultants and Gröbner
bases (polynomial systems are solved), real and imaginary parts and `conj`, and
simplification under assumptions (`refine`, `simplify`).
