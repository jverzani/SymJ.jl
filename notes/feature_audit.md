# Feature audit and roadmap: SymJ against SymPy and Giac

What SymPy and Giac provide, area by area, where SymJ stands, and what to build next.
Written 2026-10-23 against SymPy 1.14 and Giac (through Giac.jl 0.14), with SymJ 0.1.0
plus GruntzLimit (limits) and RischIntegrate (integrals); every SymPy name was checked to
exist in SymPy 1.14, and every Giac command in Giac.jl's command table. Updated through
2026-11-21 (#206). The SymPy issue review (2026-10-31 – 2026-11-18) added the feature gaps
it found, with their SymPy issue numbers; its record, including the gaps by round, is in
[sympy_review/](sympy_review/README.md).

The [capability tables](#numbers-and-the-core) come first, then the
[roadmap](#roadmap): what is deferred by decision, what is open from the SymPy issue
audit and by area, and what is done.

**Legend:** ✅ has · ◐ partial · ❌ missing · ★ planned in the first audit (all done:
series, assumptions, a basic `simplify`, summation, differential equations)

## Numbers and the core

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Exact integers and rationals of any size | `Integer`, `Rational` | built in | ✅ `SymNumber` (Int fast path, BigInt fallback) |
| Exact radicals, π, ℯ, exact function values | `sqrt`, `pi`, `E` | built in | ✅ `sqrt(Sym(8)) == 2sqrt(2)`, `sin(π/3) == sqrt(3)/2`; no `τ = 2π` constant (requested, SymPy #26454) |
| Combining and denesting radicals | `radsimp`, `sqrtdenest` | `simplify` | ✅ through `simplify`: combining (`sqrt(2)*sqrt(3) = sqrt(6)`), denesting square roots of `a + b sqrt(c)`, rationalizing denominators; ◐ no deeper denesting |
| Complex numbers: `re`, `im`, `conj`, `arg` | `re`, `im`, `conjugate`, `arg` | `re`, `im`, `conj`, `arg` | ✅ `real`, `imag`, `conj`, `angle`, `reim` (exact for constants; identities for every complex value as built; parts under assumptions through `refine`); `adjoint` for matrices of expressions |
| Symbolic infinity and undefined | `oo`, `zoo`, `nan` | `infinity`, `undef` | ◐ limits return `Sym(Inf)`; no symbolic infinity in expressions; no accumulation bounds (`AccumBounds`; SymPy #28576) |
| Numeric evaluation | `evalf(n)` | `evalf` | ✅ `N` (Float64, checked in higher precision, so cancellation is recovered), `float` (Julia's: BigFloat for big numbers), `SymJ.evalf(ex, T, prec)` (a BigFloat to `prec` bits), calling with `BigFloat`, dual numbers, … (SymPy #28537, #28514) |
| Assumptions on symbols | `symbols(positive=True)`, `ask`, `refine`, `Q` | `assume` | ✅ `@assume`, `assuming`/`@assuming` (a `SymJ.Assumptions`), per call or per block, not stored on symbols; used by `solve` (equations and systems), `limit`, `refine`, `simplify` (`sqrt(x^2) = x` for `x > 0`), `integrate`, `rsolve` and `dsolve`; expressions in parameters are generic (nonzero) by default; ◐ not integer congruences (`mod(n, 3) == 1`); no prime or composite predicate (SymPy #27446, #27441), no reasoning about unknown functions (congruence, SymPy #28154); round 5 (SymPy #23909–#25408): transcendental and algebraic predicates (#25237) |
| Parameters (named constants) | — | — | ◐ implemented, hidden until the shared-name question is settled |

## Expressions

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Canonical form as built, fast construction | ✅ | ✅ | ✅ (3–30x faster than SymEngine on the benchmarks) |
| Substitution | `subs` | `subst` | ✅ `subs`, calling with pairs or keywords |
| Compiled numeric functions | `lambdify` | — | ✅ `lambdify` (as fast as hand-written code); no `constants` keyword (substitute first; SymPy #25491) |
| Generic term interface | — | — | ✅ TermInterface |
| Constant test | `is_constant` | — | ✅ `SymJ.isconstant` (no free variables, structural; SymPy's also tests whether the value depends on them) |
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
| Rational functions | `cancel`, `together`, `apart`, `ratsimp` | `normal`, `ratnormal`, `partfrac`, `propfrac` | ✅ `cancel`, `apart`, `SymJ.numer_denom`. `cancel` is `together` followed by cancelling: one fraction in lowest terms, the denominator multiplied out (1/x + 1/(x(x − 1)) is 1/(x − 1)). `SymJ.together` is SymPy's form: the lcm of the denominators kept factored, nothing cancelled, the numerator expanded |
| Expansion | `expand`, `expand_trig`, `expand_log` | `expand`, `texpand`, `lnexpand` | ✅ `expand` (fast, multinomial); `SymJ.expand_trig` (sums and integer multiples, trigonometric and hyperbolic), `SymJ.expand_log` (under assumptions, or `force = true`) |

## Polynomials

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Polynomial type and arithmetic | `Poly` | built in | ✅ `Poly`, packed monomials |
| Division, gcd, lcm | `div`, `rem`, `gcd`, `lcm` | `quo`, `rem`, `gcd`, `lcm` | ✅ modular gcd, Zippel chosen automatically (faster than FLINT on many benchmarks) |
| Factor over ℚ, any number of variables | `factor`, `factor_list` | `factor` | ✅ Zassenhaus, van Hoeij, Wang's EEZ |
| Squarefree decomposition | `sqf_list` | `sqrfree` | ✅ `SymJ.factor_squarefree`, `SymJ.squarefree_part` |
| Factor over extensions, ℤ[i], ℤ/p | `factor(…; extension, gaussian, modulus)` | `cfactor`, factor mod p | ❌ also no algebraic number fields, irreducibility tests or `gcdex`/`invert` over them (SymPy #30845, #30613, #30196, #30158), and no polynomials over ℤ/p (#29715); `AlgebraicNumber` arithmetic (#26743), `factor` with `extension` (#25592); round 5 (SymPy #23909–#25408): exact square roots in ground domains (`is_square`, #25094), polynomial n-th roots (#25087) |
| Resultant, discriminant | `resultant`, `discriminant` | `resultant`, `discriminant` | ✅ modular or subresultant, chosen per input |
| Gröbner bases | `groebner` | `gbasis` | ✅ Buchberger with Gebauer–Möller and sugar; lex through FGLM; `reduced` |
| Roots: exact to degree 4, numeric beyond | `roots`, `nroots` | `solve`, `proot` | ✅ exactly through `solve` (cubics with three real roots in trigonometric form; `solve` is exact only, and points to `SymJ.fsolve` past degree 4), numerically through `SymJ.fsolve`; ❌ no `RootOf`/`CRootOf` objects (SymPy #26740, #25965) |
| Real root isolation and counting | `real_roots`, `count_roots`, `sturm` | `realroot`, `sturm` | ✅ `SymJ.real_roots` (correctly rounded, with multiplicity), `SymJ.root_intervals` (exact isolating intervals); ❌ counting complex roots in a region (SymPy #26202) |
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
| Numeric root finding, any equation | `nsolve` | `fsolve` | ✅ `SymJ.fsolve` (Maple's): one equation's real roots on the line or an interval, complex ones with `domain = ℂ`; a system's solution by Newton from a start, or all solutions of a polynomial system; polynomials also through `SymJ.real_roots`; Roots.jl and NonlinearSolve.jl take SymJ expressions as they are (SymPy's `nsolve`: #27799, #25897, #25895, #25447, #24727) |
| Diophantine equations | `diophantine` | `isolve` | ❌ (SymPy #29269, #26788, #26783, #26045, #25775); round 5 (SymPy #23909–#25408): #24128 |
| Recurrences | `rsolve` | `rsolve`, `seqsolve` | ✅ `SymJ.rsolve`: linear, polynomial coefficients, by Hyper (hypergeometric solutions, combined, A and B over algebraic numbers too), right sides that are sums of hypergeometric terms; d'Alembertian solutions by reduction of order (harmonic numbers, iterated sums), any right side; Liouvillian solutions by sections, piecewise in `mod(n, m)` (the double factorial); otherwise the recurrence itself, `DifferenceRoot`, evaluated at integers (Franel); ◐ no arithmetic on those |
| Ordinary differential equations | `dsolve`, `classify_ode`, `checkodesol` | `desolve` | ◐ ★ `SymJ.dsolve`: first order (linear, separable, Bernoulli, exact with integrating factors, homogeneous, Riccati), linear with constant coefficients of any order (undetermined coefficients, variation of parameters), Cauchy–Euler, second-order equations reducing to Bessel's or Airy's equation through the normal form (`I = A xᵐ + B/x²`, `A e^(λx) + B`, linear; half-integer orders elementary; Riccati equations through them; SymPy #26275, #25819), `y` missing, autonomous second order (`y'' = f(y, y')`); initial conditions (`y'(0) => 1`); linear systems with constant coefficients, of any order; series solutions (`series = n`: Taylor at an ordinary point, nonlinear too; Frobenius at a regular singular point); Laplace transforms for constant coefficients (chosen for `heaviside`/`dirac` right sides, or `method = :laplace`), systems too; every solution checked; several branches of `y'` give a solution each; Lie point symmetries for first order (`SymJ.infinitesimals`: SymPy's heuristics abaco1_simple, abaco1_product, bivariate, chi, function_sum, abaco2_similar, abaco2_unique_unknown, abaco2_unique_general, each also on the equation for x(y); dsolve uses them last, through the integrating factor 1/(η - fξ), else canonical coordinates, s - ∫ F(r) dr = C, an unevaluated definite integral when ∫ F has no closed form (`-1/y - ∫₁^(y-x) dr/(r² - sin r) = C1`); SymPy #30644, #30643); point symmetries of second-order equations (polynomial `ξ`, `η` of degree ≤ 2; the order reduced in canonical coordinates: `x²y y'' = (y - xy')²`, `y'' + y'/x + eʸ = 0`, `x⁴y'' + (xy' - y)³ = 0`) ✅; ❌ non-polynomial infinitesimals of second-order equations, symmetries of third and higher order, a second symmetry used to integrate the reduced equation (when RischIntegrate leaves `∫ dv/√(v² + 2v + 2C1)` with a constant under the radical, as for `y'' = y'/x + y³/x⁴`), and equations nonlinear in y′, Abel equations of the first and second kind with a constant invariant (✅; ❌ other Abel classes: #16387's equation, `x⁷y′ + 5x³y² − (2x² + 2)y³ = 0`, has a non-constant invariant); generalized homogeneous equations (a scaling symmetry, `u = y/xᵏ`; Chini's `y′ = y⁴ + x^(-4/3)`) and Riccati or Abel equations in `x(y)` (`y′ = 1/(x² + y)`) ✅; ❌ other nonlinear first-order classes (#25773, `x y′ + y log y = eˣ`); `x^2*u + u'' ~ 0` (SymPy #28438) now in Bessel functions; hypergeometric equations (Kummer's, Gauss's, Hermite's type, Whittaker's normal form) in `_₁F₁`/`_₂F₁`, or a polynomial solution (Legendre, Laguerre, Hermite, Chebyshev) and reduction of order; ❌ the second solution when the two hypergeometric ones coincide (`b` or `c` an integer, no polynomial solution: Legendre's equation for ν = 1/3) needs new functions, Kummer's `U` and Legendre's `P_ν`, `Q_ν` (possible addition, unexported: `SymJ.hypergeometricU`, `SymJ.legendreP`, `SymJ.legendreQ`; user, 2026-11-23); ❌ Heun, Mathieu and other equations beyond ₂F₁; differential-algebraic systems, linear with constant coefficients (unknowns without derivatives; Luenberger's shuffle, the constraints fixing constants; more initial conditions than constants when consistent; #28390, #26717, #9251); ❌ nonlinear or variable-coefficient differential-algebraic systems; round 5 (SymPy #23909–#25408): Picard iteration (#24119) |
| Partial differential equations | `pdsolve` | — | ❌ (canonical forms of second-order PDEs: SymPy #27980; separation of variables `pde_separate`: #27516, #26697; #25869) |

## Calculus

| Capability | SymPy | Giac | SymJ |
|:--|:--|:--|:--|
| Derivatives, higher and mixed | `diff` | `diff` | ✅ `diff(f, x, n)`, `diff(f, x, y, …)`; ❌ of symbolic order `n` (SymPy #29389, #29102); round 5 (SymPy #23909–#25408): with respect to an expression or function value, `diff(exp(x), exp(x))`, `diff(L, Derivative(Y(x), x))` (#25355, #25252); functions of a matrix argument (#24229) |
| Limits, one-sided, at ∞ | `limit` (Gruntz) | `limit` | ✅ `limit`, through GruntzLimit (66/66 of its benchmark; answers by case of a parameter, which neither offers; assumptions settle cases); ❌ multivariable (SymPy #30871); elliptic integrals, modified Bessel functions (#26250, #26040), `floor` of an argument tending to ∞ (#27773); round 5 (SymPy #23909–#25408): multivariable limits (#24225) |
| Series, Taylor, asymptotic, `O` | `series`, `fps`, `O` | `series`, `taylor` | ✅ Taylor, Laurent, Puiseux, at ∞, log terms (no `fps`, SymPy #29027, #26824, #26626, #26039; no coefficient of `x^n` for symbolic `n`, #27723; no `O` in expressions); one variable only: no multivariate (total-degree) Taylor series (SymPy #28583); no series of unknown functions (#26858), real branches of `abs` (#25682); asymptotic series at ±∞ with `e^g`, `cos h`, `sin h` factors (erf, Ei, Si, Ci, Bessel, Stirling; #27109, #26257, #26207, #26210, #26096); round 5 (SymPy #23909–#25408): series reversion (#24251, #24249); integrating expressions with `O` terms (#24477) |
| Padé approximants | — | `pade` | ❌ (requested, SymPy #27665) |
| Indefinite integrals | `integrate` | `integrate` | ✅ RischIntegrate (175/175 of the Apostol set, more than either), over ℝ: the variable and the parameters real (SymPy, Mathematica and Rubi work over ℂ), so an answer whose real form depends on a parameter's sign comes by cases (`∫ 1/(x² + a)`), cases equal up to expansion merged and cases where the integrand is nowhere real dropped, and `log|u|` for a log of an expression that may be negative (`dsolve` keeps `log u`); Fresnel antiderivatives (`SymJ.fresnels`, `SymJ.fresnelc`; SymPy #25093); ❌ antiderivatives in special functions beyond `erf`/`erfi` and polylogarithms: hypergeometric, Fresnel, elliptic, Owen's T (SymPy #30374, #29910, #29900, #30070, #30014, #30110; RischIntegrate's); round 5 (SymPy #23909–#25408): integrals of unevaluated `Sum`s, termwise (#25405); `atan2` integrands (#24004) |
| Definite and improper integrals | `integrate(f, (x, a, b))` | `integrate(f, x, a, b)` | ✅ split at the integrand's singular points, one-sided limits there and at ±∞, `±∞` or an error when divergent, `Cases` on parameters (`∫₀^∞ e^(-ax)`); assumptions pick antiderivative forms through RischIntegrate's sign oracle (needs a RischIntegrate with it); with no antiderivative, Beta integrals over [0, π/2] (∫ sinᵃx cosᵇx; SymPy #30664), and differentiation under the integral sign with a constant of the integrand as the parameter, checked against quadrature (#30869, #30677: ∫₀^π log((1 + cos x/2)/(1 - cos x/2))/cos x = π²/3); ❌ other Beta forms such as ∫₀¹ xᵃ(1-x)ᵇ (#29751), a parameter inserted where the integrand has none (`∫₀^π e^(cos x) cos(sin x)`, #21030), and the Leibniz rule for integrals of unknown functions (#26249); `∫₀^(π/2) log(1 + cos x/2)/cos x = 5π²/72` since its I'(a) is refined under the assumptions on a (atan(√((1-a)/(1+a))) as acos(a)/2); integral representations of special functions, digamma and ζ (#26523, #26512); round 5 (SymPy #23909–#25408): integrals giving Bessel functions, `∫₀^π e^(z cos θ) dθ` (#25194); oscillatory improper integrals, `sinc` Fourier pairs (#25137); symbolic powers `n` over [0, π/2] (#24139) |
| Numeric integration | `Integral(...).evalf()` | `romberg`, `gaussquad` | ✅ a definite integral with no elementary antiderivative stays an unevaluated `SymJ.Integral`, whose `N` is by quadrature (QuadGK; SymPy #29833) |
| Sums of series | `summation`, `Sum` | `sum` | ✅ ★ `SymJ.summation`: Faulhaber, Gosper, partial fractions (harmonic, polygamma, `ζ(2n)`), binomial theorem, hypergeometric series, Zeilberger (definite hypergeometric sums: `Σ binomial(n,k)^2`, Vandermonde; recurrences of order up to 3 by Hyper, else d'Alembertian solutions, `Σ (-1)^(k+1) binomial(n,k)/k = H_n`, else interlaced ones, `Σ (-1)^k binomial(n,k)^3` (Dixon); bounds other than the summand's own); harmonic numbers by parts (`Σ k H_k`), generalized ones `H(n, m)` (`Σ 1/k^m`, `Σ H_k/k = (H_n² + H_n^(2))/2`; not Euler sums such as `Σ H_k/k²`); infinite sums by `limit`, convergence by cases; ◐ unevaluated `Sum` when the recurrence's solutions are not d'Alembertian (`Σ binomial(n,k)^3`); ❌ bilateral sums, from -∞ (SymPy #27786); no numeric evaluation of an unevaluated `Sum` (#27765); round 5 (SymPy #23909–#25408): products of symbolic `Sum`s not merged into a multiple sum by `expand` (#25078) |
| Products | `product`, `Product` | `product` | ❌ |
| Laplace, Fourier, Mellin, Z transforms | `laplace_transform`, `inverse_laplace_transform`, `fourier_transform`, `mellin_transform` | `laplace`, `ilaplace`, `fourier`, `ifourier`, `ztrans`, `invztrans` | ◐ `SymJ.laplace_transform`, `SymJ.inverse_laplace_transform`: polynomials, exponentials, sines and cosines, `heaviside`/`dirac` shifts; inverse of rational functions times `e^(-cs)`. Not yet: Fourier (SymPy #26053), Mellin, Z (Z transforms: SymPy #29695, not in 1.14); `floor` (staircase) and held transforms (#26999), a convolution-theorem helper (#26023); round 5 (SymPy #23909–#25408): Fourier transforms (#24101, #25137) |
| Fourier series | `fourier_series` | `fourier_an`, `fourier_bn`, `fourier_cn` | ❌ (SymPy #30183, #29336); also #28755, #28724 |
| Residues | `residue` | `residue` | ✅ `SymJ.residue(f, x, x0)`, from the Laurent series; also at ∞ (SymPy #27156) |
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
| Linear equations to matrix form | `linear_eq_to_matrix` | — | ✅ `SymJ.linear_eq_to_matrix(eqs, vars)` gives `(A, b)` (SymPy #25423) |
| Eigenvalues, eigenvectors, Jordan form | `Matrix` methods `eigenvals`, `eigenvects`, `jordan_form` | `eigenvals`, `eigenvects`, `jordan` | ✅ `eigvals` (with multiplicity), `eigvecs`, `eigen` (diagonalizable), `SymJ.charpoly`, `SymJ.jordan_form` |
| LU, QR, SVD | `Matrix` methods `LUdecomposition`, `QRdecomposition` | `lu`, `qr`, `svd` | ✅ `lu` (a `LinearAlgebra.LU`), `qr` (exact Gram–Schmidt), `svd` and `svdvals` (a `LinearAlgebra.SVD`; squared singular values of degree at most 2 over the rationals) |
| Smith normal form | `smith_normal_form`, `smith_normal_decomp` (`sympy.matrices.normalforms`) | `smith`, `ismith` | ✅ `SymJ.smith_normal_form`, `SymJ.smith_normal_decomp` over ℤ and ℚ[x]; Hermite normal form `SymJ.hermite_normal_form`, `SymJ.hermite_normal_decomp` (row form) (SymPy #29139, #2860) |
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
| Special functions | `gamma`, `zeta`, `besselj`, `erf`, `polylog`, `hyper`, `meijerg`, … | many | ✅ through SpecialFunctions, HypergeometricFunctions, PolyLog (dependencies); ✅ Fresnel integrals `SymJ.fresnels`, `SymJ.fresnelc`; ❌ Meijer G; ❌ associated Legendre (`assoc_legendre`), `hyperexpand`, `besselsimp`, Owen's T (not in SymPy 1.14) (SymPy #30064, #29645, #30708, #29905, #30345; Bessel naming #29706); also Lerch Φ `lerchphi` and `expand_func` (#27276), spherical harmonics `Ynm`, `Znm` (#26906), `hyperexpand` of Meijer G (#26525), Fresnel (#26306), `Shi`, `Chi` (#26208), generalized `expint(ν, z)` (#25845); round 5 (SymPy #23909–#25408): Meijer G and `hyperexpand` (#25123, #24374), Jacobi elliptic functions `sn`, `cn`, `dn` (#24494), regularized incomplete beta `betainc_regularized` (#24186), Euler numbers and polynomials `euler` (#24156), spherical Bessel functions (#23970) |
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

## Roadmap

### Deferred by decision

| Item | Decided | Why |
|:--|:--|:--|
| Limits in SymJ itself (a Gruntz port; prototype `SymJ._glimit`, 62/66 of GruntzLimit's benchmark, as fast as cached GruntzLimit) | hold off (2026-11-19) | GruntzLimit works; a port needs symbolic exponents in the series code (`_PS`, ~95 sites) |
| Series with symbolic exponents, step 2 (`x^a + x^2`, by cases on `a`) | not started | step 1 (irrational constant exponents, a factor `x^a`) is done (#190, #195); step 2 waits on the port above |
| Parameters (`SymJ.@parameters`, named constants kept by calling) | hidden | until the question of a name shared with a variable is settled |
| The sign or content of a sum inside a product (`-(x + y)/z` against `(x + y)/(-z)`; SymPy #4596) | left (2026-11-16) | normalizing it changes how many results print |
| Rational exponents printed as `1//3`, so output pastes back (SymPy #5458) | left (2026-11-16) | tempting later, for round trips at the command line |
| `expand` idempotent where radicals merge into a sum (SymPy #5743) | dropped (2026-11-16) | every fix tried made checking an ODE system's solution 5× slower |
| An inference barrier at `show` (the precompile "exceeding 2501 frames" line) | not in (2026-11-15) | it makes Julia 1.14-beta1's precompile of SymJ run past 15 minutes |
| Julia LTS support | Julia ≥ 1.12 (2026-11-21) | RischIntegrate needs Julia 1.12 |

### Open from the SymPy issue audit

The review's [fixes table](sympy_review/README.md#fixes-by-cause) had 126 SymJ fixes with
no status, or in progress. Rechecked against main on 2026-11-21 (with RischIntegrate 1.5.0
and GruntzLimit 1.9.0): 53 were fixed by later PRs, 29 are feature gaps (now in
[Open, by area](#open-by-area)), 5 are settled by decisions (generic parameter values,
structural `==`), and 39 are still open. Each row's note in the table says what main
does now. All of the RischIntegrate and GruntzLimit fixes are closed.

Still open, by area (SymPy issue numbers):

- **Wrong answers.** `∫₀¹ 2πi e^(6πit)/(e^(2πit) - 1/2) dt` gives iπ/4, not iπ/2 (the
  integrand's jumps are not located; #13609). `x, y ∈ ℤ` is ignored by a linear system
  that has integer solutions: no Diophantine parametrization (#14170; one with none
  gives `[]`).
- **A regression.** `∫∫ sin(x² + y²)` over the disc raises a `MethodError` in
  `_fresnel_num` (`src/fresnel.jl`); it was a clean `RischNotImplemented` (#12645).
- **Definite integrals.** Parametric case trees whose branches are errors
  (`∫₀^π 1/(a²cos²x + b²sin²x)` with a, b > 0, #14394; #19046 without assumptions);
  conditions left as `(-w⁴/π⁴)^(1/4) ≶ 0` (r18-extra-4); an inner integral's cases that
  split the outer interval (#7164); no fallback for `∫₀¹ (1 - x⁴)^(-1/3)`, `∫₀^∞ x/cosh x`
  (#14154, #14043), #20171, #20043, `∫₀^∞ cos(πx/2)cos x/(1 - x²)` (and `N` fails at its
  removable pole x = 1; #18796); `∫₀^b 1/(7 + y - y⁵)` known only numerically (#5147);
  transcendental constants in RischIntegrate (#21066; a float base `N(ℯ)^(-2πix)`,
  r18-extra-6).
- **Solving.** `(7^x)^x = 2` (#17452); `angle(x + i) = 0` (#19506); `r³ = V` over ℝ with
  V < 0 left in subcases (#11379); a symmetric five-equation system whose hidden zero is
  not recognized (#19038); solutions not tidied (±√(5 ± 2√6) for ±√2 ± √3, #14244; the
  circle–line intersection, #23889).
- **Differential equations and recurrences.** `f′² = f³` and `f′ = f + b g, g′ = b f - g`
  fail their own check (#15889, #18562); `(f g)″ = 0` left as nested integrals of `g′/g`
  (#21791); `rsolve` of a cubic characteristic polynomial fails its check (#17990), and
  with symbolic coefficients (`cos θ`) needs a symbolic gcd (#18751).
- **Polynomials and linear algebra.** gcd with algebraic coefficients (`gcd(x² + 1, x + i)`,
  `gcd(4√5/7, 2√5/3)`; #18391, #14364); `cancel` with float coefficients (#13805); a
  Gröbner basis of two general cubics with 14 parameters (> 400 s; #12071);
  `jordan_form(k*[0 -c b; c 0 -a; -b a 0])` raises `DivisionByZero` (#9204).
- **Simplification and functions.** `acos(D) + 2atan((D - 1)/√(1 - D²))` for 0 < D < 1
  (#15940); `√(xy)/(xy)` (#21571); `lambertw(-π/2) = iπ/2` (#5080); `Inf + i` and
  complex-direction limits (#17701); `factorial(n) ≥ n` undecided (#13616).
- **Sums.** `Σ binomial(n, m) sin(mx) cos((n - m)x)` (#14412); `Σ 1/(√m + √(m - 1))`
  needs its denominator rationalized to telescope (#14195).
- **Other.** `piecewise` conditions with `x ∈ ℤ` (#17112); the pattern
  `(~a1*a + ~b1)/(~c1*a + ~d1) + (~a2*a + ~b2)/(~c2*a + ~d2)` (#16774); `latexify` of a
  `Poly` (#4799).

### Open, by area

In rough order of value for calculus-level use; each item's SymPy issues are in the
tables below or in the [gaps by round](sympy_review/README.md#feature-gaps-by-round).

1. **Integration.** Antiderivatives in special functions beyond `erf`, `erfi`, Fresnel
   and polylogarithms: hypergeometric, elliptic, Owen's T (RischIntegrate's; SymPy
   #30374, #29910). Definite integrals with no elementary antiderivative by contour or
   series methods (`∫₀^π e^(cos x) cos(sin x) = π`, #21030; #20171; #18796), Beta forms
   such as `∫₀¹ xᵃ(1-x)ᵇ` (#29751, #14549), polar coordinates for double integrals over a
   disc (#12645), Bessel-K results (#9461). Large parametric case trees take seconds
   (#14394). Limits of special functions the limit code does not know beyond `erf`,
   `erfc`, `erfi`, `sinint`, `cosint` at ±∞.
2. **Solving.** All solutions of trigonometric equations as image sets (#29675, #29435);
   trigonometric functions of different multiples of the angle (#30401, #19488);
   transcendental systems (`2^x + 2^y ~ 10`, #30849); inequalities with parameters
   (#30688, #30598) and transcendental or `abs` inequalities (#27589, #27573, #25983);
   mixed polynomial–trigonometric equations (`x = sin x`, #15989); `x^(yz) = x` (#10864); congruences and
   Diophantine equations (#27562, #29269); the solvable quintics (Galois group in F₂₀,
   `x^5 - 5x + 12`), which have roots in radicals (a resolvent-sextic test and the
   Dummit/Spearman–Williams formulas). Integer facts that narrow conditions (for an
   integer `x < 1000`, `x ≥ 999` is `x = 999`; #27627) and congruences as facts
   (`mod(n, 3) == 1`).
3. **Differential equations.** Bessel-type and other special-function equations
   (`x^2*u + u'' ~ 0`, #28438; #17397; `x y'' = 2y`, #7224); `f′² = g(f)` split into
   `f′ = ±√g` where the split's check fails (`f′² = f³`, #15889); Abel equations (#16387; constant invariants done);
   Lie symmetries (#30644; second order by reduction done); differential-algebraic systems (#28390, #9251; linear with constant coefficients done).
4. **Series and limits.** Formal power series (`fps`, #29027); multivariate Taylor
   series (#28583); reciprocals of sums with exponential factors (`1/(1 + e^x)`) and
   limits through asymptotic series (#26040); an exponentially small term inside a
   function's argument at ∞ (`sin(1/x + e^(-x))`, #21938); `log(log x)` at 0 (#21227);
   multivariable limits (#30871); limits of unevaluated sums (#14196).
5. **Simplification.** Trigonometric factorization (#28214); rational functions over an
   algebraic extension (#20236); denesting square roots of several surds (#21599);
   denominators with higher roots rationalized (#14184); `√(s - a)√(s + a)` with
   `s = √(a² + b²)` (#20251); a complex constant as a + bi (#11193); `asinh` of a
   constant as a logarithm (#10429); real and imaginary parts recombined (#10822);
   contiguous relations of ₂F₁ (#17702); `min`/`max` shifts (#18643).
6. **Sums and recurrences.** Arithmetic on the recurrence objects `rsolve` returns
   (holonomic closure), and `summation` returning one (it keeps `Sum` for
   `Σ binomial(n, k)^3`); Euler sums (`Σ H_k/k²`); power sums with a symbolic exponent
   (`Σ m^k`, Faulhaber; #14197); closed forms in incomplete gamma (#11643); telescoping
   detection (#14195); finite trigonometric sums (#14412); bilateral sums (#27786); numeric values of an
   unevaluated `Sum` (#27765); symbolic products ∏ (#3410).
7. **Polynomials and algebraic numbers.** Factoring over extensions, ℤ[i] and ℤ/p
   (#30845, #29715, #3298); gcd and cancel over ℚ(i) and algebraic extensions (#18391,
   #15358); `RootOf`/`CRootOf` objects (#26740); interpolation; orthogonal and cyclotomic
   polynomials; counting complex roots in a region (#26202).
8. **Linear algebra.** Symbolic Cholesky and Schur (#28477, #28255); SVD and matrix
   functions when the eigenvalues are roots of an irreducible cubic (#10355, #15362);
   permanents; rotation matrices; large parametric linear systems, whose Gauss–Jordan
   elimination spends about two minutes in multivariate gcds (#14694, 18 equations).
9. **Transforms and vector calculus.** Fourier, Mellin and Z transforms (#26053, #29695)
   and Fourier series (#30183); a table of special-function Laplace transforms; `jacobian`,
   `hessian` and the vector operators (#30186); ranges, extrema and monotonicity of a
   function (`function_range`, #30505).
10. **Functions.** Meijer G and `hyperexpand`; associated Legendre functions and spherical
    harmonics; Lerch Φ; Jacobi elliptic functions; combinatorial numbers (`fibonacci`,
    `bell`, Euler numbers).
11. **Objects SymJ does not have**, the largest gaps by issue count (some may belong in
    separate packages): symbolic sets (`ImageSet`, `Range`, `ProductSet`, symbolic
    intervals); matrix expressions (`MatrixSymbol`, block matrices, matrix calculus;
    about 70 SymPy issues); noncommutative symbols; `IndexedBase` and `KroneckerDelta`;
    Boolean logic; a symbolic infinity in expressions and `AccumBounds`; an unevaluated
    `Derivative` for users; number theory (Primes.jl has much of it).
12. **Performance.** Multivariate factoring of large products (about 50 s for six
    variables); Gröbner bases over a parametric field by FGLM (#12071); the precompile
    workload (about 80 s of inference, spread thinly).

### Done since the first audit

- The planned work (★): series (Taylor, Laurent, Puiseux, at ∞, asymptotic with
  exponential factors, #111), assumptions, `simplify` and `refine`, summation (Gosper,
  Zeilberger, Hyper, d'Alembertian and Liouvillian solutions), `SymJ.dsolve` (first-order
  and linear methods, systems, series solutions, Laplace transforms), `SymJ.rsolve`.
- The SymPy issue review, rounds 1–20 (#85–#182): see
  [sympy_review/](sympy_review/README.md).
- Defaults: generic parameter values (#154) and the default-assumption decisions 3–6
  (#155–#158); indefinite integrals over ℝ by cases on parameter signs, with `log|u|`
  (#203, #204).
- Polynomials and linear algebra: Gröbner bases over parameter fields (#88), matrix
  `sqrt`, `log`, non-integer powers and `pinv` (#98), `SymJ.linear_eq_to_matrix` (#187),
  Smith and Hermite normal forms (#188).
- Equations: `SymJ.fsolve` (#169, #170), `solve` exact only (#171).
- Calculus: unevaluated integrals valued by quadrature (#184), `SymJ.residue` (#186),
  Fresnel integrals (#189), series with non-rational exponents (#190, #195), definite
  integrals by differentiation under the integral sign (#196, #200), parametric
  definite integrals with conditions in the parameters' own terms (#198, #199, #201).
- Expressions: `SymJ.evalf` and `float` as Julia's (#183), `SymJ.together` (#194),
  `SymJ.isconstant` (#206), `==` on floats as IEEE with exact ≠ float kept (#209).
- Earlier: trigonometric and logarithmic expansion, gamma and combinatorial
  simplification, LaTeX output (a Latexify extension), resultants and Gröbner bases
  (polynomial systems are solved), real and imaginary parts and `conj`, simplification
  under assumptions.

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
