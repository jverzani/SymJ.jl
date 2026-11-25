# SymPy issue review (2026-10-31 – 2026-11-18)

A record of the review of SymPy's GitHub issues against SymJ, GruntzLimit and RischIntegrate. Each issue was classified by whether SymJ shares the problem; the problems it shared were verified, grouped by cause (a *fix*), assigned to the package that should fix it, and fixed in turn. Generated from the private tracker page's data and its status database; [issues.md](issues.md) has every issue reviewed (the page's "All issues" tab). Feature gaps are listed by round [below](#feature-gaps-by-round), and the open ones by area in the [feature audit's roadmap](../feature_audit.md#roadmap).

**Ranges:** SymPy issues #28909–#30902 (round 1, 500 issues), #27909–#28908 (round 2, 319 issues), #26909–#27908 (round 3, 379 issues), #25409–#26908 (round 4, 596 issues) #23909–#25408 (round 5, 632 issues) #22409–#23908 (round 6, 647 issues) #20909–#22408 (round 7, 659 issues) #19409–#20908 (round 8, 686 issues), #17909–#19408 (round 9, 706 issues) #16409–#17908 with #30903–#31016 (round 10, 730 issues) #14909–#16408 (round 11, 648 issues) #13409–#14908 (round 12, 728 issues), #11909–#13408 (round 13, 665 issues), #10409–#11908 (round 14, 741 issues) #8909–#10408 (round 15, 650 issues), #7409–#8908 (round 16, 727 issues) #5909–#7408 (round 17, 1394 issues) #4409–#5908 (round 18, 1499 issues) #2909–#4408 (round 19, 1342 issues) and #945–#2908 with #31017–#31179 (round 20, 197 issues), plus problems found in passing.

**Totals:** 14445 SymPy issues and 206 problems found in passing (14651 rows); 1100 shared problems; 959 fixes by cause; 41 missing features.

| | count |
|:--|--:|
| fixes fixed | 848 |
| fixes open (rechecked 2026-11-21) | 40 |
| fixes won't fix | 42 |
| fixes that are feature gaps | 29 |
| missing case | 626 |
| wrong answer | 196 |
| crash | 89 |
| hang | 48 |
| owner SymJ | 856 |
| owner RischIntegrate | 62 |
| owner GruntzLimit | 40 |
| owner other: CommonEq | 1 |

Rows per round: round 1: 514, round 2: 324, round 3: 387, round 4: 602, round 5: 646, round 6: 657, round 7: 668, round 8: 696, round 9: 717, round 10: 746, round 11: 658, round 12: 737, round 13: 679, round 14: 750, round 15: 661, round 16: 733, round 17: 1406, round 18: 1511, round 19: 1351, round 20: 208.

## Rounds at a glance

*Rows* counts the SymPy issues of the round and the problems found in passing;
*shared with SymJ* the rows whose problem SymJ shared (SymJ? = yes in
[issues.md](issues.md)); *causes* the distinct fixes those rows were grouped under.
Rounds 1–6 recorded their feature gaps in the
[feature audit](../feature_audit.md)'s tables, tagged with the round; rounds 7–20 in
[Feature gaps by round](#feature-gaps-by-round) below. The
[roadmap](../feature_audit.md#roadmap) groups the open ones by area.

| round | SymPy issues | rows | shared with SymJ | causes | gaps |
|--:|:--|--:|--:|--:|:--|
| 1 | #28909–#30902 | 514 | 68 | 52 | in the audit's tables |
| 2 | #27909–#28908 | 324 | 32 | 25 | in the audit's tables |
| 3 | #26909–#27908 | 387 | 59 | 44 | in the audit's tables |
| 4 | #25409–#26908 | 602 | 53 | 45 | in the audit's tables |
| 5 | #23909–#25408 | 646 | 71 | 62 | in the audit's tables |
| 6 | #22409–#23908 | 657 | 62 | 57 | in the audit's tables |
| 7 | #20909–#22408 | 668 | 63 | 52 | [round 7](#round-7) |
| 8 | #19409–#20908 | 696 | 59 | 53 | [round 8](#round-8) |
| 9 | #17909–#19408 | 717 | 49 | 46 | [round 9](#round-9) |
| 10 | #16409–#17908, #30903–#31016 | 746 | 63 | 49 | [round 10](#round-10) |
| 11 | #14909–#16408 | 658 | 48 | 44 | [round 11](#round-11) |
| 12 | #13409–#14908 | 737 | 75 | 60 | [round 12](#round-12) |
| 13 | #11909–#13408 | 679 | 39 | 36 | [round 13](#round-13) |
| 14 | #10409–#11908 | 750 | 55 | 51 | [round 14](#round-14) |
| 15 | #8909–#10408 | 661 | 32 | 31 | [round 15](#round-15) |
| 16 | #7409–#8908 | 733 | 47 | 42 | [round 16](#round-16) |
| 17 | #5909–#7408 | 1406 | 94 | 87 | [round 17](#round-17) |
| 18 | #4409–#5908 | 1511 | 74 | 68 | [round 18](#round-18) |
| 19 | #2909–#4408 | 1351 | 45 | 42 | [round 19](#round-19) |
| 20 | #945–#2908, #31017–#31179 | 208 | 12 | 11 | [round 20](#round-20) |

## Feature gaps by round

Missing features (not wrong answers) each round found, as recorded at the time; some
have been done since, as marked. Moved here from the feature audit, whose roadmap now
collects the open ones by area.

### Round 7

SymPy #20909–#22408. None a wrong answer.

- series: `log(log(x))` at 0, a logarithmic leading term inside log (#21227);
  irrational exponents, `exp(w^(1/π))` (#21358; done since: an exponent is a rational
  plus an irrational constant, ordered numerically, so `x + x^√2` and `1/(1 + x^π)`
  expand; symbolic exponents, `x^a + x^2`, still error); a small exponential inside an
  oscillating argument at ∞, `sin(1/x + e^(-x)) - sin(1/x)` (#21938).
- `dsolve`: Riccati equations not in the recognized forms, `u' = (x + u)²` by the
  substitution w = x + u (#21503).
- definite integrals without an elementary antiderivative by contour or series
  methods: `∫₀^π e^(cos x) cos(sin x) dx = π` (#21030).
- algebraic numbers: `atan` of a radical that is `tan(kπ/n)` not recognized, so
  `cos(25 atan(…))` stays (#21324); denesting square roots of sums of several surds,
  `√(19 + 3√2 + 6√3 - √10)` (#21599).

### Round 8

SymPy #19409–#20908. None a wrong answer.

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

### Round 9

SymPy #17909–#19408. None a wrong answer.

- solve: the 5-equation symmetric system x_i = 1/(1 + Σ others) (a zero hidden in
  nested complex radicals is taken as an identity; #19038); parameter cases e = 1,
  c = 0 of a system with a symbolic exponent (round-9 extra); nested trigonometric
  equations give the outer roots of one period only (sin(sin x) = 0; #18022).
- rsolve: irreducible cubic characteristic polynomials fail its own check (#17990);
  coefficients with cos θ reach the numeric gcd (#18751).
- dsolve: a parametric 2×2 constant-coefficient system fails its check (#18562).
- gcd and cancel over algebraic extensions (x² + 1 and x + i; #18391).
- definite integrals: ∫₀^∞ cos(πx/2)cos x/(1 - x²) (#18796); a rational function in
  parameters over ℝ splits on whether √(-g²) is real, with "known only numerically"
  cases (#19046; g, K > 0 is right).
- simplify: min/max shifts (min(0, n - 1 - c) + min(c, n - 1); #18643); -a/a for a
  sum a with radicals is not cancelled automatically (#18799); merging Or of
  relations under assumptions (n = 0 | n ≥ 1 as n ≥ 0 for integer n; #18389).
- limits of harmonic(n) (GruntzLimit has no rule for it).
- From the first pass: sets (Interval/FiniteSet/ImageSet objects), Diophantine
  equations, modular roots, Fourier transforms and series, numeric root finding
  (2^x + 3^x = 5^x), boolean simplification, Smith normal form, matrix symbols and
  block matrices, KroneckerDelta, multinomial, Bessel ODEs, B-splines, companion
  matrices, noncommutative symbols, AccumBounds, Poly over GF(p), BigFloat inputs
  (Sym(BigFloat) rounds to Float64).

### Round 10

SymPy #16409–#17908, #30903–#31016. None a wrong answer.

- solve: trigonometric systems that do not decouple (sin x sin y = 0, cos x cos y =
  1; #16618); (7^x)^x as 7^(x²) only for real x (#17452); a system of 38 solved
  equations passes the case limit, and its generic leaf is "too many cases" (#17186;
  keeping the generic side past the limit made a round-5 parametric system run on).
- piecewise conditions that are memberships (x ∈ ℤ; #17112).
- dsolve: Bessel and other special-function ODEs (#17397, #17592).
- simplify: contiguous relations of ₂F₁ (#17702); (x + y)² and (-x - y)² have
  different canonical forms (== is structural; #17768).
- infinities in the complex plane (Inf + im, directional limits; #17701); exp(Inf*k)
  with k > 0 under refine.
- patterns: a sum of two linear fractions does not match its pattern (#16774).
- GruntzLimit: binomial with two arguments, li2, impolylog (handed off).
- From the first pass: symbolic sets (ImageSet, ProductSet, Range, intersections,
  symbolic intervals), matrix expressions (MatrixSymbol, Hadamard), orthogonal and
  Euler polynomials, B-splines, Fresnel and Marcum Q functions, RisingFactorial,
  frac, noncommutative symbols, modular and discrete-log equations, PDEs, CRootOf,
  function_range/continuous_domain, change of variables in derivatives (dsubs).

### Round 11

SymPy #14909–#16408.

- dsolve: f′² = g(f) not split into f′ = ±√g (and f′ = f^(3/2)'s 4/(C - x)² is
  rejected by the numeric check on the branch where it fails; #15889); Abel ODEs
  (#16387); conditions at a singular point of the general solution are refused rather
  than taken as limits (Legendre, f(±1); #15113).
- solve: mixed polynomial–trigonometric equations (x = sin x, even over ℝ; #15989).
- simplify: the atan half-angle identity atan(√((1 - D)/(1 + D))) = acos(D)/2 (#15940).
- gcd and cancel over ℚ(i) and the Gaussian integers (#15358); matrix functions with
  irreducible cubic eigenvalues (#15362).
- Σ_{n≥2} xⁿ/n stays unevaluated (the n ≥ 1 sum is −log(1 − x)); Lerch Φ in general.
- nullspace through exp(iu) (#15585) gives its vectors in exponential, not trig, form;
  symbolic eigenvectors with eigenvalues past 3000 nodes are refused (#15076).
- numeric roots are Float64: a root below 1e-308 (x⁵ − 10⁴⁰⁰x + 1) prints as 0.0.
- RischIntegrate and GruntzLimit items handed off (round11/prompts).
- From the first pass: nonlinear and variable-coefficient ODE systems, sympy.vector
  coordinate systems, number theory (discrete log, divisor functions, diophantine),
  Legendre functions, finite-field factoring, differentiating with respect to a
  function, Boolean logic, ConditionSet, nsolve from a point, matrix symbols.

### Round 12

SymPy #13409–#14908.

- large symbolic linear systems: the 18-equation chain in 28 parameters (#14694) spent
  minutes in multivariate gcds during Gauss–Jordan elimination. Done 2026-11-25:
  fraction-free Gauss–Jordan (Bareiss) for square systems of any size, 3 s.
- parametric definite integrals whose antiderivative has the parameters in its
  denominators: ∫₀^π 1/(a²cos²x + b²sin²x) with a, b > 0 (#14394) splits into many
  cases, each taking endpoint limits (a = 2, b = 3 gives π/6 in 3 s).
- finite trigonometric sums Σ_{k=1}^{n} cos(kt), sin(kt) and Σ binomial(n, m) sin(mx)
  cos((n-m)x) (#14137, #14412): the exponential route leaves i in its closed form.
- cancel with float coefficients (#13805): no approximate gcd.
- the Beta integral folded onto [0, 1], ∫₀¹ (x^(p-1) + x^(q-1))/(1 + x)^(p+q) (#14549).
- rationalizing denominators with cube and higher roots (#14184).
- bounds such as n! ≥ n (#13616): refine knows signs only (n! > 0).

### Round 13

SymPy #11909–#13408.

- multivariate factoring of large random products (six variables, degrees 8 and 6;
  round-13 extra) takes about 50 s.
- Gröbner bases of two general bivariate cubics with fourteen parameters (#12071) take
  more than 3 minutes (FGLM over a parametric field).
- double integrals over a disc that need polar coordinates,
  ∫∫ sin(x² + y²) over x² + y² ≤ 1 (#12645).
- From the first pass: trigonometric and parametric inequalities (cos x < 1/2, K > Py;
  #13290, #13116, #12560), transcendental nonlinear systems (#13272), dsolve systems
  with variable coefficients or nonlinear ones (#13376, #13060, #11968), numeric
  evaluation of unevaluated integrals (#13371, #12646), radical forms of cos(rπ)
  (cos(π/17); #13066), rising and falling factorials (#13233), periodicity (#13207),
  function_range (#13273), symbolic derivative orders (#12484), differentiating with
  respect to a function call (#12139), Fourier and Z transforms (#12591, #12502),
  symbolic Cholesky (#12500), N to a number of digits and chop (#12680, #12209),
  gcd of radicals (#12794), integrals of floor and ceiling (#12234), the Bessel-I
  result of ∫₀^∞ e^(-(x-x₀)²)/√x (#12626), Diophantine equations (#11959, #11956,
  #11955), discrete logarithms and modular equations (#13201, #13178), and the objects
  SymJ does not have: symbolic sets (FiniteSet, ProductSet, ImageSet, ComplexRegion,
  symbolic intervals), MatrixSymbol and matrix expressions, IndexedBase and
  KroneckerDelta, noncommutative symbols, holonomic functions, sequences, formal power
  series, B-splines, PDE separation, boolean logic expressions, Q.prime.

### Round 14

SymPy #10409–#11908.

- solve: x^(yz) = x, a symbolic exponent against its base, is refused (#10864); real
  roots of r³ = V with V's sign undecided are not classified (#11379).
- eigenvalues of a symmetric 3×3 in Cardano's complex form, not the trigonometric real
  form (#11514); Aⁿ of a nilpotent matrix with symbolic n gives no cases (#10563).
- simplify: asinh of a constant rewritten as a log (#10429); re/im parts recombined
  into one complex expression (#10822); a complex constant quotient in Cartesian form,
  (1 - 2i + √3 i)/(2 - √3 + i) = -i (#11193).
- cse finds common sums, not common subproducts (#11428).
- Σ of incomplete-gamma (hypergeometric) series in closed form (#11643).
- From the first pass: holonomic functions (#11284–#11329, #11490, #11541), Diophantine
  equations and Pell's equation (#11486, #11712, #11697, #10861 and the diop issues),
  nsolve (#11768, #11860), PDE separation (#11726, #11342), Fourier transforms of
  distributions and Fourier series (#11772, #10912, #11626), integrals of spherical
  Bessel functions (#11891), ₁F₂ and its limits (#10802), more special functions
  (Struve, Kummer U, Jacobi elliptic, Riemann θ, subfactorial, Laguerre; #11579,
  #10881, #10962, #10961), rising and falling factorials (#11105, #10938), symbolic
  derivative orders and differentiation with respect to an expression or a derivative
  (#11762, #11573, #11467, #11737), partial factorization (#11869), rewriting to sums
  and to real/imaginary parts (#11377, #10897), subs on a Poly (#11447), rsolve in
  several variables (#11063), dsolve constant names (#11024), finite differences
  (#11007), factor of an integer (#10596), modular roots and discrete logarithms
  (#11502, #10886, #11512), sum-of-squares decomposition (#11332), inverses of user
  functions for solve (#10487), step-by-step solutions (#11088), and the objects SymJ
  does not have: symbolic sets (Range, ImageSet, ConditionSet, ComplexRegion, product
  sets), MatrixSymbol and matrix expressions, Indexed, noncommutative symbols, boolean
  ITE.

### Round 15

SymPy #8909–#10408.

- solve: inequalities with parameters (1.5x + 1 < a); ceil(x) - floor(x) = 1 or 0
  (all non-integers, all integers: sets SymJ cannot write).
- definite integrals: no non-elementary fallback, ∫₀^∞ e^(-√(x² + y²))/√(x² + y²) dx is
  K₀(y) (#9461; or polar coordinates).
- dsolve: differential-algebraic systems with a singular leading matrix, reducible by
  elimination (#9251).
- matrices (last): Jordan form of k·M with a parametric eigenvalue, so the Lorentz
  system's dsolve fails (#9204); svd/cond with squared singular values roots of an
  irreducible cubic (#10355).
- A Bool from Julia's a ≠ b (structural, true for distinct symbols) in a solve system
  is dropped silently as `true`; Ne(a, b) is the relation.
- Done 2026-11-14: float inequalities, #10398's conditions for real parameters, floor
  with ceil (#9118), impossible inner-integral cases (#10211), Bessel's equation
  (#9819), cse subproducts (#10228).

### Round 16

SymPy #7409–#8908.

- assumptions: linear facts with irrational coefficients ((√2 - 1)b - y < 0 from y > b > 0),
  so a definite integral keeps cases its facts settle (#8556).
- sums: Σ sin(nx) to ∞ (cases on x ∈ 2πℤ, divergent otherwise; #8329); a limit at -∞ of
  xⁿ with n an integer needs parity cases (#8634; GruntzLimit refuses).
- solve: inequalities in an integer unknown (#8221) and transcendental inequalities,
  |x|^π ≥ 0, eˣ > 1 (#8204); the 46-case answer of a system with equations in the
  parameters alone could be ~6 cases (#7621).
- integrate: cases conditioned on the integration variable itself (Σ sin(nx) closed form
  at sin(x/2) = 0, a set of points) integrated case by case (#7827).
- simplify: the Black–Scholes identity, exponential quotients with logs in the exponents
  (#7769).
- design questions: sorting symbols (`sort([y, x])`, #7413: no canonical `<` on symbols);
  error messages that show floats as the exact binary rationals they were converted to
  (#8503).
- From the first pass (feature gaps): Indexed/IndexedBase and MatrixSymbol, symbolic
  FiniteSet/ImageSet/Interval objects, Product, KroneckerDelta, Fourier transforms,
  diophantine equations, nsolve, number-theory functions (multiplicity, totient,
  subfactorial, bell, fibonacci, Catalan, Euler numbers, B-splines, Laguerre),
  boolean logic, noncommutative symbols, modular roots, time-dependent ODE systems,
  Clairaut ODEs, deep factor (inside function arguments).

### Round 17

SymPy #5909–#7408.

- dsolve: Bessel-type equations, x y'' = 2y as √x (C₁ I₁(√(8x)) + C₂ K₁(√(8x))) (#7224).
- series: a symbolic exponent of a vanishing base, log(x)^a at x = 1 as
  (x - 1)^a (1 + …), needs series with symbolic exponents (#7168; done since: such
  powers are a factor x^a of the series); Bessel functions of
  half-integer order as elementary functions (#6154).
- simplify: fractions in a tower of radicals, 1/(s + 1) + 1/(s - 1) + s with
  s = √(2 + √5), reduced over ℚ(√5)(s) (#6651); the zero of a quartic's root in a
  parameter substituted back (#6829).
- solve: equations with float exponents whose differences are integers to rounding,
  (1 + x)^0.9389 and (1 + x)^1.9389, need a numeric solver after exactifying (extra-9,
  #6363); a float residue (1e-116) in a solution where the system has x = 0 is not
  chopped (extra-8, #6785); x^y = 0 without a sign for y stays [] (y > 0 gives [0];
  #6745).
- subs: √x => y leaving a bare x (#6976): internal radical solving relies on subs not
  rewriting x + 1 inside a radical, so it was tried and reverted.
- refine: atan(sin x/cos x) outside -π/2 < x < π/2 (the + π shifts; #6495); imaginary
  assumptions (positive/negative imaginary symbols; #6106).

### Round 18

SymPy #4409–#5908.

- integrate with a float base, 2.718…^(-2πix) (RischIntegrate: more than 8
  transcendental constants; the exact exp(-2πix) works; extra-6).
- canonical form of sums inside products: -(x + y)/z is (-x - y)/z but (x + y)/(-z)
  is -(x + y)/z, and their difference does not cancel without simplify; (-x - 1)/(-y - 1)
  stays (#4596, #4592). Pulling the sign or content out of such sums would change how
  many results print; left (user decision 2026-11-16).
- printing x^(1/3): pasted back into Julia it is a Float64 exponent; printing rational
  exponents as 1//3 would round-trip (#5458). Left for now (user 2026-11-16: tempting
  later, for round trips from the command line).
- deep nesting (f(f(…f(x))) 10,000 deep) overflows the stack in printing and diff;
  SymPy raises RecursionError too (#5557).
- expand is not idempotent where radicals merge into a sum inside a product:
  expand((1 - √(1+x))√(1+x) y) is -(x + 1)y + √(x + 1)y, and a second expand gives
  -xy + √(x + 1)y - y (#5743). Every fix tried (re-expanding merged sums in the
  multiplication loop, or a second pass for small results) changed the intermediate
  forms cancel builds and made checking an ODE system's solution 5× slower, or hours on
  Julia 1.10; dropped (2026-11-16).
- N of a besseli/gamma expression cancelling in Float64: no BigFloat besseli for the
  high-precision check (#5691).

### Round 19

SymPy #2909–#4408. Recorded 2026-11-18.

- integer factorization (`factorint`, perfect powers): `factor(Sym(1030903))` stays
  1030903 = 53²·367 (#4356, #4388, #4055).
- noncommutative symbols (#4395, #4360, #4059, #3757, #3741, #3153, #3130).
- unevaluated `Integral` and `Derivative` for users (`evaluate = false`), and a change of
  variables on an Integral (#3534, #3359, #3335, #3575, #4083, #4086).
- a symbolic product ∏ (#3410, #3284); sums over the roots of a polynomial (RootSum,
  #3674); RootOf objects (#3150).
- series: symbolic exponents at a branch point, `series(w^a, w, 0, 2)` (#4047, as
  #7168; done since: a factor `w^a` of the series); a `dir` keyword (#4162); Order objects (#4279, #4266); convergence
  acceleration (#3543).
- closed forms of the alternating series Σ (-1)ᵏ log k / k = γ log 2 - log² 2/2 (#4189)
  and Σ (-1)ᵏ (ζ(k) - 1)/(k 2ᵏ) = log(3√π/4) + (γ - 1)/2 (#4188); both open in SymPy.
  Σ (2k + 1)²/(2k)! stays unevaluated too (found with the #4171 fix).
- fsolve: Newton's method from a start with a singular Jacobian gives up
  (`fsolve([x^y ~ 1, x ~ 2], [x, y], [x => 1.0, y => 1.1])` misses y = 0; #4027).
- functions: `bernoulli(n)` for users (#4193), Stieltjes constants (#4186), Kronecker
  delta (#3772), spherical harmonics (#3537), Legendre polynomials as a function
  (#3330); exact cos(π/8), cos(π/17) and other constructible values (#3804).
- helpers: `jacobian`, gradient/divergence/laplacian (#3899, #3922); `has(ex, cos)` for
  a function head (#3282); `diff(ex)` for one variable (#3649); an algebraic-number
  predicate (#4176); `as_independent` (#3784); a formula classifier (#4157).
- product-to-sum for unlike arguments, sin(x)cos(y) (#3456); factoring over extensions,
  x³ - 2 (#3298).
- LaTeX printer options: folding function brackets, inline fractional exponents, a LaTeX
  name per symbol, line wrapping (#4397, #4347, #4340, #4346, #3001).
- out of the core: symbolic optimization (#4173), L² Gram–Schmidt of functions
  (#4211), vector and tensor calculus (#4298, #4238, #4058, #3115), sets with
  cardinality (#2995), Boolean algebra (#2943), ranges attached to symbols (#3588),
  the method of brackets (#4344), a Degree constant (#4180).

### Round 20

SymPy #945–#2908, #31017–#31179. Recorded 2026-11-18.

- matrix expressions: MatrixSymbol, symbolic shapes and slicing, block matrices
  (BlockDiagMatrix, `block_collapse`) (#2778, #2748, #2452, #2459, #2460, #31178).
- number theory: a Chinese remainder function `crt` for users (#31177, #31126);
  modular equations over ℤ in solve (#31048).
- sets: Range sets and iteration over unions (#31145); symbolic intervals with a
  measure (#31030).
- transforms: the Fourier transform (#2803); residues (#2519; done since: `SymJ.residue`).
- `eliminate` (Mathematica's Eliminate; groebner by hand) (#2720); boolean constraints
  in solve's systems (#2715); symbolic Boolean algebra, Xor/Implies/Equivalent/ITE
  (#2697, #1441).
- lambdify with nested argument groups (#2790); symbol-name modifiers in display,
  xHat as x̂ (#2518); Lie-symmetry infinitesimals of ODEs (#2372).
- from the fixes: a two-sided `asech`, `asec`, `acsc` series at a branch point
  (only points inside a cut were made consistent); √(sinh² u + 1) = cosh u for real u
  in simplify (found checking Kamke 1.113); homogeneous ODEs with |x| (they are solved
  with `assume = (x ≫ 0,)` or `(x ≪ 0,)`); conditions on cos(a) and other bounded
  functions in solve's cases are not reduced (x = cos(a)√(x² + 1) has redundant cases);
  ∫₀¹ log(x^i)² and other integrals whose antiderivative jumps at infinitely many
  points are refused.

## Fixes by cause

Status from the tracker's database; the SymJ fixes without one ("open", 117) or in progress (9) were rechecked against main on 2026-11-21 and given one ("feature gap": a missing feature, listed in the [feature audit's roadmap](../feature_audit.md#open-by-area) rather than fixed here); a slug with a suffix `-rN` is a cause seen again in round N.

| fix | owner | kind | status | PR | issues | note |
|:--|:--|:--|:--|:--|:--|:--|
| `gruntz-besselk-at-zero-unsupported` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#58 | 22911 | K_ν at 0+ uses its series there (integer ν with log terms; otherwise through I_±ν). A negative argument takes the principal value K_ν(-v) = exp(-iνπ)K_ν(v) - iπI_ν(v), so x*besselk(1,x) gives 1 two-s… |
| `gruntz-bounded-oscillating-subterm` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#68 | 14503, 14298 | bounded_limit expands B^q and log(B), B = A + O with O/A -> 0 oscillating, as A^q(1+δ)^q / log A + log(1+δ) with a bounded Taylor remainder; x - sqrt(x^2+sin x) -> 0, (n^2+n)log(1+1/(n^2+cos n)) -> 1. |
| `gruntz-bounded-oscillating-subterm-r19` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#77 | 4131 | A fractional power of a sign-changing bounded factor is rated as complex, bounded in modulus only, and bounded out when the modulus tends to 0: log(2 + sqrt(atan(x)*sin(1/x))) -> log(2), sqrt(x*sin(1… |
| `gruntz-case-conditions-not-normalized` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#74 | 6843 | Case conditions and sign queries are reduced to lowest terms (exact division by one-parameter polynomial denominators) and identically zero ones dropped: with n > 0 the #6843 limit is n/2 (SymJ: assu… |
| `gruntz-case-split-on-unnormalized-zero` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#74 | r17-extra-2 | No split on n/(n+1) + 1/(n+1) - 1 (identically 0): gruntz_cases gives n == -1 => undefined, n != -1 => n/2 (1.6 s, was 13 s and 24 leaves); SymJ's plain limit gives n/2. |
| `gruntz-case-split-zero-unsolved` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#57 | r5-extra-10 | solve_zero now also solves c linear in a parameter with parameter coefficients (r = -1/x, adding x != 0) and exp(u) - 1 = 0 as u = 0; a condition still unsolved is kept in the context and the limit i… |
| `gruntz-complex-values` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#54 | 27087, 27058 | Input with im (or complex numbers, or SymJ's (-1)^(1/2)) is split into real and imaginary parts and their limits taken separately. Covers arithmetic, exp, log, powers (principal branch), sin/cos and … |
| `gruntz-complex-values-r11` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#67 | 15486 | erf/erfc of a complex argument u+iv as σ(1 - exp(-z²) erfcx(σz)) with bounded erfcx parts; limits at u -> ±oo (∫ exp(-2x²)cos(3x) over ℝ = sqrt(π/2)e^(-9/8)); finite points still unsupported. |
| `gruntz-complex-values-r14` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#71 | 10868 | acosh/asech of an argument in x take principal values (acosh(u) = i acos(u) for -1<u<1): asech(x) at Inf and acosh(x) at 0+ give iπ/2; log(x)+asech(x) at Inf is oo + iπ/2 and raises the 'direction of… |
| `gruntz-complex-values-r19` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#77, jverzani/GruntzLimit.jl#78 | r19-extra-2 | Values on branch cuts are Julia's functions at u + 0im (SymJ's N). #77 (merged, acf8b30): asin/acos, so acos(u) = -i acosh(u) for u > 1 and acos(1/(1-x)), acos(1+x) -> 0 at 0+. #78 (open): asech flip… |
| `gruntz-complex-values-r20` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#78 | r20-extra-8 | acoth on (0, 1) at its ends: GruntzLimit main (#78) gives log 2 at 1⁻ and iπ/2 at 0⁺; integrate(acoth(x), (x, 0, 1)) = log 2 - iπ/2 through SymJ (checked 2026-11-18) |
| `gruntz-complex-values-r4` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#55 | 25681 | abs(b^e) is taken as abs(b)^e (\|b^e\| = \|b\|^e for b < 0 too), so x/abs(sqrt(x^2-1)) gives -Inf at -1+ and Inf at 1-; the SymJ-side rewrite is no longer needed |
| `gruntz-complex-values-r5` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#57 | 23915 | already works on the development branch (round 3/4 complex support, unreleased since 1.3.0): (1/10 + i/10)^k and the float and ratio forms give 0; \|b\| >= 1 (not real) reports that it oscillates. Re… |
| `gruntz-complex-values-r6` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#58 | r6-extra-4 | atan of a complex argument is now i/2 (log(1 - iz) - log(1 + iz)) on the principal branch, in the complex path. atan(x - im) and atan(x + 2im) give π/2 at +oo, and x^2*(atan(x-im) - π/2 + 1/x) gives … |
| `gruntz-complex-values-r7` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#60 | 22220, 21899, 21733 | \|z\|^2 base powers cancelled, so log((x-i)/(-x-i)) -> iπ; a real log or fractional power of an argument in x that is negative near the point is retried on the principal branches ((t+1)cos(log(-t-1))… |
| `gruntz-constant-zero-sqrt-perfect-square` | GruntzLimit | wrong answer | fixed | jverzani/GruntzLimit.jl#58 | 23525 | sqrt(c p^2) with p a constant polynomial is now sqrt(c)\|p\|, with the sign of p taken from the context (Cases split on the sign of p), so m+1-sqrt(m^2+2m+1) is 0 for m > -1. #23525 gives 1/(2 sqrt(m… |
| `gruntz-constant-zero-trig-unrecognized` | GruntzLimit | wrong answer | fixed | jverzani/GruntzLimit.jl#68 | 14411 | a leading series coefficient that is numerically zero (cos(c), c = π/2 not reduced) is dropped outside cases too; 3 sec(4πx - x/3) at 3π/(24π-2) gives -oo/+oo by side. |
| `gruntz-constant-zero-trig-unrecognized-r17` | GruntzLimit | wrong answer | fixed | jverzani/GruntzLimit.jl#74 | 6364 | Parameter constants that vanish as rational functions of their atoms (cos^2 = 1 - sin^2) are recognized as 0: the #6364 limit is 1/(sin(a)^2/2 - cos(a)^2/2) = -2/cos(2a) (was 0); a = 3/10 and the 1 -… |
| `gruntz-cosint-negative-argument` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#58 | 22787 | cosint of a negative argument is now its principal value Ci(\|u\|) + iπ, through the complex path, so cosint(3x)-cosint(x) gives log(3) from both sides (two-sided, SymJ's default). cosint(x) at -oo g… |
| `gruntz-exp-large-integer-exponent` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#67 | 16222 | large integer powers of constants stay powers (exp(x) at 10^6 -> exp(1)^1000000); a product of exp/π/positive constants is never a numeric zero (exp(-x) at 1000 also failed). |
| `gruntz-float-exponent-too-large` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#58 | 23266 | A ratio of exponents with a large numerator or denominator (from floats) is kept as its own w^r leaf; the w^g rescaling is skipped for a tiny gcd; and polyPS caps its number of terms, lowering the st… |
| `gruntz-floor-of-infinite-argument` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#60 | 21310 | floor(u) with u -> ±oo is u - frac(u), frac bounded in [0,1) for the bounded-oscillation code: x*floor(1/x) -> 1 at 0 (both sides); fixed in GruntzLimit, so SymJ needs no squeeze |
| `gruntz-leadterm-nested-exp-cancellation` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#68 | 14811 | series exponents normalized (rational part split out) so w^(2-L) and w·w^(1-L) cancel, and exact all-zero series give 0; #14811 gives oo. Also fixes a wrong Inf for 2^x log(1+(2/3)^(x+1)) - (2/3)(4/3… |
| `gruntz-leadterm-power-product-cancellation` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#69 | 13332 | leadterm drops term groups whose constant coefficient is numerically zero (L/L - 1, L = log(5/6), from u^u powers); results cancel B(S)^-1 against multiples of S and combine constant roots; #13332 gi… |
| `gruntz-no-min-max` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#52 | 30586, 30568 | min/max of any number of arguments rewritten as (u + v ± \|u − v\|)/2 on input; all four examples give the expected limits. |
| `gruntz-no-mod` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#57 | 25230 | mod(u, v) = u - v*floor(u/v), rem via trunc, trunc(u) = sign(u)*floor(\|u\|); mod(x, 3) at 3- gives 3. SymEngine.jl evaluates mod(x, 3) to 0 itself before GruntzLimit sees it |
| `gruntz-parameter-coefficient-zero-not-cancelled` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#74 | 6366 | Zero coefficients like n/(n+1) + 1/(n+1) - 1 recognized, and results reduced to lowest terms: #6366 gives n/2 for dir = :+, :- and both. |
| `gruntz-real-root-negative-base` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#58 | 22824 | cbrt(u) is now the real cube root sign(u)\|u\|^(1/3), for an argument in x or a negative constant; cbrt(x) at -8 gives -2 and the #22824 ratio gives -12. GruntzLimit handles this, so SymJ needs no re… |
| `gruntz-series-blowup-nested-tanh` | GruntzLimit | hang | fixed | jverzani/GruntzLimit.jl#54 | 26991, 26990 | exp(-3x/100) is exp(x)^(-3/100), so the series were in w^(3/100) and needed ~33x the terms per nested function; the expansion variable is now w^g (g = gcd of w's exponents). The limit takes 0.5 s and… |
| `gruntz-sqrt-square-unsigned` | GruntzLimit | wrong answer | fixed | jverzani/GruntzLimit.jl#54 | 27551 | Fractional powers of even powers of atoms of x (log, sin, atan, ...) now use the atom's eventual sign, so sqrt(u^2) = \|u\|; also fixes log(-log(x)) and sqrt(-log(x)) for x < 1, and abs(sin(x))/x -> 0 |
| `gruntz-symbolic-exponent-case-split-hang` | GruntzLimit | hang | fixed | jverzani/GruntzLimit.jl#75 | 5172 | No longer hangs: the zero test merges exps (the vanishing coefficient at r = 1 - 2^c is recognized), and series coefficients over 1000 terms stop with 'too complicated'. The case tree ends in ~7 s (6… |
| `gruntz-unsupported-function` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#52 | 30641, 30581 | erf(a,b) = erf(b) − erf(a); lambertw (principal branch) supported: Taylor series at finite points, de Bruijn asymptotic series at +∞; exp(W(x)) → ∞. |
| `gruntz-unsupported-function-r10` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#66 | 17109, 16722 | binomial(a,b) handled as a gamma ratio, or a polynomial when b or a-b is an integer; li2 resolved at -oo (inversion) and at 1 from below (reflection), Taylor below 1; li0/li1/li(n,u)/polylog(n,u) for… |
| `gruntz-unsupported-function-r12` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#68 | 13681, 13571 | upper incomplete gamma(s, z) for s -> +oo (z constant > 0) and z -> +oo (s constant); gamma(n,1)/gamma(n) -> 1. impolylog (#13681) not added: needs Im Li2 of a moving complex argument; a SymJ table e… |
| `gruntz-unsupported-function-r17` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#74 | 6909 | zeta (and SymJ's _zeta_num) supported at finite limits: exact at the integers (ζ(2) = π²/6), by continuity elsewhere, 1/(u-1) + γ at the pole. Derivatives are unknown, so limits needing them fail cle… |
| `gruntz-unsupported-function-r18` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#75 | 5415 | polygamma with an order in x tending to an integer m >= 0, at a finite non-pole argument, as polygamma(m, u) + O(order - m): #5415 gives 9/4 - 2ζ(3). polygamma at positive integers is now exact for e… |
| `gruntz-unsupported-function-r7` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#60 | 21701 | besselj(ν, u) at a finite limit: (u/2)^ν/Γ(ν+1)(1+O(u²)) at 0, symbolic order too (besselj(7,z)/z^7 -> 1/645120), Taylor's series elsewhere; at ±oo still unsupported |
| `gruntz-unsupported-function-r8` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#64 | 19864 | clausen2 supported: Taylor series at finite points, u - u log\|u\| + ... at 2πk, exact at rational multiples of π, bounded oscillation at ±∞; x log 2 + Cl₂(2x-π)/2 at π/2 gives π log 2/2 |
| `limit-negative-side-symbolic-power-r17` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#74 | 7096 | x^b for an irrational constant b (√2, π) of a negative argument is retried on the principal branch: x^(-√2) at 0- gives the clean 'infinite, in a direction of the complex plane' error; x^√2 at -1 is … |
| `limit-oscillating-difference` | GruntzLimit | missing case | fixed | jverzani/GruntzLimit.jl#52 | 29567 | Oscillating sin/cos whose arguments differ by a finitely-converging amount are merged by angle addition; sin/cos(x+1/x) − sin/cos(x) → 0. x·(sin(x+1/x) − sin x) ~ cos x has no limit (the tracker's ex… |
| `definite-antiderivative-complex-on-interval` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#170 | 13792 | step LR: log(1/u) as -log(u) (equal up to a locally constant term); ∫ log(1/x)/(1 - x) is now -li2(1 - x), real on (0, 1), ∫₀¹ = π²/6 (was li2(1/x) with log(x - 1), complex there). Commit 723df9a. |
| `erf-integrand-unsupported` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#92 | 29909 | try_erf integrates A(x)·erf(ax+b) and erfi by parts (Q = ∫ A, exponentials merged into e^(∓u²)): ∫ eˣ erf x = eˣ erf x − e^(1/4) erf(x − 1/2), ∫ x eˣ erf x also |
| `erf-integrand-unsupported-r18` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#150 | 5326 | step SN: ∫ A·erf(u) (A, u rational, ∫A not rational) is now a proven NotElementary (Liouville), not 'unsupported': ∫erf(y)/y. The definite value needs SymJ's numeric fallback (definite-nonelementary-… |
| `integrate-no-variable-returns-integrand` | RischIntegrate | wrong answer | fixed | jverzani/RischIntegrate.jl#122 | r9-extra-10 | integrate(ex) with no variable (or no tuples) now raises an ArgumentError instead of returning ex; fixed in RischIntegrate, no SymJ integrate(::Sym) needed |
| `integrate-number-integrand-definite-crash` | RischIntegrate | crash | fixed | jverzani/RischIntegrate.jl#122 | 18631 | definite integrate: an antiderivative that is a plain Julia number (from 0 or 0.0) gives 0 instead of being called; fixed in RischIntegrate, no SymJ number methods needed |
| `no-polylog-for-nonelementary` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#92 | 30167 | new step XP: P(x)·R(e^(kx)) by parts when ∫ R is elementary; ∫ x/(eˣ − 1) now has an answer with Li₂ (degree-1 P; x²/(eˣ − 1) needs Li₃, still NotElementary) |
| `no-polylog-for-nonelementary-r11` | RischIntegrate | missing case | won't fix | jverzani/RischIntegrate.jl#128 | 15227 | the antiderivative needs Li3: NotElementary is correct for the indefinite integral; the value -5ζ(3)/4 belongs to a definite table in SymJ |
| `risch-antiderivative-atan-jump-real-line` | RischIntegrate | wrong answer | fixed | jverzani/RischIntegrate.jl#167 | r12-extra-2 | step KV: a radicand κ/V (κ > 0 rational) as κ^e V^(-e) where V > 0; ∫√(16/(1-x²)) is now 4 asin(x), F(1/2) - F(-1/2) = 4π/3 (was 8 atan(…/x), jump of 8π at 0). Commit 652f7b1. |
| `risch-atan-tan-complex-hang` | RischIntegrate | hang | fixed | jverzani/RischIntegrate.jl#106 | 26509 | The hang was in this package: steps US and TS2 tried about 120 substitutions over Q(i)(a,b,c,cos a,sin a), and each was refused. atan of a non-real argument has no real form on any route, so it now g… |
| `risch-check-tower-pi-growth` | RischIntegrate | hang | fixed | jverzani/RischIntegrate.jl#97 | 27957 | Step SG: not π-specific — a parameter met after the opaque radical (ℓ²)^(3/2) in SymJ's term order nested the constant field badly; plain symbols now go before opaque subterms. Indefinite ∫ dl3 now ~… |
| `risch-constant-radical-rejected` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#105 | 27374 | Radicals of π, ℯ are opaque constants; step IT takes x-free factors with π out (beside √2, π could not join the constant field), an x-free integrand with a parameter or π is f·x, and step LS flattens… |
| `risch-constant-radical-rejected-r14` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#139 | r14-extra-3 | step RK splits pi and e out of a radicand like a parameter: sqrt(PI/x) = sqrt(pi) sqrt(1/x), giving 2 sqrt(pi)/sqrt(1/x) through SymJ |
| `risch-cosint-integrand-unsupported` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#150, #153 | 5741 | the third integrand (Ci(π² z)/z beside elementary terms, digamma(1)) is a proven NotElementary (step SN on sums; unknown constant functions opaque); ∫Ci, ∫Si by parts (#150). SymJ check 2026-11-16: t… |
| `risch-ei-fallback-symbolic-exponent` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#118 | 20740 | step 8a's Ei candidates Ei(jU + k log u): x^n/log x = expinti((n + 1) log x); log² and log³ follow |
| `risch-ei-integrand-unsupported` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#147 | 7000 | step EB: ∫A·Ei(g) by parts: ∫Ei(log x) = x Ei(log x) − Ei(2 log x). No li(x) function (SymJ feature gap) unchanged |
| `risch-erf-complex-gaussian-trig` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#115 | 21362, 21041 | Step GT: xⁿ e^(αx²+βx+δ) with sin/cos of linear arguments as exponentials, by recurrence to erf of a complex shift: ∫ sin(x)e^(-x²) = (i/4)e^(-1/4)√π(erf(x+i/2) − erf(x−i/2)). Kept here, so SymJ need… |
| `risch-erf-complex-gaussian-trig-r11` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#128 | 15292 | step GT takes complex x^2 coefficients: cos(x^2)e^(-x^2), sin(x^2), cos(x^2) give erf of sqrt(1±i)x, sqrt(±i)x (indefinite); the symbolic-θ definite form and the Inf for exp(-(1-i)x^2) on (0,∞) are S… |
| `risch-erf-log-argument-missing` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#150 | 4890 | step LU after NotElementary: an exp of log x, x = e^u with the exps merged: ∫e^(-log²x) = √π e^(1/4) erf(log x - 1/2)/2 (also x^a e^(-log²x)); SymJ ∫₀^∞ = √π e^(1/4) |
| `risch-erf-rational-argument-missing` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#171 | 13634 | step XR: K exp(c + a x² + b/x²), a b > 0, by completing the squares (√\|a\| x ∓ √\|b\|/x)²: ∫ e^(-x²-1/x²) = √π/4 (e⁻² erf(x - 1/x) + e² erf(x + 1/x)); over (0, ∞) √π e⁻²/2. erfi for a, b > 0. Commit… |
| `risch-erf-rational-argument-missing-r18` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#150 | r18-extra-7 | step UR: exp of a polynomial in 1/x by u = 1/x: ∫e^(1/x²) = x e^(1/x²) - √π erfi(1/x); also e^(a/x²), e^(1/x) (Ei). e^(-x²-1/x²) still refused |
| `risch-erf-symbolic-coefficients` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#105 | 26922 | The erf kind of the parametric Risch DE takes any constant A in e^(A x² + B x + C): √π/(2√(-A)) e^E erf(√(-A)(x + h)). Branch: erf when A < 0, erfi when A > 0, the sign from the oracle (assume), else… |
| `risch-erf-symbolic-coefficients-r18` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#150 | 5547 | step GS: Gaussians k·c·z² with c symbolic or π, z = u/√c then step GT: cos(a z²), sin(a z²), e^(i a z²) (also SymJ's ℯ^w), cos(π z²), cos(2πz²/L)² |
| `risch-exp-argument-misread` | RischIntegrate | wrong answer | fixed | jverzani/RischIntegrate.jl#105 | 27586, r3-extra-4 | Exponents of exp(u) and ℯ^u (SymJ, SymEngine) are read through one helper, exp_arg (times_exp, plus lambertw, inverse_trig, polylog_parts, gamma_factors, which missed ℯ^u). verify missed it because t… |
| `risch-exp-complex-linear-argument` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#118 | 20034 | step EC: exp(α + iβ) = exp(α)exp(iβ) when the generator has no real form, through Expr (SymJ merges ℯ^z ℯ^(iz) back); fixed here, not in SymJ |
| `risch-exp-transcendental-constant` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#108 | 23949 | exp(c) the constant field cannot express (exp(exp(5)), exp(π√2)) is an opaque constant like a parameter; #23949 gives −x e^(e⁵+2) + x e^(e⁵+log x/(x+6)) − x e^(x²+e⁵); no SymJ glue needed |
| `risch-exp-transcendental-constant-r11` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#128 | 15648 | step EK: exp(p(y) + q) = exp(q) exp(p(y)) with q free of y, on the exp-of-constant refusal, through Expr (SymJ re-merges the factors); definite (0, x) also works |
| `risch-exp-transcendental-constant-r7` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#115 | 21034 | exp or a hyperbolic function of a constant with a transcendental function inside (cosh(cos 4)) is an opaque constant before its argument enters the field; asin(5) via step OC |
| `risch-float-exponent-ratio-huge-power` | RischIntegrate | hang | fixed | jverzani/RischIntegrate.jl#128 | 15498 | an exponential argument ratio with numerator or denominator over 1000 (exact floats) is refused at once; step LS then integrates the terms apart (1.5 s, was a hang) |
| `risch-float-exponent-ratio-huge-power-r12` | RischIntegrate | hang | fixed | jverzani/RischIntegrate.jl#168 | 13813 | x^2.5556 = x^(6389/2500): the substitution x = u^2500 made x/(x - 3.6^70) a degree-2500 rational integral. Radical-free summands are now integrated apart when a root index exceeds 12: 0.2 s (was a ha… |
| `risch-global-ref-flags-not-thread-safe` | RischIntegrate | wrong answer | fixed | jverzani/RischIntegrate.jl#140 | 10450, 10449 | 31 per-call flags/depths are fields of a CallState in a ScopedValue made by the outermost call; LEAF_VALUES scoped; HYP_MOD/HYP_NF Refs under a lock; 240 concurrent calls on 4 threads: 45 differed be… |
| `risch-inverse-trig-nonlinear-argument` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#108 | 24025 | elementary: t eᵗ log t atan(t eᵗ) + √(t eᵗ log t) acos(eᵗ log t); new step AM integrates by parts over a sum of arcs of tower elements, residual tested numerically |
| `risch-inverse-trig-of-constant` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#115 | 21024 | Step OC: atan, asin, acos, asinh … of a constant are opaque constants (atanh, acoth stay logs): x − atan(4) log(x + atan(4)); SymJ needs no opaque symbols |
| `risch-large-integer-power-refused` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#134 | 12446 | Step BP: K(ax+b)^n with integer \|n\| > 10000 by the power rule, base kept whole (x^100000 -> x^100001/100001, (2x)^5638654 -> (2x)^5638655/11277310); other powers beyond the cap refused with an accu… |
| `risch-li-sum-numerator-not-split` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#139 | 9367 | step LS splits SymJ's product (x^2 - 1)*log(x)^-1 term by term; through SymJ the indefinite integral is Ei(3 log x) - Ei(log x), and on (0, 1) log 3 (log 2 for x - 1) |
| `risch-log-sqrt-trig-difference-hang` | RischIntegrate | hang | fixed | jverzani/RischIntegrate.jl#147 | 6382 | step LR: log(u^e) = e log u (0<e≤1) first; the integrand is now refused (NotElementary, Clausen antiderivative) in seconds. The definite value 0 still needs SymJ's periodic rule |
| `risch-log-structure-constant-branch-lost` | RischIntegrate | wrong answer | fixed | jverzani/RischIntegrate.jl#134 | 14613 | Same fix as -r13: (cx)^m (bx^n)^p now prints both powers as written, F = (cx)^m (bx^n)^p x/(np+m+1), F'-f = 0 at c = x = -1, m = 1/2 (#14613) |
| `risch-log-structure-constant-branch-lost-r13` | RischIntegrate | wrong answer | fixed | jverzani/RischIntegrate.jl#134 | r13-extra-1 | A power u^k (k not rational) prints as the user's u^k, not (ax)^n/a^n; with such a power, logs with a parameter are parsed first (log(ax) a level, not log x + log a). x^n log(ax) = x^(n+1)(log(ax)/(n… |
| `risch-log-sum-over-denominator-not-split` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#139 | 10851, 10798 | step LS splits a quotient written as a product N*D^-1 (SymJ's form) term by term, not only N/D; #10851 gives -li2(-y)/2 + li2(y)/2 and the outer integral of #10798 gives pi^2/8 through SymJ |
| `risch-multiple-constant-square-roots` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#106 | 26502 | Did not build Q(sqrt2, sqrt3): its roots also need sqrt(1 + 2 sqrt3), and RootSums over number fields are not implemented. Step SQ instead: constant square roots with different squarefree parts are o… |
| `risch-nameof-symbol-crash` | RischIntegrate | crash | fixed | jverzani/RischIntegrate.jl#92 | 30413, 30412, extra-12 | expr_build accepts a Symbol operation and merge_exp_powers keeps unchanged subterms: loggamma, digamma, polygamma, cosint·besselj now RischNotImplemented; new step LW: ∫ W(ax+b) by x = (w eʷ − b)/a, … |
| `risch-notelementary-slow` | RischIntegrate | hang | fixed | jverzani/RischIntegrate.jl#151 | r18-extra-10 | the ~60 s was compilation (5 parameters, 2 exp levels); with one type for the constant field (CF, default in #151) 23 s on a first call, 0.1 s after |
| `risch-numeric-check-large-frequency` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#115 | 21741 | The numeric check retries with Richardson's (4D(h/2) − D(h))/3 when the central difference misses (error h²ω²/6 at ω ≈ 32000): e^(-40000ixz) passes |
| `risch-polylog-fallback-sech-squared` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#147 | 6875 | step XP takes hyperbolic functions as exps and nests for lower-degree P: ∫x²/cosh²x with Li₂; SymJ ∫₀^∞ = π²/12. x³/cosh²x (Li₃ of e^(2x)) still refused |
| `risch-power-of-power-branch-ignored` | RischIntegrate | wrong answer | fixed | jverzani/RischIntegrate.jl#147 | 6251 | step RP keeps the locally constant factor ζ = u^E/w^(E m) for odd m: ∫(x³)^(-1/3) = ((x³)^(-1/3)/x^(-1)) log x, right at every real x (principal powers) |
| `risch-radical-complex-coefficient` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#108 | 25314 | radicand with i: i treated as a parameter through Expr (step IP, via step PR), put back and checked over ℂ; both #25314 integrands answered |
| `risch-radical-of-log-unsupported` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#128 | 15640 | step LU: radicals of log x by x = e^u; sqrt(log x)/x^2 gives an erf answer (special terms kept as Expr for SymEngine's erfi) |
| `risch-radical-of-log-unsupported-r14` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#139 | 11876 | step LU reads log(x^k) as k log x (log(1/x), log(x^-1), log(sqrt x)); sqrt(log(1/x)) integrates, the definite integral on (0,1) gives sqrt(pi)/2 through SymJ |
| `risch-radical-of-log-unsupported-r19` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#157 | 4231 | Step RY extended for a square root of a tower element: conjugate denominators (also (x+y)^-1 as SymJ writes it) and a logarithmic part c log((s+y)/(s-y)). Answer 2√(x+log x) + log((x+√(x+log x))/(x-√… |
| `risch-radical-with-parameter-unsupported` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#122 | 19188 | step PR keeps its answer with a RootSum (s = d u, as for d = 1) when step PS cannot scale it away; a real form for the parametric biquadratic RootSum is still open |
| `risch-radical-with-parameter-unsupported-r12` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#169 | 14035 | param_sign read b² - 1 as 0 at the default value b = 1 and the step refused; a nonzero constant that vanishes there is now evaluated at k + 1/2. b sin t/(2(1 + (b² - 1)cos²t)^(3/2)) integrates (long … |
| `risch-radical-with-parameter-unsupported-r17` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#147 | 7109, 6926 | step PR: a radicand κ/V (κ with a parameter, deg V ≥ 2) as (V/κ)^(−e): √(a²/(a²−x²)) = a asin(x/a). #6926 already integrated on main (SymJ gives a continuous piecewise) |
| `risch-radicand-parameter-sign-branch` | RischIntegrate | wrong answer | fixed | jverzani/RischIntegrate.jl#147 | 7147 | step PR writes the answer's radical as the integrand's own √q/√κ when s or √κ is an opaque root: ∫x/√(ax²+bx+c)³ right for both signs of b²−4ac (and a<0) |
| `risch-residues-degree-gt2-over-radical-constants` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#113 | 23688 | Step CB: when constant √d comes in through step 10e and the log part over ℚ(i)(√d) has residues of degree > 2, the integrand (denominator over ℚ(i)) is split as g + √d·h and each part integrated over… |
| `risch-rootsum-term-zero-denominator` | RischIntegrate | crash | fixed | jverzani/RischIntegrate.jl#106 | 26508 | Step CI: the imaginary unit (im, I, SymJ's (-1)^(1/2)) is i of Q(i), not a parameter t (which gave the denominator (1+t^2)^4 = 0). A non-real integrand is split as u + i v and each part is put in rea… |
| `risch-sum-acos-sqrt-wrong-antiderivative` | RischIntegrate | wrong answer | fixed | jverzani/RischIntegrate.jl#147 | 7154 | step PW's sign factor σ multiplied the whole sum; now only the summands with an odd power of the radical (branches otherwise). ∫₀¹⁰ = 35π/2 |
| `risch-symbol-named-e-as-euler` | RischIntegrate | wrong answer | fixed | jverzani/RischIntegrate.jl#123 | 16951, r10-extra-3 | is_e_const: a leaf printed E is e only when it equals its term type's own exp(one(t)) (SymEngine's E, SymJ's constant); a user symbol E, and Expr's :E, stay symbols; fixed here, no renaming in SymJ |
| `risch-transcendental-constants-algebraic-field` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#108 | 24952 | √5 beside π (no parameters): the constant square root becomes opaque (step SQ), answer checked numerically; exp(-√5 s)cos(2πs) gives (2π sin − √5 cos)e^(−√5 s)/(4π²+5) in SymJ |
| `risch-try-erf-infinite-recursion` | RischIntegrate | crash | fixed | jverzani/RischIntegrate.jl#118 | 19846 | ∫ erf(u)·κe^w (w quadratic, w ± u² not constant) is now NotElementary; try_erf by parts at most two deep, no stack overflow (the definite value needs SymJ's fallback) |
| `rootsum-not-expanded-symbolic-coeffs` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#97 | 28657 | Step BN extended to x⁴ + b x² + c with b/√c rational (x = c^(1/4) u → u⁴ + γu² + 1 over ℚ): log/atan answer. SymJ ∫₀^∞ 1/(16π⁴f⁴−4π²f²+1) = 1/4; c⁴f⁴−c²f²+1 with c>0 gives π/(2c). Signs of b, c from … |
| `rootsum-not-expanded-symbolic-coeffs-r11` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#128 | 15726 | for a sum, an answer with RootSums is compared with the term-by-term answer (step LS) and the one with fewer RootSums kept: logs only here |
| `rootsum-not-expanded-symbolic-coeffs-r15` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#139 | 8954 | step BN shifts a binomial (A + x)^n + c to u = x + s first; 1/((A + x)^4 - C^4) and the issue's form give logs and atan, no RootSum |
| `sin-cos-fractional-powers` | RischIntegrate | missing case | fixed | jverzani/RischIntegrate.jl#92 | 29882 | step HS's u = sin x ₂F₁ now takes fractional exponents, after step HT's NotElementary or a radical refusal: ∫ sin^(7/2)·cos^(3/2) = ₂F₁ answer (checked numerically) |
| `abs-even-power-quotient-not-reduced` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 20108 | \|u\|^m u^(2k) = \|u\|^(m+2k) for real u: x²/\|x\| = \|x\| |
| `abs-of-complex-power-not-evaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 14274, 14260 | Fixed on main (checked 2026-11-21). |
| `abs-positive-constant-not-evaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | r14-extra-9 | abs of a real constant of known sign is ± itself at construction (abs(√2) = √2, abs(-π) = π). |
| `abs-product-not-combined` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 18840 | r9-missing-cases: \|a\|·\|b\| merged in simplify |
| `abs-quotient-common-factor-not-cancelled` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 17292 | abs factors with integer exponents merged and cancelled |
| `abs-split-one-period` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#85 | 30319 | Definite integrals of piecewise integrands split at every interior breakpoint, periodic conditions repeated over the interval. max(sin, cos) still needs solve-mixed-trig. |
| `adjoint-conj-skipped-complex-constants` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#112 | 25162 | Root cause in has(): a product's power factors (-im, x*sqrt(y)) were not seen; svdvals now (1±√5)/2 |
| `adjoint-eltype-not-inferred` | SymJ | missing case | fixed | #103 | r3-extra-3 | conj/adjoint inferred Sym |
| `algebraic-radical-ratio-not-reduced` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 18475 | radical constants with degree-1 minpoly reduced; degree bound ≤ 16 skips Cardano nests |
| `algebraic-zero-test-misses-nested-radical` | SymJ | crash | fixed | jverzani/SymJ.jl#118 | 22604 | Zero test over radicals: a constant whose value at 512 bits is below 1e-100 is 0 (nested radicals that do not cancel by form); jordan_form and the ODE systems work |
| `algebraic-zero-test-misses-nested-radical-r8` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#123 | 20427 | Numerically-zero constants (512-bit test) dropped as Poly coefficients; simplify gives 0 |
| `alpha-equivalent-integrals-not-identified` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5539 |  |
| `apart-float-coefficients` | SymJ | missing case | fixed | #110 | 26517 | float coefficients exactified, then floatified, in inverse transforms |
| `apart-float-coefficients-r5` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24481 | apart takes floats as the rationals they stand for, monic denominators, floats again |
| `assumptions-on-expressions-not-symbols` | SymJ | missing case | fixed | #89 | 29082 | Membership of ±x + k (ℤ) or c x + d (ℚ, ℝ, ℂ; c, d rational) is membership of x: refine(sin(πn); assume = (-n ∈ ℤ,)) = 0. 2n ∈ ℤ, x^2 ∈ ℝ: a clearer error. No parity facts (scope: if feasible). |
| `atan-of-radical-tan-value-not-recognized` | SymJ | missing case | fixed |  | 21324 | Fixed on main (checked 2026-11-21): simplify(cos(25 atan(…))) is 0. |
| `atan-sin-cos-quotient-not-recognized-as-tan` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6495 | Fixed for -π/2 < x < π/2; other ranges (atan(sin x/cos x) + π) not reduced. |
| `atan2-complex-args-numeric-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#161 | 7747 | atan(y, x) with complex values as -i log((x + iy)/√(x² + y²)). |
| `backslash-nonsquare-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 14548 | Fixed on main (checked 2026-11-21): least squares and least norm. |
| `base-predicates-false-on-symbols` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#86 | extra-5 | Kept by design: isodd/iseven/isinteger are false for a symbol (patterns rely on it). assume = (isodd(n),) now explains that instead of 'contradictory: 1 ~ 0'. Parity assumptions are not supported. |
| `bessel-order-symmetry-half-integer-not-expanded` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6154 | Negative integer orders; half-integer orders to elementary functions not done. |
| `bigfloat-to-float64-overflow-silent` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#120 | r7-extra-4 | Sym(::BigFloat) outside Float64's range is an ArgumentError (floats are Float64) |
| `bigint-power-inexacterror` | SymJ | crash | fixed | #92 | 28112 | integer powers beyond ~10^7 bits stay unevaluated powers |
| `bigint-power-inexacterror-r9` | SymJ | crash | fixed | jverzani/SymJ.jl#128 | 18526 | A root index past Int: the power left as it is (N gives 1.0) |
| `binomial-negative-integer-diagonal` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#116 | 23497 | binomial(n, k) for k ≤ n < 0 = (-1)^(n-k) binomial(-k-1, n-k), as binomial(n, n) = 1 (Mathematica's convention) |
| `branch-cut-rewrite-numeric-inconsistent` | SymJ | wrong answer | won't fix | jverzani/SymJ.jl#116 | 23760 | Conventions differ on branch cuts (N follows Julia's floating point; the log forms are textbook identities, which disagree even between SymPy and Mathematica for atanh); documented in rewrite's docst… |
| `calculus-ops-on-cases-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#139 | 14235, 14211 | Fixed on main (checked 2026-11-21): diff, integrate and nested summation of Cases work; a sum of terms whose sums depend on cases is refused. |
| `calculus-ops-on-cases-methoderror-r18` | SymJ | crash | fixed | jverzani/SymJ.jl#173 | r18-extra-1 | N(::Cases) added. |
| `cancel-drops-removable-singularity` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#84 | 30792, 30529, 30522, 30521 | Equations, inequalities and point solutions of systems. A free variable's exclusion in a system (x ≠ 1 with x free: #30521, #30522) is left as is, by decision. |
| `cancel-exp-complex-constant-reciprocal-not-merged` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | 6932 |  |
| `cancel-exp-multiple-generators-not-merged` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | 12787 | Powers of one base with rationally related symbolic exponents as one generator. PR open. |
| `cancel-float-coefficients-gcd-trivial` | SymJ | missing case | open |  | 13805 | Still open (checked 2026-11-21): cancel with float coefficients finds no common factor. |
| `cancel-fractional-power-common-factor` | SymJ | missing case | won't fix | jverzani/SymJ.jl#125 | 19769 | Design: cancel treats non-integer powers as independent (documented); simplify now cancels over one base's root when smaller ((x-1)/(√x-1) = √x + 1) |
| `cancel-imaginary-unit-as-generator` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | 6336 | i conjugated only when the other generators are symbols or square roots of rationals. |
| `cancel-multivariate-large-slow` | SymJ | hang | fixed | jverzani/SymJ.jl#115 | 24537 | cancel of nested fractions whose expansion would be huge: over the lcm of factored denominators, level by level (540-term numerator instead of 16M); 20 s instead of > 300 s |
| `cancel-piecewise-pieces-not-mapped` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 7022 |  |
| `cancel-radical-conj-divexact-error` | SymJ | crash | fixed | jverzani/SymJ.jl#131 | 17732, r10-extra-7 | _conj_cancel gives up when the gcd does not divide; adjugate eigenvectors for symbolic nested radicals |
| `cardano-degenerate-cases` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#86 | 30898, 30837 | Triple root (Δ0 = Δ1 = 0, once expanded) and Δ0 = 0 handled without 0/0: #30898's quintic now gives its roots. #30837 (substituting p = 0 into the generic formula) left as is, as in SymPy: it would n… |
| `cardano-symbolic-real-roots-complex-form` | SymJ | missing case | feature gap |  | 11514 | Feature gap (2026-11-21), in the feature audit's roadmap: polynomials: real roots of a symbolic cubic in trigonometric form (casus irreducibilis). |
| `cases-arithmetic-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#135 | 16182 | arithmetic and series on Cases case by case |
| `cases-limit-drops-generic-solution` | SymJ | missing case | fixed | jverzani/SymJ.jl#154 | 17186 | Fixed on main (checked 2026-11-21): one generic solution under the generic-parameter default. |
| `charpoly-uncancelled-denominator` | SymJ | wrong answer | fixed | #91 | 28793, 28212 | charpoly by Berkowitz (division-free); eigenvalue counts right |
| `coeff-non-symbol-target-silent-zero` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#164 | 7072 |  |
| `coeff-non-symbol-target-silent-zero-r18` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#172 | 4934, r18-extra-5 | Sum targets through a fresh symbol; number targets refused. |
| `coeff-non-symbol-target-silent-zero-r19` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#175 | 4167 | powers of non-symbols (e^u, (x+1)^2) through a fresh symbol |
| `coeff-power-target-variable-exponent` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#172 | 4917 |  |
| `complex-constant-cartesian-form-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | 16348 | parts of constant complex powers via exp(q log b) |
| `complex-constant-cartesian-form-missing-r14` | SymJ | missing case | feature gap |  | 11193 | Feature gap (2026-11-21), in the feature audit's roadmap: simplification: a complex constant to a + bi form. |
| `complex-constant-cartesian-form-missing-r16` | SymJ | missing case | fixed | jverzani/SymJ.jl#162 | 7689, 7582 | √(a + ci) parts; #7582's quartic-root real part not checked separately. |
| `complex-constant-parts-power-limited` | SymJ | missing case | fixed | jverzani/SymJ.jl#153 | r15-extra-10, r15-extra-8 | _const_parts: any integer power (squaring; past 8 only numeric parts), negative powers through the conjugate; abs(b^n) = abs(b)^n. |
| `complex-constant-product-not-folded` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#151 | r15-extra-5 | ishermitian/issymmetric for Sym matrices compare entries by value (difference expands or cancels to 0). Products like im*(im + 1) still stay unfolded at construction. |
| `complex-float-arithmetic-not-evaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | 12154 | Products and integer powers of complex floats folded. PR open. |
| `complex-float-functions-not-autoevaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6071 | Also exp of a complex float, and a negative float to a fractional power (principal value). |
| `complex-infinity-real-part-nan` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#160 | 8080 | _const_parts: zero parts stay zero under infinite coefficients. |
| `complex-sqrt-denesting-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 20024 | √(a + bi) with a²+b² a square and rational parts: 2 + i |
| `complex-sqrt-denesting-missing-r13` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | r13-extra-9 | √c √(a+bi) rejoined before denesting: √(8+6i) = 3 + i. PR open. |
| `complex-sqrt-denesting-missing-r19` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | 4401 |  |
| `conj-times-self-not-real` | SymJ | missing case | fixed | #89 | 29801 | real/imag of an expression real for every value (conj(w)^k w^k, abs(z)^2, real-valued functions): imag(x*conj(x)) = 0, isreal(x*conj(x)); simplify(abs(x)^2 - conj(x)*x) = 0. sqrt(e) - conj(sqrt(e)) i… |
| `conj-times-self-not-real-r10` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 16500 | f(u) f(conj u) paired for f commuting with conj |
| `const-rational-sum-of-transcendental-powers` | SymJ | missing case | fixed | #93 | r2-extra-4 | algebraic polynomials in pi or e are irrational |
| `const-zero-absolute-threshold` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#130 | 17450, r10-extra-1 | tiny is not zero: products/powers by factors, sums by rounding level at rising precision; _isreal_const/_csign likewise |
| `constant-comparison-float64-overflow` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#151 | r15-extra-9 | _decide and _csign decide at 256 bits when Float64 overflows, underflows or cannot tell (relative to the sides' size, a sum's by its terms). |
| `constant-powers-common-prime-base-not-merged` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | 15669 | rational powers of integers merged by primes |
| `constant-powers-common-prime-base-not-merged-r14` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 11544 | simplify splits rational bases, (p/q)^x = p^x q^(-x): (1/2)^x - 2^(-x) is 0. |
| `constant-powers-common-prime-base-not-merged-r17` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6370 |  |
| `constant-powers-common-prime-base-not-merged-r19` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | 3382 |  |
| `cse-drops-derivatives` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#172 | 5085 |  |
| `cse-matrix-input-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#139 | 14118 | Fixed on main (checked 2026-11-21). |
| `cse-matrix-input-methoderror-r16` | SymJ | crash | fixed | jverzani/SymJ.jl#161 | 8891 | _shape_like keeps arrays of arrays. |
| `cse-no-common-subproduct` | SymJ | missing case | fixed | jverzani/SymJ.jl#159 | 11428 | Fixed on main (checked 2026-11-21). |
| `cse-no-common-subproduct-r15` | SymJ | missing case | fixed | jverzani/SymJ.jl#159 | 10228 | _cse_subproducts!: factors two or more products share become one replacement. |
| `cse-no-common-subsum` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23297 | cse: terms two or more sums share are one subexpression (t1 => c^2/2 + d/9) |
| `cse-reciprocal-power-not-shared` | SymJ | missing case | fixed | jverzani/SymJ.jl#179 | 4203 |  |
| `cse-relation-input-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#165 | 6263 |  |
| `cse-root-of-unity-powers` | SymJ | missing case | fixed | #94 | 28051 | powers of a root of unity reduced into (-n/2, n/2] before cse counts |
| `deep-nesting-recursive-stackoverflow` | SymJ | crash | won't fix |  | r18-extra-11 | SymPy also raises RecursionError; would need non-recursive printing and diff. |
| `definite-abs-dirac-splitting` | SymJ | missing case | fixed | #106 | 27108 | conditions/heaviside/dirac of one sign settled by interval enclosures before splitting (#27283) |
| `definite-antiderivative-complex-on-interval-r17` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#164 | 7155 | Now refused as a branch cut; a correct value would need an antiderivative valid on (0, 1). |
| `definite-antiderivative-complex-on-interval-r20` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#180 | r20-extra-4 | refused (branch-cut error) rather than a value: log(x^i) jumps at x = e^-(2k+1)π, infinitely many points near 0 |
| `definite-antiderivative-simplify-blowup` | SymJ | hang | fixed | jverzani/SymJ.jl#113 | 23967 | simplify skips expansion and cancel when the expanded size estimate is huge; the definite integral takes 3 s |
| `definite-antiderivative-simplify-blowup-r14` | SymJ | hang | fixed | jverzani/SymJ.jl#146 | 11354 | The trig zero test applies only when every sin/cos argument is a rational multiple of π, so simplify of the endpoint value no longer rewrites everything in polar form. |
| `definite-beta-unit-interval-form-missing` | SymJ | missing case | feature gap |  | 14549 | Feature gap (2026-11-21), in the feature audit's roadmap: integration: Beta forms on [0, 1]. |
| `definite-beta-unit-interval-form-missing-r13` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | 12456 | B(a,b)(ψ(a) - ψ(a+b)) for a log(x) factor, ψ(b) for log(1-x). PR open. |
| `definite-case-branch-endpoint-pole` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23613 | atan(t) = ±π/2 has no root (no tan(π/2)); acos-point comparisons as conditions on C; simple poles inside refused; equal points merged: ±2π/√(C²-1), divergent for \|C\| ≤ 1 |
| `definite-case-endpoint-pole-value-finite` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#151 | r15-extra-2 | A parametric singular point equal to an end by the facts makes that end singular (the point itself as the end); a complex-infinite limit there with a real ±∞ integrand gives ±∞. The k ~ 0 cases' d/k … |
| `definite-case-split-subs-pole` | SymJ | crash | fixed | jverzani/SymJ.jl#124 | 19672 | A case where the integrand is undefined (b = 0) is an error for that case; nested denominators count as the integrand's own: b > 0 gives √(bπ) |
| `definite-case-split-subs-pole-r13` | SymJ | crash | fixed | jverzani/SymJ.jl#143 | 12568, 11992, r13-extra-6 | A fact undefined at another case's values is dropped in _build and undecided in _holds. Results still carry spurious complex singular-point cases. PR open. |
| `definite-case-tree-endpoint-limits-slow` | SymJ | hang | open |  | 14394 | Still open (checked 2026-11-21): ∫₀^π 1/(a²cos²x + b²sin²x), a, b > 0: finishes, but as Cases of errors, not π/(ab). |
| `definite-case-value-from-degenerate-parameter` | SymJ | wrong answer | won't fix |  | r13-extra-14 | Not a wrong answer: definite integrals take parameters real (definite.jl:34); for real c the case condition forces c = 0 and the value is right. The spurious √(-c²) singular-point cases remain (defin… |
| `definite-cases-error-object-as-value` | SymJ | missing case | won't fix |  | r8-extra-3 | Design: a case with no value holds its error (the case tree's convention; b = 0 in #124 too) |
| `definite-complex-limits-trig-zeros-real-only` | SymJ | crash | fixed | jverzani/SymJ.jl#180 | r20-extra-5 |  |
| `definite-complex-path-branch-cut-side` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#180 | 31017 |  |
| `definite-complex-singularity-not-discarded` | SymJ | missing case | fixed | #107 | r4-extra-6 | a Cases with the same value in every case collapses to it (acos(A) split for A > 1 gone) |
| `definite-complex-singularity-not-discarded-r20` | SymJ | crash | fixed | jverzani/SymJ.jl#181 | 2537, r20-extra-1 |  |
| `definite-convergence-needs-param-sign` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#120 | 21054 | A limit failing on an undecided sign: a case per sign of the parameter; a limit not found is 'cannot find', not divergence: a > 0 gives 1/(a²+1). (GruntzLimit cannot split on sign(a²) given a ≠ 0.) |
| `definite-convergence-needs-param-sign-r11` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | 15432 | gamma moment divergence is its own case |
| `definite-degenerate-subst-abs-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#115 | 24518 | SymJ #115: degenerate sums of squares by their real components, equal alternatives merged; RischIntegrate #111 (trig products) does the degenerate integrals: all five cases right with both main branc… |
| `definite-dirac-product-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5588 |  |
| `definite-dirac-sifting-unknown-function` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6968 |  |
| `definite-divergent-log-endpoint-unsigned` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23337 | log(u) = log(-u) + iπ at an end: a real infinite part beside finite iπ terms is the value: ∫₁³ dx/(x-3) = -∞ |
| `definite-endpoint-pole-substituted` | SymJ | wrong answer | fixed | #101 | 27185 | trig functions at a pole raise DivisionByZero; the end is then a one-sided limit |
| `definite-equal-limits-not-zero` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5098 |  |
| `definite-floor-ceil-integrand-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#163 | 8480, 8420 |  |
| `definite-fourier-rational-residues-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | 6136, 5909 | Residues for simple poles; multiple poles not yet. |
| `definite-ftc-on-branch-cut` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#138 | 14147 | Fixed on main (checked 2026-11-21): an antiderivative off on the interval (a branch cut) is refused, not used; no value by quadrature yet. |
| `definite-gaussian-imaginary-parameter-coefficient` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | 6123 |  |
| `definite-gaussian-trig-factor-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 14573 | Fixed on main (checked 2026-11-21). |
| `definite-hyper2f1-endpoint-limit` | SymJ | missing case | fixed | #109 | 26477, 26068 | 2F1 terms by their endpoint values (continuity; Gauss at 1); Beta integral works; #26068 depends on RischIntegrate giving a 2F1 antiderivative |
| `definite-hyper2f1-endpoint-limit-r5` | SymJ | missing case | fixed | jverzani/SymJ.jl#115 | 24879 | limit through ₂F₁: Taylor at z₀, the expansion at z = 1, the 1/z connection formula at -∞, truncation checked by agreement; definite integrals take the ₂F₁ form where its radicand is positive: ∫₋₁^∞ … |
| `definite-impossible-case-not-refuted-under-assumptions` | SymJ | missing case | fixed | jverzani/SymJ.jl#159 | 10211 | _isnonzero proves a + ib ≠ 0 for real a, b with b ≠ 0; y ≠ √(-h²) for real y, h > 0. |
| `definite-impossible-case-not-refuted-under-assumptions-r16` | SymJ | missing case | won't fix |  | 8556 | Feature audit 19: linear facts with irrational coefficients. |
| `definite-impossible-case-not-refuted-under-assumptions-r19` | SymJ | missing case | fixed |  | r19-extra-9 | not reproduced on main after round-19 work |
| `definite-inner-cases-unordered-endpoints` | SymJ | missing case | open |  | 7164 | Still open (checked 2026-11-21): an inner integral's Cases on x's sign split the outer interval: refused even with a, b > 0. |
| `definite-integral-parameter-no-cases` | SymJ | wrong answer | fixed | #101 | 27280, 27073 | a case per parameter factor of the value's denominator not in the integrand's |
| `definite-jump-points-unlocated-ignored` | SymJ | wrong answer | open | jverzani/SymJ.jl#138 | 13609 | Still open (checked 2026-11-21): the e^(2πit)/(e^(2πit) - 1/2) case is refused; ∫₀¹ 2πi e^(6πit)/(e^(2πit) - 1/2) still gives iπ/4 (iπ/2). |
| `definite-known-table-missing-zeta` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24431 | ∫₀^∞ xᵖ/(e^(ax) ∓ 1): Γ(p+1)ζ(p+1)/a^(p+1), with (1 - 2^-p) for +1; 2ζ(3), 3ζ(3)/2 |
| `definite-nonelementary-no-fallback` | SymJ | missing case | fixed | #89 | 30869, 30677, 30664 | #30664 done: Beta integrals over [0, π/2] when no antiderivative is found (∫ 1/√(sin x) = Γ(1/4)√π/(2Γ(3/4))). #30869/#30677 (differentiation under the integral sign) deferred to a standalone PR, by … |
| `definite-nonelementary-no-fallback-r11` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | 16396, 16084 | King's rule; log sin/cos table; clausen2(k pi) = 0 |
| `definite-nonelementary-no-fallback-r12` | SymJ | missing case | open |  | 14154, 14127, 14065, 14043 | Still open (checked 2026-11-21): ∫₀¹ (1 - x⁴)^(-1/3): RischNotImplemented; ∫₀^∞ x/cosh x: the limit fails (impolylog) with no fallback; the others fixed. |
| `definite-nonelementary-no-fallback-r13` | SymJ | missing case | open |  | 12645 | Needs polar coordinates; not addressed; in notes/feature_audit.md (item 16). Still open (checked 2026-11-21): ∫∫ sin(x² + y²) over the disc now crashes: MethodError in _fresnel_num (fresnel.jl:26); was RischNotImplemented. |
| `definite-nonelementary-no-fallback-r14` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 11845 | ∫₀ᵇ e^(c - a xᵖ) dx as a lower incomplete gamma; also ∫₀¹ xᵐ log(1/x)ˢ = Γ(s+1)/(m+1)^(s+1) (which covers #11876's sqrt(log(1/x)) and log(1/x)^(a-1)). |
| `definite-nonelementary-no-fallback-r15` | SymJ | missing case | feature gap |  | 9461 | Feature gap (2026-11-21), in the feature audit's roadmap: integration: polar coordinates or Bessel K₀ for a double integral over the quadrant. |
| `definite-nonelementary-no-fallback-r17` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | 6348 | exp(ix)/(1 + x²) by residues: π/e. |
| `definite-nonelementary-no-fallback-r6` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 22487 | Fallback over (-∞,∞)/(0,∞): sums term by term, Fresnel sin/cos(ax²), parts for n(x)/x² with n(0)=0: π - √(2π) |
| `definite-nonelementary-no-fallback-r7` | SymJ | missing case | feature gap |  | 21030 | Feature gap (2026-11-21), in the feature audit's roadmap: integration: contour or series methods (∫₀^π e^(cos x) cos(sin x)); the value by quadrature since #184. |
| `definite-nonelementary-no-fallback-r8` | SymJ | missing case | open | jverzani/SymJ.jl#125 | 20831, 20685, 20171, 20043, r8-extra-5 | Divergence by the limit at ∞ (Inf) and Γ moments ∫₀^∞ xᵖlogᵏx e^(-ax^q) (Γ(4/3), Γ''(2)); #20171 (contour) and the finite-interval fallback open (feature_audit 11) Still open (checked 2026-11-21): partly fixed; #20171 and #20043 still NotElementary or RischNotImplemented. |
| `definite-nonelementary-no-fallback-r9` | SymJ | missing case | open |  | 18796 | Still open (checked 2026-11-21): ∫₀^∞ cos(πx/2)cos(x)/(1 - x²) unevaluated, and N fails at the removable pole x = 1. |
| `definite-odd-part-symmetry-unused` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23246 | Over (-a, a): odd terms continuous there dropped before integrating: π |
| `definite-odd-part-symmetry-unused-r7` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 21091 | Odd integrand over (-∞, ∞), no singular points, absolutely convergent at ±∞: 0 |
| `definite-odd-part-symmetry-unused-r9` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 19175 | r9-missing-cases: even part used over the whole line |
| `definite-parameter-case-split-divideerror` | SymJ | crash | open |  | 19046 | The DivideError no longer occurs on main; the answer is a poor case split on √(-g²) ∈ ℝ with 'known only numerically' branches (g, K > 0 fine): a missing case Still open (checked 2026-11-21): no DivideError; without assumptions, cases on √(-g²) ∈ ℝ whose branches are errors; g, K > 0 gives the answer. |
| `definite-pole-case-conditions-unnormalized` | SymJ | missing case | open |  | r18-extra-4 | Still open (checked 2026-11-21): the conditions stay as (-w⁴/π⁴)^(1/4) ≶ 0 and the like; with w > 0 assumed, 3√2w/8. |
| `definite-result-not-refined-under-assumptions` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | r5-extra-5 | Each case's value refined under its assumptions: n ≠ 1 gives 0 |
| `definite-singularity-search-transcendental` | SymJ | missing case | fixed | #106 | 27139 | _zeros_in returns none when an enclosure excludes 0 (#27139 now NotElementary) |
| `definite-singularity-search-transcendental-r4` | SymJ | wrong answer | fixed | #107 | 26283 | split also at the antiderivative's singular points (tan(t/2) jump at π) |
| `definite-termwise-split-false-divergence` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 17058 | Frullani rule; no ∞ - ∞ claim from opposite divergent terms |
| `derivative-at-point-display-ambiguous` | SymJ | missing case | fixed | #105 | r3-extra-8 | Shown as Derivative(f(t, x), t)(t => a); differentiated in all arguments. |
| `diff-abs-assumes-real-argument` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#172 | 4757 |  |
| `diff-integral-bound-variable-in-limits` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5086 |  |
| `diff-max-min-unevaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 22726 | ∂max(a,b)/∂a = heaviside(a - b), ∂min(a,b)/∂a = heaviside(b - a) |
| `diff-negative-order-silent` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#138 | 13873, r12-extra-3 | Fixed on main (checked 2026-11-21). |
| `diff-of-sum-not-moved-inside` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5312 |  |
| `dirac-argument-not-normalized` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | 12257 | δ(c u) = δ(u)/\|c\|, u with leading coefficient 1. PR open. |
| `dirac-nonlinear-argument` | SymJ | missing case | fixed | #89 | 29001 | diracdelta of a polynomial: sum over simple real roots; transcendental arguments still unevaluated |
| `dirac-times-vanishing-factor-not-reduced` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 13996 | Fixed on main (checked 2026-11-21). |
| `division-by-zero-ordinary-inputs` | SymJ | crash | fixed | jverzani/SymJ.jl#113 | 23952, 23942, 23914 | A sum with 0^(c - k) is summed in the other direction (k => a + b - k): 1 - 0^(n+1) |
| `dsolve-abel-ode-unsupported` | SymJ | missing case | feature gap |  | 16387 | Feature gap (2026-11-21), in the feature audit's roadmap: ODEs: Abel equations. |
| `dsolve-bessel-equation-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#132, jverzani/SymJ.jl#174 | 17397 | common polynomial factor divided out (#132); round 18 (#4414): y'' + 2y'/x + y = 0 through the normal form y = w v, constant invariant (#174). Bessel ODEs in general remain a feature gap. |
| `dsolve-bessel-equation-unsupported-r15` | SymJ | missing case | fixed | jverzani/SymJ.jl#159 | 9819 | x²y'' + xy' + (α²x² - ν²)y = 0 gives C1 J_ν(αx) + C2 Y_ν(αx). |
| `dsolve-bessel-equation-unsupported-r17` | SymJ | missing case | won't fix |  | 7224 | Feature gap: notes/feature_audit.md item 20. |
| `dsolve-bessel-equation-unsupported-r18` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 4414 | Fixed on main (checked 2026-11-21): u'' + 2u'/x + u = 0 solved. |
| `dsolve-bvp-two-value-conditions-only-crash` | SymJ | crash | fixed | jverzani/SymJ.jl#124 | 20547 | Constant solution only when the value conditions agree (only() crashed); C'(L) => 0 at a symbolic point accepted |
| `dsolve-constant-names-ignore-function-names` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | 15056 | constants avoid function names |
| `dsolve-derivative-power-not-solved` | SymJ | missing case | open |  | 15889 | Still open (checked 2026-11-21): f'² = f³: "the solution found does not satisfy". |
| `dsolve-derivative-power-not-solved-r12` | SymJ | missing case | fixed |  | 14760 | Fixed on main (checked 2026-11-21): u'² = g(u) split into u' = ±√g. |
| `dsolve-derivative-power-not-solved-r17` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | 7093 |  |
| `dsolve-exact-linear-ode-not-recognized` | SymJ | missing case | fixed | #110 | 25479 | exact linear ODEs integrated once per order when Euler fails |
| `dsolve-exact-linear-ode-not-recognized-r7` | SymJ | missing case | open |  | 21791 | Still open (checked 2026-11-21): (f g)'' = 0 solved, but as nested integrals of g'/g, not (C1 + C2 x)/g. |
| `dsolve-factored-ode-branch-dropped` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#151 | 10043 | Without initial conditions, an ODE with y in two or more factors is solved factor by factor (a vector of solutions). |
| `dsolve-homogeneous-radical-not-recognized` | SymJ | missing case | fixed | jverzani/SymJ.jl#182 | r20-extra-7 | with assume = (x ≫ 0,) or (x ≪ 0,); without one only \ Fixed on main (checked 2026-11-21): with x > 0 or x < 0 assumed. |x\| makes it homogeneous |
| `dsolve-ics-at-singular-point-ignored` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#134 | 15113 | singular initial point an error; float Gauss-Jordan rejects non-finite entries |
| `dsolve-ics-multiple-solutions-error` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 16880 | candidates checked near the initial point |
| `dsolve-linear-system-cubic-eigenvalues` | SymJ | hang | fixed | jverzani/SymJ.jl#121 | 21867 | Parametric eigenvalues of degree ≥ 3: clear error; numeric: eigenvalues in a + ib form, conjugates matched by value, N not shadowed: real solution (17 s) |
| `dsolve-linear-system-parameter-verify-fails` | SymJ | missing case | open |  | 18562 | Still open (checked 2026-11-21): f' = f + b g, g' = b f - g: "the solution found does not satisfy"; b = 2 solves. |
| `dsolve-misses-singular-solution` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#86 | 29827 | ICs at an equilibrium (y(a) = c, derivatives 0, y ≡ c solves) give the constant: both #29827 cases. Also found: with y(1) = 2 the implicit answer has Integral(1/log(2), 2): an indefinite Integral's v… |
| `dsolve-nonzero-factor-not-divided-out` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5279 |  |
| `dsolve-other-function-derivative-taken-as-unknown` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#138 | 13691 | Fixed on main (checked 2026-11-21). |
| `dsolve-reduction-of-order-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23351 | 2nd-order linear, variable coefficients: a homogeneous solution guessed (polynomial deg ≤ 4, e^±x), then reduction of order |
| `dsolve-riccati-rational-particular-solution-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#182 | r20-extra-6 | Fixed on main (checked 2026-11-21). |
| `dsolve-riccati-substitution-missing` | SymJ | missing case | fixed |  | 21503 | Fixed on main (checked 2026-11-21): u' = (x + u)² gives tan(x + C) - x. |
| `dsolve-separable-radical-product-not-split` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 4866 | Now an implicit solution. |
| `dsolve-special-riccati-integrable-case` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5348 |  |
| `dsolve-system-singular-leading-matrix-dae` | SymJ | missing case | feature gap |  | 9251 | Feature gap (2026-11-21), in the feature audit's roadmap: ODEs: differential-algebraic systems. |
| `dsolve-unevaluated-integral-ic-subs` | SymJ | wrong answer | fixed | #101 | 27479 | unevaluated antiderivatives taken from the conditions' point; empty integrals are 0 |
| `eigvals-float-entries-radical-formula-inaccurate` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#116 | 22685 | Floats with symbols and eigenvalues of degree ≥ 3: an ArgumentError saying to substitute numbers first (the radical formula has no accuracy in floating point) |
| `eigvals-triangular-charpoly-expanded-then-factored` | SymJ | hang | fixed | jverzani/SymJ.jl#135 | 15395 | triangular (5x5 up): the diagonal |
| `eigvecs-adjugate-quartic-radical-blowup` | SymJ | hang | fixed | jverzani/SymJ.jl#135 | 15076 | oversized symbolic eigenvalues: an error, not a hang |
| `eigvecs-float-eigenvalues-empty-nullspace` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#112 | r5-extra-12 | A matrix of numbers with a float goes to LAPACK for eigvals/eigvecs/eigen (exact order kept); the 4x4 no longer errors |
| `eigvecs-numeric-eigenvalues-exact-matrix-refused` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 20351 | Float eigenvalues of an exact matrix: numeric eigenvectors for them, exact for exact ones (0, 1) |
| `ellipe-ellipk-at-one-not-evaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | r18-extra-8 |  |
| `empty-symbol-name-accepted` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | 4101 |  |
| `eq-real-imag-parts-not-split` | SymJ | missing case | fixed | #103 | 27759 | real and imaginary parts of real-valued sides |
| `eq-real-imag-parts-not-split-r17` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#164 | 6821 |  |
| `eq-real-imag-parts-not-split-r4` | SymJ | missing case | fixed | #110 | 25783 | call-form error with keywords is an ArgumentError, not a MethodError |
| `eq-sum-of-signed-terms-not-split` | SymJ | missing case | fixed | #103 | 27755 | terms of one sign summing to 0: each 0 |
| `erf-erfc-complement-rule` | SymJ | missing case | fixed | #93 | 28900, 28894 | simplify tries erfc(u) = 1 - erf(u) |
| `error-message-float-shown-as-exact-rational` | SymJ | missing case | won't fix |  | r16-extra-4 | Design question (feature audit 19). |
| `even-power-base-sign-not-canonical` | SymJ | missing case | won't fix |  | 17768 | Won't fix (2026-11-21): == is structural; (x + y)² and (-x - y)² stay distinct. |
| `exactfloats-nan-invalid-rational` | SymJ | crash | fixed | jverzani/SymJ.jl#143 | 12414 | Non-finite floats left alone. PR open. |
| `exp-i-pi-not-evaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#115 | 24082, r5-extra-4 | Numeric roots of a charpoly with roots of unity guessed exactly (r or r² rational) and verified exactly modulo Φ_N or in radicals: DFT5 gives ±√5, ±i√5 (√5 twice), DFT3 ±√3, i√3 |
| `exp-log-not-rewritten-as-power` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 21360 | simplify candidate exp(c log u) = u^c: 0 |
| `exp-log-not-rewritten-as-power-r19` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | r19-extra-8 |  |
| `exp-log-not-rewritten-as-power-r9` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | r9-extra-2 | r9-missing-cases: exp(n log x) → x^n in antiderivatives |
| `expand-not-idempotent-radical-product` | SymJ | missing case | won't fix |  | 5743 | Fixes tried in #174 (re-expanding merged sums in the multiplication loop, small products only, a final pass for small results) all changed cancel's intermediate forms and made checking a 3x3 ODE syst… |
| `expand-trig-float-multiple-angle` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6098 |  |
| `fact-on-real-imag-marks-symbol-real` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#134 | 14978, r11-extra-1 | order facts make real only symbols outside function calls |
| `factor-exp-common-base-generators` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 21511 | factor: e^(c u), c rational, as powers of t = e^(u/q): (y - e^(x/7))(y - e^(x/2)), e^(2x) - 1 too; not a simplify candidate |
| `factor-exp-common-base-generators-r8` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 19776 | Integer bases as powers of their smallest root (4^x = (2^x)²): (2^x - 2)(2^x + 1) |
| `factor-exp-constant-as-polynomial-in-e` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24998 | Factors free of the symbols combined into one expanded constant: x*(1 - exp(-10)) |
| `factor-exp-huge-exponent-as-polynomial` | SymJ | hang | fixed | jverzani/SymJ.jl#139 | 14300 | Fixed on main (checked 2026-11-21). |
| `factor-multivariate-slow` | SymJ | hang | fixed |  | r13-extra-12 | Fixed on main (checked 2026-11-21): an expanded product of five factors in 6 variables factors in 2.5 s. |
| `factor-radical-common-factor-not-extracted` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5149 | Under assumptions only: new factor(ex; assume) splits radicals where valid (x real here); over C the split is invalid. |
| `factor-zero-argumenterror` | SymJ | crash | fixed | jverzani/SymJ.jl#146, jverzani/SymJ.jl#173 | 11151 | factor(0) and factor_squarefree(0) are 0 (#146). Round 18 (#5153): a sum that cancels to 0 only after its exponentials share a base is 0 too (#173). |
| `factor-zero-argumenterror-r18` | SymJ | crash | fixed | jverzani/SymJ.jl#173 | 5153 | Fixed on main (checked 2026-11-21). |
| `factorial-nonpositive-integer-not-pole` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | r9-extra-4 | r9-missing-cases: factorial of negative integer raises DivisionByZero |
| `ferrari-degenerate-cases` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#130 | 17563 | q = 0, Δ0 = 0, S = 0 are cases for symbolic coefficients |
| `float-coercion-enters-function-arguments` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#116 | 22524 | _floatify leaves unknown functions' arguments alone: -4.0*f(2) |
| `float-double-root-spurious-imag` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#86 | 30346 | With a real unknown (assumptions), a float root whose imaginary part is ≤ cbrt(eps)·\|root\| is taken as real with a warning; near-duplicates merged. #30346 gives [264.0]. |
| `float-exact-not-merged` | SymJ | missing case | fixed | #103 | 27796, 27706 | integer float exponents made exact in simplify |
| `float-exponent-difference-roundoff` | SymJ | missing case | won't fix |  | r17-extra-9 | Needs a numeric solver; notes/feature_audit.md item 20. |
| `float-exponent-not-merged` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#86, jverzani/SymJ.jl#174 | 29566, extra-3 | Float exponent prints as a float (x^0.5), so x^0.5 − sqrt(x) shows as -sqrt(x) + x^0.5. Float exponents not converted to rationals (by decision, as in SymPy). Round 18 (#5822): simplify tries short b… |
| `float-exponent-not-merged-r18` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5822 | Fixed on main (checked 2026-11-21). |
| `float-matrix-elimination-no-tolerance` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#130 | 17218 | largest pivot and rounding tolerance for float entries |
| `float-matrix-elimination-no-tolerance-r14` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#145 | 11549, 10780 | Float numeric matrices go to Julia's LinearAlgebra (det, inv, lu, \) with pivoting by magnitude. |
| `float-one-coefficient-hidden` | SymJ | wrong answer | fixed | #101 | 27198, r3-extra-6 | a float coefficient 1.0 is printed (1.0*k) |
| `float-power-overflow-silent` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#145 | 11004 | Float * and ^ outside Float64's range raise an ArgumentError instead of giving 0 or Inf. |
| `float-rational-comparison-not-exact` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#160 | 8245 | SymNumber </== compare a float with an exact number exactly; _compare decides two numbers directly. |
| `float-roots-imag-residue-not-chopped` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#151 | r15-extra-1 | _floatify drops a rounding-level imaginary part when _isreal_const confirms the constant is real. |
| `floatify-imaginary-unit-in-product` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 14371 | Fixed on main (checked 2026-11-21). |
| `floor-ceil-exact-integer-value-not-detected` | SymJ | missing case | fixed | jverzani/SymJ.jl#123 | 20865 | floor/ceil/trunc of a constant within rounding of an integer: simplify decides it exactly: 0 |
| `floor-ceil-negation-identity` | SymJ | missing case | fixed | #89 | 30300 | floor(-x) = -ceil(x), likewise ceil/trunc/round |
| `floor-ceil-nested-not-folded` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 11207 | floor/ceil/round/trunc of a rounding function return it: floor(ceil(x)) = ceil(x). |
| `floor-evaluated-in-float64` | SymJ | wrong answer | fixed | #107 | 26368 | floor/ceil/round/sign decided in BigFloat near a jump |
| `floor-frac-bounds-unknown` | SymJ | missing case | fixed | #103 | 27564 | signs by cases through floor/ceil (u - w, w in [0,1)) |
| `floor-huge-constant-insufficient-precision` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#164 | 6167 |  |
| `floor-integer-summand-not-extracted` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 18689 | r9-missing-cases: floor/ceil pull out integer-valued terms |
| `free-symbols-no-array-method` | SymJ | crash | fixed | jverzani/SymJ.jl#165 | 7204 |  |
| `fsolve-misses-roots-steep-exponential` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#175 | r19-extra-6 | noise candidates no longer exhaust maxzeros; whole-line σ ≥ \|f\| at the middle |
| `fsolve-newton-complex-iterate-inexacterror` | SymJ | crash | fixed | jverzani/SymJ.jl#177 | 4027 |  |
| `fsolve-positive-dimensional-solution-set` | SymJ | missing case | fixed | jverzani/SymJ.jl#179 | 4031 | a clear error pointing to solve |
| `function-of-own-inverse-not-reduced` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5346 |  |
| `gamma-nonpositive-integer-unevaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | r6-extra-1 | gamma, digamma, trigamma, polygamma at 0, -1, -2, …: DivisionByZero, as tan(π/2) (limits give 1/gamma(x) → 0) |
| `gather-ignores-denominators` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | 15394 | gather inside denominators |
| `gaussian-complex-coefficient-taken-nonpositive` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#134 | r11-extra-10 | complex a decided on Re a |
| `gaussian-integral-with-parameter` | SymJ | missing case | fixed | #93 | 28485 | definite Gaussian rule for a > 0 by cases (Inf for a <= 0) |
| `gcd-algebraic-coefficients-unsupported` | SymJ | missing case | open |  | 18391 | Still open (checked 2026-11-21): gcd(x^2 + 1, x + im), gcd(x^2 - 2, x - √2) refused; cancel works. |
| `gcd-algebraic-coefficients-unsupported-r11` | SymJ | missing case | feature gap |  | 15358 | Feature gap (2026-11-21), in the feature audit's roadmap: polynomials: gcd over ℤ[i] and ℚ(i). |
| `gcd-algebraic-coefficients-unsupported-r12` | SymJ | missing case | open |  | 14364, r12-extra-4 | Still open (checked 2026-11-21): gcd of rationals fixed (#172); gcd(4√5/7, 2√5/3) refused. |
| `gcd-lcm-no-int-method` | SymJ | crash | fixed | jverzani/SymJ.jl#87 | 28998 | gcd/lcm of a Sym and a Julia number (gcd(2x, 0) = 2x, lcm(2x, 3) = 6x) |
| `gcd-lcm-rational-constants-field-convention` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#172 | r18-extra-12 |  |
| `gcd-lcm-varargs-no-method` | SymJ | crash | fixed | jverzani/SymJ.jl#173 | 5135 |  |
| `groebner-fglm-parametric-slow` | SymJ | hang | open |  | 12071 | Still open (checked 2026-11-21): two general cubics with 14 parameters: no basis after 400 s. |
| `groebner-lex-parametric-slow` | SymJ | hang | fixed | jverzani/SymJ.jl#88 | 30609 | Lex with parameters now goes through grevlex + FGLM, as without parameters; FGLM's linear algebra over Q(params) with RatFun. #30609: 0.34 s (was > 400 s); agrees with the numeric basis at b = 2, c =… |
| `groebner-unreduced-zero-leading-coeff` | SymJ | crash | fixed | jverzani/SymJ.jl#87 | 29750 | solve: when the basis over Q(sqrt(3)) has a vanishing leading coefficient, each algebraic constant becomes a new variable (r^2 - 3) — only then, since 10^(-189/20) would give r^20. groebner itself th… |
| `gruntz-unsupported-function-r14` | SymJ | missing case | fixed | jverzani/GruntzLimit.jl#71, jverzani/GruntzLimit.jl#72, jverzani/SymJ.jl#149 | 11776, 11465, 10804 | GruntzLimit #71 and #72 (merged, released in 1.6.0): Airy functions for a positive argument, at -oo, airyaiprime/airybiprime. SymJ #149: Airy at finite points by substitution, sphericalbesselj in ele… |
| `gruntz-unsupported-function-r3` | SymJ | missing case | fixed | #106 | 27143, 27036, 26916 | unknown functions of x raise ArgumentError; x-free ones treated as constants (#27143) |
| `harmonic-noninteger-numeric-evaluation` | SymJ | missing case | fixed | jverzani/SymJ.jl#162 | 7443 |  |
| `harmonic-zero-general-order-not-evaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 14185 | Fixed on main (checked 2026-11-21). |
| `has-misses-numeric-coefficients` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#134 | 16054 | has sees a sum's constant/coefficients and a product's coefficient |
| `huge-integer-power-not-known-integer` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 10607 | Constants built from integers by +, × and positive integer powers are integers; sin and tan of π times one are 0. The test in sin/tan is structural only (the assumptions version made precompiling han… |
| `hurwitz-zeta-digamma-rational-values-not-evaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | r11-extra-9 | zeta(s, a) and digamma at integers and half-integers |
| `hurwitz-zeta-digamma-rational-values-not-evaluated-r13` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | r13-extra-13 | ψ and ζ(s,a) at negative half-integers; ψ(1) kept by design. PR open. |
| `hurwitz-zeta-negative-a-specialfunctions-convention` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#145 | r14-extra-1 | ζ(s, a) for a ≤ 0 numerically by the Hurwitz shift; ζ(0, a) = 1/2 - a exactly. |
| `hyperbolic-of-inverse-hyperbolic-not-reduced` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | 4136 |  |
| `ilaplace-dirac-times-heaviside` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#160 | 8227 | Shifted dirac terms without heaviside. |
| `ilaplace-exp-positive-shift-not-refused` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#112 | 24595 | exp(c s) with c > 0 (or a symbol, taken positive, as in exp(a s)) is an ArgumentError, unless assumptions say otherwise |
| `ilaplace-float-repeated-pole-tolerance` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#138 | 14326 | Fixed on main (checked 2026-11-21). |
| `imag-of-rational-in-real-symbols-not-expanded` | SymJ | missing case | fixed | #110 | 26102 | negative powers via conj/\|w\|^2; integrand refined under assumptions |
| `imaginary-unit-powers-not-canonical` | SymJ | missing case | fixed | #110 | 26140 | powers of i reduced in constant simplify |
| `imaginary-unit-symbolic-power-not-exp-form` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 16797 | constant^symbolic as exp in rewrite and a simplify candidate |
| `imaginary-unit-symbolic-power-not-exp-form-r12` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 14262 | Fixed on main (checked 2026-11-21). |
| `inequality-algebraic-not-supported` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 21173 | Inequalities with radicals, log, exp, cbrt: breakpoints (roots, domain ends, poles), each piece decided numerically: (2√5/5, 1) |
| `inequality-complex-coefficients-sampled-as-real` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#142 | 13105 | Complex coefficients: n/d = (P + iQ)/\|d\|² for real x; holds only where Q = 0. PR open. |
| `inequality-cubic-root-not-recognized-real` | SymJ | crash | fixed | #108 | 26069, r4-extra-5 | Cardano uses the real cube root of a negative radicand; _rv tolerates rounding imag |
| `inequality-evaluates-uncancelled-relation` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#84 | 30584 | The cancelled relation is evaluated: a relation constant once cancelled is all reals or none. |
| `inequality-float64-breakpoint-sampling` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#151 | 10268 | _set_from sorts breakpoints past 2^50 (or too close for Float64) at 256 bits and samples exact rationals; membership compares at 256 bits when Float64 ties. |
| `inequality-infinite-bound-divisionbyzero` | SymJ | crash | fixed | jverzani/SymJ.jl#131 | 16582, r10-extra-10 | a finite side against ±∞ decided |
| `inequality-integer-unknown-unsupported` | SymJ | missing case | won't fix |  | 8221 | Feature audit 19. |
| `inequality-inverse-trig-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#123 | 19977 | Inequalities in atan/acot of rational, log, exp arguments, with acot's jump: atan(log x) > 0 is (1, ∞) |
| `inequality-root-endpoint-decided-numerically` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#160 | 8518 | Numerator roots in the set exactly for ≤/≥. |
| `inequality-system-inexact-no-retry` | SymJ | missing case | fixed | jverzani/SymJ.jl#160 | 8040, 7982, r16-extra-6 | One-unknown systems use solve's float retry. |
| `inequality-transcendental-unsupported` | SymJ | missing case | won't fix |  | 8204 | Feature audit 19. |
| `inf-power-float-result` | SymJ | wrong answer | fixed | #101 | r3-extra-2 | exact base to an infinite power that is 0 gives exact 0 |
| `inf-power-float-result-r12` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#138 | r12-extra-8 | Fixed on main (checked 2026-11-21): 1^Inf stays unevaluated, not 1.0. |
| `inf-power-float-result-r15` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#151 | r15-extra-6 | _floatpow: an exact constant to ±∞ is exact 0 where \|b\| < 1 (> 1 for -∞), \|b\| = 1 kept as a power. |
| `infinity-arithmetic-unevaluated` | SymJ | missing case | fixed | #103 | 27605 | real terms plus ±Inf is ±Inf in refine |
| `infinity-arithmetic-unevaluated-r10` | SymJ | missing case | open |  | 17701, r10-extra-13 | Still open (checked 2026-11-21): Inf + im, Inf*im + 1/2 unevaluated; complex-direction limits refused. |
| `infinity-times-sum-distributed` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#175 | 3802, 3738 | ±∞ times a constant sum: ±∞ by its sign; a symbolic sum kept as a product |
| `integer-power-base-not-reduced` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24090 | n^(u p/q) with n a perfect q-th power reduced at construction: 256^(x/2) = 16^x |
| `integer-radicand-not-split-into-prime-powers` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | 7046 | Prime factors below 10⁴ by trial division, in simplify. |
| `integer-sequence-limits-unsupported` | SymJ | missing case | fixed | #97 | 28313 | n ∈ ℤ at ±∞: sin/cos(πv) drop v's integer polynomial part with sign (-1)^p |
| `integral-rootsums-not-sym` | SymJ | crash | fixed | jverzani/RischIntegrate.jl#92, jverzani/SymJ.jl#87, jverzani/SymJ.jl#173 | extra-9 | RischIntegrate part: IntegralWithRootSums prints '- clausen2(...)' instead of '+ -1 * clausen2(...)'; diff/Sym() of an IntegralWithRootSums is SymJ's. SymJ part (dev, committed, not pushed, 2026-11-0… |
| `integral-rootsums-not-sym-r18` | SymJ | crash | open | jverzani/SymJ.jl#173 | 5231, 5147 | Still open (checked 2026-11-21): #5231 fixed; #5147: ∫₀^b 1/(7 + y - y⁵) known only numerically. |
| `integrate-abs-real-log-branch` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 22533 | Piecewise antiderivatives: log(u) of a linear u negative on the interval as log(-u): -log(-x) for x < 0 |
| `integrate-abs-unsolved-sign-condition-dropped` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#123 | 20755 | acot (and acoth, asec, acsc, …) discontinuous for signs; breakpoints from the inequality set (acot x < 0 for x < 0): π/2 + log 2 |
| `integrate-cases-condition-on-integration-variable` | SymJ | missing case | won't fix |  | r16-extra-3 | Feature audit 19. |
| `integrate-cases-not-mapped` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23610 | integrate(Cases, x): case by case with each case's conditions as assumptions (Inf and errors kept) |
| `integrate-derivative-shifted-argument-stackoverflow` | SymJ | crash | fixed | jverzani/SymJ.jl#165 | 6407 |  |
| `integrate-dirac-derivative-boundary-term-dropped` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#130 | 17383 | integration-by-parts δ terms kept |
| `integrate-dirac-times-piecewise-no-piece` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 19962 | dirac beside pieces: their value at its point (not a breakpoint): 1 |
| `integrate-equality-assumption-not-substituted` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#164 | 6985 |  |
| `integrate-float-input-exactified` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | r8-extra-4 | Indefinite integral of a float integrand: the answer's numbers floats |
| `integrate-float-input-exactified-r17` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | r17-extra-7 |  |
| `integrate-float-input-exactified-r9` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | r9-extra-1, r9-extra-3 | r9-missing-cases: definite results floated for inexact integrands |
| `integrate-non-symbol-variable-silent` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#172 | 4459 |  |
| `integrate-piecewise-dirac-piece-rejected` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 20872 | A dirac inside a piece's value taken out by linearity: 100 (the Kronecker sum stays) |
| `integrate-piecewise-multivariable-condition` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 22023 | Conditions left in another variable kept as a piecewise of it (constant in x); also a simplify closure bug that returned nothing: either order gives 1/2 |
| `integrate-piecewise-parameter-boundary-piece-nan` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#175 | r19-extra-1 | point piece 0·log 0 replaced by the limit; piecewise in parameters only no longer crashes |
| `integrate-poly-not-converted` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | 4196 |  |
| `integrate-symseries-not-converted` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 4831 |  |
| `integrate-unknown-function-raises` | SymJ | missing case | fixed | #109 | 26116, 25457, r4-extra-3 | terms with unknown functions held as Integral, rest integrated (indefinite and definite) |
| `integrate-unknown-function-raises-r12` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 14796 | Fixed on main (checked 2026-11-21). |
| `integrate-unknown-function-total-derivative-not-recognized` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6958 | Products of the unknown functions (and f²/2) are tried. |
| `interval-unbounded-end-contains-infinity` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#151 | 10326 | ±∞ are no members of real intervals. |
| `inverse-function-composition-unsimplified` | SymJ | missing case | fixed | #93 | 28376 | refine: atanh(tanh), asinh(sinh) for real v; asin/acos/atan/acosh compositions on principal ranges |
| `inverse-hyperbolic-exact-values-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | 3906, 3663 |  |
| `inverse-laplace-cancel-cubic-hang` | SymJ | hang | fixed | #108 | 26560, 25796 | no cancel over cube roots in _lap_neat (result still a long radical form) |
| `inverse-laplace-float-e-base` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | r5-extra-11 | _floatify keeps exp as exp (#114); the #25352 float cubic now by numeric partial fractions with simple poles: 0.2 s instead of 68 s, same values (r5-remaining) |
| `inverse-of-periodic-numeric-argument-not-reduced` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 17300, 17208 | atan(tan c), asin(sin c), acos(cos c), log(exp z) reduced for constants |
| `inverse-trig-exact-values-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5342 |  |
| `inverse-trig-half-angle-identity-not-simplified` | SymJ | missing case | open |  | 15940 | Still open (checked 2026-11-21): acos(D) + 2atan((D - 1)/√(1 - D²)) not simplified to 0 for 0 < D < 1 (#200 did the acos(a)/2 form, not this one). |
| `inverse-trig-pole-unevaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 17143, r10-extra-8 | poles of acsc, asec, acsch, asech, acoth, atanh, atan raise DivisionByZero |
| `isless-symbols-no-canonical-order` | SymJ | missing case | won't fix |  | 7413 | Design question (feature audit 19). |
| `isposdef-hermitian-type-sym` | SymJ | crash | fixed | jverzani/SymJ.jl#124 | r8-extra-1 | hermitian_type/hermitian for Sym; isposdef by Sylvester's criterion, signs decided exactly (undecided: error) |
| `isrational-transcendental-at-rational-undecided` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6314 |  |
| `isreal-negative-base-rational-power` | SymJ | missing case | fixed | #89 | 29335 | negative base to a real non-integer power is not real; signs of sin/cos over a bounded interval |
| `iterated-integral-outer-ranges-not-assumed` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 14096 | Fixed on main (checked 2026-11-21). |
| `jordan-form-parametric-eigenvalue-nullspace-missed` | SymJ | missing case | open |  | 9204 | Still open (checked 2026-11-21): jordan_form(k*[0 -c b; c 0 -a; -b a 0]): DivisionByZero. |
| `jordan-form-trig-eigenvalue-nullspace` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24942 | Simple eigenvalues in cos/acos (Cardano) or higher radicals: eigenvector from an adjugate column, in eigvecs, eigen and jordan_form |
| `lambdify-function-arguments` | SymJ | missing case | fixed | #103 | 27006 | function calls as variables, replaced by fresh symbols |
| `lambdify-integer-only-numeric-functions` | SymJ | crash | fixed | jverzani/SymJ.jl#113 | 24246 | factorial/binomial at non-integers through Γ (N and lambdify); an integer base to a negative integer power is a float. Float64 still overflows on 365!; BigFloat gives 0.58856 |
| `lambdify-relation-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#146 | 11248 | lambdify of a Relation returns a Bool-valued function. |
| `lambertw-exact-value-missing` | SymJ | missing case | open |  | 5080 | Still open (checked 2026-11-21): lambertw(-π/2) not evaluated to iπ/2. |
| `lambertw-identity-not-simplified` | SymJ | missing case | fixed | #89 | 29737 | simplify: W(u)^k exp(k W(u) + r) -> u^k exp(r) |
| `lambertw-numeric-fallback-real-below-branch` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#172 | 5673, 5215 |  |
| `laplace-constant-base-power` | SymJ | missing case | fixed | #93 | 28136 | b^t = exp(t log b) |
| `laplace-dirac-derivatives-unsupported` | SymJ | missing case | fixed | #110 | 26775 | forward δ^(k); polynomial inverses give dirac(t, k) |
| `laplace-dirac-removable-singularity-divzero` | SymJ | crash | fixed | jverzani/SymJ.jl#113 | 25293 | dirac(t - a)·g(t): g's limit at a for a removable singularity; gives 4π e^(-s) |
| `laplace-integral-bound-variable-same-as-t` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 17094 | bound variable t allowed |
| `laplace-no-derivative-rule` | SymJ | missing case | fixed | #93 | 28720 | held LaplaceTransform(u(t), t, s); derivative rule forward and back |
| `laplace-no-integral-rule` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24545 | L{∫₀ᵗ h} = H(s)/s; a held ∫ u(t) dt is taken from 0 |
| `laplace-piecewise-not-rewritten-heaviside` | SymJ | missing case | fixed | jverzani/SymJ.jl#153 | 10354 | Piecewise of t with half-line conditions rewritten as heaviside steps before transforming. |
| `laplace-sqrt-argument-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24309 | sin(a√t), sin(a√t)/√t, cos(a√t)/√t in the table |
| `laplace-table-gaps-t-division-symbolic-power` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 21372 | e^(αt) sin(at)/t ↔ atan(a/(s-α)); t^n ↔ Γ(n+1)/s^(n+1) for symbolic n, gamma ratios simplified: tⁿe^(-at) |
| `laplace-unevaluated-sum-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | 7177 |  |
| `laplace-unknown-function-symbolic-factor-dropped` | SymJ | wrong answer | fixed | #107 | 25815 | t-free factors kept |
| `latex-conj-not-overline` | SymJ | missing case | fixed | jverzani/SymJ.jl#182 | r20-extra-3 | Fixed on main (checked 2026-11-21). |
| `latex-derivative-not-rendered` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 4438 |  |
| `latex-derivative-power-not-parenthesized` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#175 | 3179 |  |
| `latex-float-exponent-not-standard-form` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 4429 |  |
| `latex-greek-symbol-names-not-converted` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 4576 |  |
| `latex-negative-fraction-sign-in-numerator` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | 3973 |  |
| `latex-special-function-names-as-greek` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#175 | 3911 |  |
| `latex-unevaluated-sum-not-sigma` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 4767 |  |
| `latexify-poly-no-recipe` | SymJ | missing case | open |  | 4799 | Still open (checked 2026-11-21): latexify(Poly) prints Poly(…) as text. |
| `limit-abs-negative-base-power` | SymJ | missing case | fixed | #110 | 26513 | \|u^r\| = \|u\|^r rewritten before the limit |
| `limit-complex-parameter-taken-real` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#134 | 16242 | declared-complex parameters refused in limit and definite integrals |
| `limit-float-input-exactified` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | r14-extra-2 | limit of float input floats its answer again: 0.36 + 1/k at ∞ is 0.36. |
| `limit-harmonic-not-rewritten-as-digamma` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 14204, 14125 | Fixed on main (checked 2026-11-21). |
| `limit-integral-rootsums-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#113 | r5-extra-14 | An indefinite answer whose only root sums are ₂F₁ terms is a Sym (x*₂F₁(…)); limit then gives GruntzLimit's clear 'unsupported ₂F₁' (the definite-hyper2f1 item) |
| `limit-negative-side-symbolic-power` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24067 | c b^e with e free of x and positive by the assumptions, b → 0: 0 (SymJ side; GruntzLimit's log error is theirs) |
| `limit-negative-side-symbolic-power-r16` | SymJ | missing case | won't fix |  | 8634 | Feature audit 19 (parity cases). |
| `limit-operation-name-mismatch` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 19379, 18992, 18969 | r9-missing-cases: factorial/binomial/lambertw/harmonic/li/heaviside/dirac passed to GruntzLimit by name |
| `limit-symbolic-point-generic` | SymJ | missing case | won't fix |  | 21029 | By design: a limit at a symbolic point gives the generic answer |
| `limit-unevaluated-sum-unsupported` | SymJ | missing case | feature gap |  | 14196 | Feature gap (2026-11-21), in the feature audit's roadmap: sums: limits of unevaluated sums (needs symbolic Faulhaber); refused cleanly (#140). |
| `linalg-float-symbol-mix-unsupported` | SymJ | missing case | fixed | #110 | 26271 | floats with symbols exactified (rationalize), then floatified |
| `linsolve-complex-float-coefficients-refused` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 21445 | Numeric systems with complex floats solved in ComplexF64 (largest pivot, rounding residues 0) |
| `linsolve-complex-float-coefficients-refused-r16` | SymJ | missing case | fixed | jverzani/SymJ.jl#162 | 8639 |  |
| `linsolve-surd-products-zero-pivot` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#116 | 22852 | Linear systems with radical coefficients over the field they generate (_AlgF, radicals merged); Cramer's determinant tested with radicals merged: no 0/0 solutions |
| `log-abs-power-not-expanded` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | 15854 | log(\|u\|^q) = q log\|u\| |
| `log-base-one-not-undefined` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#172 | 5100 |  |
| `log-base-relations-not-combined` | SymJ | missing case | fixed | #110 | 25860 | log quotients and prime logs in simplify |
| `log-gaussian-rational-not-split` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | r19-extra-4 |  |
| `log-integer-base-not-evaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 20953 | log(b, x) for positive rationals: the integer k with b^k = x: 3, and the sum 31/168 |
| `log-minus-one-not-reduced` | SymJ | missing case | fixed | #89 | extra-4 | log of a negative number or ±im gives its principal value; solve(exp(x) ~ -1) = [im*pi] |
| `log-ratio-of-integer-powers-not-reduced` | SymJ | missing case | fixed | jverzani/SymJ.jl#123 | r8-extra-7 | solve uses the exact two-argument log when it is an integer: 2^x = 8 gives 3 (over ℤ too) |
| `log-ratio-of-integer-powers-not-reduced-r17` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | r17-extra-6 |  |
| `log-unit-complex-algebraic-not-evaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 10428 | log of a complex constant of modulus 1 at a rational multiple of π is i·θ, checked exactly: log(1/2 + i√3/2) = iπ/3. |
| `log-zero-not-neg-infinity` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#116 | 22612 | log(0) = -∞ (so 0*log(0) is NaN); cancel leaves expressions with infinities as they are |
| `matrix-exp-complex-eigenvalues-not-realified` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 14330 | Fixed on main (checked 2026-11-21). |
| `matrix-exp-quartic-radical-eigenvalues-hang` | SymJ | hang | fixed | jverzani/SymJ.jl#124 | 19565 | exp/sin/cos of a matrix with large cubic/quartic radical eigenvalues: ArgumentError at once (as the SVD); no closed form attempted |
| `matrix-exp-scalar-multiple-sqrt-square` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | r5-extra-8 | f(aB) for a constant B through B's Jordan form, the eigenvalues aλ: no sqrt(s^2) |
| `matrix-function-irreducible-cubic-eigenvalues-refused` | SymJ | missing case | feature gap |  | 15362 | Feature gap (2026-11-21), in the feature audit's roadmap: linear algebra: matrix functions with eigenvalues roots of an irreducible cubic (an explicit error). |
| `matrix-symbolic-power-nilpotent-no-cases` | SymJ | missing case | feature gap |  | 10563 | Feature gap (2026-11-21), in the feature audit's roadmap: linear algebra: symbolic matrix powers by cases (n = 0). |
| `max-min-common-factor-not-extracted` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24578 | simplify takes common positive factors, common terms and a common sign out of max/min |
| `min-max-absorption-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | 12638 | min(r, max(…, s, …)) dropped when s ≥ r (numeric difference); max likewise. PR open. |
| `min-max-shift-not-combined` | SymJ | missing case | feature gap |  | 18643 | Feature gap (2026-11-21), in the feature audit's roadmap: simplification: min/max shifts. |
| `minpoly-exp-i-pi-rational-not-root-of-unity` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 19744 | e^(iqπ) as cos + i sin first: x² - x + 1 |
| `minpoly-trig-casus-irreducibilis` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23677 | cos(acos(c)/m + rπ) with m r integer: Chebyshev T_m(t) = ±c, resultant over c's minimal polynomial |
| `minpoly-zero-constant-factor-undecided` | SymJ | missing case | fixed | jverzani/SymJ.jl#163 | 8354 |  |
| `mixed-partials-not-canonically-ordered` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 11107 | Nested Derivative variables in one canonical order: ∂²f/∂x∂y - ∂²f/∂y∂x = 0. |
| `mod-infinite-modulus-float` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#112 | r5-extra-9 | Sym(Inf) is a float, so mod(3, Inf) = 3.0 matches Julia (not a bug); the real wrong answer, mod(sqrt(2), Inf) = sqrt(2) + NaN, now evaluates as a float |
| `mod-integer-coefficient-not-reduced` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | 15493 | refine reduces integer coefficients mod m |
| `mod-multiple-of-modulus-not-dropped` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | 12592 | Integer multiples of a symbolic modulus dropped; integer-valued ones under refine. PR open. |
| `mod-multiple-of-modulus-not-dropped-r16` | SymJ | missing case | fixed | jverzani/SymJ.jl#162 | 8531 |  |
| `mod-nested-same-modulus-not-reduced` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 13543 | Fixed on main (checked 2026-11-21). |
| `multiarg-unknown-function-chain-rule` | SymJ | missing case | fixed | #89 | extra-11 | chain rule for unknown functions of several arguments via Derivative per slot |
| `n-bigfloat-special-function-errorexception` | SymJ | crash | fixed | jverzani/SymJ.jl#139 | 14716 | Fixed on main (checked 2026-11-21). |
| `n-complex-infinity-arithmetic-nan` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#172 | 5412 |  |
| `n-float64-cancellation` | SymJ | wrong answer | fixed | #91 | 28602, 28579, r2-extra-3 | N checks its Float64 value in BigFloat at rising precision (to 2048 bits) and uses that where Float64 loses it |
| `n-float64-cancellation-r12` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#138 | 14637 | Fixed on main (checked 2026-11-21). |
| `n-float64-cancellation-r18` | SymJ | wrong answer | won't fix |  | 5691 | SpecialFunctions' besseli has no BigFloat method, so N's high-precision re-check cannot run; a known gap. |
| `n-high-precision-imag-residue-not-chopped` | SymJ | missing case | fixed | jverzani/SymJ.jl#162 | r16-extra-2 |  |
| `n-huge-argument-reduction-hang` | SymJ | hang | fixed | jverzani/SymJ.jl#113 | 24279 | N of sin/cos/… of a BigFloat beyond 2^(2^20) throws SymJ.HugeArgumentError instead of reducing for hours |
| `nan-handling` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#85 | 30343, 29296, extra-1 | NaN prints (printer used sign(), which throws for NaN); diff(NaN) is NaN, not 0. |
| `nan-handling-r15` | SymJ | missing case | fixed | jverzani/SymJ.jl#153 | 9385 | NaN absorbs in max/min. |
| `nan-handling-r17` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6635 | NaN absorbs in products and sums with symbols. |
| `negated-quotient-not-auto-cancelled` | SymJ | missing case | won't fix |  | 18799 | Would change the canonical form of sums inside products; left by user decision 2026-11-16 (feature_audit 23). |
| `negated-quotient-not-auto-cancelled-r18` | SymJ | missing case | won't fix |  | 4596 | Won't fix (2026-11-21): the sign of a sum in a product is left (deferred 2026-11-16). |
| `nested-floor-division-unsimplified` | SymJ | missing case | fixed | #103 | 27400 | floor(floor(u)/m) = floor(u/m), ceil alike |
| `no-half-angle-simplify` | SymJ | missing case | fixed | #89 | 29371 | simplify tries half-angle forms of sqrt(1 ± cos) radicands (rejected if a root of a trig square remains) |
| `nsimplify-integer-valued-float-kept` | SymJ | missing case | fixed | jverzani/SymJ.jl#162 | 8239 |  |
| `nsimplify-loose-tol-rejects-simple-fractions` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6527 |  |
| `nsimplify-quartic-algebraic-missed` | SymJ | missing case | fixed | jverzani/SymJ.jl#179 | 4143 |  |
| `numberset-contains-infinity` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#160 | 8203 | Non-finite numbers in no NumberSet. |
| `numberset-no-issubset` | SymJ | crash | fixed | jverzani/SymJ.jl#135 | 16268 | issubset on number sets |
| `numer-denom-linear-nesting` | SymJ | crash | fixed | #92 | 27957, r2-extra-2 | flat common denominator for sums with many denominators; simplify skips factor when the expansion would be large |
| `numeric-comparison-float64-precision` | SymJ | wrong answer | fixed | #91 | 28280 | comparisons of constants by the checked sign of a - b; an error when undecidable |
| `numeric-comparison-float64-precision-r8` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#123 | 20897 | BigFloat evaluation computed log(6) etc. in Float64 (Int arguments); now at the working precision; a constant at a jump decided by simplify: sign 0 |
| `numeric-derivative-at-point-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 11802 | N of a derivative at a point by central differences in 256-bit arithmetic (ζ'(4), ζ'''(4)). |
| `numeric-eval-complex-fallback-at-singularity` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#130 | 16847, r10-extra-11, r10-extra-9 | float poles raise DivisionByZero; N of an infinity real when it can be |
| `numeric-eval-unsupported-argument-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#131 | r10-extra-12 | no numeric method: ArgumentError |
| `numeric-nan-reported-as-pole` | SymJ | missing case | fixed | jverzani/SymJ.jl#138 | 13499, r12-extra-9 | Fixed on main (checked 2026-11-21). |
| `numeric-roots-float64-overflow-huge-coefficients` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#134 | 15355 | BigFloat Aberth for out-of-range or wide coefficients |
| `numeric-roots-precision-loss` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#86 | 29140 | Float solutions of polynomial systems polished by Newton (Gauss–Newton) on the numerators in 256-bit precision; kept as found if no convergence. #29140's physical H = 9.458052e-14, H·OH = 1e-14. |
| `numeric-roots-precision-loss-r11` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#134 | r11-extra-2 | Newton polishing at 256 bits, exact conjugate pairs; (u)^n = c through its base |
| `odd-symmetry-at-branch-point` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#85 | 30588 | acot/acoth take a minus sign out of a nonzero number only; acot(-x) stays, so subs(acot(-x), x=>0) is π/2. N(acoth(0)) is still NaN+NaN*im. |
| `ode-conjugate-root-pairing-syntactic` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | 15594, 15520 | constant conjugates paired by value; exact conjugate float roots (#134) |
| `ode-series-float-coefficients-verify-fails` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 22525 | Series check: a float coefficient that is a rounding residue counts as 0 |
| `one-power-complex-infinity-evaluated-to-one` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#160 | 8582 | 1^(non-finite) kept as a power. |
| `pattern-optional-coefficient-missing-power-term` | SymJ | missing case | fixed | jverzani/SymJ.jl#179 | 3208 |  |
| `pattern-power-not-split-into-factors` | SymJ | missing case | fixed | jverzani/SymJ.jl#179 | 3778 |  |
| `pattern-reciprocal-slot-no-match` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 4559 |  |
| `pattern-scaled-slot-terms-not-collected` | SymJ | missing case | fixed | jverzani/SymJ.jl#179 | 3883 |  |
| `pattern-sum-of-linear-fractions-no-match` | SymJ | missing case | open |  | 16774 | Still open (checked 2026-11-21): (~a1*a + ~b1)/(~c1*a + ~d1) + (~a2*a + ~b2)/(~c2*a + ~d2) does not match; ~u/~v + ~w/~q does. |
| `perfect-root-huge-index-oom` | SymJ | crash | fixed | jverzani/SymJ.jl#121 | 21860, 21762, 21647 | _iroot: a root index ≥ the bit length gives 0 or 1 (no 2^k); _introot_pow tries only j up to the radicand's bit length: huge rational powers instant |
| `piecewise-antiderivative-continuity-complex` | SymJ | wrong answer | fixed | #91 | 28190 | continuity from the antiderivative's value where the one-sided limit fails |
| `piecewise-antiderivative-inf-constant` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#86 | 30163 | No continuity matching at singular breakpoints (or where the limit fails); fixed pieces on an interval integrated directly; zeros of \|u\| found as u's. ∫₁² 1/\|sin x\| = 1.04761, the z-version no lo… |
| `piecewise-cond-linear-test-subs-divzero` | SymJ | crash | fixed | jverzani/SymJ.jl#152 | 9567 | The linearity test's substitutions of 0 and 1 are both guarded against poles; abs(1/(x - 1)) > 1 gives (0, 1) ∪ (1, 2). |
| `piecewise-cond-linear-test-subs-divzero-r16` | SymJ | missing case | fixed | jverzani/SymJ.jl#161 | 7840 | ArgumentError (float power out of range) treated like a division by zero in the linearity test. |
| `piecewise-condition-or-and-ops-refused` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | r6-extra-5 | piecewise takes & and \| conditions: an Or is one piece per alternative, And of Ors distributed |
| `piecewise-cut-points-unresolved` | SymJ | missing case | fixed | #106 | 27358, 27283, 27234 | linear shifts, quadratic conditions, term-by-term; periodic antiderivative; non-periodic infinitely-many breakpoints raise |
| `piecewise-cut-points-unresolved-r17` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | 7158, 6900, 6416 |  |
| `piecewise-cut-points-unresolved-r6` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23761 | Breakpoints linear in one parameter: its line cut where they meet, integral per interval: the triangle function |
| `piecewise-diff-ignores-jumps` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#112 | 24127 | diff of a piecewise adds J*dirac(x - p) at each jump J (breakpoints known, not periodic), as for heaviside |
| `piecewise-diff-point-piece-differentiated` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#145 | 11402 | A point piece's derivative is the limit of its neighbour's derivative. |
| `piecewise-equal-branches-not-merged` | SymJ | missing case | fixed | #105 | 27627 | A piece at a point merges with the piece taking over there. The integer-assumption narrowing (x ≥ 999 to x = 999) is not done; noted in notes/feature_audit.md. |
| `piecewise-equality-condition-not-substituted` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 16690 | equality conditions substituted into piece values |
| `piecewise-float-breakpoints` | SymJ | missing case | fixed | #89 | 29064 | piecewise conditions that evaluate to 1.0 are treated as true |
| `piecewise-fold-condition-poles-not-breakpoints` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#142 | 12095 | Poles of a condition are breakpoints; a pole point takes the next interval's formula. PR open. |
| `piecewise-fold-relation-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#146 | 11120 | piecewise_fold maps over a Relation's sides. |
| `piecewise-infeasible-conditions-not-pruned` | SymJ | missing case | fixed | #93 | 28323 | refine: contradictory And conditions are false |
| `piecewise-infeasible-conditions-not-pruned-r17` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | r17-extra-4 |  |
| `piecewise-infeasible-conditions-not-pruned-r4` | SymJ | missing case | fixed | #110 | 25930 | pieces with contradictory And conditions dropped |
| `piecewise-infeasible-conditions-not-pruned-r7` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 22066 | A linear form known ≠ 0 that the linear relations force to 0 makes the case impossible: 3 cases |
| `piecewise-integrate-condition-free-of-variable` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | 6747 |  |
| `piecewise-limit-needs-param-sign` | SymJ | missing case | fixed | #89 | 30249 | limit splits a piecewise whose conditions are free of x and returns a piecewise of limits |
| `piecewise-limit-needs-param-sign-r5` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24042 | The side sign of a piecewise argument at a symbolic point uses the assumptions: gives L |
| `piecewise-membership-condition-refused` | SymJ | missing case | open |  | 17112 | Still open (checked 2026-11-21): piecewise with x ∈ ℤ in a condition refused. |
| `piecewise-no-otherwise-returns-nothing` | SymJ | missing case | fixed | #91 | 28469, r2-extra-1 | intervals with no piece are left undefined (NaN inside, omitted at the end) |
| `piecewise-rewrite-misses-abs-in-product` | SymJ | missing case | fixed | jverzani/SymJ.jl#162 | 8849 |  |
| `piecewise-shadowed-branch-not-pruned` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 22148 | simplify folds piecewise functions in each condition symbol under the assumptions: 0 |
| `piecewise-sign-split-not-recognized-as-abs` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 19167 | r9-missing-cases: piecewise ±u split recognized as abs |
| `piecewise-subs-evaluates-unused-branch` | SymJ | crash | fixed | jverzani/SymJ.jl#118 | 23295 | subs into a piecewise piece by piece: a false condition's value is not substituted, the first true one ends: 0 |
| `poly-bound-variable-treated-free` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5602 |  |
| `poly-coefficient-type-mismatch` | SymJ | crash | fixed | jverzani/SymJ.jl#87 | extra-10 | gcd/lcm of Polys with symbolic coefficients (sqrt(2)) through cancel, monic |
| `poly-coefficients-not-zero-tested` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#123 | 20389 | Poly{Sym} drops coefficients that expand to 0: degree 0 |
| `poly-diff-higher-order-missing` | SymJ | crash | fixed | jverzani/SymJ.jl#124 | 19590 | diff(::Poly, x, n) |
| `poly-div-scalar-missing` | SymJ | crash | fixed | jverzani/SymJ.jl#87 | 29671 | Poly / c and Poly // c for a constant; by a polynomial in the generators: ArgumentError pointing to divrem/divexact |
| `poly-divrem-mixed-coefficient-stackoverflow` | SymJ | crash | fixed | jverzani/SymJ.jl#118 | 23828 | divrem of Polys with different coefficient types: both as Poly{Sym} over the generators of either (no Base fallback recursion) |
| `poly-generator-power-of-exponential-not-recognized` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5853 |  |
| `poly-generator-validation` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#86 | 29080, 29058 | Poly(ex, gens...) refuses a coefficient containing an expression generator's symbols (#29080). Numeric generators like sqrt(2) (#29058) left as is, as in SymPy. |
| `poly-missing-array-scalar-methods` | SymJ | crash | fixed | jverzani/SymJ.jl#124 | r8-extra-2 | conj, real, imag, adjoint (coefficientwise) and transpose for Poly: [p, p]' * [p, p] works |
| `poly-zero-no-generators-methoderror` | SymJ | crash | fixed | #108 | 25899 | empty generator list typed Vector{Sym} |
| `poly-zero-power-empty-reduce` | SymJ | crash | fixed | jverzani/SymJ.jl#143 | 12114 | 0^k returns 0; the 7-unknown system now solves (families checked numerically). PR open. |
| `polylog-numeric-noninteger-order` | SymJ | crash | fixed | jverzani/SymJ.jl#113 | r5-extra-13 | Numeric li for real non-integer order: Σzᵏ/kˢ for \|z\| ≤ 1/2, the μ = -log z expansion near 1 (\|μ\| < 4), ζ(s) or ∞ at 1 |
| `polylog-special-values-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 7132 |  |
| `power-log-base-change-not-simplified` | SymJ | missing case | fixed | #110 | 26006 | refine: b^(k log a/log b) = a^k for a, b > 0 |
| `power-log-base-change-not-simplified-r14` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 10410 | b^e with logs in e becomes exp(e log b), logs combined: 4^(3iπ/(2 log 2)) = -1. |
| `power-log-base-change-not-simplified-r17` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6887 |  |
| `power-of-product-integer-part-distributed` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5728 |  |
| `powsimp-common-exponent-quotient-not-combined` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 22546 | simplify_powers joins b^(-p) with the group of p (b > 0): ((a+b)/(a+2b))^(n+1) |
| `printing-rational-exponent-not-julia-roundtrip` | SymJ | missing case | won't fix |  | 5458 | Left for now by user decision 2026-11-16 (may print 1//3 later for round trips); feature_audit 23. |
| `qr-factorization-ldiv-missing` | SymJ | crash | fixed | jverzani/SymJ.jl#177 | 4006 |  |
| `radical-denesting-cube-roots` | SymJ | missing case | fixed | #106 | 27082 | cube roots of quadratic surds denested; simplify stays size-based |
| `radical-denesting-cube-roots-r4` | SymJ | missing case | fixed | #110 | 26857 | cube-root denesting with a non-cube norm |
| `radical-denesting-not-applied-in-subterms` | SymJ | missing case | fixed | #110 | 25915 | radicals then simplify_rational candidate |
| `radical-denominator-higher-roots-not-rationalized` | SymJ | missing case | feature gap |  | 14184 | Feature gap (2026-11-21), in the feature audit's roadmap: simplification: denominators with higher roots rationalized. |
| `radical-root-substitution-zero-not-recognized` | SymJ | missing case | won't fix |  | 6829 | notes/feature_audit.md item 20. |
| `radical-tower-construction-slow` | SymJ | hang | fixed | jverzani/SymJ.jl#113 | 24565 | Construction was instant; show re-printed subterms to sort them (exponential). Sort keys and degrees memoized per show |
| `rank-cases-minor-gcd-slow` | SymJ | hang | fixed | jverzani/SymJ.jl#165 | r17-extra-3 | Also a wrong answer: plain 22 where rank is 16/21/20 at x = 0, ±2, ±4; now Cases in 10 s. |
| `rank-cases-zero-pivot-not-split` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#164 | 7146 |  |
| `rank-generic-no-cases` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24613 | rank(A; generic = false): Cases where the r×r minors all vanish (k == 11 => 2); generic = true (default) as before |
| `rank-misses-zero-pivot-hyperbolic-identity` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#116 | 23531 | rref/rank/nullspace: hyperbolic functions as exponentials and exp(c u) as powers of w = exp(u/N), so their relations are seen: rank 2, one null vector |
| `rank-misses-zero-pivot-trig-identity` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#134 | 15585 | related trig arguments through exp(iu) in _exp_powers |
| `rank-misses-zero-pivot-trig-identity-r14` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#145 | 11905 | lu tests trig/hyperbolic pivots with simplify, so sin² + cos² - 1 is not a pivot. |
| `rank-misses-zero-pivot-trig-identity-r15` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#151 | 10120 | rref (rank, nullspace), inv's adjugate determinant, Cramer's determinant and Gauss-Jordan's nonzero test see zeros by trig/hyperbolic identities. |
| `rational-power-of-power-not-merged` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | r5-extra-6 | (b^p)^q = b^(pq) at construction for b a positive number or constant, p rational |
| `real-imag-acos-out-of-range` | SymJ | missing case | fixed | #93 | 28541 | real/imag of acos, asin of real arguments outside [-1, 1] |
| `real-imag-complex-constant-factor-not-split` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5166 |  |
| `real-imag-of-log-complex-not-evaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 17948 | r9-missing-cases: real/imag of log of complex constant = log\|z\|, angle(z) |
| `real-imag-real-factor-not-extracted` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6428 |  |
| `refine-angle-not-reduced-to-atan` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | 15945 | angle with Re > 0 is atan(Im/Re) |
| `refine-angle-not-reduced-to-atan-r13` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | 11988 | angle(x + iy) = atan(y, x) for real parts. PR open. |
| `refine-asin-of-real-argument-real-unknown` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 10598 | asin/acos of sin or cos of a real argument is real. |
| `refine-dirac-nonzero-argument` | SymJ | missing case | fixed | #110 | 26691 | dirac(u) = 0 when the sign of u is decided |
| `refine-exp-imaginary-periodicity` | SymJ | missing case | fixed | #89 | 29196 | simplify(...; assume) reduces exp(i*pi*w) by periodicity for integer w |
| `refine-factorial-binomial-sign-unknown` | SymJ | missing case | open | jverzani/SymJ.jl#140 | 13617, 13616 | Still open (checked 2026-11-21): signs of binomial and factorial fixed; factorial(n) ≥ n not decided. |
| `refine-factorial-lower-bound-one-unknown` | SymJ | missing case | fixed | jverzani/SymJ.jl#153 | 9115 | _sign_by_split: n! as 1 + w, w ≥ 0, for integer n ≥ 0. |
| `refine-floor-ceil-argument-bounds-unused` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 21429 | floor/ceil of an argument with decided bounds (n in -2:2): -1, 1 |
| `refine-infinity-times-signed-symbol` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6106 | Imaginary assumptions: feature gap. |
| `refine-interval-bounds-ignore-imaginary-unit` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#123 | 20453 | Signs by extrema: none for an expression with im or a non-real value (real parts were taken): imag(exp(iφ)) = sin φ, abs 1 |
| `refine-max-min-not-decided` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 16926 | numeric differences and assumptions drop dominated arguments |
| `refine-membership-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#128 | 19145 | refine of a membership via _holds; e strictly between consecutive integers is not an integer: false |
| `refine-mod-factorial-divisibility` | SymJ | missing case | fixed | jverzani/SymJ.jl#162 | 8677 |  |
| `refine-mod-integer-range-unused` | SymJ | missing case | fixed | jverzani/SymJ.jl#162 | 8486 |  |
| `refine-mod-nonnegative-unknown` | SymJ | missing case | fixed | jverzani/SymJ.jl#162 | 8581 |  |
| `refine-nan-comparison-crash` | SymJ | crash | fixed | #103 | 27440 | NaN has no sign mask; comparisons false as in IEEE |
| `refine-ne-assumption-unused` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24183 | x ≠ 0 in a tuple is Julia's Bool true and was dropped; a true fact is now an ArgumentError pointing to x ≶ 0 or @assume(x ≠ 0), both of which already work |
| `refine-negative-base-rational-power` | SymJ | missing case | fixed | #110 | 25724 | (c y^k)^q with y < 0 rewritten |
| `refine-negative-base-rational-power-r19` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | 4362, 2920 |  |
| `refine-or-relation-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#128 | 18389 | refine of FactAnd/FactOr: parts decided, combined (merging n = 0 \| n ≥ 1 into n ≥ 0 still open) |
| `refine-piecewise-ne-condition-undecided` | SymJ | missing case | fixed | #110 | 25522 | u ≠ v facts decide nonzero differences |
| `refine-real-implies-finite-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 17113 | a real side against ±∞ decided |
| `refine-real-part-coth-csch-sech-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 5998 |  |
| `refine-relation-skips-reduce-fact` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 14815 | Fixed on main (checked 2026-11-21). |
| `refine-trig-range-not-used-for-sign` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 21359, 21351 | Signs by splitting: sin/cos of a real argument as w ∈ [-1, 1]: \|cos x + 1\| = cos x + 1 |
| `relation-infinite-sides-subtracted-nan` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#127 | 19245 | A relation with an infinite side compared as it is (no ∞ - ∞): Inf ~ Inf true, Inf ≦ Inf true |
| `relation-nonzero-quotient-undecided-complex` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 11740 | A quotient with a nonzero constant numerator is never 0: simplify(1/x ~ 0) is false. |
| `relation-with-poly-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#131 | r10-extra-14 | relations with a Poly side |
| `relational-abs-max-not-decided` | SymJ | missing case | fixed | #103 | 27364 | signs by cases through abs, max, min |
| `replace-piecewise-pattern-keyerror` | SymJ | crash | fixed | jverzani/SymJ.jl#113 | 25328 | Showing a piecewise whose condition is not a relation (a pattern slot) prints it as is |
| `resultant-content-gcd-recursive-prs-variable-order` | SymJ | hang | fixed | jverzani/SymJ.jl#143 | 13126 | Quick constant-gcd test (a variable one side lacks, a numeric coefficient in it): >100 s to 0.03 s. PR open. |
| `resultant-float-multivariate-refused` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 16646 | exact over the floats' rationals, floats after |
| `reversed-sum-limits-inconsistent` | SymJ | wrong answer | fixed | #101 | r3-extra-1 | numeric reversed bounds by Karr's convention, as the closed forms; Sym(Inf) bounds as Inf |
| `rewrite-complex-to-polar-exp` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24522 | rewrite(z, exp) for a complex constant whose argument is a rational multiple of π: exp(iπ/3) |
| `rewrite-symbol-target-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#143 | r13-extra-4 | Non-rule arguments: ArgumentError pointing to rewrite(ex, f) / piecewise(ex). PR open. |
| `rewrite-tan-half-angle-not-quotient` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 16499 | sin/cos ratios as tan first |
| `risch-float-base-power-constants-overflow` | SymJ | missing case | open |  | r18-extra-6 | Feature gap (feature_audit 23): a float base 2.718...^(-2*pi*i*x) in RischIntegrate; exact exp works. Still open (checked 2026-11-21): a float base N(ℯ)^(-2πix): more than 8 transcendental constants; exact ℯ works. |
| `risch-transcendental-constants-algebraic-field-r7` | SymJ | missing case | open | jverzani/SymJ.jl#122 | 22075, 21066 | x-free non-numeric factors pulled out before integrating: #22075 gives a² + 1; #21066 not rechecked Still open (checked 2026-11-21): #22075 fixed; #21066's second integral: RischNotImplemented (transcendental constants). |
| `round-complex-numeric-unevaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 21513 | floor/ceil/trunc/round of a complex constant componentwise: 3 + 1.0im |
| `rref-exp-powers-sum-exponents-independent` | SymJ | crash | fixed | jverzani/SymJ.jl#131 | 17130 | exponents split into rational multiples of their atoms |
| `rref-float-pivot-not-normalized` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#175 | 4227 |  |
| `rsolve-check-zero-test-unexpanded-quotient` | SymJ | missing case | fixed |  | 6042 | Fixed on main (checked 2026-11-21). |
| `rsolve-constant-coefficient-slow-hyper` | SymJ | hang | fixed | #103 | 27902 | characteristic roots; held recurrence when they are large radicals |
| `rsolve-cubic-characteristic-roots-check-fails` | SymJ | missing case | open |  | 17990 | Still open (checked 2026-11-21): rsolve of a cubic characteristic: "the solution found does not satisfy the recurrence". |
| `rsolve-singular-leading-coefficient-hypergeometric-missed` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5507 |  |
| `rsolve-symbolic-coefficient-gcd` | SymJ | missing case | fixed | #106 | 27901 | n-free denominators not cleared (#27901) |
| `rsolve-symbolic-coefficient-gcd-r9` | SymJ | missing case | open |  | 18751 | Still open (checked 2026-11-21): rsolve with symbolic coefficients (cos θ): gcd requires numeric coefficients. |
| `series-asymptotic-special-functions` | SymJ | missing case | won't fix | #106 | 27109 | atan at ±∞ done (#106); erf at ∞ stays an error; e^g·series form planned as feature work (feature_audit.md, Gaps item 9) |
| `series-asymptotic-special-functions-r6` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23563 | Airy Ai, Bi asymptotic series at +∞ (e^∓ζ pieces); limit(airyai) at ∞ stays GruntzLimit's |
| `series-bessel-noninteger-order-taylor` | SymJ | wrong answer | fixed | #107 | 25678 | x^ν·(series in x²) for J, I, Y, K of non-integer order; integer Y, K raise |
| `series-branch-point-puiseux` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 25175 | acos, asin, acosh, atanh, acoth at ±1 through square roots and logs: Puiseux series |
| `series-expint-at-zero-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | 6725 |  |
| `series-expint-log-singularity` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | 6057 | The series at 0; the limit part is GruntzLimit's (unsupported expint). |
| `series-float-exponent-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6726 |  |
| `series-floor-ceil-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | 4156 |  |
| `series-generic-functions-truncate-or-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#139 | r12-extra-7 | Fixed on main (checked 2026-11-21): expand, factor and cancel of a series still a plain MethodError. |
| `series-infinity-exp-small-term-in-argument` | SymJ | missing case | feature gap |  | 21938 | Feature gap (2026-11-21), in the feature audit's roadmap: series and limits: an exponentially small term inside a function's argument at ∞. |
| `series-integer-power-full-expansion` | SymJ | hang | fixed | jverzani/SymJ.jl#165 | 6753 |  |
| `series-inverse-hyperbolic-log-singularity` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 10885 | asech and acsch at 0 through their log forms: log(2) - log(x) ∓ x²/4 - 3x⁴/32 + … |
| `series-inverse-trig-branch-point-wrong-branch` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#175 | 4004 | acos/asin at 1 from the right: -i·acosh, Julia's value above the cut (user decision: Julia semantics) |
| `series-inverse-trig-on-cut-wrong-side` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#180 | r20-extra-10 |  |
| `series-irrational-exponent-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#195 | 21358 | Fixed on main (checked 2026-11-21). |
| `series-lambertw-at-zero-removable-singularity` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 7259 |  |
| `series-laurent-gamma-digamma-pole` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23635 | Γ, ψ, ψ⁽ᵏ⁾ at a pole -m through the recurrences: Laurent series; γ = -digamma(1); GruntzLimit main gives y*gamma(y) → 1 |
| `series-log-log-leading-term` | SymJ | missing case | feature gap |  | 21227 | Feature gap (2026-11-21), in the feature audit's roadmap: series and limits: log(log x) at 0. |
| `series-piecewise-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#179 | 4315 |  |
| `series-polylog-branch-point-term-missing` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#112 | 25159 | li at z→1: Γ(1-s)μ^(s-1) + Σζ(s-k)(-μ)^k/k! (log term for integer s); li at 0 as Σzᵏ/kˢ (errored); zeta at negative integers, gamma at negative half-integers evaluated |
| `series-printing-reciprocal-term-unit-numerator` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | r19-extra-5 |  |
| `series-sym-infinity-point-not-recognized` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#130 | 16777, r10-extra-2 | Sym(±Inf) is the point at infinity |
| `series-symbolic-exponent-vanishing-base` | SymJ | missing case | won't fix |  | 7168, 7129, 6236, 6235 | Series with symbolic exponents: notes/feature_audit.md item 20. |
| `series-symbolic-exponent-vanishing-base-r19` | SymJ | missing case | won't fix |  | 4047 | feature_audit item 24 (symbolic exponents in series; as #7168) |
| `series-unknown-function-derivative-misleading-error` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | r14-extra-8 | Series of an unknown function's derivative: an error that says so. |
| `series-unknown-function-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | 3978 |  |
| `sign-complex-constant-not-evaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5763 |  |
| `sign-of-transcendental-unresolved` | SymJ | missing case | fixed | #96 | 29342 | sign from the range on an interval (critical points, ends) for one-symbol, non-periodic, smooth expressions; limit drops cases whose conditions are false |
| `sign-of-transcendental-unresolved-r7` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 21811 | A sum grouped by a symbol of known sign: h(cos A + 1) + p sin A > 0; the issue's expression is p |
| `sign-positive-factor-not-extracted` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 4708 |  |
| `simplify-abs-sign-product` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 20402, 19484 | simplify: \|u\| sign(u) = u (any u) |
| `simplify-algebraic-extension-rational-not-reduced` | SymJ | missing case | feature gap |  | 20236 | Feature gap (2026-11-21), in the feature audit's roadmap: simplification: rational functions over an algebraic extension. |
| `simplify-asinh-constant-not-rewritten-as-log` | SymJ | missing case | feature gap |  | 10429 | Feature gap (2026-11-21), in the feature audit's roadmap: simplification: asinh of a constant as a log, so powers combine. |
| `simplify-cases-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#113 | r5-extra-1 | simplify(::Cases) simplifies each case's value |
| `simplify-equation-common-terms` | SymJ | missing case | fixed | #93 | 28666 | simplify(relation) drops terms common to both sides |
| `simplify-equation-nonzero-factor-not-cancelled` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 21108 | simplify(::Relation) drops factors known nonzero from = and ≠ |
| `simplify-exp-hyperbolic-rewrite-wins` | SymJ | missing case | fixed | #110 | 26329 | expansion candidates with a common denominator |
| `simplify-exp-quotient-exponents-not-combined` | SymJ | missing case | won't fix |  | 7769 | Feature audit 19. |
| `simplify-exp-rewrite-radical-reciprocal-not-rationalized` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | 2958 |  |
| `simplify-fact-combination-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#118 | r6-extra-6 | simplify of FactAnd/FactOr: relations simplified, facts the others imply (And) or that imply another (Or) dropped: x ≪ 1 |
| `simplify-hyperbolic-compound-argument` | SymJ | missing case | fixed | #110 | r4-extra-1 | hyperbolic functions with exp powers rewritten to exp |
| `simplify-hypergeometric-contiguous-relations-missing` | SymJ | missing case | feature gap |  | 17702 | Feature gap (2026-11-21), in the feature audit's roadmap: simplification: contiguous relations of ₂F₁. |
| `simplify-identically-zero-denominator` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#120 | 21340 | simplify_trig: a denominator 0 by the identities (before cancel) leaves the expression as it is |
| `simplify-inverse-trig-complementary-sum` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 11446 | c asin(u) + c acos(u) = cπ/2 among simplify's candidates. |
| `simplify-inverse-trig-constant-sum-zero` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24948 | A constant with inverse trig near qπ (small denominator) is qπ when cos(e - qπ) simplifies to 1: 0 |
| `simplify-inverse-trig-constant-sum-zero-r10` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 16730 | inverse hyperbolic constants through logs |
| `simplify-large-power-radical-slow` | SymJ | hang | fixed | jverzani/SymJ.jl#113 | r5-extra-2 | A constant with large powers of sums of radicals is expanded first: 15.9 s to 0.16 s |
| `simplify-log-quotient-real-not-expanded` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23569 | log(n/d) = log n - log d for n > 0, d ≥ 0 (sums over their denominator): 0 for real x |
| `simplify-negated-quotient-sign-not-normalized` | SymJ | missing case | won't fix |  | 4592 | Same canonical-form question as negated-quotient-not-auto-cancelled; left by user decision 2026-11-16. |
| `simplify-not-idempotent-sqrt-product` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | r5-extra-7 | Constant simplify repeats until it stops shrinking (up to 3 passes): 1 in one call |
| `simplify-powers-relation-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#128 | r9-extra-11 | simplify_* on a relation, each side; an array of Polys times a Poly entry by entry |
| `simplify-radical-fractions-not-combined` | SymJ | missing case | won't fix |  | 6651 | Radical towers: notes/feature_audit.md item 20. |
| `simplify-rational-content-not-cancelled-in-product` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 17233 | sum contents that cancel the coefficient |
| `simplify-re-im-recombination` | SymJ | missing case | fixed | #110 | 25564 | real(a) + i imag(a) recombined |
| `simplify-re-im-recombination-r14` | SymJ | missing case | feature gap |  | 10822 | Feature gap (2026-11-21), in the feature audit's roadmap: simplification: real and imaginary parts recombined. |
| `simplify-re-im-recombination-r9` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 18271 | r9-missing-cases: conj recombined as real − i·imag candidate |
| `simplify-reciprocal-trig-power-not-rewritten` | SymJ | missing case | fixed | jverzani/SymJ.jl#166, jverzani/SymJ.jl#174 | 6303 | Round 18 (#4621): 1/cos(x)^2 - 1 to tan(x)^2, numerator and denominator reduced each its own way (#174). |
| `simplify-reciprocal-trig-power-not-rewritten-r18` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 4621 | Fixed on main (checked 2026-11-21). |
| `simplify-relation-not-collected` | SymJ | missing case | fixed | #103 | 27754, 27623 | one side without numeric content when shorter |
| `simplify-relation-trivial-not-decided` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | r13-extra-10 | Constant difference after simplify decided. PR open. |
| `simplify-size-ignores-integer-digits` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23210 | simplify's size counts long exact numbers (one more per 4 digits past 6): 30 log(5) stays |
| `simplify-sqrt-conjugate-quotient-not-cancelled` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23641 | simplify candidate: denominator linear in a symbolic √ times its conjugate, cancel with half powers as variables: √c/√(cx²+bx+a) |
| `simplify-trig-constant-sum-inside-function-arg` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | 15129 | constant function arguments simplified first |
| `simplify-trig-missed-in-larger-product` | SymJ | missing case | fixed | #110 | 26835 | circular factors of one argument simplified apart |
| `simplify-trig-missed-in-larger-product-r9` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 18231 | second simplify pass on results ≤ 40 nodes |
| `simplify-trig-missed-in-larger-sum` | SymJ | missing case | fixed | #103 | 27356 | a sum's trig terms simplified apart |
| `simplify-trig-multiple-angle-not-collapsed` | SymJ | missing case | fixed |  | 20091 | Fixed on main (checked 2026-11-21). |
| `simplify-trig-rational-pi-shift-power-sum` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6995 |  |
| `simplify-trig-rational-pi-sum-zero` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | r6-extra-9 | Constants in sin/cos of rational multiples of π: roots of unity reduced mod Φ_N (im included): 0 |
| `simplify-zero-numerator-not-reexpanded` | SymJ | missing case | fixed | jverzani/SymJ.jl#182 | r20-extra-9 | a simplify candidate through s = u^(1/q), reduced by s^q = u; expand unchanged Fixed on main (checked 2026-11-21). |
| `sinc-not-rewritten` | SymJ | missing case | fixed | #89 | 30545 | sinc rewritten to sin(pi u)/(pi u) in limits and series |
| `solve-abs-branch-identity-reported-global` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#120 | 22259, r7-extra-1 | An \|u\| branch that is an identity: the set where its sign holds, less undefined points: (-∞, 0], (0, ∞), [1, ∞) |
| `solve-abs-not-inverted` | SymJ | missing case | fixed | #93 | 28742 | real unknown: split by sign of the argument, by cases; complex x stays an error (a curve) |
| `solve-abs-not-inverted-r4` | SymJ | missing case | fixed | #110 | 25779 | \|b^q\| = r inverted for real b, q |
| `solve-abs-not-inverted-r5` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24189 | abs in systems: sign cases, unknowns real, numeric solutions kept in their case (single equations already worked with x real) |
| `solve-abs-quotient-singular-point-kept` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#130 | 17438 | abs split drops roots where the equation is undefined |
| `solve-abs-system-drops-parametric-solutions` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#142 | 13306 | Families/parametric solutions refused with an error; one equation in one unknown gives the scalar Cases. PR open. |
| `solve-abs-zero-treated-as-curve` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | r11-extra-3 | \|u\| = 0 gives u = 0 over C |
| `solve-angle-not-inverted` | SymJ | missing case | open |  | 19506 | Still open (checked 2026-11-21): solve(angle(x + im) ~ 0, x): cannot invert angle. |
| `solve-atan2-not-inverted` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6977 |  |
| `solve-backsub-cancel-radical-blowup` | SymJ | hang | fixed | jverzani/SymJ.jl#143 | 12345 | A basis element c x + c0 with numeric c gives x directly, unexpanded: minutes to 2.5 s. PR open. |
| `solve-bool-equation-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#165 | r17-extra-11 |  |
| `solve-call-unknown-other-occurrences-of-argument` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#142 | 12129 | expand_trig in terms of the unknown call; refused when the variable remains (sin(2a) for sin(a)). PR open. |
| `solve-case-split-inequality-cancel-blowup` | SymJ | hang | fixed | jverzani/SymJ.jl#165 | 6818 | Refused (ArgumentError) in ~26 s: conditions on Cardano roots too large to split on. Solving it would need the cubic's discriminant as the condition. |
| `solve-cases-error-object-as-value` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | r13-extra-5 | The capped case's error now says it is left unsolved and how to settle it; still an error value. PR open. |
| `solve-cases-numeric-eval-free-parameter` | SymJ | missing case | fixed | jverzani/SymJ.jl#173 | 4671 | Same cause as solve-weierstrass-numeric-check-on-parameter (root scale). |
| `solve-composite-generator-not-substituted` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24418 | x only through h (s² + x² in log(s² + x²)): solved for v = h, then h = v; four Lambert W roots |
| `solve-conj-not-inverted` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 20305 | conj(t) = r: t = conj(r) |
| `solve-cot-not-inverted` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 19115 | r9-missing-cases: cot inverted with periods |
| `solve-csc-sec-not-inverted` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 17536 | csc/sec through sin/cos of 1/r |
| `solve-cube-root-not-isolated-and-cubed` | SymJ | missing case | fixed | jverzani/SymJ.jl#163 | 8789 |  |
| `solve-cyclotomic-roots-numeric` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 25192 | A cyclotomic factor Φ_N: roots cos(2πk/N) + i sin(2πk/N); x^5 + 1 too, before Ferrari |
| `solve-cyclotomic-roots-numeric-r7` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | r7-extra-9 | An irreducible factor equal to cos(2π/N)'s minimal polynomial (from Φ_N): roots cos(2πk/N); the 11×11 eigenvalues exact |
| `solve-degenerate-irrational-coefficient-missed` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#112 | 24344, r5-extra-3 | Degenerate-case factors may have constant coefficients (a - √3 b, a - π); the issue's equation gets its √3 a0/8 + r2(2-Z) = 0 case |
| `solve-degenerate-nonlinear-factor-missed` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#116 | 23459 | Cases on the roots of a²+ab+b² (per-root alternatives, #115); solving for derivatives keeps their names inside Cases (a Cases is a vector: mapped by case) |
| `solve-degenerate-nonpolynomial-coefficient-missed` | SymJ | wrong answer | won't fix | jverzani/SymJ.jl#120 | r7-extra-5 | Won't fix (2026-11-21): generic parameter values by default (2026-11-14). |u\| coefficients give u's factors as cases (\|x\| y = 0: x = 0); sin(x) y stays generic (infinitely many zeros) |
| `solve-degenerate-nonpolynomial-coefficient-missed-r11` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | r11-extra-4 | y^x = z: cases y = 1, y = 0 |
| `solve-exp-linear-arguments-not-common-generator` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 20199 | Exponents with constant parts: b^(x+1) = b·b^x in the common generator |
| `solve-exp-sum-common-ratio-substitution` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24138 | Σ cᵢ bᵢ^(αx+β) with prime exponents on a line: polynomial in t = r^x: log((1+√5)/2)/log(3/2) |
| `solve-exp-sum-common-ratio-substitution-r8` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 20709 | Constant bases that are small powers of one another (reciprocals); small integer powers recognised: ±2 |
| `solve-float-exponent-not-inverted` | SymJ | missing case | fixed | #93 | 28772 | principal root for float or symbolic exponents |
| `solve-float-exponent-not-inverted-r19` | SymJ | crash | fixed | jverzani/SymJ.jl#177 | 4260 | float exponents with other symbols via shortest rationals |
| `solve-float-exponent-not-inverted-r6` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23284, 22584 | Floats where exact methods need exact numbers: shortest rationals, solved exactly, results floats: [77064.73] |
| `solve-float-quadratic-radicand-scaled` | SymJ | missing case | fixed | jverzani/SymJ.jl#179 | 2993 |  |
| `solve-float-result-not-refloated-in-arguments` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | r14-extra-6 | Floating an e^u result floats u's numbers: exp(3.33…t). |
| `solve-float-root-integer-domain-dropped` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#116 | 23731 | A float with an integer value is an integer for the assumptions: [1.0], [-2.0, 2.0] |
| `solve-float-system-exactified` | SymJ | wrong answer | fixed | #107 | 26223 | Gauss-Jordan float tolerance and largest pivot |
| `solve-float-system-exactified-r5` | SymJ | missing case | fixed | jverzani/SymJ.jl#115 | 24193 | Floats among parameters as opaque constants, put back after (conditions of them decided by value); results floats again; a degenerate factor with several roots gives a case per root instead of 'not s… |
| `solve-float-system-residue-not-chopped` | SymJ | wrong answer | won't fix |  | r17-extra-8 | Minor; notes/feature_audit.md item 20. |
| `solve-floor-ceil-mixed-not-supported` | SymJ | missing case | fixed | jverzani/SymJ.jl#159 | 9118 | floor(x) with ceil(x): non-integer x via ceil = floor + 1 (intervals without integer ends), integer x as points. ceil(x) - floor(x) = 1 or 0 still refused (infinite point sets). |
| `solve-floor-not-inverted` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 25001 | c floor(u) + d = 0 (or ceil) gives the solution set, as an inequality: [-1, 1) |
| `solve-floor-not-inverted-r16` | SymJ | missing case | fixed | jverzani/SymJ.jl#163 | 8472 |  |
| `solve-generic-default-cases-not-dropped` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | r17-extra-1 | Ferrari's degenerate points now generic (decision 2); 0.4 s. |
| `solve-generic-default-cases-not-dropped-r18` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | r18-extra-2, r18-extra-3 |  |
| `solve-heaviside-not-inverted` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 18939 | r9-missing-cases: r = 1/2 gives 0; r = 0 or 1 errors (holds on a set) |
| `solve-hyperbolic-not-rewritten-as-exp` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 19489 | Several hyperbolic functions (and exponentials): as exponentials, roots checked: log(2)/8 |
| `solve-identity-ignores-excluded-points` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | r6-extra-8 | The every-value error names the points where the equation is not defined: except x = 1 |
| `solve-implicit-ode-log-high-degree-hang` | SymJ | hang | fixed | jverzani/SymJ.jl#118 | 22462 | Equations in logs: coefficients over their rational gcd (72/5 log(a) - 72/5 log(b) is log(a/b) times 72/5), not a degree-72 polynomial: the ODE in 0.06 s |
| `solve-inexact-retry-infinite-recursion` | SymJ | hang | fixed | jverzani/SymJ.jl#152 | 9037 | The exact-numbers retry stops when it changes nothing; an inequality with an Inf/NaN coefficient gets the same ArgumentError as an equation. Inexact inequalities like 1.5x > 1 remain unsupported (mis… |
| `solve-inexact-symbolic-coefficients` | SymJ | missing case | fixed | #103 | 27712, 27528 | floats solved as the rationals they stand for, solutions' numbers floats again |
| `solve-inexact-symbolic-coefficients-r4` | SymJ | wrong answer | fixed | #107 | 26820 | float equations with parameters split as exact rationals |
| `solve-infinite-coefficient-misreported-complex` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 11015 | An Inf or NaN coefficient gets its own error message, not 'complex floating point coefficients'. |
| `solve-integer-domain-condition-not-decided` | SymJ | missing case | fixed | jverzani/SymJ.jl#153 | 9152 | _holds_ decides a constant's membership through `in`. |
| `solve-invdigamma-not-used` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6605 |  |
| `solve-lambertw-log-form-unsupported` | SymJ | missing case | fixed | #110 | 26461 | log(x) with x via x = e^y; k-th power Lambert W; numeric W fallback |
| `solve-lambertw-log-form-unsupported-r6` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23253 | A polynomial of degree ≥ 2 in log(…) solved for it first, then each log = root (Lambert W) |
| `solve-lambertw-not-inverted` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | 12783 | W(t) = r: t = r e^r; none for real r < -1. PR open. |
| `solve-linear-pythagorean-denominator-unsimplified` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23140 | Linear solutions: denominators reduced by simplify_trig when smaller (cos² + sin² = 1) |
| `solve-log-branch-spurious-roots` | SymJ | wrong answer | won't fix |  | r9-extra-8 | Design: each root is valid for some real parameters (root 2 for x < 0, root 3 for x > 0, root 1 at C1 = 5.3, x = -0.61); roots valid somewhere are kept without conditions |
| `solve-log-symbolic-base-not-inverted` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6176 |  |
| `solve-logical-combination-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#124 | 20080, r8-extra-6 | solve(FactAnd/FactOr, x): the parts' sets intersected or joined: [0, 1] |
| `solve-max-min-not-rewritten-piecewise` | SymJ | missing case | fixed | jverzani/SymJ.jl#153 | 10158 | An equation with max/min of the unknown is solved as piecewise. |
| `solve-mixed-polynomial-trig-unsupported` | SymJ | missing case | feature gap |  | 15989 | Feature gap (2026-11-21), in the feature audit's roadmap: solving: mixed polynomial–trigonometric equations (x = sin x); fsolve finds the roots. |
| `solve-mixed-trig-not-rewritten` | SymJ | missing case | fixed | #89 | 29745 | solve with sin and cos of one argument via t = tan(w/2), solutions in [0, 2pi) verified |
| `solve-negative-base-exponential-roots-dropped` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#123 | 19816 | b^t = r, b < 0, t real: integer t only, \|b\|^t = \|r\| with the parity's sign: [3] |
| `solve-nested-sqrt-generator-not-flattened` | SymJ | missing case | fixed | jverzani/SymJ.jl#162 | 8591 |  |
| `solve-nested-trig-outer-roots-one-period` | SymJ | missing case | feature gap |  | r9-extra-9 | Feature gap (2026-11-21), in the feature audit's roadmap: solving: all solutions (sin(sin x) = 0 takes outer roots 0, π, 2π); for the solveset-like function. |
| `solve-nonpolynomial-system-refused` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24868, 24269 | √R in a system: a new unknown w, w² = R, principal root kept; the Lagrange system gives its 2 points |
| `solve-nonpolynomial-system-refused-r10` | SymJ | missing case | feature gap |  | 16618 | Feature gap (2026-11-21), in the feature audit's roadmap: solving: transcendental systems. |
| `solve-nonpolynomial-system-refused-r6` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23662, 23318, 22723, 22717, 22413 | log(A) = log(B) as A = B (A ≠ 0); elimination of linear unknowns; polynomial equations first then the rest; floats via exact rationals |
| `solve-nonpolynomial-system-refused-r7` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 21908 | E^u factors dropped from product equations: (-2, 0), (0, 0) |
| `solve-nonpolynomial-system-refused-r8` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 19587 | c₁Πb^e = c₂Πb^e in positive numbers: logs, prime powers apart: m = -4, n = 1 |
| `solve-nonpolynomial-system-refused-r9` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 17941 | r9-missing-cases: systems polynomial in exponential generators |
| `solve-nonpolynomial-system-triangular` | SymJ | missing case | fixed | #93 | 28064 | non-polynomial systems solved one unknown at a time when an equation has a single unknown |
| `solve-parametric-cubic-slow` | SymJ | hang | fixed | #108 | 25920 | numeric nonzero test before cancel; cube-root roots skip parameter-factor search; 220 s -> 0.6 s |
| `solve-parametric-linear-system-hang` | SymJ | hang | fixed | jverzani/SymJ.jl#121 | 21694 | Degenerate cases whose parameter values are large radicals are left unsupported (an error leaf), and √(-w²) rewriting only for small radicands (a regression from #119): Cases in 0.6 s |
| `solve-parametric-polynomial-system-hang` | SymJ | hang | fixed | jverzani/SymJ.jl#118 | 23604 | Parametric systems checked at sample values first: an irreducible eliminant of degree > 4 (not decomposable), or 3-4 with >= 6 parameters, is a clear ArgumentError (0.7 s); also a size guard in Buchb… |
| `solve-piecewise-inequality-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 19718 | Inequalities in piecewise, abs, max, min, piece by piece: (10, ∞) |
| `solve-piecewise-nested-not-folded` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 18359 | r9-missing-cases: nested piecewise folded before solving |
| `solve-piecewise-unsupported` | SymJ | missing case | fixed | #103 | 27037 | piece by piece, roots kept where the condition holds |
| `solve-piecewise-zero-interval-dropped` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#120 | 21079 | A piece that holds identically gives its interval: a SymSet with the point roots, (-∞, -1) ∪ {0} |
| `solve-power-equals-one-base-cases` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | 15731 | base 1, exponent 0 (b nonzero), base -1 at even e |
| `solve-power-of-exponential-not-merged` | SymJ | missing case | open |  | 17452 | Still open (checked 2026-11-21): (7^x)^x = 2: cannot invert. |
| `solve-power-of-linear-expanded` | SymJ | hang | fixed | jverzani/SymJ.jl#131 | 17867, r10-extra-5 | large powers of polynomials solved by their bases; _solve_floor checks first |
| `solve-radical-parameter-case-dropped` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#130 | r10-extra-4 | a residual c·u^p vanishing for some parameters is a case |
| `solve-radical-parameter-case-dropped-r15` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#155 | 10398 | Decision 3 merged (#155): conditions for real parameters; complex parameters keep the root (the user's choice). #10398's own residual has no \|u\| to split, so it stays unconditional. |
| `solve-radical-spurious-root` | SymJ | wrong answer | fixed | #101 | 27381 | with parameters left, candidates checked at sampled parameter values; dropped when the equation fails at all |
| `solve-radical-spurious-root-r14` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#145 | 11534 | Root checks split on the sign of real parameters under \|u\|; satisfiability sampled within the assumption sets. |
| `solve-radical-spurious-root-r20` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#180 | r20-extra-2 | A·B ≤ 0 condition per root with real parameters; conditions on cos(a) leave some redundant cases |
| `solve-radical-spurious-root-r4` | SymJ | wrong answer | won't fix |  | r4-extra-2 | design: a root valid for some parameter values is kept (as sqrt(x) = a gives a^2); conditions on y would need Cases on ranges |
| `solve-radical-squared-identity-unchecked` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#160 | 8716, 8545 | Squared identity: √t = -A/B as a set (-A/B ≥ 0) for real x, refused over ℂ. |
| `solve-radical-system-hang` | SymJ | hang | fixed | jverzani/SymJ.jl#128 | 17960 | Buchberger gives up after 120 s with a clear error (reductions with growing coefficients); passed through the radical path |
| `solve-rational-exponents-drops-roots` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#86 | 29180 | x = u^q kept unexpanded (the expanded Cardano power lost the root numerically); Aberth starts from Fujiwara's bound (degree 181 gave NaN). All three #29180 cases now solve, x^(3/2)+x^(1/5)=3 included… |
| `solve-real-domain-root-membership-undecided` | SymJ | missing case | open |  | 11379 | Still open (checked 2026-11-21): r³ = V over ℝ, V < 0: four subcases on undecided memberships instead of -(-V)^(1/3). |
| `solve-reciprocal-trig-poles-dropped` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#142 | 12527 | tan/cot/sec/csc in a denominator: solved in sin and cos. PR open. |
| `solve-root-at-parametric-pole-no-case` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#123 | r8-extra-10 | A denominator's value at a root with radicals of a parameter solved for it: the case a = 0 => [2] |
| `solve-root-check-cancel-radical-blowup` | SymJ | hang | fixed | jverzani/SymJ.jl#139 | 13625 | Fixed on main (checked 2026-11-21). |
| `solve-roots-unsimplified-duplicates` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | r6-extra-2 | Roots through a power substitution simplified when smaller and deduplicated: [-1] |
| `solve-several-exp-bases-not-logged` | SymJ | missing case | fixed | #103 | 27330 | products of constant-base powers by logarithms, constant logs into primes |
| `solve-several-exp-bases-not-logged-r10` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 17604 | two-term difference as two sides |
| `solve-several-exp-bases-not-logged-r17` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6819 | Symbolic bases known positive under the assumptions. |
| `solve-several-radicals-not-supported` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 21167, 21026 | Radical systems with cube roots (real cube root or principal root checked): [0], [1, 2, 10] |
| `solve-sides-cancel-before-domain-check` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#151 | 9522 | Each side's domain checked, not lhs - rhs after cancellation. |
| `solve-sides-cancel-before-domain-check-r16` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#160 | 8715 | Each side's poles removed from inequality sets. |
| `solve-sign-not-inverted` | SymJ | missing case | fixed | #103 | 27322 | sign(u) = 0 is u = 0 |
| `solve-small-polynomial-system-slow` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | r6-extra-3 | _extract_power memoized for big integers (trial division repeated in cancel): 41 s → 1 s |
| `solve-solutions-not-simplified` | SymJ | missing case | open | jverzani/SymJ.jl#119 | 23889 | The a² + b² = 0 case: √(-a²) as ±i a in case substitutions, lifted solve retried with the new variable first; forms not tidied further Still open (checked 2026-11-21): a² + b² = 0 case done; the circle–line solutions are not tidied. |
| `solve-solutions-not-simplified-r12` | SymJ | missing case | open |  | 14244 | Still open (checked 2026-11-21): x⁴ - 10x² + 1: ±√(5 ± 2√6), not denested to ±√2 ± √3. |
| `solve-spurious-root-pi` | SymJ | wrong answer | fixed | #91 | 28213 | x = pi from tan(x/2) kept only where the equation provably holds |
| `solve-sum-of-squares-not-split` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24207 | Σ cᵢ uᵢ^(2k) = 0 with x real: one uᵢ solved, roots kept where all vanish: [5] |
| `solve-symbolic-exponent-power-vs-base-refused` | SymJ | missing case | feature gap |  | 10864 | Feature gap (2026-11-21), in the feature audit's roadmap: solving: x^(yz) = x. |
| `solve-symbolic-exponent-roots-dropped` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#120 | 21014 | x^n = r with n an assumed integer: refused with the roots exp(2πik/n) named, not the one principal root (a non-integer symbolic exponent keeps the principal root, by design) |
| `solve-symbolic-exponent-roots-dropped-r13` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#142 | 12774 | User decision: Cases on a (a = 0: every x or none), principal root for a ≠ 0; other roots documented, not listed. PR open. |
| `solve-system-abs-complex-treated-real` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#127 | r9-extra-6 | Systems over ℂ: abs/real/imag/conj through real and imaginary parts (a family); assume = (I ∈ ℝ, …) for systems narrows to ℝ (±2) and filters solutions |
| `solve-system-bool-equation-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#146 | 11659 | true among a system's equations is dropped; false makes the system unsolvable. |
| `solve-system-domain-check-skips-parametric` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#142 | 12932, 12546 | Families narrowed by the unknowns' sets (one free unknown) and integer domains (affine, residues over one period). PR open. |
| `solve-system-drops-degenerate-components` | SymJ | wrong answer | fixed | #101 | 27458, r3-extra-7 | components where a solution denominator or a level's common factor vanishes are solved too (memoized, cap 256) |
| `solve-system-drops-parameter-cases` | SymJ | missing case | fixed | #90 | 30224 | solve returns Cases where a coefficient with parameters can vanish, for single equations and systems (a*x + b*y = 0, x^2 + y^2 = 1 keeps a = 0) |
| `solve-system-drops-parameter-cases-r8` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#123 | 19564 | An equation free of the unknowns is a condition on the parameters, factor by factor: k = 5, k = -3 cases |
| `solve-system-drops-parameter-cases-r9` | SymJ | missing case | won't fix |  | r9-extra-5 | Won't fix (2026-11-21): generic parameter values by default (2026-11-14). |
| `solve-system-factors-not-split` | SymJ | hang | fixed | jverzani/SymJ.jl#113 | 24500 | A polynomial equation that factors in the unknowns splits the system before the Gröbner basis: 1.2 s |
| `solve-system-hidden-zero-equation-identity-error` | SymJ | missing case | open |  | 19038 | Still open (checked 2026-11-21): a five-equation symmetric system: a hidden zero not recognized. |
| `solve-system-integer-domain-ignored` | SymJ | wrong answer | open | jverzani/SymJ.jl#138 | 14170 | Still open (checked 2026-11-21): no integer solutions: fixed; x, y ∈ ℤ with solutions: no Diophantine parametrization. |
| `solve-system-ne-not-checked-after-equalities` | SymJ | missing case | fixed | jverzani/SymJ.jl#153 | 10408 | Equations solved first, ≠ checked on each solution; undecided ≠ keeps the error. Julia's Bool a ≠ b in a system is still dropped as true (feature audit 18). |
| `solve-system-nonvanishing-factor-not-removed` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 4509 |  |
| `solve-system-numeric-backsub-unverified` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#138 | 14799 | Fixed on main (checked 2026-11-21): solve is exact only; fsolve's 22 roots have residuals ~1e-15. |
| `solve-system-parameter-cases-redundant-infeasible` | SymJ | missing case | won't fix |  | r16-extra-1 | Feature audit 19. |
| `solve-system-re-im-of-unknowns-returns-empty` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#127 | 19349 | Unknowns as a + bi, real and imaginary parts of each equation: 1 + 2i; x y = 220, re x = 2 a family |
| `solve-system-redundant-components` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23110 | Solutions within a family (free unknowns taking their values) dropped: 2 families |
| `solve-system-variable-order-blowup` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23322 | Unknowns with an equation of their own last in the lex order: x = ±2√2 ± √3 |
| `solve-tiny-float-coefficient-dropped` | SymJ | wrong answer | fixed | #107 | 25792 | _holds_somewhere uses exact sample values |
| `solve-transcendental-returns-empty` | SymJ | wrong answer | fixed | #107 | 26253 | x outside sin/cos in the tan(x/2) path raises 'not supported' |
| `solve-transcendental-returns-empty-r17` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#164 | 6613 | Refused (must be polynomial) instead of []; solving it would be a feature. |
| `solve-trig-half-angle-roots-unsimplified` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 19102 | r9-missing-cases: 2atan(r) tidied to atan(q) (+π) |
| `solve-trig-interval-assumption-no-periods` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#123 | 20798 | sin/cos/tan with a linear argument and x in a bounded assumed interval: every root there (four for cos 2x = 1/2 on [0, 2π]) |
| `solve-trig-interval-assumption-no-periods-r10` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#130 | 17579 | nonlinear argument: its range (ends, critical points); unbounded is an error |
| `solve-trig-multiple-angles-not-reduced` | SymJ | missing case | feature gap |  | 19488 | Feature gap (2026-11-21), in the feature audit's roadmap: solving: different multiples of one angle (tan(x/6) = cos(x/2)). |
| `solve-trig-multiple-angles-not-reduced-r9` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 18315 | r9-missing-cases: multiple angles reduced to one argument, then Weierstrass |
| `solve-trig-parameter-cases-unchecked` | SymJ | missing case | fixed | #91, #95 | r2-extra-5 | excluded points by cases (#91); b + c = 0 case with x = π, and leading coefficients in one function of x (#95) |
| `solve-weierstrass-complex-coefficient-atan-pole` | SymJ | crash | fixed | jverzani/SymJ.jl#146 | 10809 | The tan(x/2) substitution skips t = ±i; cos x + i sin x = 1 gives [0]. |
| `solve-weierstrass-numeric-check-on-parameter` | SymJ | crash | fixed | jverzani/SymJ.jl#173 | 4467 |  |
| `solve-zero-base-power-not-inverted` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 20205 | b^e = 0 through b = 0, roots checked: [1] |
| `solve-zero-base-power-not-inverted-r14` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#145 | 11536 | 0^t = 1 gives t = 0. |
| `solve-zero-base-power-not-inverted-r17` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6745 | Under y > 0; without assumptions still []. |
| `special-function-derivative-missing` | SymJ | missing case | fixed | #103 | 26958 | spherical Bessel j and y derivatives |
| `special-function-derivative-missing-r17` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6007 | factorial's derivative; digamma recurrence in simplify not done. |
| `special-function-derivative-missing-r5` | SymJ | missing case | fixed | jverzani/SymJ.jl#114 | 24318 | diff of real, imag, conj (the part of u') and angle ((Re u Im u' - Im u Re u')/\|u\|²) in a real variable |
| `special-function-infinity-evaluated-as-float` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6602 | Any function at ±∞ whose limit is an integer or pπ/q. |
| `sqrt-denest-multiple-surds` | SymJ | missing case | feature gap |  | 21599 | Feature gap (2026-11-21), in the feature audit's roadmap: simplification: denesting square roots of several surds. |
| `sqrt-negative-algebraic-constant-not-extracted` | SymJ | missing case | fixed | jverzani/SymJ.jl#136 | r11-extra-7 | sqrt of a negative real constant sum is i sqrt(-u) |
| `sqrt-negative-transcendental-constant-not-extracted` | SymJ | missing case | fixed | jverzani/SymJ.jl#178 | r19-extra-7 |  |
| `sqrt-over-radicand-product-not-reduced` | SymJ | missing case | open |  | 21571 | Still open (checked 2026-11-21): √(xy)/(xy) not reduced to 1/√(xy). |
| `sqrt-product-nested-radicand-not-combined` | SymJ | missing case | feature gap |  | 20251 | Feature gap (2026-11-21), in the feature audit's roadmap: simplification: √(s - a)√(s + a) with s = √(a² + b²) to b. |
| `sqrt-product-positive-factor-not-split` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 17634 | simplify candidate: positive factors out of radicals |
| `sqrt-square-factor-above-trial-bound` | SymJ | missing case | fixed | #109 | 26903, 26317, r4-extra-4 | trial division up to 2^17 for q-th power factors |
| `stdlib-linalg-numeric-only` | SymJ | crash | fixed | #98 | extra-6, extra-7 | sqrt, log, non-integer powers through the Jordan form; pinv from the rank factorization (Matrix{Sym}, no SymMatrix) |
| `stdlib-linalg-numeric-only-r12` | SymJ | crash | fixed | jverzani/SymJ.jl#139 | 13985 | Fixed on main (checked 2026-11-21). |
| `stdlib-linalg-numeric-only-r13` | SymJ | crash | fixed | jverzani/SymJ.jl#143 | r13-extra-8 | cholesky for Matrix{Sym}. PR open. |
| `subs-captures-bound-variable` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#86 | 30788 | Sum, Product and Integral(f, v, a, b) keep v bound; renamed to v′ when a substituted value contains it. |
| `subs-captures-bound-variable-r12` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#138 | 13795 | Fixed on main (checked 2026-11-21). |
| `subs-exact-power-blowup` | SymJ | hang | fixed | #103 | r3-extra-5 | coefficient's huge power stays a power; N rechecks overflow in BigFloat |
| `subs-exact-power-blowup-r14` | SymJ | hang | fixed | jverzani/SymJ.jl#146 | r14-extra-4 | Huge fractional powers stay powers; N evaluates them through logs. |
| `subs-function-call-inside-derivative` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#120 | r7-extra-3 | f(x) => y (a symbol free of x) keeps f'(x) inside derivatives; f(x) => an expression in x still differentiates it |
| `subs-function-replacement-unsupported` | SymJ | crash | fixed | jverzani/SymJ.jl#118 | r6-extra-10 | subs(ex, f => g) with f an unknown function: each call f(args) by g(args), held derivatives then taken: sin(g(x)) |
| `subs-indefinite-integral-variable-replaced` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#164 | 6928 |  |
| `subs-number-target-hits-stored-exponents` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#120 | r7-extra-8 | Number keys leave numeric exponents alone (structure, as coefficients): subs(x*y, 1 => 2) = x*y |
| `subs-power-target-not-matched-in-other-base` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 10829 | For a target b^u with integer b, powers (b^m)^w are written b^(mw) first: subs(4^x - 3·2^x + 2, 2^x => y) = y² - 3y + 2. |
| `subs-power-target-not-matched-in-other-base-r17` | SymJ | missing case | won't fix |  | 6976 | Tried and reverted: radical solving relies on it; notes/feature_audit.md item 20. |
| `subs-product-with-numeric-coefficient-not-matched` | SymJ | missing case | fixed | jverzani/SymJ.jl#125, jverzani/SymJ.jl#174 | 19422 | Keys c·t matched in products (integer multiples) and as sum terms: yz, z + 4 (#125). Round 18 (#5217): a key c·M, M a power or product, as M => v/c: subs(2x^2 + y, 2x^2 => s) = s + y (#174). |
| `subs-product-with-numeric-coefficient-not-matched-r10` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 17225 | scaled sum factors match product keys |
| `subs-product-with-numeric-coefficient-not-matched-r18` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 5217 | Fixed on main (checked 2026-11-21). |
| `subs-reciprocal-target-inconsistent` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#151 | 10175 | A bare base is not rewritten as a negative power of the target (x with 1/x => w stays x, as in SymPy). |
| `subs-relation-no-method` | SymJ | crash | fixed | jverzani/SymJ.jl#128 | 18412 | subs(::Relation) on both sides |
| `subs-relation-no-method-r10` | SymJ | crash | fixed | jverzani/SymJ.jl#131 | 17530, 16627 | subs on FactAnd, FactOr, Membership |
| `subs-result-not-canonical` | SymJ | missing case | won't fix |  | 20616 | Won't fix (2026-11-21): == is structural; the substituted form is not expanded (as SymPy). |
| `subs-symseries-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#128 | 19120 | subs(::SymSeries): parameters into coefficients; x => c x^m about 0 the composed series (order × m); otherwise a clear error |
| `summation-abs-not-refined-on-range` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | r8-extra-8 | Summand refined on the range: \|x + 9\| = x + 9, divergent: Inf |
| `summation-alternating-divergent-not-detected` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#127 | 17942 | (-1)^k R(k) paired only when R → 0 (deg num < deg den); otherwise 'the series diverges' |
| `summation-alternating-log-series-missing` | SymJ | missing case | won't fix |  | 4189 | feature_audit item 24 (open in SymPy too) |
| `summation-alternating-rational-unevaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 21148 | (-1)^k R(k) summed in pairs R(a+2j) - R(a+2j+1): π/(2 sinh π) + 1/2 |
| `summation-binomial-infinite-upper-not-truncated` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6347 |  |
| `summation-binomial-theorem-subs-divisionbyzero` | SymJ | crash | fixed | jverzani/SymJ.jl#173 | 5790 |  |
| `summation-binomial-trig-unevaluated` | SymJ | missing case | open |  | 14412 | Still open (checked 2026-11-21): Σ binomial(n, m) sin(mx) cos((n-m)x) stays a Sum (2^(n-1) sin(nx)). |
| `summation-binomial-two-variable-form-unsimplified` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 14142 | Fixed on main (checked 2026-11-21). |
| `summation-cancelled-factorial-start-term` | SymJ | missing case | fixed | jverzani/SymJ.jl#162 | 8254 |  |
| `summation-cases-sum-of-terms-unsupported` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | 6173 | Decided by the sign of c e^u + d; cases that stay undecided are still refused. |
| `summation-cos-full-period-unevaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#119 | 23576 | sin/cos in k summed through exponentials (kept when im leaves): 0 |
| `summation-difference-split-divergent-nan` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#116 | 23119 | Divergent parts of a sum (∞ - ∞): the summand as a whole through its antidifference or rational sum: 1 and 1/10 |
| `summation-difference-split-divergent-nan-r16` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#160 | 8251 | NaN over an infinite range: divergence test. |
| `summation-difference-split-divergent-nan-r7` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#120 | 21434, 21315, 21151 | An infinite part beside others: the whole summand decided (Σ c/k - log(1+c/k) = cγ + logΓ(1+c); divergence by its limit or a c/k term; else unevaluated); 1 - 1/k is Inf, not NaN |
| `summation-divergent-nonrational-unevaluated` | SymJ | missing case | fixed | jverzani/SymJ.jl#140 | 14579, 14112, 14111 | Fixed on main (checked 2026-11-21). |
| `summation-drops-small-n-cases` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#86 | 29731 | Binomial theorem with a polynomial factor p(k) (falling-factorial basis); x = -1 gives piecewise values at n = 0..deg p. #29731's sum is piecewise(n ~ 1 => 1, true => 0). |
| `summation-drops-small-n-cases-r16` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#160 | 8077 | _binomial_theorem keeps the cofactor as a product. |
| `summation-exp-series-shifted-start` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 20464 | A series from a small start left unevaluated: from 0 less the first terms: e² - 1 |
| `summation-expand-retry-infinite-recursion` | SymJ | crash | fixed | jverzani/SymJ.jl#177 | 4171 | retry keeps its terms; term-by-term through summation; no reduction that puts a pole on the first terms |
| `summation-factorial-ratio-not-simplified` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | r14-extra-5 | Summands with factorial/gamma go through simplify_gamma first when that is smaller. |
| `summation-geometric-ratio-not-recognized` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#85 | extra-8 | A constant ratio r gives f = f(0) r^k directly; (-1)^n x^n and x^(2n) now sum by cases on \|r\| < 1. GruntzLimit's side (oscillation message) not needed. |
| `summation-geometric-ratio-one-case-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 11823 | Equalities among the assumptions are substituted first: Σ aᵏ with a ~ 1 is the number of terms. |
| `summation-gosper-exp-term-gcd-refused` | SymJ | missing case | fixed | jverzani/SymJ.jl#174 | 4599, 4598 |  |
| `summation-hypergeometric-series-closed-form-missing` | SymJ | missing case | feature gap |  | 11643 | Feature gap (2026-11-21), in the feature audit's roadmap: sums: closed forms in incomplete gamma. |
| `summation-hypergeometric-series-no-convergence-cases` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#138 | 14213, r12-extra-1 | Fixed on main (checked 2026-11-21): Cases on convergence. |
| `summation-infinite-range-pole-term-unchecked` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#138 | r12-extra-5 | Fixed on main (checked 2026-11-21). |
| `summation-negative-infinite-bound` | SymJ | missing case | fixed | #103 | 27893, 27824, 27256 | k -> -k; over all integers split at 0 |
| `summation-noninteger-bound-closed-form` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#123 | 19745 | Mathematica's convention: k = a, a+1, … ≤ b; a non-integer span ends at a + floor(b - a): 1 |
| `summation-piecewise-condition-holds-on-range` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 19687 | Summand refined on a ≤ k ≤ b, k an integer: n/2 |
| `summation-piecewise-point-condition-not-split` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | 18388 | r9-missing-cases: point conditions split off the sum |
| `summation-piecewise-point-condition-not-split-r10` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 16532 | in a product, ≠ first, symbolic point by assumptions |
| `summation-polylog-form-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#167 | 6396 |  |
| `summation-product-summand-not-expanded` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 5953 |  |
| `summation-rational-infinite-unevaluated` | SymJ | missing case | fixed | #89 | 29670 | rational infinite sums via roots and digamma; fixed wrong value of sum 1/((k+1)(2k+1)) (was log 4 - 1) |
| `summation-shift-to-zero-divisionbyzero` | SymJ | crash | fixed | jverzani/SymJ.jl#135 | 14909 | retry at 0 catches DivisionByZero; -log(1-z) series in the table |
| `summation-shifted-limits-not-combined` | SymJ | missing case | fixed | jverzani/SymJ.jl#125 | 19685 | simplify combines Sums with nearby limits: u(n) |
| `summation-summand-zero-not-recognized` | SymJ | missing case | fixed | jverzani/SymJ.jl#122 | 22260 | A summand that simplify_trig makes 0 sums to 0 |
| `summation-symbolic-power-faulhaber-missing` | SymJ | missing case | feature gap |  | 14197 | Feature gap (2026-11-21), in the feature audit's roadmap: sums: Σ m^k with symbolic k. |
| `summation-telescoping-not-detected` | SymJ | missing case | open | jverzani/SymJ.jl#140 | 14195 | Still open (checked 2026-11-21): differences telescope; Σ 1/(√m + √(m - 1)) needs its denominator rationalized first. |
| `summation-trig-finite-closed-form-missing` | SymJ | missing case | fixed |  | 14137 | Fixed on main (checked 2026-11-21): 1/2 + Σ cos(kt) by cases on sin(t/2). |
| `summation-trig-finite-closed-form-missing-r13` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | 12107 | Σ sin/cos(ck + d) closed form, Cases on sin(c/2) = 0. PR open. |
| `summation-trig-infinite-oscillating-cases` | SymJ | missing case | won't fix |  | 8329 | Feature audit 19. |
| `summation-trig-over-k-series-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#149 | 11700 | Σ sin(kθ)/k = (π - θ)/2 and Σ cos(kθ)/k = -log(2 sin(θ/2)) for 0 < θ < 2π. |
| `summation-zero-base-power-summand` | SymJ | missing case | fixed | jverzani/SymJ.jl#162 | 8414 |  |
| `summation-zeta-series-missing` | SymJ | missing case | won't fix |  | 4188 | feature_audit item 24 (open in SymPy too) |
| `summation-zeta-symbolic-exponent-missing` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | r13-extra-3 | Σ 1/k^p = ζ(p) less first terms for Re p > 1, Inf otherwise. PR open. |
| `svd-empty-matrix-reduce-error` | SymJ | crash | fixed | jverzani/SymJ.jl#146 | 10768 | svd/svdvals of an empty matrix are empty. |
| `svd-float-entries-refused-degree-gt2` | SymJ | missing case | fixed | jverzani/SymJ.jl#132 | 17333, 17332 | numeric SVD for float matrices |
| `svd-irreducible-cubic-singular-values-refused` | SymJ | missing case | feature gap |  | 10355 | Feature gap (2026-11-21), in the feature audit's roadmap: linear algebra: SVD when the singular values are roots of an irreducible cubic. |
| `sym-hash-inconsistent-with-isequal` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#172 | 5720, 5072 |  |
| `sym-irrational-constant-methoderror` | SymJ | crash | fixed | jverzani/SymJ.jl#131 | r10-extra-15 | golden and catalan exact; others an ArgumentError |
| `symbolic-det-radical-blowup` | SymJ | hang | fixed | #91 | 28606, 28507 | fixed by Berkowitz charpoly (no division): 3x3 and 5x5 eigvals of square-root matrices under 0.3 s; eigvecs of #28212 too |
| `symbolic-det-radical-blowup-r5` | SymJ | hang | fixed | jverzani/SymJ.jl#113 | 24057 | Eigenvectors for simple eigenvalues with cube or higher roots from an adjugate column of A - tI at t = λ: 2.7 s, residual 0 |
| `symbolic-det-radical-blowup-r6` | SymJ | hang | fixed | jverzani/SymJ.jl#118 | 23408 | det: im and radicals as generators for Bareiss (a polynomial identity), rows cleared of symbolic denominators (det(DA)/det(D)): 0.002 s |
| `symbolic-linsolve-gauss-jordan-blowup` | SymJ | hang | fixed | jverzani/SymJ.jl#113 | 25302, 25032, 25009, 24679, 24135 | Cramer/Bareiss for \, inv, square solve (n ≤ 8; rows cleared of denominators; im/√ atoms as generators); quick monomial-gcd test; D factored once for large quotients. #25302 0.6 s, #24679 0.001 s, #2… |
| `symbolic-linsolve-gauss-jordan-blowup-r11` | SymJ | hang | fixed | jverzani/SymJ.jl#135 | 16245, 16125 | constant zero tests by value; large components uncancelled; radicals as symbols for IVP constants |
| `symbolic-linsolve-gauss-jordan-blowup-r12` | SymJ | hang | fixed |  | 14694 | Fixed on main (checked 2026-11-21): #14694's 18-equation chain finishes in 117 s (was > 500 s). 2026-11-25: 3 s by fraction-free Gauss–Jordan (one Bareiss elimination for all unknowns, no size cap); `inv` and `\` of a sparse 12 × 12 went from minutes to under a second. |
| `symfunction-builtin-name-not-undefined` | SymJ | missing case | fixed | jverzani/SymJ.jl#144 | 12716 | Docstring: a function named like a built-in is the built-in (name-clash question still deferred). PR open. |
| `symfunction-builtin-name-not-undefined-r14` | SymJ | missing case | won't fix |  | 11841 | By design: a SymFunction named like a built-in (gamma) is the user's function; use another name for an undefined one. |
| `symset-no-issubset` | SymJ | crash | fixed | jverzani/SymJ.jl#87 | 29730, 29728 | ⊆ for SymSet and SymInterval |
| `symset-numberset-equality` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#160 | 8553 | SymSet == ℝ when the whole line. |
| `trig-exact-values-missing-pi-12-pi-5` | SymJ | missing case | fixed | jverzani/SymJ.jl#166 | 6048 | π/8 values not autoevaluated (as in SymPy). |
| `trig-imaginary-argument-not-hyperbolic` | SymJ | missing case | fixed | jverzani/SymJ.jl#129 | r9-extra-7 | sin/cos/tan(i·u) → hyperbolic in simplify, when smaller (not automatic: compile hang) |
| `trig-multiple-angle-expansion-unbounded` | SymJ | hang | fixed | #103 | 27293 | simplify does not expand multiples beyond 32 |
| `trig-of-inverse-trig-multiple-angle-unreduced` | SymJ | missing case | fixed | #110 | 26373, 25685 | trig of inverse trig algebraic; integer multiples of π in sums |
| `trig-of-inverse-trig-sum-not-simplified` | SymJ | missing case | won't fix |  | 24328 | cos(acos(4/5) + acos(8/17)) is not expanded at construction, by design (no trig expansion of sums there); simplify gives -13/85 |
| `trig-pythagorean-numeric-arguments` | SymJ | missing case | fixed | #110 | 25705 | trig identities tried on constants |
| `zero-base-power-exponents-combined` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#164 | 6861 |  |
| `zero-base-power-exponents-combined-r18` | SymJ | wrong answer | fixed | jverzani/SymJ.jl#172 | 5804 | The (0, π/2) case fixed; the issue's (0, 2π) case is still a GruntzError (limit-negative-side-symbolic-power). |
| `zero-base-power-unevaluated` | SymJ | missing case | fixed | jverzani/RischIntegrate.jl#92 | 30756 | RischIntegrate part only: the refusal now says 'log of the constant 0 is undefined (a power with base 0)', not 'negative constant 0'; simplifying 0^(x−1) itself is SymJ's |
| `zero-denominator-divideerror` | SymJ | crash | fixed | jverzani/SymJ.jl#85 | extra-2 | Decision (a): keep throwing, now SymJ.DivisionByZero with a message instead of a bare DivideError; internal pole checks accept both. |
| `eq-matrix-fallback-returns-nothing` | other: CommonEq | missing case | open |  | r17-extra-12 |  |

## Missing features

| feature | owner | status | issues |
|:--|:--|:--|:--|
| `erf-integrands` | RischIntegrate | won't fix | 30110 |
| `fresnel-antiderivatives` | RischIntegrate | won't fix | 30070, 30014 |
| `nonelementary-antiderivatives` | RischIntegrate | won't fix | 30374, 29910, 29900 |
| `algebraic-field-domains` | SymJ | won't fix | 30845, 30613, 30196, 30158 |
| `associated-legendre` | SymJ | won't fix | 29645 |
| `bessel-functions` | SymJ | won't fix | 29706 |
| `besselsimp` | SymJ | won't fix | 29905 |
| `beta-function-integrals` | SymJ | won't fix | 29751 |
| `eml-form` | SymJ | won't fix | 29688 |
| `euler-lagrange` | SymJ | won't fix | 29123 |
| `finite-field-polys` | SymJ | won't fix | 29715 |
| `formal-power-series` | SymJ | won't fix | 29027 |
| `fourier-series` | SymJ | won't fix | 30183, 29336 |
| `fresnel-functions` | SymJ | won't fix | 30064 |
| `function-range` | SymJ | won't fix | 30505, 30499, 30386, 30293 |
| `hyperexpand` | SymJ | won't fix | 30708 |
| `imageset` | SymJ | won't fix | 29675, 29435, 29419, 29418, 29416 |
| `jacobian-hessian` | SymJ | won't fix | 30186 |
| `lie-symmetries-ode` | SymJ | won't fix | 30644, 30643 |
| `matrix-calculus` | SymJ | won't fix | 30902, 30527, 29355 |
| `matrix-expressions` | SymJ | won't fix | 30900, 30860, 30278, 30272, 29837, 29758, 29292, 28916 |
| `monotonicity-predicates` | SymJ | won't fix | 30582 |
| `multivariable-limit` | SymJ | won't fix | 30871 |
| `owens-t-function` | SymJ | won't fix | 30345 |
| `parametric-inequalities` | SymJ | won't fix | 30688, 30598 |
| `smith-normal-form` | SymJ | won't fix | 29139 |
| `solve-floor-sign` | SymJ | won't fix | 29026, 29007, 28970 |
| `special-function-rewrites` | SymJ | won't fix | 30352, 30349 |
| `symbolic-order-derivative` | SymJ | won't fix | 29389, 29102 |
| `symbolic-pinv` | SymJ | fixed | 30759 |
| `transcendental-systems` | SymJ | won't fix | 30849, 29315 |
| `trig-sum-equations` | SymJ | won't fix | 30401 |
| `unevaluated-integral` | SymJ | won't fix | 29833 |
| `z-transform` | SymJ | won't fix | 29695 |
| `array-expressions` | other: array-expressions package | won't fix | 30387, 30271, 30270, 30269, 30268 |
| `linear-programming` | other: linear-programming package | won't fix | 29382, 29368 |
| `sat-solver` | other: logic/SAT package | won't fix | 30861, 30729, 30719, 30714 |
| `diophantine` | other: number-theory package | won't fix | 29269 |
| `prime-counting` | other: number-theory package | won't fix | 29125 |
| `quaternions` | other: quaternion package | won't fix | 30560 |
| `vector-calculus` | other: vector-calculus package | won't fix | 29107 |
