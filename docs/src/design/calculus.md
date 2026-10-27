# [Design notes: calculus](@id design-calculus)

## Piecewise functions

- **`SymJ.piecewise(c₁ => v₁, …)`** is an expression, not a `Cases`: a function node
  `piecewise(c₁, v₁, c₂, v₂, …)` whose conditions are nodes too, `Lt(a, b)` and its kin
  and `And(…)`. A relation between constants evaluates to 1 or 0, so substitution picks
  the piece; a condition linear in one symbol is solved for it (`1 - x < 0` as `x > 1`).
  As built, false pieces are dropped, a true one ends the list, and a piece with the
  value of the true one after it goes; with none left the value is NaN.
- **Taken piece by piece**: `diff`; `N` and `lambdify` evaluate only the piece taken (an
  `if` chain), so `piecewise(x < 0 => log(-x), true => sqrt(x))` is real throughout;
  `refine` decides conditions, `heaviside`, `sign` and `abs` under assumptions;
  `limit` replaces each piecewise, `heaviside` and `dirac` by its value on the side
  approached (the sign of a condition from the limit of `g`, or of `g/t^j` when that
  is 0), the two sides compared for a two-sided limit.
- **`heaviside(x, h0 = 1/2)`**, **`dirac(x, k = 0)`**: `diff(heaviside) = dirac`,
  `diff(sign(x)) = 2 dirac(x)`. `piecewise(ex)` writes `heaviside`, `sign`, `abs`,
  `max`, `min` as pieces (real arguments); `piecewise(cs)` turns cases into one.
- **Integration** (in the RischIntegrate extension, a method of `integrate` for `Sym`):
  a `dirac(m x + c, k)` term by `∫ g δ⁽ᵏ⁾(x - a) = (-1)^k g⁽ᵏ⁾(a) heaviside(x - a)`; the
  rest between consecutive breakpoints (sorted real numbers, or a single one with
  parameters), the integrand read at a point inside each interval and the constants
  chosen for continuity. `abs`, `sign`, `max` and `min`, which RischIntegrate does not
  take, are rewritten as pieces when it fails on them.
- **Folding** (`SymJ.piecewise_fold(e, x)`): `heaviside`, `sign`, `abs`, `max`, `min`
  rewritten as pieces, the breakpoints (roots of the conditions) found and ordered, as
  numbers or by the assumptions, the parameters taken to be real; an order they do not
  decide, and a condition they do not decide at a point, are split into cases (the
  machinery of `Cases`: `a < b`, `a = b`, `a > b`, the equal points merged), then
  the value on each open interval (at a point inside) and at each breakpoint (exactly),
  each condition decided by the sign of its difference under the assumptions. Pieces
  are merged: an interval into the piece before it when the formulas agree, a
  breakpoint when the formula before it gives its value there (else it starts a piece
  that the next interval joins when that formula gives the value); each condition is
  by the piece's right end (`x < b`, `x ≤ b`, `true`). Without `x`, the one symbol of
  the conditions. `simplify` tries it for two or more piecewise functions, kept when
  smaller: `∫₀¹ |x - a| dx`, a difference of two piecewise functions from the
  integrator, becomes three pieces in `a`.
- **Printing**: `piecewise(x ≪ 0 => -x, true => x)`; in LaTeX a `cases` environment,
  `heaviside` and `dirac` as θ and δ.

## Summation

- **`SymJ.summation(f, k, a, b)`**, a chain of methods as in SymPy: linearity; direct
  addition for a short numeric range; an antidifference `F` (`F(k+1) - F(k) = f(k)`)
  giving `F(b+1) - F(a)`, by Bernoulli polynomials for polynomials and Gosper's
  algorithm for hypergeometric terms (the term ratio from polynomials, `r^k`,
  factorials, binomials, `gamma`; the Gosper form through `resultant` and `gcd`; the
  polynomial `x(k)` from its degree bound and a linear system); partial fractions for
  rational functions without one (harmonic numbers, polygamma, exact `ζ(2n)` by
  Bernoulli numbers); known sums (the binomial theorem; `₀F₀`, `₁F₀`, `₀F₁`
  hypergeometric series: `exp`, `cos`, `sin`, `cosh`, `sinh`).
- **Zeilberger's algorithm** for a definite sum whose upper bound involves `n` and
  whose summand has no antidifference: creative telescoping, Gosper's equation with
  the coefficients `c_j(n)` of a recurrence `Σ_j c_j(n) S(n+j) = 0` as further
  unknowns of the same linear system. A first-order recurrence is solved as a product
  (factorials, half-integer shifts through `(2n)!/(4^n n!)`, presented as binomials);
  a recurrence of order 2 or 3 by its hypergeometric solutions (Hyper, below), combined
  to match the sum's first values. Each closed form is checked against direct sums for
  several `n` (and integer values of other parameters); a failed check leaves the sum
  unevaluated. With the summand's own bounds the recurrence is homogeneous, and that is
  tried first. Otherwise (an upper bound `s n + t`, `s` a positive integer, a lower one
  free of `n`) the telescoper `Σ_j c_j F(n+j, k) = G(n, k+1) - G(n, k)`, `G = R F` with
  `R` from Gosper's solution, summed over `k = a..b(n)`, gives the right side
  `G(n, b+1) - G(n, a) + Σ_j c_j Σ_{k=b(n)+1}^{b(n+j)} F(n+j, k)`; `G(n, b+1)` is taken
  through the ratio `F(n, b+1)/F(n, b)` as a rational function, so a zero of `F` there
  meets a pole of `R`. The inhomogeneous recurrence is made homogeneous as in `rsolve`
  (below) and solved by Hyper: `Σ_{k=0}^{n-1} binomial(n, k) = 2^n - 1`,
  `Σ_{k=0}^{n} binomial(n, k)/(k + 1) = (2^{n+1} - 1)/(n + 1)`,
  `Σ_{k=0}^{n} binomial(2n, k) = (4^n + binomial(2n, n))/2`.
- **Right sides.** `rsolve` takes `Σ_i p_i y(n+i) = g` with `g` a sum of hypergeometric
  terms: grouped into classes (a rational ratio between members), each class
  annihilated in turn by `L ↦ ρ_den(n) L(n+1) - ρ_num(n) L(n)`, which maps every other
  class `c` to `c (ρ_den ρ_c - ρ_num)` (rational multiples, so still one term each), one
  order higher per class, then Hyper on the homogeneous recurrence. Powers with shifted
  exponents are merged in the answer (`8·2^(n - 3)` as `2^n`).
- **Hyper** (Petkovšek), for `Σ_i p_i(n) y(n+i) = 0`: a hypergeometric solution has
  `y(n+1)/y(n) = Z A(n)/B(n) C(n+1)/C(n)` with `A` a monic factor of `p_0(n)`, `B` of
  `p_d(n-d+1)`; for each pair, `Z` a nonzero root of the leading coefficients of
  `P_i = p_i Π_{j<i} A(n+j) Π_{j≥i} B(n+j)`, and `C` a polynomial solution of
  `Σ Z^i P_i C(n+i) = 0`, its degree bounded by the indicial polynomial of the
  recurrence written in differences (`E = 1 + Δ`), its coefficients from a nullspace.
  When that finds too few solutions, it is run again with the coefficients' irreducible
  factors of higher degree split over their roots (radicals), so `A = n - sqrt(2)` can
  divide `n^2 - 2`; the values used to fit the solution are simplified first, as they
  have Γ at algebraic shifts.
  `SymJ.rsolve(eq, y(n), inits...)` reads the
  recurrence, raises a hypergeometric right side `g` to a homogeneous recurrence one
  order higher (`L(y)(n+1) - ρ(n) L(y)(n) = 0`, `ρ = g(n+1)/g(n)`), fits the solutions
  to the initial values (or to `y(0), …, y(d-1)` themselves) and checks the result
  against values from the recurrence.
- **d'Alembertian solutions** (Abramov–Petkovšek), when Hyper's solutions do not span:
  reduction of order. A hypergeometric solution `h` of the homogeneous recurrence, with
  `y = h z` and `z(n+i) = z(n) + Σ_{j<i} u(n+j)` (`u = Δz`), gives
  `Σ_j s_j u(n+j) = g/h`, `s_j = Σ_{i>j} p_i h(n+i)/h(n)`, one order lower (the
  coefficient of `z(n)` is `L(h)/h = 0`). At order 0, `u = g/s_0`, so the right side
  need not be hypergeometric. Each level adds a solution `h Σ_{m=n0}^{n-1} u(m)`; `n0`
  is past the integer zeros and poles of every level's ratios and leading coefficient.
  The constants are fitted with the sums unevaluated (added up at numbers, so
  `polygamma` need not evaluate), then the sums are summed, innermost first, and the
  answer is checked in the unevaluated form: `y(n+1) - y(n) = 1/(n+1)` gives
  `harmonic(n)`, `y(n+1) = (n+1) y(n) + 1` gives `n! (1 + Σ 1/(k+1)!)`. Zeilberger's
  recurrences go the same way when Hyper falls short:
  `Σ_{k=1}^{n} (-1)^(k+1) binomial(n, k)/k = harmonic(n)`. A level without a
  hypergeometric solution ends it, so the solution is not d'Alembertian: the Franel
  recurrence of `Σ binomial(n, k)^3` is irreducible of order 2, and that sum stays
  unevaluated (Liouvillian solutions, Hendriks–Singer, would not help either).
- **Harmonic numbers**, `Σ_{k=a}^{b} g(k) H(k + c)` with `c` an integer and `g` with an
  antidifference `G`, by parts: `G(b+1) H(b+c+1) - G(a) H(a+c) - Σ G(k+1)/(k+c+1)`.
  A polynomial `G` is shifted by a constant to vanish at `-c`, so the last sum is of a
  polynomial; the answer is written with `H(b+c+1)` or `H(b+c)`, whichever is smaller:
  `Σ_{k=1}^{n} k H(k) = n(n+1) H(n)/2 - n(n-1)/4`. The sum stays unevaluated when the
  last sum does not close (`Σ 2^k H(k)`), and for `Σ H(k)/(k+1)`, which needs
  generalized harmonic numbers.
- **Infinite sums** are limits of the antidifference (`limit`). A polynomial times
  `r^k` converges exactly when `-1 < r < 1`, decided by cases on a parameter.
- **Otherwise** the sum stays unevaluated, `Sum(f, k, a, b)`, which adds up once its
  bounds become numbers. The name is `summation`, not a method of `Base.sum`, since the
  call is different (an index and bounds, not an iterable).

## Differential equations

- **`SymJ.dsolve(eq, y(x), ics...)`** reads the equation as `F(x, Y₀, …, Yₙ) = 0`, the
  derivatives of `y(x)` replaced by placeholder symbols, and tries the methods in turn:
  a missing `y` (solve for `v = y⁽ʲ⁾`, then integrate), linear (`F` affine in the `Yₖ`:
  the integrating factor at first order; constant coefficients by the characteristic
  roots, with `x^j` for a multiple root and `e^(αx) cos βx`, `e^(αx) sin βx` for a
  conjugate pair; Cauchy–Euler by `x = eᵗ`, which turns `xᵏ y⁽ᵏ⁾` into
  `θ(θ - 1)…(θ - k + 1) y`), then first order `y' = f(x, y)`: separable (`f` factors as
  `g(x)h(y)`), Bernoulli, exact (with an integrating factor in `x` or `y` alone),
  homogeneous (`f(x, vx)` free of `x`).
- **Riccati equations**, `y' = q₀ + q₁ y + q₂ y²`, the last first-order method:
  `y = -u'/(q₂ u)` gives `u'' - (q₁ + q₂'/q₂) u' + q₀ q₂ u = 0`, solved when it has
  constant coefficients or is Cauchy–Euler, and `u = C u₁ + u₂` gives `y` with one
  constant. Otherwise a particular solution `y₁ = c xᵏ` (`|k| ≤ 2`) or `c e^(kx)`
  (`|k| ≤ 2`), with `c` constant, and `y = y₁ + 1/v`, `v` from a linear equation. A
  first-order method whose answer has an unevaluated integral is kept only if no later
  one does better: `y' = y² + 4` is separable, but without RischIntegrate the
  integral is left undone, while as a Riccati equation it has a closed form in `sin`
  and `cos`.
- **Autonomous equations**, second order without `x`: with `p = y'` a function of
  `s = y`, `y'' = p dp/ds`, a first-order equation for `p(s)` solved by the methods
  above; then `∫ ds/p(s) = x + C₂`, solved for `y` when that gives one solution and
  implicit otherwise. When `p` has several branches (`p = ±sqrt(q(s))`, as for
  `y'' = -sin y`) each gives a solution, and `dsolve` returns a vector. With initial
  conditions `y(a) = v`, `y'(a) = w`, `C₁` is fixed first from the relation between
  `p` and `s` at `(v, w)`, the branches through `(v, w)` kept, and only then
  `∫ ds/p` taken, which is often elementary when the general one is not
  (`y'' = 2y³`, `y(0) = 1`, `y'(0) = 1`: `C₁ = 0`, `y = 1/(1 - x)`); without a closed
  form it is the definite `∫ᵥʸ ds/p = x - a`, an `Integral(f, y, v, y(x))` node printed
  `Integral(f, (y, v, y(x)))` as `integrate` takes a definite integral; `y` is bound
  there (for `subs` and `free_symbols`), and `diff` uses the fundamental theorem. `sqrt(exp(2y))` is `exp(y)`, `y` being
  real. An implicit solution of any order is checked, and its constants
  fixed by initial conditions at one point, through implicit differentiation: `y' =
  -φₓ/φ_y` and each higher derivative by `d/dx + y' d/dy`. When the explicit form gives
  a nonlinear system for the constants (`y = e^(C₁(C₂ + x))`), the implicit one is used.
- **Particular solutions.** Undetermined coefficients when every term of the right side
  is `c xᵐ e^(αx) cos(βx)` (or `sin`); the trial form is multiplied by `xˢ`, `s` the
  multiplicity of `α + iβ` as a root, and the coefficients come from a linear system.
  Otherwise variation of parameters, with the Wronskian.
- **Solutions** are explicit `y = S` or implicit `φ(x, y) = C`; an implicit one is solved
  for `y` when that gives one solution (with an initial condition, the branch through
  the point). Each is checked by substitution (implicitly, `y' = -φₓ/φ_y`),
  symbolically and else numerically at a few points. Constants are absorbed where one
  appears once: `k C` is `C`, and `e^C` in an exponent is `C`.
- **Derivatives as values.** `Derivative(f, v)` binds `v`: substituting a value gives
  the derivative at a point, `Derivative(f, v, a)`, and substituting into `f` evaluates
  it (`Derivative(y(x), x)` with `y(x) => sin(x)` is `cos(x)`). `y'` on a `SymFunction`
  builds these, so initial conditions read `y'(0) => 1`, and a derivative at a point
  prints the same way (`y''(0)`, not `Derivative(Derivative(y(x), x), x, 0)`).
- **Series solutions** (`series = n`, about `x0`). At an ordinary point, the equation
  solved for `y⁽ᵏ⁾ = f`: the first `k` Taylor coefficients are constants (or initial
  values), and each next `cₘ` is the coefficient of `tᵐ⁻ᵏ` in the series of
  `f(x₀ + t, y, …)` with the coefficients so far, which is all it depends on. This is
  the series code doing the work, so nonlinear equations need nothing more. A pole or a
  fractional power in that series marks a singular point. There, for a second-order
  linear homogeneous equation with `t p(t)`, `t² q(t)` analytic, Frobenius: the
  indicial equation, the recurrence for `aₘ(r)` with `r` a symbol (rational functions,
  cancelled), and for the second solution `∂/∂r` of `tʳ Σ aₘ(r) tᵐ` (equal roots) or of
  `(r - r₂) tʳ Σ aₘ(r) tᵐ` (roots an integer apart) at `r₂`, which is where the `log t`
  term comes from. Each series is checked: its residual is `O(tⁿ⁻ᵏ)` (Frobenius: over
  `a₂/t²`, `O(t^(r₂ + n))`).
- **Laplace transforms** (`src/laplace.jl`). Forward, term by term: `dirac(t - a) g` is
  `e^(-as) g(a)` and `heaviside(t - a) g` is `e^(-as) L{g(t + a)}`; anything else has its
  `sin`, `cos`, `sinh`, `cosh` written as exponentials, so each term is `c tⁿ e^(αt)`
  (`α` possibly complex), with transform `c n!/(s - α)ⁿ⁺¹`; the sum is cancelled, which
  makes it real again. Inverse: the terms grouped by their `e^(-cs)` factor (a
  `heaviside(t - c)` shift each), each group split by `apart` into `N/qᵏ` with `q`
  irreducible, and inverted by the residues of `N e^(st)/qᵏ` at the roots of `q`, then
  made real as for systems. A quadratic's roots are `α ± iβ`, `β = sqrt(c/a - α²)` with
  the parameters taken to be positive (so `s² + w²` gives `sin(wt)/w`, not exponentials
  of `sqrt(-w²)`); either square root gives the same function, a conjugate pair being
  even in `β`. `dsolve` uses them for a linear equation (or system) with
  constant coefficients when asked (`method = :laplace`) or when the right side has
  `heaviside` or `dirac`, where the other methods leave integrals undone:
  `Y(s) = (L{g} + Σₖ aₖ Σⱼ sᵏ⁻¹⁻ʲ y⁽ʲ⁾(a)) / Σₖ aₖ sᵏ`, with `t` shifted to the point `a`
  of the initial conditions; for a system `Z(s) = (sI - M)⁻¹ (Z(a) + L{F})`.
- **Systems.** `SymJ.dsolve([eq₁, …], [y₁(t), …], ics...)`, linear with constant
  coefficients. The highest derivative of each unknown is solved for (the matrix of
  their coefficients must be invertible: the system in normal form), and the system made
  first order, `Z' = M Z + F(t)`, with `Z` the unknowns and their lower derivatives. With
  `M = P J P⁻¹` (Jordan form), the columns of `P exp(J t)` are solutions; a conjugate
  pair of complex columns is replaced by its real and imaginary parts (`i` rewritten out
  of denominators first, and `e^((a + ib)t)` as `e^(at)(cos bt + i sin bt)`), giving a
  real fundamental matrix `Φ`. The general solution is `Φ C`, so a defective `M` gives
  `t^k e^(λt)` terms. A particular solution is `P exp(J t) ∫ exp(-J t) P⁻¹ F dt` when
  the integrand, with `sin`, `cos`, `sinh`, `cosh` as exponentials, is a sum of
  `t^k e^(at)` (integrated in closed form, then made real); otherwise `Φ ∫ Φ⁻¹ F dt`
  through the integrator. Parameters in `M` are taken to be real.
- **Integrals** go through an integrator hook that the RischIntegrate extension sets;
  without it a small integrator does polynomials, powers of a linear function, and
  `exp`, `sin`, `cos`, `sinh`, `cosh` of linear arguments, and anything else is an unevaluated `Integral(f, x)` whose
  derivative is `f`.
