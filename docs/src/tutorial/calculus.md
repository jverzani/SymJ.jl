# Calculus

Limits, derivatives, series, sums, integrals and differential equations.

```@setup tut
using SymJ
using RischIntegrate
@syms x y z a b c t f()
```

## Limits

`limit(f, x => x0)` is the limit of `f` as `x` tends to `x0`, a
number, an expression or `Inf`. It is computed by Gruntz's algorithm,
through GruntzLimit. An infinite limit is returned as `Sym(Inf)` or `Sym(-Inf)`:

```@example tut
limit(sin(x) / x, x => 0)
```

```@example tut
limit((1 + 1 / x)^x, x => Inf), limit(x * ((1 + 1 / x)^x - ℯ), x => Inf)
```

Using exact, symbolic numbers is important in many limits.

```@example tut
limit((sin(x) - sin(Sym(π) / 6)) / (x - Sym(π) / 6), x => Sym(π) / 6)
```

The limit is two-sided unless `dir = :+` or `dir = :-` asks for a one-sided one:

```@example tut
limit(x * exp(1 / x), x => 0; dir = :-), limit(x^x, x => 0; dir = :+)
```

All of `-`, `:-`, `"-"`, `"+"`, `:+`, and `+` can be used to specify a one-sided direction.


The algorithm handles cancellations that numerical evaluation cannot follow:

```@example tut
limit(exp(x) * (exp(1 / x - exp(-x)) - exp(1 / x)), x => Inf)
```

Other symbols may appear:

```@example tut
limit((1 + a / x)^(b * x), x => Inf)
```

The other symbols are taken to be real numbers. When the limit depends
on their values, `limit` gives it for each case:

```@example tut
limit(x^a, x => Inf)
```

Assumptions, as in the section on solving, settle cases:

```@example tut
(
limit(x^a, x => Inf; assume = @assume(a > 1)),
limit((a^x + b^x)^(1 / x), x => Inf; assume = @assume(a > b > 0))
)
```

## Derivatives

`diff(f, x)` differentiates with respect to `x`, and `diff(f, x, n)` gives the `n`-th derivative:

```@example tut
diff(sin(x)^2 * exp(x), x)
```

```@example tut
diff(x^x, x), diff(atan(x), x, 3)
```

Mixed partial derivatives take the variables in order: `diff(f, x, y)` differentiates with respect to `x`, then `y`:

```@example tut
(
diff(x^3 * y^2, x, y),
diff(sin(x * y), x, y),
diff(x^2 * y * z, x, x, y, z)
)
```

The gradient, Jacobian and Hessian (not exported):

```@example tut
@syms r θ
(
SymJ.gradient(x^2 * y + sin(y), [x, y]),
SymJ.jacobian([r * cos(θ), r * sin(θ)], [r, θ]),
SymJ.hessian(x^3 + x * y^2, [x, y])
)
```

and the divergence, curl and Laplacian, in Cartesian coordinates:

```@example tut
(
SymJ.divergence([x^2, x * y, y * z], [x, y, z]),
SymJ.curl([-y, x, 0], [x, y, z]),
SymJ.laplacian(x^2 * y + y^3, [x, y])
)
```

A symbolic function gives an unevaluated derivative:

```@example tut
diff(f(x^2), x)
```

For implicit differentiation, write `y` as a function `u(x)`, differentiate, and solve for the derivative:

```@example tut
@syms u()
F = x^2 + u(x)^2 - 1
solve(diff(F, x) ~ 0, diff(u(x), x))
```

### Example: critical points

The critical points of a function are the zeros of its derivative, and the second derivative classifies them:

```@example tut
h = x^3 - 6x^2 + 9x + 1
cps = solve(diff(h, x) ~ 0, x)
```

```@example tut
h2 = diff(h, x, 2)
[(c, h(c), h2(c)) for c in cps]
```

So 1 is a local maximum and 3 a local minimum.

More complicated equations for finding the critical points are also possible.

```@example tut
solve(diff((a / x)^x, x) ~ 0, x), solve(diff(x^2 * exp(-x), x) ~ 0, x)
```


### Example: tangent lines

The tangent line at `c` is `h(c) + h'(c)(x - c)`:

```@example tut
c₀ = 2
expand(h(c₀) + diff(h, x)(c₀) * (x - c₀))
```

## Integrals

Antiderivatives come from `integrate`, from RischIntegrate, an
implementation of the Risch algorithm. It decides whether an
elementary antiderivative exists, and finds it when it can.

```@example tut
integrate(x^2 * exp(x), x)
```

```@example tut
integrate(1 / (x^3 + 1), x)
```

```@example tut
integrate(1 / (x^3 - x), x), integrate((2x + 3) / (x^2 + 1), x)
```

SymJ integrates over the real numbers: the variable of integration and the parameters
are taken to be real. This differs from Rubi, Mathematica and SymPy, which take them to be
complex and give antiderivatives valid over ℂ (on principal branches). Over ℝ, an
antiderivative whose real form depends on the sign of a parameter comes by cases, each
case real (`assume` picks one):

```@example tut
integrate(1 / (x^2 + a), x)
```

```@example tut
integrate(1 / (x^2 + a), x; assume = @assume(a > 0))
```

For the same reason a logarithm of an expression that may be negative has its absolute
value: `∫ 1/x dx` is `log|x|`, an antiderivative on both sides of 0 (the others give
`log(x)`, complex for `x < 0`). An argument known to be positive keeps the plain log:

```@example tut
integrate(1 / x, x), integrate(1 / x, x; assume = @assume(x > 0)), integrate(2x / (x^2 + 1), x)
```

```@example tut
integrate(sin(x)^3, x), integrate(sqrt(1 - x^2), x)
```

```@example tut
integrate(x * log(x), x), integrate(1 / (x * log(x)), x)
```

Piecewise functions, `heaviside`, `abs` and `sign` are integrated interval by interval, the antiderivative made continuous; a `dirac` picks out a value:

```@example tut
H = SymJ.heaviside
integrate(x * H(x - 1), x), integrate(cos(x) * SymJ.dirac(x - 1), x)
```

### Definite integrals

A definite integral is the difference of an antiderivative's values at the end points. Calling the antiderivative with exact numbers keeps the answer exact:

```@example tut
F = integrate(x^2, x)
F(3) - F(0)
```

```@example tut
integrate(abs(x - 1), (x, 0, 3))
```

`integrate(f, (x, a, b))` splits the interval at the points where `f` is singular and takes one-sided limits there and at infinite limits, so improper integrals come out right; one that diverges is `±∞`, or an error when it has no value:

```@example tut
(integrate(1 / (1 + x^2), (x, 0, Inf)),
integrate(exp(-x^2), (x, -Inf, Inf)), integrate(log(x), (x, 0, 1)),
integrate(1 / x^2, (x, -1, 1)))
```

Some antiderivatives are known only numerically: RischIntegrate's
`₂F₁` terms and sums over the roots of a polynomial, which it can
evaluate but SymJ cannot hold. A definite integral through one, over
numeric limits where `f` is not singular, is a float:
`integrate(x / sqrt(1 - x^3), (x, 0, 1//2))` is `0.12832…`.

An integrand proven to have no elementary antiderivative, over numeric limits where the
integral converges, gives the integral unevaluated; `N` (or `float`, or calling it) finds its
value by numeric quadrature (QuadGK's), and `SymJ.evalf` to a `BigFloat` precision:

```@example tut
ex = integrate(exp(x^3), (x, 0, 1))
ex, N(ex), SymJ.evalf(ex, BigFloat, 100)
```

```@example tut
@syms y
ex = integrate(exp(-x^2 * y^2), (x, 0, 1), (y, 0, 1))
ex, N(ex)
```

Before keeping such an integral, SymJ tries differentiation under the integral sign
(Feynman's trick). A constant `c` of the integrand becomes a parameter `a`, and the
integral is found as `I(c) = I(a₀) + ∫ I'(a) da` from a base point `a₀` (0 or 1) where it
is known, with `I'(a) = ∫ ∂f/∂a dx`. Here `I'(a) = 2π/√(1 - a²)`:

```@example tut
integrate(log((1 + cos(x)/2) / (1 - cos(x)/2)) / cos(x), (x, 0, π))
```

A value found this way is kept only when it agrees with the integral's numeric value, so
a constant for which the exchange of derivative and integral is not valid gives no answer,
never a wrong one.

`SymJ.Integral(f, (x, a, b))` (or `SymJ.Integral(f, x)`) writes an integral without
integrating it; `integrate(ex)`, with no variable, carries out the integrals in `ex`,
keeping those it cannot do:

```@example tut
I = SymJ.Integral(exp(-x^2), (x, 0, Inf))
I, N(I), integrate(I)
```

`principal_value = true` asks for the Cauchy principal value instead:
a singular point inside is approached from both sides at once,
lim_{ε→0⁺} (∫ₐ^{s-ε} + ∫_{s+ε}^b) f, and an integral over (-∞, ∞) is
the limit of ∫ from -R to R. An integral that converges keeps its value; a
singular end point is not crossed, so ∫₀¹ 1/x is still `∞`:

```@example tut
pv(f, ab) = integrate(f, ab; principal_value = true)
(pv(1 / x, (x, -1, 2)),
pv(1 / x^3, (x, -1, 2)),
pv(tan(x), (x, 0, π)),
pv(x / (1 + x^2), (x, -Inf, Inf)))
```

With a parameter the answer may depend on its sign; the cases are
listed unless the assumptions decide them:

```@example tut
integrate(exp(-a * x), (x, 0, Inf))
```

```@example tut
integrate(exp(-a * x), (x, 0, Inf); assume = @assume(a > 0))
```

Where the formula divides by a parameter, the parameter takes a
generic value (nonzero, as for `solve`); the special value is a case
only when assumed:

```@example tut
@syms c
integrate(exp(c * x), (x, 0, 1))
```

```@example tut
integrate(exp(c * x), (x, 0, 1); assume = @assume(c == 0))
```

An iterated integral lists its intervals innermost first. The inner
integral may have cases on the outer variable (here on the sign of
`x`, for `y^m`); the outer interval keeps the ones it allows, or is
split where they change:

```@example tut
@syms n m
integrate(x^n * y^m, (y, 0, x), (x, 0, 1); assume = @assume(n > 0, m > 0))
```

### Implicit assumptions for integration

The variable of integration, by default, is assumed to be real, as are
the parameters, so the `conj` that `dot` puts
on symbols drops out. A line integral, the work of ⟨-y, x⟩ around the
unit circle:

```@example tut
using LinearAlgebra: dot
r = [cos(t), sin(t)]
integrate(dot([-r[2], r[1]], diff.(r, t)), (t, 0, 2Sym(π)))
```

```@example tut
G = integrate(sin(x), x)
G(x => Sym(π)) - G(x => 0)
```

## Differential equations

`SymJ.dsolve` solves an ordinary differential equation for an unknown
function. The function is made with `@syms y()`, and its derivatives
are written `y'(x)`, `y''(x)` (or `diff(y(x), x)`). The answer is an
equation, with constants `C1, C2, …`:

```@example tut
@syms y()
SymJ.dsolve(y'(x) + 2y(x) ~ exp(x), y(x))
```

!!! note "The `'` operator"
    The prime goes on the function, before the call. After the call, `'`
    is Julia's adjoint of the value `y(x)`, its complex conjugate, and not
    a derivative. That is, these differ: `y'(x)`, `y(x)'` are quite different.

First-order equations that are linear, separable, Bernoulli, exact or
homogeneous are recognized:

```@example tut
(
SymJ.dsolve(y'(x) ~ x * y(x)^2, y(x)),
SymJ.dsolve(y'(x) ~ (x + y(x)) / x, y(x))
)
```

A Riccati equation, `y' = q₀ + q₁ y + q₂ y²`, is made linear (of
second order), or solved from a particular solution:

```@example tut
SymJ.dsolve(y'(x) ~ y(x)^2 - 2 / x^2, y(x))
```

An equation nonlinear in `y'` gives its general solutions and its
singular solutions; Clairaut's equation gives a family of lines and
their envelope, and d'Alembert's equation a solution in the parameter
`p = y'`:

```@example tut
(
SymJ.dsolve(y(x) ~ x * y'(x) + y'(x)^2, y(x)),
SymJ.dsolve(y(x) ~ 2x * y'(x) + y'(x)^2, y(x))
)
```

When `y` cannot be solved for, the solution is implicit:

```@example tut
SymJ.dsolve(2x * y(x) + (x^2 + 3y(x)^2) * y'(x) ~ 0, y(x))
```

Linear equations with constant coefficients, of any order, are solved
through the roots of the characteristic polynomial; a right side is
handled by undetermined coefficients, or by variation of parameters:

```@example tut
(
SymJ.dsolve(y''(x) + 2y'(x) + 5y(x) ~ 0, y(x)),
SymJ.dsolve(y''(x) + y(x) ~ sin(x), y(x))
)
```

```@example tut
SymJ.dsolve(y''(x) + y(x) ~ tan(x), y(x))
```

Cauchy–Euler equations, with coefficients `x^k` on `y⁽ᵏ⁾`, are solved too:

```@example tut
SymJ.dsolve(x^2 * y''(x) - 2y(x) ~ 0, y(x))
```

Equations that reduce to Bessel's or Airy's equation have solutions in
those functions, or in elementary ones for Bessel functions of
half-integer order:

```@example tut
(
SymJ.dsolve(x * y''(x) + y'(x) + y(x) ~ 0, y(x)),
SymJ.dsolve(y''(x) - x * y(x) ~ 0, y(x)),
SymJ.dsolve(y''(x) + 3 / (2x) * y'(x) + y(x) / (4x^3) ~ 0, y(x))
)
```

An equation of second order without `x`, `y'' = f(y, y')`, is reduced
to first order with `p = y'` as a function of `y`:

```@example tut
SymJ.dsolve(y(x) * y''(x) ~ y'(x)^2, y(x))
```

When `y'` has two branches, `±sqrt(…)`, there is a solution for each;
initial conditions pick one:

```@example tut
SymJ.dsolve(y''(x) ~ 2y(x)^3, y(x), y(0) => 1, y'(0) => 1)
```

Initial conditions, one for each order, fix the constants:

```@example tut
SymJ.dsolve(y''(x) + y(x) ~ 0, y(x), y(0) => 0, y'(0) => 1)
```

```@example tut
SymJ.dsolve(y'(x) ~ -x / y(x), y(x), y(0) => 2)
```

With `series = n` the solution is a series, here about `0` (Airy's
equation), with error `O(x^n)`:

```@example tut
SymJ.dsolve(y''(x) - x * y(x) ~ 0, y(x); series = 8)
```

At a regular singular point it is a Frobenius series, which can have a
`log` term (Bessel's equation of order 0):

```@example tut
SymJ.dsolve(x^2 * y''(x) + x * y'(x) + x^2 * y(x) ~ 0, y(x); series = 6)
```

A right side with a step or an impulse (`SymJ.heaviside`,
`SymJ.dirac`) is handled by Laplace transforms:

```@example tut
SymJ.dsolve(y''(x) + y(x) ~ SymJ.heaviside(x - 1), y(x), y(0) => 0, y'(0) => 0)
```

The transforms themselves are `SymJ.laplace_transform` and
`SymJ.inverse_laplace_transform`:

```@example tut
@syms s
(
SymJ.laplace_transform(x * sin(x), x, s),
SymJ.inverse_laplace_transform(exp(-2s) / (s^2 + 1), s, x)
)
```

A system of linear equations with constant coefficients is given as
vectors of equations and unknowns. The answer is a vector of
equations:

```@example tut
@syms t u() v()
SymJ.dsolve([u'(t) ~ u(t) - v(t), v'(t) ~ 5u(t) - 3v(t)], [u(t), v(t)])
```

Higher derivatives, right sides and initial conditions (one for each
unknown and each order below its highest) are allowed:

```@example tut
SymJ.dsolve([u''(t) + v(t) ~ 0, v'(t) ~ u(t) + t], [u(t), v(t)], u(0) => 1, u'(0) => 0, v(0) => 0)
```


## Series

`series(ex, x, x0, n)` expands `ex` about `x0` with an error of order `(x - x0)^n`; `x0` defaults to 0 and `n` to 6:

```@example tut
series(sin(x), x)
```

```@example tut
series(exp(x), x, 1, 3)
```

Laurent series (negative powers), Puiseux series (fractional powers), logarithmic terms and expansions at infinity (in powers of `1/x`) are included:

```@example tut
(
series(cot(x), x, 0, 4),
series(sqrt(x + x^2), x, 0, 3)
)
```

```@example tut
(
series(x^x, x, 0, 3),
series(sqrt(x^2 + x), x, Inf, 3)
)
```

At `±Inf` an expansion may need a factor that is not a power of `x`:
terms `exp(g)*(series)`, `cos(h)*(series)` and `sin(h)*(series)`, from
`exp`, `sin` and `cos` of growing arguments and from the asymptotic
expansions of SpecialFunctions' `erf`, `erfc`, `erfi`, `erfcx`,
`expinti`, `expint`, `sinint`, `cosint`, the Bessel functions and
`gamma` (Stirling's series; the special functions need to be brought
in from SpecialFunctions to be used directly). Factors that cancel do
cancel:

```@example tut
series(sin(x) / x + exp(-x) * (1 + 1/x), x, Inf, 4)
```


```@example tut
using SpecialFunctions: erf, erfc
series(erf(x), x, Inf, 4)
```
Exponentially small terms are kept: `erf(x)` at `Inf` is not just `1`.


```@example tut
series(exp(x^2) * erfc(x), x, Inf, 4)
```


A series is an object that keeps its error term. Series about the same
point can be added, multiplied and divided, with the order tracked:

```@example tut
s, c = series(sin(x), x), series(cos(x), x)
s / c
```

Calling `Sym` on a series drops the `O` term:

```@example tut
Sym(s)
```

## Sums

`SymJ.summation(f, (k, a, b))` (not exported) sums `f` over `k = a, …, b`, in closed form where it finds one; `b` may be `Inf`:

```@example tut
@syms k n r
(
SymJ.summation(k^2, (k, 1, n)),
SymJ.summation(r^k, (k, 0, n)),
SymJ.summation(1 / (k * (k + 1)), (k, 1, n))
)
```

```@example tut
(
SymJ.summation(1 / k^2, (k, 1, Inf)),
SymJ.summation(x^k / factorial(k), (k, 0, Inf))
)
```

Definite sums of binomials, which have no antidifference, come from Zeilberger's
algorithm:

```@example tut
@syms m
(
SymJ.summation(binomial(n, k)^2, (k, 0, n)),
SymJ.summation(binomial(n, k) * binomial(m, n - k), (k, 0, n))
)
```

When Zeilberger's recurrence is of higher order, its hypergeometric solutions (found by Petkovšek's algorithm Hyper) may still give the sum; here Binet's formula for the Fibonacci numbers:

```@example tut
SymJ.summation(binomial(n - k, k), (k, 0, n))
```

Bounds other than the summand's own give an inhomogeneous recurrence, solved as well:

```@example tut
(
SymJ.summation(binomial(n, k), (k, 0, n - 1)),
SymJ.summation(binomial(2n, k), (k, 0, n))
)
```

The same algorithm solves linear recurrences, `SymJ.rsolve` (not exported), from initial values or in terms of them:

```@example tut
@syms w()
(
SymJ.rsolve(w(n + 2) ~ w(n + 1) + w(n), w(n), w(0) => 0, w(1) => 1),
SymJ.rsolve(w(n + 1) ~ 2w(n) + 1, w(n))
)
```

When a recurrence has too few hypergeometric solutions, its order is
reduced through one of them; the solutions are then iterated sums of
hypergeometric terms (harmonic numbers, or a `Sum` that does not
close), and the right side may be any function of `n`:

```@example tut
(
SymJ.rsolve(w(n + 1) - w(n) ~ 1 / (n + 1), w(n), w(0) => 0),
SymJ.rsolve(w(n + 1) ~ (n + 1) * w(n) + 1, w(n), w(0) => 1)
)

```

```@example tut
SymJ.summation((-1)^(k + 1) * binomial(n, k) / k, (k, 1, n))
```

Solutions that are hypergeometric only on every second (or `m`-th) term are found too, as piecewise functions of `mod(n, m)`:

```@example tut
SymJ.summation((-1)^k * binomial(n, k)^2, (k, 0, n))
```

With no closed form, `rsolve` returns the recurrence itself, which is evaluated at integers by unrolling it:

```@example tut
F = SymJ.rsolve((n + 2)^2 * w(n + 2) - (7n^2 + 21n + 16) * w(n + 1) - 8(n + 1)^2 * w(n) ~ 0, w(n), w(0) => 1, w(1) => 2)
```

```@example tut
[subs(F, n => m) for m in 0:6]
```

Harmonic numbers times a term that has an antidifference are summed by parts:

```@example tut
SymJ.summation(k * SymJ.harmonic(k), (k, 1, n))
```

`SymJ.harmonic(n, m)` is the generalized harmonic number `1 + 1/2^m + ... + 1/n^m`:

```@example tut
(
SymJ.summation(1 / k^2, (k, 1, n)),
SymJ.summation(SymJ.harmonic(k) / k, (k, 1, n))
)
```

Whether a series converges can depend on a parameter; the answer then
has cases, which assumptions settle:

```@example tut
SymJ.summation(r^k, (k, 0, Inf))
```

```@example tut
SymJ.summation(k * r^k, (k, 0, Inf); assume = @assume(0 < r < 1))
```
