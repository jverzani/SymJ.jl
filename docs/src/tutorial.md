# Tutorial

This page works through the symbolic side of calculus with SymJ: making variables and
expressions, evaluating and substituting, polynomials and rational functions, solving
equations, and limits, derivatives and integrals. It follows the outline of the
[Giac page of Calculus with Julia](https://github.com/jverzani/CalculusWithJuliaNotes.jl/blob/main/quarto/alternatives/giac.qmd).

SymJ is one of several ways to do symbolic mathematics in Julia:

- [SymPy.jl](https://github.com/JuliaPy/SymPy.jl) connects to the Python library SymPy,
  a full-featured computer algebra system, which needs a Python installation.
- [Giac.jl](https://github.com/s-celles/Giac.jl) wraps Giac, a mature computer algebra
  system written in C++.
- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) is a numerically
  oriented system written in Julia, the basis of ModelingToolkit.
- [SymEngine.jl](https://github.com/symengine/SymEngine.jl) wraps SymEngine, a fast
  symbolic manipulation library written in C++.

SymJ is written in Julia, aims to be fast, and keeps numbers exact. Limits are computed
by [GruntzLimit](https://github.com/jverzani/GruntzLimit.jl) (a dependency) and integrals
by [RischIntegrate](https://github.com/jverzani/RischIntegrate.jl), two packages that work
on the expressions of several symbolic systems, SymJ's among them.

!!! note "Julia versions"
    SymJ and GruntzLimit support Julia 1.10 and later. RischIntegrate needs Julia 1.12 or
    later, so the section on integrals needs that version.

```@example tut
using SymJ
using RischIntegrate
```

A few ideas run through the package:

- Variables are made with `@syms`, and generic Julia code (`+`, `sin`, `^`, …) builds
  expressions from them.
- Expressions are kept in a canonical form as they are built.
- Numbers stay exact unless a floating point number is involved.
- Expressions are callable: with numbers they evaluate, with pairs they substitute.

## Symbolic variables

The `@syms` macro makes variables and binds them to the same names:

```@example tut
@syms x y
```

`Sym(:z)` makes one without binding a name, and `symbols` makes several from a string:

```@example tut
z = Sym(:z)
a, b, c = symbols("a b c")
```

A symbolic function, with no rules of its own, is made with `@syms f()`:

```@example tut
@syms f()
f(x^2)
```

## Symbolic numbers

When a symbolic value meets a Julia integer or rational, the number becomes an exact
symbolic number:

```@example tut
2x, x / 3, (1 // 2) * x
```

`Sym` makes an exact number directly. Exact numbers keep radicals and π exact, and
simplify standard values:

```@example tut
Sym(1 // 3), sqrt(Sym(8)), cos(Sym(π) / 6), Sym(2)^100
```

### Watch out for early conversion

Julia evaluates `1/2` before SymJ sees it, giving the float `0.5`. A float in an
expression makes the result a float, so `x * (1/2)` is not `x/2`:

```@example tut
x / 2 == x * (1 // 2), x / 2 == x * (1 / 2)
```

The same holds for `sqrt(2)` or `pi * 2` evaluated before meeting a symbol. Use `1//2`,
`Sym(2)`, or let a symbolic value come first:

```@example tut
sqrt(2) * x, sqrt(Sym(2)) * x, 2 * (pi * x)
```

## Expressions

Expressions are built by applying functions to variables. They are kept in a canonical
form: like terms are collected, numbers are combined, and products of powers merge.

```@example tut
ex = x + 2x - y + x * y * x
```

```@example tut
exp(x) * exp(2x), (x * y)^3, sin(-x)
```

Functions with no numeric argument stay unevaluated, while exact values at standard
points are simplified:

```@example tut
sin(x), sin(Sym(π) / 4), log(Sym(1)), sin(2)
```

An expression is a tree: an operation applied to arguments. SymJ implements
TermInterface, so `operation` and `arguments` take it apart:

```@example tut
using SymJ: operation, arguments
F = 9 // 5 * (y - 32)
operation(F), arguments(F)
```

## Substitution

`subs` replaces symbols, or whole subexpressions, by other values:

```@example tut
ex = a * x^2 + b * x + c
subs(ex, x => y - 1)
```

```@example tut
subs(ex, a => 1, b => -2, c => 1)
```

Calling an expression with pairs does the same:

```@example tut
ex(x => 2)
```

## Numerical evaluation

Calling an expression with a number evaluates it with Julia's arithmetic on that number,
so the result is a Julia number of the matching type:

```@example tut
u = x^2 - 2x + 1
u(3), u(3.0), u(1 // 3), u(big(2.5))
```

Because the arithmetic is generic, automatic differentiation and other number types work
through it. An expression can also be broadcast over a collection:

```@example tut
u.(0:4)
```

With several variables, keywords give their values by name. When every variable has a
numeric value the result is a Julia number, and otherwise
the expression with the given values substituted:

```@example tut
ex = a * x^2 + b * x + c
ex(; x = 2, a = 1, b = -3, c = 2), ex(; x = 2)
```

Names that don't appear are ignored, so one set of values can be used for several
expressions:

```@example tut
vals = (; x = 1.5, y = 2.0)
[e(; vals...) for e in (x + y, x * y, x^y)]
```

`N` gives the floating point value of a constant expression, and calling it with no
arguments gives a Julia number:

```@example tut
N(sqrt(Sym(2)) + Sym(π)), (Sym(3) // 4)()
```

Calling interprets the expression each time. For many evaluations, `lambdify` compiles
it into a Julia function, as fast as one written by hand:

```@example tut
g = lambdify(sin(x)^2 + exp(-x / 2), x)
g(1.0)
```

## Plotting

Since an expression in one variable is callable, it can be passed where a function is
expected, for instance to a plotting package:

```julia
using Plots
plot(sin(x)^2 * exp(-x / 4), 0, 4pi)
```

## Polynomials

### Construction

A polynomial expression can be expanded with `expand`:

```@example tut
p = expand((x - 1) * (x - 2) * (x + 3))
```

`Poly` converts an expression into a polynomial in given generators. Its coefficients
and exponents are available, and `Sym` converts it back:

```@example tut
P = Poly(p, x)
coefficients(P), exponent_vectors(P), degree(P)
```

Other symbols become coefficients:

```@example tut
Q = Poly(a * x^2 + b * x + c, x)
coefficients(Q)
```

### Evaluation and shifts

A polynomial is evaluated by calling it; a shift is a substitution:

```@example tut
p(2), expand(p(x => x + 1))
```

### Factoring

`factor` factors over the rationals, in any number of variables:

```@example tut
factor(p), factor(x^4 + 4), factor(x^6 - y^6)
```

`factor_squarefree` only separates repeated factors, which is cheaper:

```@example tut
factor_squarefree(expand((x - 1)^3 * (x + 2)^2))
```

The greatest common divisor of polynomials:

```@example tut
gcd(expand((x - 1)^2 * (x + y)), expand((x - 1) * (x + y)^2))
```

### Resultants

The resultant of two polynomials in `x` is zero exactly when they have a common root.
It is a polynomial in the other variables, so it eliminates `x` from a system:

```@example tut
resultant(x^2 + y^2 - 1, x - y, x)
```

The discriminant is zero exactly when a polynomial has a repeated root:

```@example tut
@syms a b c
discriminant(a*x^2 + b*x + c, x), discriminant(x^3 + a*x + b, x)
```

### Gröbner bases

A Gröbner basis of a set of polynomials generates the same ideal and has the same
common zeros. In the lex order, the default, it is triangular: the last polynomials
involve only the last variables.

```@example tut
groebner([x^2 + y^2 - 1, x*y - 1//4], x, y)
```

The basis is `[1]` exactly when the polynomials have no common zero:

```@example tut
groebner([x^2 + y^2 - 1, x^2 + y^2 - 4], x, y)
```

`reduced` divides by the basis. The remainder is zero exactly when a polynomial is in
the ideal:

```@example tut
G = groebner([x^2 + y^2 - 1, x - y], x, y)
reduced(x^2 - y^2, G, x, y), reduced(x^3, G, x, y)
```

### Roots

`solve` finds the roots of a polynomial. Degrees up to four are solved exactly, with
real roots of cubics given in trigonometric form (no complex numbers needed):

```@example tut
solve(p ~ 0, x)
```

```@example tut
solve(x^3 - 3x + 1 ~ 0, x)
```

```@example tut
solve(x^4 - 10x^2 + 1 ~ 0, x)
```

From degree five, roots are computed numerically:

```@example tut
solve(x^5 - x + 1 ~ 0, x)
```

## Rational functions

`cancel` reduces a rational expression to lowest terms:

```@example tut
r = (x^3 - 1) / (x^2 - 1)
cancel(r)
```

```@example tut
cancel(1 / x + 1 / (x + 1))
```

`numer_denom` gives the numerator and denominator without reducing:

```@example tut
numer_denom(r)
```

### Partial fractions

`apart` decomposes a rational function into a polynomial part and simple fractions:

```@example tut
apart((x^3 + 1) / (x^2 - 1), x)
```

```@example tut
apart(1 / (x^4 - 1), x)
```

Other symbols are allowed in the coefficients:

```@example tut
apart(1 / ((x - a) * (x - b)), x)
```

## Simplification

Expressions are kept in a canonical form as they are built, but that form does not
cancel fractions, combine radicals or apply identities that need facts. `simplify`
tries several transformations and keeps the smallest result:

```@example tut
simplify((x^2 - 1) / (x - 1)), simplify((x + 1)^2 - x^2), simplify(sqrt(Sym(2)) * sqrt(Sym(3)))
```

```@example tut
simplify(x / (x + 1) + 1 / (x + 1)), simplify(x^2 + 2x + 1)
```

Identities such as `sqrt(x^2) = x` hold only for some values. Symbols are complex unless
a fact says otherwise, so these need assumptions:

```@example tut
simplify(sqrt(x^2)), simplify(sqrt(x^2); assume = @assume(x > 0)), simplify(sqrt(x^2); assume = @assume(x in Reals))
```

Trigonometric identities are recognised, including multiple angles and sums of angles:

```@example tut
simplify(1 - cos(x)^2), simplify(2sin(x) * cos(x)), simplify(sin(x) * cos(y) + cos(x) * sin(y))
```

```@example tut
simplify(sin(3x) - 3sin(x) + 4sin(x)^3), simplify(1 / (1 + tan(x)^2))
```

Hyperbolic functions are treated the same way, and exponentials combine with both:

```@example tut
simplify(cosh(x)^2 - sinh(x)^2), simplify(exp(x) + exp(-x)), simplify(cos(x) + im * sin(x))
```

```@example tut
simplify(1 - tanh(x)^2), simplify(1 + tan(x)^2), simplify(cosh(x) / sinh(x))
```

Logarithms combine and radicals denest; with symbols, where the assumptions allow:

```@example tut
simplify(log(Sym(2)) + log(Sym(3))), simplify(sqrt(3 + 2sqrt(Sym(2)))), simplify(1 / (1 + sqrt(Sym(2))))
```

```@example tut
simplify(log(x) + log(y); assume = @assume(x > 0, y > 0)), simplify(sqrt(x^2 + 2x + 1); assume = @assume(x > -1))
```

`real`, `imag`, `conj` and `angle` evaluate constants exactly; for symbols, which are
complex, they take apart what holds for every value, and `refine` does the rest under
assumptions:

```@example tut
real(Sym(2) + 3im), conj(x + im), refine(real(exp(im * x)); assume = @assume(x in Reals))
```

`refine` applies only those identities, without trying anything else:

```@example tut
@syms n
refine(log(exp(x)); assume = @assume(x in Reals)), refine(sin(x + 2n*Sym(π)); assume = @assume(n in Integers))
```

The transformations `simplify` tries are also available on their own, unexported:
`SymJ.simplify_rational`, `SymJ.simplify_powers`, `SymJ.simplify_trig`,
`SymJ.simplify_hyperbolic`, `SymJ.simplify_log` and `SymJ.simplify_radicals`.

## Solving equations

An equation is written with `~`, which binds more loosely than everything else, so no
parentheses are needed:

```@example tut
solve(x^2 + 2x ~ 3, x)
```

Equations in which `x` appears through one function are solved by solving for that
function first:

```@example tut
solve(exp(2x) - 3exp(x) + 2 ~ 0, x)
```

```@example tut
solve(sqrt(x + 1) ~ x - 1, x)
```

Several logarithms of the unknown are combined; roots where a logarithm would be complex
(here `x = -2`) are dropped:

```@example tut
solve(log(x) + log(x - 1) ~ log(Sym(6)), x)
```

Some equations need the Lambert W function:

```@example tut
solve(x * exp(x) ~ 1, x)
```

Symbols other than the unknown are allowed:

```@example tut
solve(a * x^2 + b * x + c ~ 0, x)
```

### Inequalities

`<` and its relatives compare numbers and return `true` or `false`, so inequalities use
other operators: `≪` (typed `\ll`), `≦` (`\leqq`), `≫` (`\gg`), `≧` (`\geqq`) and `≶`
(`\lessgtr`). The solution is a union of intervals:

```@example tut
solve(x^2 - 1 ≫ 0, x)
```

```@example tut
solve((x - 1) / (x + 2) ≧ 0, x)
```

### Systems

A system of equations gives a vector of solutions, each a `Dict` from the unknowns to
their values. A linear system has at most one:

```@example tut
solve([x + y ~ 3, x - y ~ 1], [x, y])
```

```@example tut
solve([a * x + y ~ 1, x - y ~ b], [x, y])
```

Polynomial systems are solved through a lex Gröbner basis, which is triangular:

```@example tut
solve([x^2 + y^2 ~ 1, x ~ y], [x, y])
```

```@example tut
solve([x*y ~ 1, x + y ~ 3], [x, y])
```

With infinitely many solutions, the unknowns left free map to themselves:

```@example tut
solve([x^2 + y^2 ~ 1], [x, y])
```

A system of linear inequalities gives a description of the region, each unknown between
bounds in the ones after it:

```@example tut
solve([x ≧ 0, y ≧ 0, x + y ≦ 1], [x, y])
```

### Assumptions

Symbols are complex numbers unless a fact says otherwise. Facts are relations, such as
`x > 0`, or membership in one of `Integers ⊂ Rationals ⊂ Reals ⊂ Complexes` (also
written `ℤ ⊂ ℚ ⊂ ℝ ⊂ ℂ`), such as `n in Integers`. The `@assume` macro writes them with
Julia's operators, including `&&`, `||`, `!` and chained comparisons. They are given to
`solve` with the `assume` keyword, and only the solutions they allow are kept:

```@example tut
solve(x^2 ~ 4, x; assume = @assume(x > 0))
```

```@example tut
solve(x^2 + 1 ~ 0, x; assume = @assume(x in Reals))
```

Without the macro, relations use the symbolic operators: `x ≫ 0` (or `Gt(x, 0)`), and
`n ∈ ℤ` for membership.

```@example tut
solve(x^2 ~ 2, x; assume = (x ∈ ℚ,))
```

`@assuming` gives facts for a block of code (as does the function `assuming`, with a
`do` block):

```@example tut
@assuming x > 0 begin
    solve(x^3 ~ 8, x), solve(x^2 - 5x + 6 ~ 0, x; assume = @assume(x > 2))
end
```

When the answer depends on something the facts leave open, the result lists the cases:

```@example tut
solve(x^2 ~ a, x; assume = @assume(x in Reals, a in Reals))
```

Without `a ∈ ℝ`, `a` may be complex, so the cases are a fact and its negation:

```@example tut
cs = solve(x^2 ~ a, x; assume = @assume(x > 0))
```

Calling the cases with a value gives the answer for it:

```@example tut
cs(a => 4), cs(a => -1)
```

Contradictory facts are an error (`@assume(x > 1 && x < 0)`), and the `Assumptions`
value `@assume` makes can be reused: `A = @assume(x > 0, n in Integers)`, then
`assume = A`.

## Matrices

Matrices of expressions are Julia arrays of `Sym`. With `LinearAlgebra`, `det`, `inv`,
`lu` and `\` work exactly, pivoting on the first nonzero entry rather than comparing
magnitudes:

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

Row reduction, rank, null spaces and eigenvalues (with multiplicity) and eigenvectors:

```@example tut
rref(Sym[1 2 3; 4 5 6; 7 8 9]), nullspace(Sym[1 2 3; 4 5 6; 7 8 9])
```

```@example tut
eigvals(Sym[1 2; 3 4]), eigvecs([a b; b a])
```

`A - I` works as for numbers; a symbolic multiple of the identity is `λ*I(n)`:

```@example tut
@syms λ
det([x 1; 1 x] - λ * I(2))
```

## Limits

`limit(f, x => x0)` is the limit of `f` as `x` tends to `x0`, a number, an expression or
`Inf`. It is computed by Gruntz's algorithm, through GruntzLimit. An infinite limit is
`Sym(Inf)` or `Sym(-Inf)`:

```@example tut
limit(sin(x) / x, x => 0)
```

```@example tut
limit((1 + 1 / x)^x, x => Inf), limit(x * ((1 + 1 / x)^x - ℯ), x => Inf)
```

```@example tut
limit((sin(x) - sin(Sym(π) / 6)) / (x - Sym(π) / 6), x => Sym(π) / 6)
```

The limit is two-sided unless `dir = :+` or `dir = :-` asks for a one-sided one:

```@example tut
limit(x * exp(1 / x), x => 0; dir = :-), limit(x^x, x => 0; dir = :+)
```

The algorithm handles cancellations that numerical evaluation cannot follow:

```@example tut
limit(exp(x) * (exp(1 / x - exp(-x)) - exp(1 / x)), x => Inf)
```

Other symbols may appear:

```@example tut
limit((1 + a / x)^(b * x), x => Inf)
```

The other symbols are taken to be real numbers. When the limit depends on their values,
`limit` gives it for each case:

```@example tut
limit(x^a, x => Inf)
```

Assumptions, as in the section on solving, settle cases:

```@example tut
limit(x^a, x => Inf; assume = @assume(a > 1)), limit((a^x + b^x)^(1 / x), x => Inf; assume = @assume(a > b > 0))
```

## Derivatives

`diff(f, x)` differentiates with respect to `x`, and `diff(f, x, n)` gives the `n`-th
derivative:

```@example tut
diff(sin(x)^2 * exp(x), x)
```

```@example tut
diff(x^x, x), diff(atan(x), x, 3)
```

Mixed partial derivatives take the variables in order: `diff(f, x, y)` differentiates
with respect to `x`, then `y`:

```@example tut
diff(x^3 * y^2, x, y), diff(sin(x * y), x, y), diff(x^2 * y * z, x, x, y, z)
```

A symbolic function gives an unevaluated derivative:

```@example tut
diff(f(x^2), x)
```

### Critical points

The critical points of a function are the zeros of its derivative, and the second
derivative classifies them:

```@example tut
h = x^3 - 6x^2 + 9x + 1
cps = solve(diff(h, x) ~ 0, x)
```

```@example tut
h2 = diff(h, x, 2)
[(c, h(c), h2(c)) for c in cps]
```

So 1 is a local maximum and 3 a local minimum.

### Tangent lines

The tangent line at `c` is `h(c) + h'(c)(x - c)`:

```@example tut
c₀ = 2
expand(h(c₀) + diff(h, x)(c₀) * (x - c₀))
```

## Series

`series(ex, x, x0, n)` expands `ex` about `x0` with an error of order `(x - x0)^n`;
`x0` defaults to 0 and `n` to 6:

```@example tut
series(sin(x), x)
```

```@example tut
series(exp(x), x, 1, 3)
```

Laurent series (negative powers), Puiseux series (fractional powers), logarithmic terms
and expansions at infinity (in powers of `1/x`) are included:

```@example tut
series(cot(x), x, 0, 4), series(sqrt(x + x^2), x, 0, 3)
```

```@example tut
series(x^x, x, 0, 3), series(sqrt(x^2 + x), x, Inf, 3)
```

A series is an object that keeps its error term. Series about the same point can be
added, multiplied and divided, with the order tracked, and `Sym` drops the `O` term:

```@example tut
s, c = series(sin(x), x), series(cos(x), x)
s / c
```

```@example tut
Sym(s)
```

## Integrals

Antiderivatives come from `integrate`, from RischIntegrate, an implementation of the
Risch algorithm. It decides whether an elementary antiderivative exists, and finds it:

```@example tut
integrate(x^2 * exp(x), x)
```

```@example tut
integrate(1 / (x^3 + 1), x)
```

```@example tut
integrate(1 / (x^3 - x), x), integrate((2x + 3) / (x^2 + 1), x)
```

```@example tut
integrate(sin(x)^3, x), integrate(sqrt(1 - x^2), x)
```

```@example tut
integrate(x * log(x), x), integrate(1 / (x * log(x)), x)
```

### Definite integrals

A definite integral is the difference of an antiderivative's values at the end points.
Calling the antiderivative with exact numbers keeps the answer exact:

```@example tut
F = integrate(x^2, x)
F(3) - F(0)
```

```@example tut
G = integrate(sin(x), x)
G(x => Sym(π)) - G(x => 0)
```
