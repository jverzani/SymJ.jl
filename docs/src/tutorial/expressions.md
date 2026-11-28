# Symbolic expressions

This tutorial works through the use of SymJ to explore concepts from calculus symbolically.
It covers making variables and expressions, evaluating and substituting,
polynomials and rational functions, solving equations, and the core concepts of limits, derivatives and integrals. It originally mirrored the outline of the [Giac page of
Calculus with
Julia](https://github.com/jverzani/CalculusWithJuliaNotes.jl/blob/main/quarto/alternatives/giac.qmd).

SymJ is one of several ways to do symbolic mathematics in Julia. Other packages include:

- [SymPy](https://github.com/JuliaPy/SymPy.jl) connects to the
  Python library SymPy, a full-featured computer algebra system, which
  needs a Python installation.

- [Giac](https://github.com/s-celles/Giac.jl) wraps Giac, a mature
  computer algebra system written in C++.

- [Symbolics](https://github.com/JuliaSymbolics/Symbolics.jl) is a
  numerically oriented system written in Julia, the basis of
  ModelingToolkit.

- [SymEngine](https://github.com/symengine/SymEngine.jl) wraps
  SymEngine, a fast symbolic manipulation library written in C++.

SymJ is written in Julia, aims to be fast, and keeps numbers
exact. Limits are computed by
[GruntzLimit](https://github.com/jverzani/GruntzLimit.jl) (a
dependency) and integrals by
[RischIntegrate](https://github.com/jverzani/RischIntegrate.jl), two
packages that work on the expressions of several symbolic systems,
SymJ's among them.

!!! note "Typeset output in notebooks"
    With [Latexify](https://github.com/korsbo/Latexify.jl) loaded
    (`using Latexify`), SymJ's results also display as LaTeX, as on
    tutorial pages. Notebooks that render LaTeX, such as Jupyter
    (through IJulia) and Pluto, then show them as typeset mathematics;
    `latexify(ex)` gives the LaTeX itself.

Once installed, the SymJ package loads in the typical manner. Below we
also load RischIntegrate, which is called on to perform integration.

```@example tut
using SymJ
using RischIntegrate
```

A few ideas run through the package:

- Symbolic variables are made with `@syms` or `Sym`, and generic Julia code (`+`, `sin`, `^`, …) builds symbolic expressions from them.

- Expressions are kept in a canonical form as they are built.

- Numbers stay exact unless a floating point number is involved.

- Expressions are callable: depending on how they are called, a Julia number or a symbolic expression will be returned.

## Symbolic variables

The `@syms` macro makes variables and binds them to the same names:

```@example tut
@syms x y
```

`Sym(:z)` or `Sym("z")` makes one without binding a name when called with a symbol or string.

```@example tut
z = Sym(:z)
```

(`Sym` is also used to create symbolic numbers.)

A symbolic function is made with `@syms f()`:

```@example tut
@syms f()
f(x^2)
```

(`SymFunction` is also used for this task.)

An indexed family of variables is an array of symbols named by their
indices (the array itself is indexed from 1):

```@example tut
@syms xs[0:3] A[1:2, 1:2]
xs[1], A[2,1]
```

On each page of this tutorial, the following variables and symbolic functions are assumed to be present:

```@example tut
@syms x y z a b c t f()
```


## Symbolic numbers

Symbolic numbers can be created by promotion or construction.

When a symbolic value meets a Julia integer or rational, the number becomes an exact symbolic number:

```@example tut
2x, x / 3, (1 // 2) * x
```

`Sym` makes an exact number directly. Exact numbers keep radicals and π exact, and simplify standard values:

```@example tut
Sym(1 // 3), sqrt(Sym(8)), cos(Sym(π) / 6), Sym(2)^100
```

### Watch out for early conversion

Julia evaluates `1/2` before SymJ sees it, giving the float `0.5`. A float in an
expression makes the result a float, so `x * (1/2)` is not `x/2`:

```@example tut
x / 2 == x * (1 // 2), x / 2 == x * (1 / 2)
```

The same holds for `sqrt(2)` or `pi * 2` evaluated before meeting a symbol. Use `1//2`, `Sym(2)`, or let a symbolic value come first:

```@example tut
sqrt(2) * x, sqrt(Sym(2)) * x, 2 * (pi * x)
```

A common idiom is to use `one(x)` or `zero(x)` instead of `1` or `0`.

## Expressions

Expressions are formed with combinations of symbolic numbers, variables, and deferred function calls.

An expression is a tree: an operation applied to arguments. SymJ implements
TermInterface, so `operation` and `arguments` take it apart; `isexpr` checks if there is an expression tree (or a leaf):

```@example tut
using SymJ: operation, arguments, isexpr
F = 9 // 5 * (y - 32)
operation(F), arguments(F), isexpr(F), isexpr(Sym(9//5))
```

Expressions are built by applying functions to variables. Expressions are kept in a canonical form: like terms are collected, numbers are combined, and products of powers merge.

```@example tut
ex = x + 2x - y + x * y * x
```

```@example tut
(
exp(x) * exp(2x),
(x * y)^3,
sin(-x)
)
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

A power, sum or product is also found where it occurs algebraically, not only as written:
`exp(-x)` in `exp(x)` (as `1/w`), `x*log(x)` in `x^2*log(x)^2`, `a + b` in `2a + 2b + c`:

```@example tut
@syms w
(
subs(exp(x - exp(-x)), exp(-x) => w),
subs(x^2 * log(x)^2, x * log(x) => w),
subs(2a + 2b + c, a + b => w)
)
```

Calling an equation (written with `~`, as in the section on solving equations below)
substitutes into both sides; `lhs` and `rhs` give the sides:

```@example tut
eq = y ~ ex
eq(x => z + 1)
```

```@example tut
lhs(eq), rhs(eq)
```

## Numerical evaluation of symbolic constants


SymJ has both symbolic numbers and symbolic expressions which are numeric constants (like `Sym(2)` and `sqrt(Sym(2))`).

On construction, values at standard points are simplified, returning a symbolic number instead of a constant symbolic expression. For example:

```@example tut
sin(Sym(π) / 4), log(Sym(1)), sin(Sym(2))
```

To work with these symbolic values in the rest of the Julia ecosystem, there are numerous ways to produce related values.

In SymJ, symbolic expressions are subtypes of `Function` (that is, `Sym <: Function`) and have an associated call syntax. What is returned differs by what is called and how. For a symbolic number or numeric constant, the expression is evaluated to a Julia object:

```@example tut
ex, ex1 = Sym(2)^8, sqrt(Sym(2))
ex(), ex1()
```

The type should match what would result if the expression were built non-symbolically.


The `unwrap_const` function is a named functional equivalent.

```@example tut
ex = Sym(pi)^2/4
SymJ.unwrap_const(ex)
```

When an expression is non-constant (it contains a free symbol), both 0-argument evaluation and `unwrap_const` simply return the expression.

The conversion of a symbolic constant to a numeric Julia value can be controlled in different ways, as this table summarizes:

| call | gives | note |
|:-----|:------|:-----|
| `ex()`, `SymJ.unwrap_const(ex)` | the exact Julia value | `Int`, `Rational`, `π`, …; a float where Julia's arithmetic gives one (`sqrt(2)`) |
| `N(ex)` | a `Float64` (`ComplexF64` when not real) | checked in higher precision |
| `float(ex)` | a `Float64` or a `BigFloat` (or their complex) | as Julia's `float`: `BigFloat` when `ex` holds a big integer or rational, or is past `Float64`'s range |
| `SymJ.evalf(ex, T, prec)` | a `T` (`Complex{T}` when not real) | for `BigFloat` type, use `prec` bits |
| `T(ex)` (`Int`, `Rational`, `Float64`, `BigFloat`, …) | a `T` | an `InexactError` when the value is not a `T` |

`N` is a function which returns a `Float64` value for a numeric expression; `float` returns an `AbstractFloat` value, possibly a `BigFloat`. All but `T(ex)` return an expression with a free variable unchanged, so `SymJ.unwrap_const` and `float` can be mapped over a mix of constants and expressions:

```@example tut
map(SymJ.unwrap_const, [Sym(3) / 4, sqrt(Sym(2)), x^2 + x]), float(Sym(pi) / 4)
```

```@example tut
N(Sym(10)^400), float(Sym(10)^400), SymJ.evalf(sqrt(Sym(2)), BigFloat, 100)
```

The value from `N`, `float` and `SymJ.evalf` is computed with more precision where `Float64` (or the requested precision) loses it, so digits lost to cancellation are recovered:

```@example tut
ex = sqrt(Sym(10)^20 + 1) - Sym(10)^10
sqrt(1e20 + 1) - 1e10, N(ex), SymJ.evalf(ex, BigFloat, 64)
```

Julia's number types convert a value that fits them: `Int(Sym(3))`, `Rational(Sym(7)/2)`, `BigFloat(Sym(pi))` at the current precision. Rounding works as for numbers, with a rounding mode, `digits` or `sigdigits` (exactly: a rational), or to an integer type:

```@example tut
(
round(Sym(pi); digits = 2), floor(Int, Sym(7) / 2),
round(Sym(5) / 2, RoundNearestTiesAway)
)
```

`isreal`, `isfinite`, `isinf` and `isnan` answer `true` only when it
is known: `isreal(x)` is `false` for a symbol, which is complex unless
an enclosing `@assuming` block (described in "Solving") says otherwise (as
`isinteger(x)` is `false`). A real expression is finite.

### Rationalizing numbers

The `rationalize` function replaces the floating-point numbers in an expression by nearby rationals:

```@example tut
rationalize(0.5x + 1.25)
```

`SymJ.nsimplify` recognizes floats as exact numbers: rationals, algebraic numbers of low degree in radicals, and rational combinations of given constants:

```@example tut
(
SymJ.nsimplify(0.5x + 1.4142135623730951),
SymJ.nsimplify(1.618033988749895),
SymJ.nsimplify(2.0943951023931953; constants = [pi])
)
```

## Evaluation of expressions with variables

For expressions with a single free symbol, calling that expression with a number evaluates it with Julia's arithmetic on that number, so the result is a Julia number of the matching type:

```@example tut
u = x^2 - 2x + 1
u(3), u(3.0), u(1 // 3), u(big(2.5))
```


Because the arithmetic is generic, automatic differentiation and other
number types work through it.

!!! note "ForwardDiff and compile time"
    A package that adds methods to Base's number functions, as ForwardDiff does,
    invalidates some of SymJ's compiled code, and the first symbolic operations after
    loading it compile again, once per session: on Julia 1.12 a first `sin(x)` took 7 s
    and a first `diff` 15–25 s. Loading ForwardDiff before SymJ (`import ForwardDiff`, then
    `using SymJ`) avoids it for `sin(x)`, not for `diff`.

An expression can also be broadcast over a collection:

```@example tut
u.(0:4)
```

Calling interprets the expression each time. For many evaluations, `lambdify` compiles it into a Julia function, as fast as one written by hand:

```@example tut
g = lambdify(sin(x)^2 + exp(-x / 2), x)
g(1.0)
```


With several variables, the use of positional arguments is ambiguous, so keyword arguments are used to assign values by name. When every variable has a numeric value, the result is a Julia number.

```@example tut
ex = a * x^2 + b * x + c
ex(; x = 2, a = 1, b = -3, c = 2)
```

Names that don't appear are ignored, so one set of values can be used
for several expressions:

```@example tut
vals = (; x = 1.5, y = 2.0)
[e(; vals...) for e in (sqrt(x), x * y, x^y)]
```


When the number of specified values that match is *smaller* than the number of free variables, the values are substituted and a symbolic value returned:

```@example tut
ex = a * x^2 + b * x + c
ex(; x = 2)
```

Substitution is *always* performed when replacements are called using pairs notation:

```@example tut
ex(a => 1, b=>-3, c=>2)
```
The use of pairs is a shortcut notation for `subs`.

## Plotting

Since an expression in zero or one variable is callable, it can be passed where a function is expected, for instance to a plotting package:

```julia
using Plots
plot(sin(x)^2 * exp(-x / 4), 0, 4pi)
```

For plot recipes requiring a multivariate function, a specification of variables might be needed:


```julia
ex = x^2 - 2y^2
xs = ys = range(-3, 3, 100)
contour(xs, ys, (x,y) -> ex(; x, y))
```

## Manipulating symbolic expressions

There are several functions to extract and manipulate parts of expressions.

 The `free_symbols` method returns a vector of all free symbols, empty if the expression is constant. (`SymJ.isconstant` tests for this without allocating.)

```@example tut
ex = x^2 - 2y^2
free_symbols(ex), SymJ.isconstant(ex)
```

`SymJ.gather` groups the terms of an expression by powers of a variable (SymPy's `collect`); `SymJ.coeff` picks out one coefficient:

```@example tut
@syms a b c
ex = a * x^2 + b * x^2 + a * x - b * x + c
```

```@example tut
SymJ.gather(ex, x), SymJ.coeff(ex, x^2)
```

For evaluation, `SymJ.horner` writes a polynomial in nested form, and `SymJ.cse` names the subexpressions that occur more than once:

```@example tut
(
SymJ.horner(3x^4 - 2x^3 + x - 7, x),
SymJ.cse((x + y)^2 + sin(x + y))
)
```

However, the most commonly used manipulations involve expansion and simplification.

## Simplification

Expressions are kept in a canonical form as they are built, but that form does not cancel fractions, combine radicals or apply identities that need facts. `simplify` tries several transformations and keeps the smallest result:

```@example tut
(simplify((x^2 - 1) / (x - 1)),
simplify((x + 1)^2 - x^2),
simplify(sqrt(Sym(2)) * sqrt(Sym(3))))
```

```@example tut
simplify(x / (x + 1) + 1 / (x + 1)), simplify(x^2 + 2x + 1)
```

Trigonometric identities are recognised, including multiple angles and sums of angles:

```@example tut
(simplify(1 - cos(x)^2),
simplify(2sin(x) * cos(x)),
simplify(sin(x) * cos(y) + cos(x) * sin(y)))
```

```@example tut
simplify(sin(3x) - 3sin(x) + 4sin(x)^3), simplify(1 / (1 + tan(x)^2))
```

Hyperbolic functions are treated the same way, and exponentials combine with both:

```@example tut
(simplify(cosh(x)^2 - sinh(x)^2),
simplify(exp(x) + exp(-x)),
simplify(cos(x) + im * sin(x)))
```

```@example tut
(simplify(1 - tanh(x)^2),
simplify(1 + tan(x)^2),
simplify(cosh(x) / sinh(x)))
```

Logarithms combine and radicals denest; with symbols, where the assumptions allow:

```@example tut
(
simplify(log(Sym(2)) + log(Sym(3))),
simplify(sqrt(3 + 2sqrt(Sym(2)))),
simplify(1 / (1 + sqrt(Sym(2))))
)
```

### Assumptions

Identities such as `sqrt(x^2) = x` hold only for some values.
In most uses, symbolic variables are complex unless the assumptions say otherwise, so without them only the identities that hold for complex numbers apply. Assumptions can be given to a function call to change this behavior. Facts are specified through assumptions, which are passed to the `assume` keyword argument. In the following, we use the `@assume` macro which affords a more familiar syntax (`x > 0` evaluates to a `Bool`, but in the macro is the same as the symbolic statement `x ≫ 0`).

```@example tut
(
simplify(sqrt(x^2)),
simplify(sqrt(x^2); assume = @assume(x > 0)),     # or assume = (x ≫ 0)
simplify(sqrt(x^2); assume = @assume(x in Reals)) # or assume = (x ∈ ℝ)
)
```

The Unicode shortcut `ℝ` (`\bbR<tab>`) can be used instead of the named domain, as can shortcuts for other elements of `SymJ.NumberSet`: `Integers, Rationals, and Complexes`.



```@example tut
(
simplify(log(x) + log(y); assume = @assume(x > 0, y > 0)),
simplify(sqrt(x^2 + 2x + 1); assume = @assume(x > -1))
)
```

`real`, `imag`, `conj` and `angle` evaluate constants exactly; for
symbols, which are complex, they take apart what holds for every
value, and `refine` does the rest under assumptions:

```@example tut
real(Sym(2) + 3im), conj(x + im), refine(real(exp(im * x)); assume = @assume(x in Reals))
```

`refine` applies only those identities, without trying anything else:

```@example tut
@syms n
(
refine(log(exp(x)); assume = @assume(x in Reals)),
refine(sin(x + 2n*Sym(π)); assume = @assume(n in Integers))
)
```

### More fine-grained simplification

The transformations `simplify` tries are also available on their own,
unexported: `SymJ.simplify_rational`, `SymJ.simplify_powers`,
`SymJ.simplify_trig`, `SymJ.simplify_hyperbolic`, `SymJ.simplify_log`,
`SymJ.simplify_radicals` and `SymJ.simplify_gamma`. `SymJ.radsimp` writes an
expression with a denominator free of radicals, whatever the size:

```@example tut
(
SymJ.radsimp(1 / (1 + sqrt(Sym(2)) + sqrt(Sym(3)))),
SymJ.radsimp(1 / (1 + Sym(2)^(Sym(1)//3))),
SymJ.radsimp(1 / (sqrt(x) + 1))
)
```

`SymJ.simplify_gamma` combines
factorials, binomials and gamma functions whose arguments differ by
integers, and writes factorial ratios as binomials:

```@example tut
@syms k
(
SymJ.simplify_gamma(factorial(n + 2) / factorial(n)),
SymJ.simplify_gamma(binomial(n + 1, k) / binomial(n, k)),
SymJ.simplify_gamma(factorial(n) / (factorial(k) * factorial(n - k)))
)
```

With SpecialFunctions' `gamma` it also applies the reflection formula
(`gamma(x)gamma(1 - x) = π/sin(πx)`) and the duplication formula.


Expansion goes the other way. `SymJ.expand_trig` writes trigonometric
and hyperbolic functions of sums and integer multiples through
functions of the terms:

```@example tut
(
SymJ.expand_trig(sin(x + y)),
SymJ.expand_trig(cos(3x)),
SymJ.expand_trig(tan(2x))
)
```

`SymJ.expand_log` splits logarithms of products and powers where the assumptions allow it; `force = true` takes every base as positive:

```@example tut
(
SymJ.expand_log(log(x^2 * y); assume = @assume(x > 0, y > 0)),
SymJ.expand_log(log(x * y)),
SymJ.expand_log(log(x / y); force = true)
)
```

## Rewriting parts of an expression

`SymJ.rewrite(ex, target)` writes the functions of an expression
through another one: trigonometric and hyperbolic functions through
`exp`, inverse ones through `log`, `tan` and `sec` through `sin` and
`cos`, exponentials with imaginary parts through `cos` and `sin`,
factorials through `gamma` and back. `fs => target` rewrites only the
functions `fs`:

```@example tut
(
SymJ.rewrite(sin(x), exp),
SymJ.rewrite(exp(im * x), cos),
SymJ.rewrite(asinh(x), log),
SymJ.rewrite(tan(x) + cos(x), tan => sin)
)
```

### Pattern matching

`SymJ.pattern` (not exported) makes a pattern from a quoted
expression, with `~x` marking a slot (the syntax comes from the
SymbolicUtils package). `match` gives the first match as a `Dict`,
`eachmatch` every match, and `occursin` whether a part of an
expression matches:

```@example tut
p = SymJ.pattern(:(~!a * sin(~x) + ~~r))
match(p, 3sin(2x) + x + 1)
```

In a sum or product a slot takes a group of terms (factors), so there
can be several matches; `~!a` is 1 when there is no factor for it:

```@example tut
@syms a b c
eachmatch(SymJ.pattern(:(~x + ~y)), a + b + c)
```

```@example tut
match(p, sin(x)), occursin(SymJ.pattern(:(cos(~x)^2)), 1 + cos(x)^2)
```

A head slot matches a function and binds it:

```@example tut
match(SymJ.pattern(:((~f)(~x)^2)), cos(x)^2)
```

A rule rewrites what a pattern matches. `SymJ.@rule lhs => rhs` takes
the right side as Julia code, with `~x` the matched value;
`SymJ.rewrite` applies rules to every part of an expression until
nothing changes, and `replace` applies them once:

```@example tut
r = SymJ.@rule log(~x * ~y) => log(~x) + log(~y)
SymJ.rewrite(log(a * b * c), r), replace(log(a * b * c), r)
```

```@example tut
pyth = SymJ.@rule sin(~x)^2 + cos(~x)^2 + ~~r => 1 + ~r
SymJ.rewrite(3 + sin(2x)^2 + cos(2x)^2 + x, pyth)
```

`replace` also takes `old => new` pairs, applied in one pass from the
leaves up: a quoted pattern and template, a predicate and a function
of the subexpression, or a subexpression and its value:

```@example tut
replace(sin(a) + sin(b)^2, :(sin(~u)) => :(cos(~u))),
replace(x^2 + 3x + 5, (e -> SymJ.isnumber(e) && e > 2) => (e -> e - 1))
```

A predicate on a slot, or a right side that gives `nothing`, limits
where a rule applies:

```@example tut
d = SymJ.@rule (~x)^(~n::isinteger) => ~n * (~x)^(~n - 1)
d(x^3), d(x^a)
```


## Polynomials

Polynomial expressions have several methods. Some are defined on symbolic expressions that are polynomial and others only for expressions of the `Poly` type, which is used to more efficiently represent sparse, multivariate polynomials.

For example, a polynomial expression in factored form is not expanded by default, but can be expanded with `expand`:

```@example tut
p = expand((x - 1) * (x - 2) * (x + 3))
```
Related is `factor`, which factors over the rationals, in any number of variables:

```@example tut
(
factor(p),
factor(x^4 + 4),
factor(x^6 - y^6)
)
```


Assumptions may change how an expression factors:

```@example tut
ex = sqrt(1 - x^2) + sqrt(1 - x)
(
factor(ex),
factor(ex; assume = (x ∈ ℝ,))
)
```

(Without facts, symbols are complex, and
$\sqrt{ab} = \sqrt{a}\sqrt{b}$ can fail.)

`SymJ.factor_squarefree` only separates repeated factors, which is cheaper:

```@example tut
SymJ.factor_squarefree(expand((x - 1)^3 * (x + 2)^2))
```

Other polynomial operations include the least common multiple (`lcm`) and the greatest common divisor (`gcd`) of polynomials:

```@example tut
gcd(expand((x - 1)^2 * (x + y)), expand((x - 1) * (x + y)^2))
```


`SymJ.decompose` writes a polynomial as a composition of polynomials
that do not decompose further, `p = g₁(g₂(⋯))`. `solve` uses it: an
irreducible polynomial of degree above 4 whose components have degree
at most 4 has its roots found through them, exactly:

```@example tut
q = expand((x^2 + x)^3 - 2)
SymJ.decompose(q, x), solve(q ~ 0, x)[1:2]
```

(The `solve` function is introduced later in this tutorial.)

`SymJ.minimal_polynomial` gives the irreducible polynomial with
integer coefficients that an algebraic number is a root of:

```@example tut
(
SymJ.minimal_polynomial(sqrt(Sym(2)) + sqrt(Sym(3)), x),
SymJ.minimal_polynomial(cos(Sym(pi) / 7), x)
)
```

### Construction of Poly objects

`Poly` converts an expression into a sparse, multivariate polynomial in given generators. Its coefficients and exponents are available; `Sym` converts it back to an expression.

```@example tut
P = Poly(p, x)
coefficients(P), exponent_vectors(P), degree(P)
```

Other symbols become coefficients:

```@example tut
Q = Poly(a * x^2 + b * x + c, x)
coefficients(Q)
```

Some evaluations are more performant using the `Poly` type, such as those involving the expansion of high powers of a polynomial.

### Evaluation and shifts

A polynomial is evaluated by calling it; a shift is a substitution:

```@example tut
P(2), expand(P(x => x + 1))
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

The `numerator` and `denominator` generics return the top and bottom of a rational expression. The function `SymJ.numer_denom` gives the numerator and denominator without reducing:

```@example tut
SymJ.numer_denom(r)
```

### Partial fractions

`apart` decomposes a rational function into a polynomial part and
simple fractions:

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

### Resultants

The resultant of two polynomials in `x` is zero exactly when they have
a common root. It is a polynomial in the other variables, so it
eliminates `x` from a system:

```@example tut
resultant(x^2 + y^2 - 1, x - y, x)
```

The discriminant is zero exactly when a polynomial has a repeated root:

```@example tut
@syms a b c
discriminant(a*x^2 + b*x + c, x), discriminant(x^3 + a*x + b, x)
```

### Gröbner bases

A Gröbner basis of a set of polynomials generates the same ideal and
has the same common zeros. In the lex order, the default, it is
triangular: the last polynomials involve only the last variables.

```@example tut
groebner([x^2 + y^2 - 1, x*y - 1//4], x, y)
```

The basis is `[1]` exactly when the polynomials have no common zero:

```@example tut
groebner([x^2 + y^2 - 1, x^2 + y^2 - 4], x, y)
```

`reduced` divides by the basis. The remainder is zero exactly when a
polynomial is in the ideal:

```@example tut
G = groebner([x^2 + y^2 - 1, x - y], x, y)
reduced(x^2 - y^2, G, x, y), reduced(x^3, G, x, y)
```

## Piecewise functions

`SymJ.piecewise(c₁ => v₁, c₂ => v₂, …)` (not exported) is the value of
the first piece whose condition holds; a condition is a relation, a
tuple of relations that all hold, or `true`. It is an expression like
any other: substitution picks a piece, and `diff`, `limit`, `refine`
and `lambdify` take it piece by piece:

```@example tut
pw = SymJ.piecewise(x ≪ 0 => -x, (0 ≦ x, x ≪ 1) => x^2, true => 1)
pw, pw(x => -2), diff(pw, x)
```

```@example tut
limit(pw, x => 0), limit(pw, x => 1; dir = :-)
```

The unit step `SymJ.heaviside` (½ at 0) and the delta `SymJ.dirac` are
its derivative pair; `SymJ.piecewise(ex)` writes `heaviside`, `sign`,
`abs`, `max` and `min` as pieces:

```@example tut
H = SymJ.heaviside
diff(H(x - 1), x), SymJ.piecewise(abs(x - 1))
```

`SymJ.piecewise_fold` combines piecewise functions in a sum, product
or other expression into one, merging neighbouring pieces that agree;
`simplify` does it when the result is smaller (a definite integral
with a parameter gives such sums):

```@example tut
SymJ.piecewise_fold(abs(x) + abs(x - 1), x)
```

With parameters whose order is not known, the answer is by cases:

```@example tut
@syms a b
SymJ.piecewise_fold(abs(x - a) + abs(x - b), x)
```
