# SymJ

!!! warning "Work in progress"
    SymJ is under active development. Its interface may change, including in ways that break existing code, before a stable release.

SymJ is a symbolic manipulation kernel written in pure Julia. It
builds expressions in a canonical form as you go, keeps numbers exact,
and is fast: construction, expansion, differentiation and substitution
are within a small factor of SymEngine, and usually
faster. Expressions are callable, so they can stand in for functions
in numeric code, and they implement
[TermInterface](https://github.com/JuliaSymbolics/TermInterface.jl),
so generic term-rewriting code, such as
[GruntzLimit](https://github.com/jverzani/GruntzLimit.jl) for limits
and [RischIntegrate](https://github.com/jverzani/RischIntegrate.jl)
for integrals, works on them.

- The tutorial works through the symbolic side of calculus with SymJ, in four pages: [Expressions](tutorial/expressions.md), [Solving](tutorial/solving.md), [Calculus](tutorial/calculus.md) and [Matrices](tutorial/matrices.md).

- [Implementation](@ref) discusses the design choices: one node type, exact numbers, packed polynomials, and the gcd and factoring algorithms.

- The [Design notes](@ref) record each decision in detail, with measurements.

- The [Reference](@ref) lists the functions. Most functions are not exported and need qualifying (`SymJ.name`).

## General usage

```@example usage
using SymJ
```

(The package loads relatively quickly, though it can take a while to precompile.)

### Symbolic variables and numbers

Variables are made with `@syms`, or with `Sym` from a name. Numbers
stay exact: `Sym` makes an exact number, and arithmetic with symbols
keeps integers and rationals exact.

```@example usage
@syms x y z
p = (x + 1)^2 / 3
```

```@example usage
sqrt(Sym(8)), sin(Sym(π) / 3), Sym(2)^100
```

A floating point number in an operation makes the result a float, so
use `1//2`, not `1/2`, for an exact half:

```@example usage
x / 2 == x * (1 // 2), x * (1 / 2)
```

A useful idiom is to use `one(x)/2`, for some variable `x` say, in
contrast to `1/2`.

### Expressions are built in canonical form

Sums collect like terms and products collect like bases as they are
built, with no separate simplification step:

```@example usage
x + 2x - y + x * y * x, exp(x) * exp(y), (x * y)^2
```

### Calling an expression

The basic `Sym` type is a subtype of `Function`. `Sym` values are callable, though with some idiosyncratic conventions.

* A constant expression called with no arguments gives its value as a Julia number.

* An expression of just one variable may be called using the first
positional argument and returns a Julia number. It is evaluated with
Julia's own arithmetic, so it works in numeric code (with dual
numbers, `BigFloat`, intervals, …). (Cf. [`lambdify`](@ref) for more performant evaluation.)

* For expressions with multiple variables, keyword arguments can be used to specify the values. If all the symbolic variables are substituted, a Julia number is returned. If some variables are left unspecified, a symbolic value is returned with the specified variables substituted through [`subs`](@ref).

* Use pairs to substitute in values, always leaving a symbolic value.

```@example usage
(Sym(2)^10)()
```

```@example usage
u = x^2 - 2x
u(3), u(3.0), u(1 // 2)
```

```@example usage
ex = x^2 * y + sin(x)
ex(; x = 2, y = 3), ex(; x = 2)
```

```@example usage
ex(x => 2), subs(ex, x => y + 1)
```




For speed in a loop, `lambdify` compiles an expression into a Julia function:

```@example usage
ex = sin(x)^2 + 3x^5 - exp(-x / 2)
f = lambdify(ex, x)
ex(0.5), f(0.5)      # the ex call is much costlier
```

For symbolic constants, a call without a value will return a Julia number (as will `SymJ.unwrap_const`) with a matching type; the `N` function gives a `Float64`:

```@example usage
ex = one(x)/2
ex(), N(ex)
```

```@example usage
ex = Sym(10)^309
N(ex), ex() |> log10
```

The unexported `nsimplify` function is a partial inverse to `N`; it will identify some numeric values as symbolic constants.

```@example usage
SymJ.nsimplify(1.4142135623730951)
```

### Algebra

Several methods for algebraic manipulation, such as `expand`, `factor`, `simplify`, `cancel`, and `apart`, are exported.

```@example usage
expand((x + y + 1)^3)
```

```@example usage
factor(x^4 - 1), factor(x^2 * y^2 - 1)
```

```@example usage
ex = (x^2 - 1)/(x-1)
simplify(ex)
```

```@example usage
cancel((x^2 - 1) / (x + 1)), apart(1 / (x^2 - 1), x)
```

#### Polynomials


For polynomial expressions, `gcd` works:

```@example usage
gcd(x^2 - 1, x^2 + 2x + 1)
```

The `Poly` constructor is used to create sparse, multivariate, symbolic polynomials. These have some performance gains over a simple symbolic expression and some additional methods.

```@example usage
q = Poly(x^2 * y + 3x - 1, x, y)
degree(q), degree(q, y), nterms(q)
```

### Equations and inequalities

An equation is written with `~`; `solve` (from CommonSolve) finds its
solutions. An inequality uses `≪` (typed `\ll<tab>`), `≦` (`\leqq<tab>`), `≫` (`\gg<tab>`), `≧` (`\geqq<tab>`), or `≶` (`\lessgtr<tab>`) (from CommonEq);
base comparison operators, e.g., `<`, are kept for comparing numbers and return a `Bool`.

```@example usage
solve(x^2 - 5x + 6 ~ 0, x), solve(x^3 - 3x + 1 ~ 0, x)
```

```@example usage
solve(x^2 ≫ 1, x)
```

```@example usage
solve([x + y ~ 3, x - y ~ 1], [x, y])
```

Assumptions may be specified to `solve` and many other methods; assumptions are not attached to a variable.

### Calculus

The basic operations `limit`, `diff`, and `integrate` of calculus are implemented in `SymJ` or a related package.

For derivatives, `diff` differentiates an expression with respect to a specified variable; `diff(f, x, n)` gives the `n`-th derivative.

```@example usage
diff(x^x, x), diff(sin(x)^2, x, 2)
```

The notation allows for mixed derivatives:

```@example usage
diff(x^2 * y, x, y)
```

Limits are computed by GruntzLimit (loaded by default) and
integrals by RischIntegrate, two companion packages that work on
SymJ expressions.

!!! note "Julia versions"
    SymJ supports Julia 1.10 and later. RischIntegrate needs Julia 1.12
    or later, so the integrals here need that version.

```@example usage
limit(sin(x) / x, x => 0), limit((1 + 1 / x)^x, x => Inf)
```

```@example usage
using RischIntegrate
integrate(x * sin(x), x), integrate(1 / (1 + x^2), (x, -1, 1))
```
