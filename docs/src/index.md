# SymJ

!!! warning "Work in progress"
    SymJ is under active development. Its interface may change, including in ways that
    break existing code, before a stable release.

SymJ is a symbolic manipulation kernel written in pure Julia. It builds expressions in a
canonical form as you go, keeps numbers exact, and is fast: construction, expansion,
differentiation and substitution are within a small factor of SymEngine, and usually
faster. Expressions are callable, so they can stand in for functions in numeric code, and
they implement [TermInterface](https://github.com/JuliaSymbolics/TermInterface.jl), so
generic term-rewriting code, such as
[GruntzLimit](https://github.com/jverzani/GruntzLimit.jl) for limits and
[RischIntegrate](https://github.com/jverzani/RischIntegrate.jl) for integrals, works on
them.

- The [Tutorial](@ref) works through the symbolic side of calculus with SymJ.
- [Implementation](@ref) discusses the design choices: one node type, exact numbers,
  packed polynomials, and the gcd and factoring algorithms.
- The [Design notes](@ref) record each decision in detail, with measurements.
- The [Reference](@ref) lists the functions.

## General usage

```@example usage
using SymJ
```

### Symbolic variables and numbers

Variables are made with `@syms`, or with `Sym` from a name. Numbers stay exact: `Sym`
makes an exact number, and arithmetic with symbols keeps integers and rationals exact.

```@example usage
@syms x y z
p = (x + 1)^2 / 3
```

```@example usage
sqrt(Sym(8)), sin(Sym(π) / 3), Sym(2)^100
```

A floating point number in an operation makes the result a float, so use `1//2`, not
`1/2`, for an exact half:

```@example usage
x / 2 == x * (1 // 2), x * (1 / 2)
```

### Expressions are built in canonical form

Sums collect like terms and products collect like bases as they are built, with no
separate simplification step:

```@example usage
x + 2x - y + x * y * x, exp(x) * exp(y), (x * y)^2
```

### Calling an expression

An expression is a `Function`. Called with a number, it evaluates with Julia's own
arithmetic, so it works in numeric code (with dual numbers, `BigFloat`, intervals, …):

```@example usage
u = x^2 - 2x
u(3), u(3.0), u(1 // 2)
```

Called with pairs, it substitutes; `subs` does the same:

```@example usage
ex = x^2 * y + sin(x)
ex(x => 2), subs(ex, x => y + 1)
```

A constant expression called with no arguments gives its value, and `N` gives a float:

```@example usage
(Sym(2)^10)(), N(sqrt(Sym(2)))
```

For speed in a loop, `lambdify` compiles an expression into a Julia function:

```@example usage
f = lambdify(sin(x)^2 + 3x^5 - exp(-x / 2), x)
f(0.5)
```

### Algebra

```@example usage
expand((x + y + 1)^3)
```

```@example usage
factor(x^4 - 1), factor(x^2 * y^2 - 1)
```

```@example usage
cancel((x^2 - 1) / (x + 1)), apart(1 / (x^2 - 1), x)
```

`gcd` works on polynomial expressions, and `Poly` gives the polynomial itself:

```@example usage
gcd(x^2 - 1, x^2 + 2x + 1)
```

```@example usage
q = Poly(x^2 * y + 3x - 1, x, y)
degree(q), degree(q, y), nterms(q)
```

### Equations and inequalities

An equation is written with `~`; `solve` (from CommonSolve) finds its solutions. An
inequality uses `≪` (typed `\ll`), `≦`, `≫`, `≧` or `≶`, since `<` is kept for comparing
numbers and returns a `Bool`.

```@example usage
solve(x^2 - 5x + 6 ~ 0, x), solve(x^3 - 3x + 1 ~ 0, x)
```

```@example usage
solve(x^2 ≫ 1, x)
```

```@example usage
solve([x + y ~ 3, x - y ~ 1], [x, y])
```

### Calculus

`diff` differentiates; `diff(f, x, n)` gives the `n`-th derivative.

```@example usage
diff(x^x, x), diff(sin(x)^2, x, 2)
```

Limits and integrals come from two companion packages, which work on SymJ expressions:

!!! note "Julia versions"
    SymJ and GruntzLimit support Julia 1.10 and later. RischIntegrate needs Julia 1.12 or
    later, so the integrals here need that version.

```@example usage
using GruntzLimit, RischIntegrate
gruntz_limit(sin(x) / x, x, 0), gruntz_limit((1 + 1 / x)^x, x, Inf)
```

```@example usage
integrate(x * sin(x), x), integrate(1 / (1 + x^2), x)
```

### Parameters

A parameter is a named constant, as opposed to a variable. Calling an expression
substitutes for its variable and keeps the parameters:

```@example usage
@parameters a b
g = a * x^2 + b * x
g(2), solve(g ~ 0, x)
```
