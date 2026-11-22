# Solving

Solving equations: roots of polynomials, equations and inequalities in one unknown, and systems.

```@setup tut
using SymJ
using RischIntegrate
@syms x y z a b c t f()
```

SymJ provides solvers for different types of equations, the most general being `solve`, described here along with a few other methods. In the Calculus section we also introduce `rsolve` and `dsolve`.

Equations are specified with `~` (as is common in Julia) or `Eq`. This
infix operator binds more loosely than everything else, so no
parentheses are typically needed. Equations are held internally as a
`Relation`. The left- and right-hand sides of a relation are returned
by `lhs` and `rhs` respectively.

## Roots of polynomials

We will see more uses later, but `solve` finds the roots of a
polynomial. Degrees up to four are solved exactly, with real roots of
cubics given in trigonometric form (no complex numbers needed):

```@example tut
p = expand((x^2 + x)^3 - 2)
solve(p ~ 0, x)
```

```@example tut
solve(x^3 - 3x + 1 ~ 0, x)
```

```@example tut
solve(x^4 - 10x^2 + 1 ~ 0, x)
```

(Neither of the last two examples can be factored with `factor`, which provides factorization into irreducibles over the rationals.)

From degree five, it is possible there are no roots in
radicals. Those of special forms, such as `x^5 ~ 2`, are still solved
exactly by `solve`.

When exact roots aren't identifiable, `SymJ.fsolve` finds the
roots numerically, the real ones by default (with `domain = ℂ` all of them):

```@example tut
(
SymJ.fsolve(x^5 - x + 1 ~ 0, x),
SymJ.fsolve(x^5 - x + 1 ~ 0, x; domain = ℂ)
)
```

The real roots of a polynomial with rational coefficients can be found
numerically, with certainty, when the closed form is
unwieldy. `SymJ.real_roots` gives them correctly rounded, each
repeated by its multiplicity, and `SymJ.root_intervals` exact
isolating intervals:

```@example tut
(
SymJ.real_roots(x^3 - 2x - 5, x),
SymJ.real_roots((x - 1)^2 * (x^2 - 2), x)
)
```

```@example tut
SymJ.root_intervals(x^2 - 2, x; width = 1 // 100)
```

## Solving equations

Non-polynomial equations can also be solved.

```@example tut
solve(a*sin(x) + b ~ 0, x)
```

Answers to `solve(eqn, var)` are a vector of `Sym` values, empty when there are no identified solutions, or `Cases` when the answer depends on the values of the parameters.

Equations in which `x` appears through one function are solved by solving for that function first:

```@example tut
solve(exp(2x) - 3exp(x) + 2 ~ 0, x)
```

```@example tut
solve(sqrt(x + 1) ~ x - 1, x)
```

Several logarithms of the unknown are combined; roots where a logarithm would be complex (here `x = -2`) are dropped:

```@example tut
solve(log(x) + log(x - 1) ~ log(Sym(6)), x)
```

A product is zero where one of its factors is, so a product (or a sum
with a common factor) is solved factor by factor. A factor that is
never zero, such as `(a/x)^x`, gives no roots, and a root where
another factor is undefined is dropped (`x*log(x)` none at `0`).

```@example tut
solve(x * log(x) ~ 0, x)
```

Some equations need the Lambert W function, from the
[LambertW](https://github.com/jlapeyre/LambertW.jl) package:

```@example tut
solve(x * exp(x) ~ 1, x)
```

### Generic values

Symbols other than the unknown are allowed.

Such symbols, or parameters, take *generic values*: an expression in
them alone is taken to be nonzero, so the quadratic formula comes out
as it is usually written:

```@example tut
solve(a * x^2 + b * x + c ~ 0, x)
```

An assumption that makes a coefficient vanish gives that case:

```@example tut
solve(a * x^2 + b * x + c ~ 0, x; assume = (a ~ 0))
```

A symbolic exponent gives the principal root only. `x^n = 1` has `n` roots for a
positive integer `n`, but their number depends on `n`, so with `n` symbolic the answer is the principal one:

```@example tut
@syms n
solve(x^n ~ 1, x; assume = (n in Integers, n ≫ 0))
```

With a number for `n`, all of them are identified:

```@example tut
solve(x^3 ~ 1, x)
```

### Inequalities

`<` and its relatives compare numbers and return `true` or `false`, so inequalities use other operators: `≪` (typed `\ll`), `≦` (`\leqq`), `≫` (`\gg`), `≧` (`\geqq`) and `≶` (`\lessgtr`).

The solution to an inequality is a union of intervals:

```@example tut
solve(x^2 - 1 ≫ 0, x)
```

```@example tut
solve((x - 1) / (x + 2) ≧ 0, x)
```

With parameters, which an order makes real, the set comes by cases on them:

```@example tut
solve(x^2 ≦ a, x)
```


### More on assumptions

Symbols are complex numbers unless a fact says otherwise. Facts are relations, such as `x > 0`, or membership in one of `Integers ⊂ Rationals ⊂ Reals ⊂ Complexes` (also written `ℤ ⊂ ℚ ⊂ ℝ ⊂ ℂ`), such as `n in Integers`. Facts are passed to the `assume` keyword, and only the solutions they allow are kept.

The `@assume` macro lets facts be written with Julia's familiar Boolean operators, including `&&`, `||`, `!` and chained comparisons.

```@example tut
solve(x^2 ~ 4, x; assume = @assume(x > 0))
```

```@example tut
solve(x^2 + 1 ~ 0, x; assume = @assume(x in Reals))
```

Without the macro, relations use the symbolic operators: `x ≫ 0` (or
`Gt(x, 0)`), and `n ∈ ℤ` for membership.

```@example tut
solve(x^2 ~ 2, x; assume = (x ∈ ℚ,))
```

`@assuming` gives facts for a block of code (as does the function
`assuming`, with a `do` block):

```@example tut
@assuming x > 0 begin
    solve(x^3 ~ 8, x), solve(x^2 - 5x + 6 ~ 0, x; assume = @assume(x > 2))
end
```

The facts are read by the operations that take `assume` (`refine`,
`simplify`, `solve`, `limit`, …), not by arithmetic: an expression is
built the same way inside a block as outside it, so it never depends
on the facts in force. `(a^4)^(3/4)` is not `a^3` for every `a` (it is
`1` at `a = -1` and at `a = im`), so `16(a^4)^(3/4)/(9a^2)` keeps
`(a^4)^(3/4)` even with `a > 0` assumed; `simplify` or `refine` in the
block applies the facts, and `refine = true` refines the block's value
with them:

```@example tut
@assuming a > 0 begin
    e = 16(a^4)^(3//4) / (9a^2)
    e, simplify(e)
end
```

```@example tut
@assuming refine = true a > 0 begin
    16(a^4)^(3//4) / (9a^2)
end
```

#### Cases

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

Contradictory facts are an error (`@assume(x > 1 && x < 0)`), and the
`Assumptions` value `@assume` makes can be reused: `A = @assume(x > 0,
n in Integers)`, then `assume = A`.


As with polynomials, for an equation `solve` cannot solve exactly, `SymJ.fsolve` (with name
and functionality mirroring Maple's function) finds its roots numerically: the
real ones, on the whole line or an interval, and with `domain = ℂ`
complex ones too.

```@example tut
(
SymJ.fsolve(cos(x) ~ x, x),
SymJ.fsolve(log(x) + sin(x) ~ 1, x, (0, 4))
)
```

### External zero finding

The packages for numerical root finding also take SymJ expressions as
they are, since an expression is callable (`lambdify` makes a faster
function): `Roots.find_zero(cos(x) - x, 0.5)`, or with Newton's method
and the exact derivative, `find_zero((lambdify(f, x), lambdify(diff(f,
x), x)), x0, Roots.Newton())`.

### Systems

A system of equations gives a vector of solutions, each a `Dict` from the unknowns to their values. A linear system's answer is a single solution (possibly in terms of free unknowns), yet a vector is returned:

```@example tut
solve([x + y ~ 3, x - y ~ 1], [x, y])
```

With parameters, a generic value is assumed—the pivots are taken to be nonzero (here `a + 1`; with `a ~ -1` assumed the system has no solution):

```@example tut
eqns, vars = [a * x + y ~ 1, x - y ~ b], [x, y]
solve(eqns, vars)
```

Polynomial systems are solved through a lex Gröbner basis, which is triangular:

```@example tut
solve([x^2 + y^2 ~ 1, x ~ y], [x, y])
```

```@example tut
solve([x*y ~ 1, x + y ~ 3], [x, y])
```

Unknowns that appear only through exponentials, through one function, or through `sin` and `cos` of one argument are solved for those first:

```@example tut
solve([x + y ~ 4, 2^x + 2^y ~ 10], [x, y])
```

```@example tut
solve([sin(x) + cos(y) ~ 1, sin(x) - cos(y) ~ 0], [x, y])
```

With infinitely many solutions, the unknowns left free map to themselves:

```@example tut
solve([x^2 + y^2 ~ 1], [x, y])
```

A system of linear inequalities gives a description of the region, each unknown between bounds in the ones after it:

```@example tut
solve([x ≧ 0, y ≧ 0, x + y ≦ 1], [x, y])
```
