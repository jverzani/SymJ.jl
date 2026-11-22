# Scope: a `solveset`-like function beside `solve`

Status: a first version built (2026-11-22, `src/solveset.jl`); decisions below. In it:
finite sets, image sets (several bound variables, conditions on them), ConditionSet;
periodic equations through `solve`'s one period; exp, bᵘ, log, sin, cos, tan, cot, sec,
csc, sinh, cosh, tanh inverted with every branch, nested; `in`, `∩` with an interval
set, `collect`. Added the same day: parameters over ℝ, real, by cases (equal neighbouring
cases merged into one interval); Lambert's W with every branch over ℂ (checked at
k = -2…2 without parameters) and W₀, W₋₁ over ℝ; LaTeX through Latexify. Not yet:
membership in a family not linear in its variable (`lambertw(1, n)`), periodic equations
with parameters through more than one trigonometric function (a ConditionSet).

## Why

`solve` gives finite lists. For a periodic or multivalued inverse it picks one period or
the principal branch: `solve(sin(x) ~ 1//2, x)` is `[π/6, 5π/6]`, `exp(x) ~ 2` is
`[log(2)]`, `x^n ~ 1` is `[1]`, and `sin(sin(x)) ~ 0` keeps sin x = 0, π, 2π but not
-π, -2π (the roadmap's "all solutions"; SymPy #29675, #29435). A separate function
gives every solution as a set, so `solve`'s answers and their users stay as they are.

## The function

`SymJ.solveset(eq, x; domain = ℝ)` (`domain = ℂ` for complex solutions), not exported at
first. It returns a set:

| Equation | Result |
|:--|:--|
| `x^2 ~ 4` | `{-2, 2}` |
| `sin(x) ~ 1//2` | `{π/6 + 2πn : n ∈ ℤ} ∪ {5π/6 + 2πn : n ∈ ℤ}` |
| `tan(x) ~ 1` | `{π/4 + πn : n ∈ ℤ}` |
| `exp(x) ~ 2` | `{log(2) + 2πin : n ∈ ℤ}` (ℂ); `{log(2)}` (ℝ) |
| `x^(1//3) ~ 2` | `{8}` (principal), `x^3 ~ 8`: `{2, 2ω, 2ω²}` |
| `x*exp(x) ~ a` | `{W_k(a) : k ∈ ℤ}` (ℂ) |
| `x ~ sin(x)` | `ConditionSet(x, x - sin(x) = 0, ℝ)` (not solved: says so, no error) |
| `x^2 ≪ 1` | `(-1, 1)` (the existing `SymSet`) |

`domain = ℝ`, the default, keeps the real members (an image set's real members: those n
with a real value, e.g. of `log(2) + 2πin` only n = 0). With an interval domain, an image
set over ℤ is cut down to its finitely many members there. `solve` is left as it is.

## What it needs

1. **Set objects:** a finite set; an *image set* `{f(n) : n ∈ ℤ}` with `n` bound (so no
   free parameter is introduced, unlike `C1` in `solve`'s integer families); unions;
   `ConditionSet`; the empty set; the existing `SymSet` for real intervals. `in`, `∪`,
   `isempty`, display, LaTeX.
2. **Inverses with all branches:** sin, cos, tan, cot, sec, csc (and their period
   through the argument, `sin(2x + 1)`); sinh, cosh, tanh over ℂ; exp (`+ 2πin`); log;
   rational powers (the q roots); `lambertw` (branches k).
3. **Composition:** `sin(sin(x)) ~ 0` as `sin(x) ∈ {πn}`, each member inverted again,
   giving a union indexed by two integers (`{asin(πm) + 2πn} ∪ …`, real only for m = 0
   over ℝ).
4. **Intersection with an interval** for an image set linear in n: the n with
   `lo ≤ a + bn ≤ hi`, enumerated (a finite set).
5. **Reuse:** polynomial and rational equations, radicals, inequalities and parameter
   cases come from `solve`; `Cases` holds sets as it holds `solve`'s answers.

## Not in a first version

Intersections of two image sets (`{πn/2} ∩ {πn/3}`, a lattice problem), systems,
symbolic periods (`sin(a x)` with `a` a parameter needs `a ≠ 0`, generic, but its period
is symbolic), and simplifying unions whose members overlap (`{πn} ∪ {π/2 + πn}` as
`{πn/2}`).

## Decisions (2026-11-21)

1. The name is `solveset`.
2. The default domain is ℝ; `domain = ℂ` gives the complex solutions.
3. `solve` is left as it is (its one-period answers for sin and cos are not rebuilt on
   `solveset`).
