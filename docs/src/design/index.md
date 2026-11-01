# Design notes

SymJ is a symbolic manipulation kernel written in pure Julia. It aims for SymEngine-class
speed with a small, explicit core that later layers (simplification, solving, matrices)
build on.

The notes are in three pages: this one (the core: the node type, numbers, functions,
evaluation, relations, pattern matching, open questions), [Algebra](@ref design-algebra)
(polynomials, solving, simplification, linear algebra, assumptions) and
[Calculus](@ref design-calculus) (piecewise functions, summation, differential equations).

## Goals

- **Speed comparable to SymEngine.** Construction, canonicalization, expansion,
  differentiation and substitution should be within a small factor of SymEngine, and
  ideally faster.
- **Callable objects.** Symbolic objects subtype `Function` (via `AbstractSymbolic`).
  Calling is meant to make expressions drop-in replacements for functions in numeric
  code:
  - `ex(x => 1, y => 2)` substitutes and returns a `Sym`.
  - `ex(a)` with a symbolic `a` substitutes into the expression's single free symbol.
  - `ex(a)` with any other `a` evaluates with Julia arithmetic on `a`, returning a Julia
    value: if `u = x^2` then `u(3) === 9`, `u(3.0) === 9.0` and `u(1//2) === 1//4`.
    Because the arithmetic is generic, dual numbers, intervals and `BigFloat` work.
    Constants such as `sqrt(2)` and `π` are computed in the argument's floating point
    type.
  - `ex()` gives the value of a constant expression as a Julia number (`Sym(2)() === 2`,
    `Sym(π)() === π`), and returns `ex` unchanged if it has a free symbol. It
    deliberately does not throw: a Julia `throw` costs about 40 μs (the backtrace),
    and generic code reading leaves pays that on every symbol. This matches what generic
    term code expects from a callable numeric leaf; GruntzLimit, for example, reads
    numbers this way.
  - An expression with several free symbols must be called with pairs or keywords. A
    single value would leave implicit which symbol receives it.
  - `ex(; x = 2, y = 3)` matches keywords to variables and parameters by name. When all
    of them get non-symbolic values it evaluates with Julia arithmetic, returning a Julia
    number; otherwise it substitutes. Names not in `ex` are ignored (one set of keywords
    can serve several expressions), and keywords aren't mixed with pairs or a positional
    value. A name shared by a variable and a parameter is an error.
  - A `SymFunction` called on arguments builds an unevaluated application (`f(x)`).

  Evaluation interprets the tree: about 25 ns for `x^2` and 2.4 μs for
  `sin(x)^2 + 3x^5 - 2x + exp(-x/2)`. For speed, use `lambdify` (below).
- **Scalar-like collections.** Like `Number`, a `Sym` or `Poly` iterates once (yielding
  itself), has `length` 1 and `size` `()`, and supports `x[]` and `x[1]`. Generic code
  that iterates or indexes a scalar therefore works unchanged.
- **Exact numbers.** `sqrt(Sym(8))` is `2*sqrt(2)`, `sin(Sym(π)/3)` is `sqrt(3)/2`, and
  `Sym(2)^100` is exact.
- **Canonical form as you go.** Sums collect like terms and products collect like bases
  at construction time. No separate `simplify` pass is needed for those identities.
- **TermInterface v2.** SymJ implements the contract so that generic rewriting tools
  work on its expressions.
- **Commutative scalars.** `+` and `*` are commutative. Noncommutative objects such as
  symbolic matrices are future work and will need their own node kinds.
- **Independent of SymbolicUtils.** The representation choices below are SymJ's own.

## One concrete node type

```julia
abstract type AbstractSymbolic <: Function end

mutable struct Sym <: AbstractSymbolic
    const kind::Kind              # NUMBER SYMBOL CONSTANT ADD MUL POW FN
    const hash::UInt              # structural hash, computed once
    const name::Symbol            # symbol / constant / function name
    const num::SymNumber          # value, constant term, or coefficient
    const args::Vector{Sym}       # terms / bases / [base, exp] / fn args
    const exps::Vector{Sym}       # MUL: exponent of each base
    const coefs::Vector{SymNumber}# ADD: coefficient of each term
end
```

This follows the same idea as SymEngine's `TypeID`: one concrete type with a tag, rather
than a Julia type per node kind. As a result:

- Child containers are concretely typed (`Vector{Sym}`), so traversals never box and the
  compiler can infer everything.
- Hot paths such as `add`, `mul` and `pow` branch on a byte instead of dispatching
  dynamically.
- `===` is pointer identity (the struct is `mutable` with all-`const` fields). That gives
  an O(1) fast path for equality and avoids copying a wide struct when vectors are
  stored or passed.
- Every node stores its structural hash. `==` compares hashes first, and `Dict{Sym,…}`
  never rehashes a subtree.

Empty child vectors are shared constants (`NOSYMS`, `NONUMS`), so leaves do not allocate
them. Integers in −64:64 are cached nodes.

**Nodes are immutable by convention.** Child vectors are shared between nodes, for
example when an `ADD` changes one coefficient and reuses its `args`. Nothing may mutate
them after construction.

### Representation by kind

| kind | meaning | `num` | children |
|:--|:--|:--|:--|
| `NUMBER` | exact or float scalar | value | — |
| `SYMBOL` | variable `x` | — | — |
| `CONSTANT` | `π`, `ℯ` | — | — |
| `ADD` | `c + Σ cᵢ tᵢ` | `c` | `args = tᵢ`, `coefs = cᵢ` |
| `MUL` | `c · Π bᵢ^eᵢ` | `c` | `args = bᵢ`, `exps = eᵢ` (symbolic) |
| `POW` | `b^e` | — | `args = [b, e]` |
| `FN` | `f(a₁,…,aₙ)` unevaluated | — | `args`, `name = f` |

This is SymEngine's dictionary representation (`Add`: term → coefficient, `Mul`:
base → exponent) laid out as sorted parallel vectors instead of hash maps:

- **Canonical storage order** is `_cmp`. It compares hashes first and falls back to
  structure on a collision. The order is total and deterministic across sessions,
  because symbol hashes come from the name string. It has no mathematical meaning.
- Combining two sums or two products is a **linear merge** of sorted vectors. Adding a
  single term is a binary search plus one copy. Sums of many terms (`sum`, `expand`,
  `diff`, `subs`) collect into a `Dict` and sort once.
- Equality is elementwise comparison of the vectors.

Display order (higher degree first, then lexical) is computed only by the printer and by
`sorted_arguments`.

### Canonical-form invariants

- **ADD.** No numeric or `ADD` terms. Terms carry no coefficient (the coefficient lives
  in `coefs`). Coefficients are nonzero. The node has at least two terms, or one term
  plus a nonzero constant.
- **MUL.** No `MUL` bases. Bases are distinct and exponents are nonzero. A `NUMBER` base
  appears only with a non-integer exponent, as in `2^(1/2)`. The node is never
  coefficient 1 times a single factor (that is a `POW` or the base itself). It is never
  a number times a single sum, because numbers distribute: `2(x+y) → 2x+2y`. This keeps
  `a - a == 0` true for sums.
- **POW.** The exponent is never 0 or 1. Integer powers distribute over products
  (`(xy)² → x²y²`) and multiply through powers (`(x^a)^n → x^(an)`). Non-integer outer
  powers are left alone, because `(x²)^(1/2) ≠ x` in general.
- **`exp(x)` is `ℯ^x`** (as in SymEngine), so `exp(x)·exp(y) → exp(x+y)` comes free from
  the product rules. `exp(log(z)) → z`.

The canonicalizing constructors are `add`, `mul` and `pow`. The `_mk*_node` functions
build raw nodes. The `_mk*` helpers (`_mkadd`, `_mkmul`, `_mkpow`) assemble nodes from
parts that are already canonical and handle the degenerate cases.

## Numbers: `SymNumber`

```julia
struct SymNumber <: Real
    tag::UInt8      # SMALL (Int n//d), BIG (Rational{BigInt}), FLOAT (Float64)
    n::Int; d::Int; f::Float64
    big::Union{Nothing,Rational{BigInt}}
end
```

- **Exact by default.** Arithmetic on small values uses `Int` with overflow checks
  (`add_with_overflow`, `mul_with_overflow`, gcd cross-cancellation). It promotes to
  `BigInt` on overflow and demotes when a result fits again. Most coefficient arithmetic
  never allocates.
- **Floats are contagious.** A float in an operation makes the result a float.
- **Equality.** `==` is numeric (`1 == 1.0`). `isequal` and `hash` are structural, so the
  exact `1` and the float `1.0` are different expression nodes.
- **Exact radicals.** A number to a rational power is canonicalized as `n^(p/q) =
  c·m^(r/q)` with `m` free of `q`-th powers (trial division up to 1000, then a
  perfect-power check). For example `sqrt(8) → 2·sqrt(2)`, `4^(1/4) → sqrt(2)`,
  `2^(-1/2) → sqrt(2)/2`.
  - Negative bases split off `(-1)^e`, with `e` reduced into (−1, 1]. So `sqrt(-1)` is
    `(-1)^(1/2)`, which prints as `im` and squares to −1.
  - Rational bases split into numerator and denominator.
- **Symbolic numbers.** `Sym(2)`, `Sym(1//3)`, `Sym(π)` are the user-facing exact
  numbers, and `sqrt(Sym(2))` stays exact. `N(ex)` evaluates to `Float64` or
  `ComplexF64`.

## Functions

Unevaluated applications are `FN` nodes identified by name. Behaviour lives in a
registry:

```julia
register_function!(:sin; julia = sin,                 # numeric evaluation, TermInterface.operation
                         eval  = args -> ...,         # simplification, or `nothing`
                         diff  = (args, i) -> ...)    # ∂/∂argᵢ
```

- **Built-in functions** (`sin`, `cos`, `tan`, inverse and hyperbolic trig, `abs`,
  `log`) extend `Base` on `Sym`. Their rules handle exact special values (`sin(kπ/6)`,
  `log(1)`, …), odd and even symmetry (via `_could_extract_minus`, so `sin(-x) →
  -sin(x)`), and inverse compositions.
- **Float arguments** evaluate numerically. Exact arguments stay symbolic: `sin(2)` is
  unevaluated.
- **Coverage of Base's mathematical functions.**
  - *Unevaluated functions with rules:*
    - the reciprocal trig and hyperbolic functions (`sec csc cot sech csch coth`) and
      their inverses;
    - `cbrt` (the real cube root), `sinc`/`cosc`;
    - `sign floor ceil trunc round`, `mod rem`, `max min` (n-ary, flattened,
      canonically ordered), `factorial`, and `binomial` (a polynomial for a literal `k`).
    - `atan(y, x)` (internally `atan2`, so that rules for `atan`'s one argument never
      meet it): exact for constants of known sign, placed in its quadrant by `refine`.
    - Rounding takes Base's rounding modes, `digits` and `sigdigits` (an exact result),
      and integer types (`floor(Int, x)`).
    - Exact values come from the function they are defined through, accepted only when
      that function disappears: `sec(π/3) = 1/cos(π/3) = 2` and `sec(acos(x)) = 1/x`,
      while poles such as `sec(π/2)` stay unevaluated.
    - Piecewise constant functions evaluate constant arguments numerically only away
      from their jumps. Their derivatives are left unevaluated.
  - *Rewrites into existing forms:* the degree variants (`sind`, `asind`, …), `sinpi`,
    `cospi`, `tanpi`, `log2`, `log10`, `log1p`, `exp2`, `exp10`, `expm1`, `hypot`,
    `abs2`, `fourthroot`, `cis`, `cispi`, `deg2rad`, `rad2deg`, `div`/`fld`/`cld`,
    `mod2pi`, `clamp`, `minmax`, `muladd`/`fma` (so `evalpoly` works), `\`, and
    `numerator`/`denominator` (via `numer_denom`).
  - *Complex parts:* `real`, `imag`, `conj`, `angle`, `reim`, and `adjoint` and
    `transpose` for arrays (see Complex numbers below).
  - *Predicates and conversions:* `isreal`, `isfinite`, `isinf`, `isnan` are `true`
    only when known (from the expression or an enclosing `@assuming` block);
    `Int`, `Integer`, `Rational`, `BigFloat` and the other number types convert
    values that fit them, else an `InexactError`.
  - *Not yet:* the bit and integer-representation functions.
- **Special-function packages.** SymJ depends on **SpecialFunctions**,
  **HypergeometricFunctions** (`_₁F₁`, `_₂F₁`, `_₃F₂`, and `pFq` of those sizes) and
  **PolyLog** (`li0`–`li6`, `li(n, z)`, `reli…`); import a package to name its functions.
  - Calling a package function with a `Sym` argument, in any position, builds the
    unevaluated function: `sinint(x)`, `besselj(n, x)`.
  - The package's own function is registered as the numeric method. It runs only when
    asked for (`N`, calling with numbers, `lambdify`) or when an argument is a float.
  - Exact special values are given where standard: `gamma(7/2) = 15√π/8`,
    `zeta(2) = π²/6`, `li2(1) = π²/6`.
  - Derivative rules are given where standard: Bessel recurrences,
    `d/dz ₂F₁ = ab/c ₂F₁(a+1, b+1, c+1, z)`, `d/dz Li_n = Li_{n-1}/z`.
  - A partial derivative with no closed form, such as with respect to an order or a
    parameter, is left as an unevaluated `Derivative`. A derivative rule may return
    `nothing` for this.
  - The symbolic methods are defined at the extension's top level, so they're
    precompiled. The rules are registered in its `__init__`, because SymJ's rule table
    can't be modified during the extension's precompilation.
- **Integer arguments stay integers.** An exact integer passed directly to a function
  stays an `Int` when evaluated, as in Julia (`besselj(2, x)`). Some parameters must be
  integers: `polygamma`'s order, `besselh`'s kind.
- **`N` uses real methods for real arguments.** Many special functions have no complex
  method, so `N` evaluates with real numbers when all arguments are real. It falls back
  to complex only when the real call raises a `DomainError`.
- **`log(b, x)`** is stored unevaluated as a two-argument `FN`, with `log(b, b) = 1`,
  `log(b, 1) = 0` and derivatives in both arguments.
- **User functions** (`SymFunction(:f)` or `@syms f()`) have no rules. Their derivatives
  are `Derivative(f(x), x)`.

The registry is a `Dict` lookup per application, and the rule closures are called
dynamically. Function application is not on the innermost hot path, but if it becomes a
bottleneck the fix is a function-id field and a rule table indexed by integer.

## TermInterface

- **Calls.** `ADD`, `MUL`, `POW` and `FN` satisfy `isexpr` and `iscall`. `head` is
  `operation`, and `children` is `arguments`.
- **`operation`.** Returns `+`, `*`, `^`, the registered Julia function (`sin`, …), or a
  `SymFunction`.
- **`arguments`.** For `ADD` and `MUL` it is the list of summands or factors, with the
  numeric term or coefficient included. For `POW` and `FN` it is the node's own vector,
  which callers must treat as read-only.
- **`sorted_arguments`.** Returns the arguments in display order.
- **`maketerm(Sym, op, args, metadata)`.** Rebuilds through the canonicalizing
  constructors. `metadata` is ignored, and nodes deliberately carry no metadata field.

## Compiled evaluation: `lambdify`

`lambdify(ex, x, y)` compiles `ex` into a Julia function of `x` and `y`, in that order.
Without variables, `ex` may have at most one free symbol. The result is a
RuntimeGeneratedFunction, so it can be called right away, from the same function that
created it, with no world-age problem.

- **Same semantics as calling.**
  - Exact constants are Julia integers and rationals, converted to the arguments'
    float type when there is one: `F = _liftfloat(promote_type(...))` is a compile-time
    constant, so for `Float64` the constants fold and `BigFloat` keeps its precision.
  - Integer powers are literal powers (`x^5`), and negative ones use `inv`.
  - `ℯ^u` is emitted as `exp(u)`.
  - A negative base with a fractional power gives the principal complex value.
- **Shared subexpressions.** Nodes used more than once become local variables. Nodes
  carry hashes, so counting uses is one pass over the DAG. Arrays and tuples of
  expressions compile to one function returning the same shape, sharing work between
  entries, which suits gradients and Jacobians.
- **Speed.** `sin(x)^2 + 3x^5 - 2x + exp(-x/2) + sqrt(2)x/3` takes 9.3 ns per call,
  against 9.2 ns hand-written, and is inferred as `Float64`. Evaluating it at 1,000
  points takes 7.49 μs against 7.48 μs.
- **When it pays.** `bench/newton.jl` runs Newton's method on `x^5 - x - 1` from 1,
  which takes 7 iterations and 14 evaluations. Calling the expressions takes 3.6 μs per
  solve, the `lambdify`'d functions 86 ns (42x faster; hand-written: 89 ns).
  Compiling costs about 2.9 ms per new expression, the first call compiling the
  generated method, so it pays off after roughly 800 solves or 11,600 evaluations. The
  compiled code is cached by expression, so compiling the same expression again is
  free.
- **Generic arithmetic.** Arguments may be any number type: `ForwardDiff` derivatives
  and gradients work through both `lambdify` and plain calls (tested).
- `lambdify_expr` returns the generated code.

## Relations

`==`, `<` and `≤` keep returning `Bool`. Structural equality is what `Dict`s, `in` and
generic code need, and a symbolic result in `if x < 0` fails confusingly. Equations and
inequalities use other operators:

- **Equations:** `lhs ~ rhs`, the convention of Symbolics.jl, ModelingToolkit and
  SimpleExpressions. Binary `~` is free in Base, and its low precedence means
  `x^2 + 1 ~ 3` needs no parentheses.
- **Everything else** uses CommonEq.jl's operators and names, defined in SymJ without
  depending on that package:

  | relation | infix | typed as | function |
  |:--|:--|:--|:--|
  | `=` | `~`, `⩵` | `\Equal` | `Eq` |
  | `<` | `≪` | `\ll` | `Lt` |
  | `≤` | `≦` | `\leqq` | `Le` |
  | `>` | `≫` | `\gg` | `Gt` |
  | `≥` | `≧` | `\geqq` | `Ge` |
  | `≠` | `≶`, `≷` | `\lessgtr`, `\gtrless` | `Ne` |

- **Representation.** A `Relation` (an operation plus `lhs` and `rhs`) is its own small
  type, not a `Sym` node, because relations aren't algebraic quantities.
- **TermInterface.** `operation` is `~`, `Lt`, …
- **Calling.** It substitutes into both sides, like calling an expression. When no
  symbols remain, the result is a `Bool`: exact numbers compare exactly, other values
  numerically. A relation whose two sides are identical is decided outright.
- **Precedence caveat.** The infix inequalities have comparison precedence, so Julia
  chains them like `<`: `a ≪ b == c` means `(a ≪ b) && (b == c)`. Parenthesize:
  `(a ≪ b) == c`. `~` has the lowest precedence and doesn't chain.

## Parameters (internal)

!!! note "Internal, not exported"
    Parameters are hidden for now: `@parameters`, `parameter`, `isparameter` and
    `free_parameters` are not exported or documented for users (they are reachable as
    `SymJ.parameter` and so on), so that the current approach does not become
    permanent before the shared-name question below is settled.

`parameter(:a)` or `@parameters a b` makes a symbolic parameter: a named constant of
unknown value, as opposed to a variable. It is a node kind, `PARAMETER`, beside
`SYMBOL`, not a new Julia type. Expressions stay one concrete type, so containers
stay concretely typed and dispatch can't distinguish expressions that contain
parameters; `isparameter` and `free_parameters` answer that instead. The parameter
`a` and the variable `a` hash differently and are distinct.

**Shared names (undecided).** Because the two are distinct nodes, a parameter and a
variable can have the same name. Their sum prints as `a + a` and is not `2a`, which is
confusing. Calling with keywords raises an error for such a name; nothing else checks.
Prohibiting it needs a session-wide record of parameter names, checked whenever either
kind is created: `@syms a` after `@parameters a` (or the reverse) would be an error. The
costs are that a name can't change role within a session, and that everything run in
one process (the test suite, the documentation build, other packages creating
variables by name) shares that record. The decision is deferred until it is clearer how
much parameters are used. The lighter alternative is to keep them separate and print
parameters distinctly. Either way the current behavior will almost certainly change,
and the user-facing documentation says so.

| | variables | parameters |
|:--|:--|:--|
| `free_symbols` / `free_parameters` | listed by `free_symbols` | listed by `free_parameters` |
| `ex(v)` | substituted, or evaluated if no parameters | kept: `(a*x^2)(3) == 9a` |
| `lambdify(ex)` | first argument | after it, by name: `(x, a, b) -> …` |
| printing order | by degree | degree 0, so `a*x^2 + b*x + c` |
| `diff`, `subs`, `expand`, `Poly`, `cancel`, `factor`, `gcd` | as symbols | as symbols |
| constant tests (`sign`, `floor`, `mod`, `max`) | not constant | not constant |

**Cost.** `bench/parameters.jl` times 22 core workloads, recording the minimum time of
each. One run of identical code differs from another by up to ±13% per workload on
this machine, and a baseline recorded earlier drifts with load. So the comparison put
the previous commit in a git worktree and ran the old and new versions alternately,
three rounds each, taking each workload's best time.

- The geometric mean ratio over the suite was 0.995.
- The scan for parameters before a numeric call costs about 4% (0.2 μs) on a 40-node
  expression.
- `free_symbols` became 1.8x faster, because loops replaced the closure-based recursion.

## Performance snapshot

Run with `julia --project=bench bench/bench.jl` (Apple Silicon, Julia 1.13, SymEngine.jl
0.13). Times are medians. The ratio is SymJ time / SymEngine time, so lower is better.

| benchmark | SymJ | SymEngine | ratio |
|:--|--:|--:|--:|
| `x + y` | 32 ns | 98 ns | 0.33 |
| build Σ i·xⁱ·y²⁰⁰⁻ⁱ term by term | 216 μs | 593 μs | 0.36 |
| build Π (x+i)ⁱ, i ≤ 30 | 5.1 μs | 14.7 μs | 0.35 |
| `expand((x+y+z+1)^15)` | 135 μs | 534 μs | 0.25 |
| `expand((a+x)^100)` | 26 μs | 29 μs | 0.91 |
| `diff` of expanded poly · sin | 38 μs | 52 μs | 0.72 |
| `subs` into 286-term poly | 152 μs | 170 μs | 0.89 |
| `expand(e*(e+w))`, `e=(x+y+z+w)^15` | 5.1 ms | 415 ms | 0.01 |

The `expand(e*(e+w))` row takes 5.1 ms in total. Of that, the `Poly` product (816 × 817
terms, giving 6 272) takes 1.8 ms, `Poly → Sym` takes 3.2 ms, and `Sym → Poly` takes
1.5 ms. SymEngine.jl's `expand` works on the generic tree, so this compares like with
like at the user level. SymEngine's own sparse polynomial classes are not exposed by
SymEngine.jl. SymEngine.jl's numbers also include the cost of the C wrapper (a `ccall`
and a finalizer per operation), which dominates the small benchmarks.

## Pattern matching

The design, and the decisions behind it, are in `notes/pattern_matching.md`; stages 1–4,
the matcher, rules, their speed and their first use in `simplify`, are done.

- **Patterns are expressions.** `SymJ.pattern(:(…))` reads a quoted expression with
  SymJ's arithmetic, a slot `~x` being a private function node `_slot(x, flavor)` and a
  head slot `(~f)(…)` a node `_hslot(f, …)`. So a pattern is in canonical form and means
  what that form means: `~x + ~x` is `2*~x`, `exp(log(~x))` is `~x`. A slot that
  vanishes (`~x - ~x`) is an error.
- **Sums and products are associative and commutative.** Their elements are the terms
  (the constant last) or the factors (the coefficient last). A pattern's parts without
  slots are found by hash; its other parts each match one element, tried only against
  elements with the same head; its slots take groups, `~x` one or more and `~!x`, `~~x`
  any number, bound to the group's sum or product. A group is a subsequence of the
  subject's terms, already in canonical order, so its sum is built directly.
- **Numbers divide into slots.** A product's coefficient over the pattern's is one more
  element (`2*~x` against `6a` gives `x => 3a`), and a scaled slot in a sum takes a group
  divided by its scale (`~x - ~y` against `a - 2b` gives `y => 2b`).
- **The first match.** Slots already bound take their value's elements, then defaults
  (offered the number first, else taking none), scaled slots (negative terms first),
  slots with predicates and plain slots, each in the order written; the last takes what
  is left. Ambiguous splits (`~x + ~y` against `a + b`) follow the stored order of the
  terms, not the printed one, which would cost a string comparison.
- **The search** is backtracking with the continuation as a linked list of concrete
  goal records, not closures, so it is fully typed.

- **Rules.** `@rule lhs => rhs` reads `lhs` as a pattern, with its predicates (and any
  `$v`) evaluated where the macro is called, and makes `rhs` a function of the bindings,
  `~x` becoming the bound value. A rule's search goes on past a match whose right side
  gives `nothing`. `SymJ.rule(lhs => rhs)` does the same at run time, `rhs` a quoted
  template (read like a pattern, with the slots filled in) or a function. A right side
  naming a slot the left side lacks is an error when the rule is made.
- **Rewriting.** `SymJ.rewrite` walks the parts of an expression (the terms of a sum, the
  factors of a product, the arguments of a function), bottom up or top down, and at each
  tries the rules in order; the parts are rebuilt by the canonical constructors. It
  repeats until a pass changes nothing, up to `maxsteps` passes (with a warning).
  `replace(ex, rules...)` is one pass. The rules are indexed by the head of their left
  side (a function's name, any function, anything), so a node is tried only against
  rules that could match it, and a quick check of head, arity and the number of terms
  rejects most of those before any search.
- **Speed.** On ACPM's 62 benchmark cases the first match takes a median 0.78 μs and at
  most 3.6 μs (ACPM's associative-commutative matcher: 13.1 and 628 μs). A 66-rule set
  rewrites a 30-term expression in 0.43 ms.
- **Use in `simplify`.** For an expression with trigonometric or hyperbolic functions,
  `simplify` adds three candidates: a small rule set (`src/simprules.jl`) for the
  identities the trig normal form does not reach (`(1 - cos(2x))/2 = sin(x)^2`,
  `sin(x)/(1 + cos(x)) = tan(x/2)`), `expand_trig`, and `expand_trig` then
  `simplify_rational`. Each is kept only when it is smaller, so `simplify` stays
  size-based. A rule for `cos(2u)` asks with a predicate that the 2 be written:
  otherwise coefficient division would match `cos(x)` as `cos(2 * (x/2))`.

## Not yet done / open questions

- **Parameters and variables sharing a name.** Prohibit it, or print them distinctly?
  See Parameters above.
- **Complex numbers.** `im` is the exact power `(-1)^(1/2)`, and `Sym(z::Complex)`
  builds `a + b*im`. `real`, `imag`, `conj` and `angle` apply, as they are built, only
  the identities that hold for every complex value: constants are evaluated exactly,
  `real` and `imag` are linear, and `conj` distributes over sums, products, integer
  powers and functions without a branch cut (not `log` or fractional powers). Under
  assumptions, `refine` takes expressions apart into real and imaginary parts
  (`real(exp(im*x)) = cos(x)` for real `x`).
- **Assumptions.** The design is settled (see Assumptions above); the implementation is
  under way. Integer congruences (`mod(n, 3) ~ 1`) are to follow the core.
- **Relations.** `==` is structural and returns `Bool`. Symbolic `Eq`, `<`, etc. need
  their own constructors.
- **Combining radicals.** Construction does not combine `sqrt(2)*sqrt(3)` into
  `sqrt(6)`; `simplify` does. `72^(1/3)` prints as `2*9^(1/3)` rather than
  `2*3^(2/3)`. Prime-base normalization would fix the latter.
- **Polynomials, next steps.**
  - Zippel: transposed Vandermonde solves per degree when the leading coefficient is a
    single term (now one dense solve).
  - Brown recursion: done so far are bivariate images in dense buffers and in-place
    Newton updates. Mod-p reduction is a few percent of the time, so Montgomery or
    Barrett reduction isn't worth it.
  - A univariate dense representation.
  - Faster LLL for van Hoeij: floating-point LLL on truncated traces (van
    Hoeij–Novocin). The integer LLL is now the main cost of S₆'s 0.76 s.
  - Graded monomial orders.
  - Evaluation by Horner's scheme.
  - Threads for large products.
- **Hash-consing.** Interning nodes would make equality pointer comparison and save
  memory, at the cost of a global, thread-safe table.
- **Further features:** series, `simplify`, `solve`,
  and symbolic matrices (noncommutative `MUL`).
