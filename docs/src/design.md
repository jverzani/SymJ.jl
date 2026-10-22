# Design

SymJ is a symbolic manipulation kernel written in pure Julia. It aims for SymEngine-class
speed with a small, explicit core that later layers (simplification, solving, matrices)
build on.

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
    `Sym(π)() === π`), and throws if `ex` has a free symbol. This matches what generic
    term code expects from a callable numeric leaf; GruntzLimit, for example, reads
    numbers this way.
  - An expression with several free symbols must be called with pairs. A single value
    would leave implicit which symbol receives it.
  - A `SymFunction` called on arguments builds an unevaluated application (`f(x)`).

  Evaluation interprets the tree: about 25 ns for `x^2` and 2.4 μs for
  `sin(x)^2 + 3x^5 - 2x + exp(-x/2)`. A future `lambdify` will compile expressions for
  speed.
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

## Polynomials: `Poly{C}`

`Poly{C}` is a sparse multivariate polynomial. Expressions convert with
`Poly(ex, x, y, …)` and convert back with `Sym(p)`.

```julia
struct Poly{C} <: AbstractSymbolic
    gens::Vector{Sym}       # generators, most significant first
    bits::Int               # 64 ÷ length(gens) bits per exponent
    monos::Vector{UInt64}   # packed exponent vectors, strictly decreasing
    coefs::Vector{C}        # nonzero coefficients
end
```

- **Packed monomials.** Each exponent vector is packed into one `UInt64`, with the first
  generator in the most significant field. Multiplying monomials is integer addition
  (degree checks rule out carries), and lex order is integer order. So terms are sorted
  parallel vectors, addition is a linear merge, `diff` is a subtraction, and changing
  generators means repacking. The all-ones field is reserved, so each exponent is at most
  `2^bits − 2`, and too high a degree throws an `OverflowError`. For example, 4
  generators allow degree 65 534 each, and 16 allow 14.
- **Products.** Products accumulate in a purpose-built open-addressing table keyed by the
  packed monomial (multiplicative hash, linear probing), then one sort. A generic `Dict`
  was about 3x slower.
- **Coefficients.**
  - The coefficient type is `C = SymNumber` for numeric polynomials and `Sym` for
    symbolic ones (`Poly(a*x^2 + b*x, x)`). The arithmetic is generic in `C`.
  - **Integer fast path.** When every coefficient of both factors is a machine integer
    and `max|a|·max|b|·min(#p,#q) < 2^125`, products accumulate in `Int128`. That
    removes about 98% of the coefficient cost. Otherwise the code falls back to exact
    `SymNumber` arithmetic, which overflows into `BigInt` as needed.
- **Conversion.**
  - `Sym → Poly` packs monomial terms directly, without building intermediate
    polynomials.
  - `Poly → Sym` builds canonical nodes directly when the generators are symbols: bases
    visited in `_cmp` order are already sorted, and distinct monomials give distinct
    keys.
- **`expand` uses `Poly`** whenever the expression is a polynomial with numeric
  coefficients in at most 16 symbols. Anything else (functions, radicals, symbolic
  exponents, or an exponent overflow) goes through the generic tree expansion. The
  generic expansion still calls `expand`, so polynomial subexpressions such as the
  `(x+1)^2` in `sin(x)*(x+1)^2` take the fast path.

### Division

`divrem(p, q)` is multivariate division with respect to lex order. For one generator it
is ordinary long division. `div`, `rem`, `divides` and `divexact` are built on it.

- **Algorithm.** The dividend lives in the same open-addressing table used for products,
  with a max-heap of its monomials. Each step pops the largest monomial `m`.
  - If `LM(q)` divides `m`, the term `t = m/LM(q)` goes to the quotient and
    `-t·(q − LT(q))` is added into the table.
  - Otherwise the term goes to the remainder.
- **Why it is cheap.** Everything added is smaller than `m`, so no monomial is revisited,
  and the quotient and remainder come out already sorted. The work is about
  `#quotient × #q` table operations; no intermediate polynomial is ever formed.
- **Coefficients.** Division is exact in the coefficient field (`SymNumber` rationals, or
  `Sym`), so dividing integer polynomials can produce rational coefficients.
- **Overflow.** Lex division can raise the degrees of the later generators, so each
  quotient term is checked against the exponent limit.

### GCD

`gcd(p, q)` (and `lcm`) works on polynomials with exact numeric coefficients, and on
expressions that are such polynomials.

- **Normalization.**
  - Integer inputs give the integer gcd, including the gcd of the contents, with a
    positive leading coefficient.
  - Rational inputs are scaled to integers first, and the result is made monic.
  - Floating point coefficients are rejected.

**The algorithm is Brown's dense modular GCD**, in the form of Geddes–Czapor–Labahn
MGCD/PGCD (`src/modgcd.jl`).

- **Over ℤ.**
  - Compute images of the gcd modulo primes below 2^31, so products fit in an `Int`.
  - Normalize each image so its leading coefficient is `γ = gcd(lc(a), lc(b))`.
  - A prime whose image has a larger leading monomial is unlucky and is discarded. A
    smaller one means every earlier prime was unlucky, so the work so far is discarded
    instead.
  - Combine images by Chinese remaindering in the symmetric range. Then take the
    primitive part and accept it only if it divides both inputs.
  - Trial division therefore makes the answer certain, not probabilistic. It runs only
    when the reconstruction looks finished: unchanged by the last prime, or every
    coefficient at least 8 bits below the modulus. A wrong reconstruction has residues
    spread across the whole range, so this avoids failed divisions.
- **Over ℤ_p.**
  - Eliminate the last active variable `x_k` by evaluating it at points, recurse, and
    rebuild `x_k` by Newton interpolation. The base case is univariate Euclid on dense
    vectors.
  - Leading coefficients in `ℤ_p[x_k]` are handled by scaling each image to `g(α)`, with
    `g` the gcd of the leading coefficients.
  - Unlucky points give images that are too large, and are detected the same way as
    unlucky primes.
  - The number of points is `min(deg_k a, deg_k b) + deg g + 1`.
- **Data layout.** Mod-p polynomials reuse the packed lex monomials. Because `x_k` is the
  least significant active field, evaluating it keeps terms sorted. Terms sharing a
  monomial in the other variables form contiguous runs, which are exactly the
  coefficients in `ℤ_p[x_k]`. Evaluation, contents and interpolation are therefore
  linear scans.

**The recursive primitive PRS** (`_gcd`) remains as the reference implementation. The
tests check the two against each other.

**Performance.** `bench/gcd.jl` compares against Nemo, which wraps FLINT. The inputs are
`gcd(g·a, g·b)` with `g = (x+y+z+w+1)^k`:

| case | modular | PRS | Nemo (FLINT) |
|:--|--:|--:|--:|
| `k = 10`, linear cofactors | 3.4 ms | 6.4 ms | 0.61 ms |
| `k = 4`, degree-4 cofactors | 0.42 ms | 2.1 ms | 0.34 ms |
| `k = 4`, degree-6 cofactors | 0.64 ms | 18 ms | 0.47 ms |
| `k = 6`, degree-6 cofactors | 1.4 ms | 66 ms | 0.79 ms |
| `k = 8`, degree-6 cofactors | 3.4 ms | 803 ms | 1.4–2.3 ms |
| big coefficients (`3^20`, `5^15`, `7^10`), cubed | 0.85 ms | 33 ms | 0.73 ms |

The modular algorithm is within 1.2–2.4x of FLINT, except for a high-degree gcd with tiny
cofactors (first row). That case is where Brown's dense method is weakest: it needs about
`(deg+1)^(nvars−1)` univariate gcds however sparse the answer is. Zippel's sparse
interpolation addresses it: it learns the support of the gcd from one dense image, then
solves for coefficients with far fewer evaluations.

### Rational expressions

`numer_denom(ex)` writes `ex` as a single fraction without expanding anything.
`cancel(ex)` then reduces it to `p/q`:

- **Polynomials.** The numerator and denominator are converted to `Poly`, and their gcd
  is divided out. The denominator gets a positive leading coefficient, and the integer
  content is shared out.
- **Generators.** These are the expression's polynomial atoms: symbols, constants,
  function applications, and non-integer powers.
  - Treating them as independent variables is always sound, because a polynomial
    identity in them still holds after substituting them back. So
    `cancel((sin(x)^2 - 1)/(sin(x) + 1)) == sin(x) - 1`.
  - Relations among atoms are not used: `(x - 1)/(sqrt(x) - 1)` stays as it is.
- **Floats.** Expressions with floating point coefficients are returned unchanged.

## Assumptions

Assumptions (real, positive, integer) will be passed per call, e.g.
`simplify(ex; assume = …)`, not stored on symbols. This keeps canonical form
context-free: `x` is always the same node with the same hash, and construction never
consults assumptions. The trade-off is that identities such as `sqrt(x^2) = x` apply
only inside the operations that take assumptions, never during construction.

## Performance snapshot

Run with `julia --project=bench bench/bench.jl` (Apple Silicon, Julia 1.13, SymEngine.jl
0.13). Times are medians. The ratio is SymJ time / SymEngine time, so lower is better.

| benchmark | SymJ | SymEngine | ratio |
|:--|--:|--:|--:|
| `x + y` | 32 ns | 98 ns | 0.33 |
| build Σ i·xⁱ·y²⁰⁰⁻ⁱ term by term | 216 μs | 593 μs | 0.36 |
| build Π (x+i)ⁱ, i ≤ 30 | 5.1 μs | 14.7 μs | 0.35 |
| `expand((x+y+z+1)^15)` | 240 μs | 433 μs | 0.56 |
| `diff` of expanded poly · sin | 38 μs | 52 μs | 0.72 |
| `subs` into 286-term poly | 152 μs | 170 μs | 0.89 |
| `expand(e*(e+w))`, `e=(x+y+z+w)^15` | 5.1 ms | 415 ms | 0.01 |

The `expand(e*(e+w))` row takes 5.1 ms in total. Of that, the `Poly` product (816 × 817
terms, giving 6 272) takes 1.8 ms, `Poly → Sym` takes 3.2 ms, and `Sym → Poly` takes
1.5 ms. SymEngine.jl's `expand` works on the generic tree, so this compares like with
like at the user level. SymEngine's own sparse polynomial classes are not exposed by
SymEngine.jl. SymEngine.jl's numbers also include the cost of the C wrapper (a `ccall`
and a finalizer per operation), which dominates the small benchmarks.

## Not yet done / open questions

- **Complex numbers.** `im` is the exact power `(-1)^(1/2)`. That covers the basics, but
  there is no `re`/`im` decomposition or `conj`, and `Sym(z::Complex)` builds `a + b*im`.
- **Assumptions.** These will be passed per call, as described above. Neither the API
  nor the rules exist yet.
- **Relations.** `==` is structural and returns `Bool`. Symbolic `Eq`, `<`, etc. need
  their own constructors.
- **Combining radicals.** `sqrt(2)*sqrt(3)` is not combined into `sqrt(6)`, and
  `72^(1/3)` prints as `2*9^(1/3)` rather than `2*3^(2/3)`. Prime-base normalization
  would fix the latter.
- **Polynomials, next steps.**
  - Zippel's sparse modular GCD for sparse, high-degree gcds.
  - Faster mod-p arithmetic (Montgomery or Barrett reduction) and fewer small allocations in the Brown recursion.
  - A univariate dense representation.
  - Squarefree decomposition and factorization.
  - `cancel` for rational expressions, built on `gcd`.
  - Graded monomial orders.
  - Wider monomials (`UInt128` or exponent vectors) for many generators or high degree.
  - `Poly^k` by the multinomial theorem rather than repeated multiplication.
  - Evaluation by Horner's scheme.
  - Threads for large products.
- **Hash-consing.** Interning nodes would make equality pointer comparison and save
  memory, at the cost of a global, thread-safe table.
- **Further features:** `lambdify`/code generation, series, `simplify`, `solve`,
  and symbolic matrices (noncommutative `MUL`).
