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
    `Sym(π)() === π`), and returns `ex` unchanged if it has a free symbol. It
    deliberately does not throw: a Julia `throw` costs about 40 μs (the backtrace),
    and generic code reading leaves pays that on every symbol. This matches what generic
    term code expects from a callable numeric leaf; GruntzLimit, for example, reads
    numbers this way.
  - An expression with several free symbols must be called with pairs. A single value
    would leave implicit which symbol receives it.
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
  - *Not yet:* `real`, `imag`, `conj`, `angle`, `reim`, which need to know whether
    symbols are real (assumptions), and the bit and integer-representation functions.
- **Package extensions.** These load automatically when the package is loaded with
  SymJ: **SpecialFunctions**, **HypergeometricFunctions** (`_₁F₁`, `_₂F₁`, `_₃F₂`,
  and `pFq` of those sizes) and **PolyLog** (`li0`–`li6`, `li(n, z)`, `reli…`).
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
- **Powers.** `p^k` uses the multinomial theorem when the base is sparse: each
  composition of `k` is computed once, its monomial is `Σ aᵢ·mᵢ` in packed integer
  arithmetic, and the multinomial coefficients are updated incrementally. Integer
  coefficients are computed in `Int128` when `(Σ|cᵢ|)^k < 2^125` guarantees no
  overflow, or else in `BigInt`, and converted to `SymNumber` once per term. (Doing
  the arithmetic in `SymNumber`, whose big values are rationals, re-normalized a gcd
  on every operation.) When the number of compositions is far above the degree bound
  on the result size, many compositions share a monomial, as in a dense univariate
  base, and repeated multiplication is cheaper. `(a+x)^100` went from 867 μs and 56k
  allocations to 26 μs.
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

### Squarefree factorization

`factor_squarefree(p)` returns `(unit, [f₁ => m₁, f₂ => m₂, …])` with
`p = unit·Π fᵢ^mᵢ`. The guarantees:
- each `fᵢ` is squarefree, primitive, and has a positive leading coefficient;
- the `fᵢ` are pairwise coprime;
- the multiplicities are distinct and increasing.

`squarefree_part(p)` is `Π fᵢ`.

- **Algorithm.** A primitive `f` splits into its content with respect to its most
  significant generator `xᵢ` (a polynomial in the later generators) and its primitive
  part, all of whose irreducible factors involve `xᵢ`.
  - **Primitive part:** decomposed by Yun's algorithm, a few gcds with `∂f/∂xᵢ`
    using the modular gcd.
  - **Content:** decomposed recursively in the remaining generators.
  - **Merging:** the two sets of factors are coprime, so factors of equal multiplicity
    are multiplied.
- **Speed.** 1.5 ms for a 612-term polynomial in 4 variables, 0.09 ms for a degree-55
  univariate one.
- **Expressions.** `factor_squarefree(ex)` and `squarefree_part(ex)` treat atoms the way
  `cancel` does. Rational expressions are cancelled first, then numerator and
  denominator are factored: `(x^2 - 1)/(x^2 + 2x + 1)^2 → (x - 1)/(x + 1)^3`.

### Factorization

`factor(p)` returns `(unit, [f₁ => m₁, …])`. The factors are irreducible over ℚ,
primitive, and have positive leading coefficients. They're sorted by multiplicity,
then by degree. The input is first reduced to squarefree primitive parts.

**Univariate: Zassenhaus.**
- **Choosing a prime.** Pick a small odd prime that keeps the polynomial squarefree and
  doesn't divide the leading coefficient. Of up to five candidates, keep the one with
  the fewest modular factors, since that keeps recombination cheap.
- **Factoring mod p.** Distinct-degree factorization, then Cantor–Zassenhaus to split
  equal-degree parts, with a fixed xorshift generator so results are reproducible.
- **Hensel lifting.** The two-factor lift is linear, one power of `p` per step; several
  factors are lifted through a balanced split. It runs to `p^k > 2|lc|·2ⁿ‖f‖₂`, past
  the Landau–Mignotte bound.
- **Recombination.** Try subsets of the lifted factors, smallest first. Each candidate
  `lc·Π mod p^k`, taken in the symmetric range, is tested by exact division over ℤ.
- **van Hoeij recombination** (`src/vanhoeij.jl`). The number of subsets grows
  exponentially with the number of modular factors, so from 9 factors up a lattice
  replaces the search:
  - **Traces.** For each lifted factor `gᵢ`, the power sums of its roots,
    `tᵢⱼ = ℓʲ·Pⱼ(gᵢ) mod pᵃ` (Newton's identities, `ℓ` the leading coefficient). For a
    true factor, the sum over its pieces is an integer bounded by `n(|ℓ|R)ʲ`, where `R`
    is Fujiwara's root bound. Cauchy's bound was far too loose for polynomials with
    large coefficients and small roots, and inflated every number in the lattice.
  - **Lattice.** Rows `[eᵢ | tᵢⱼ/Dⱼ]` and `[0 | pᵃ/Dⱼ]`, with the trace columns scaled
    down so that true factors give short vectors. LLL uses exact rational Gram–Schmidt.
  - **Reading off factors.** Basis vectors whose Gram–Schmidt norm exceeds the proven
    bound are dropped. If one vector is left, `f` is irreducible. Otherwise the reduced
    echelon form should be a 0/1 partition, each part confirmed by trial division.
    Until it is, the number of traces doubles, and the precision grows with it.
  - **Fallback.** Subset search remains the fallback.
  - **Speed.**

    | Swinnerton-Dyer polynomial | modular factors | subset search | van Hoeij |
    |:--|--:|--:|--:|
    | S₅, degree 32 | 16 | 3.2 s | 0.85 s |
    | S₆, degree 64 | 32 | out of reach | 6.1 s |

    These are irreducible, but split into many factors modulo every prime.

**Multivariate: Wang's EEZ approach** (`src/eez.jl`).
- **Setup.** The main variable `x` is the active variable of least degree. The content
  in `x`, a polynomial in the other variables, is factored on its own.
- **Evaluation point.** An integer point `α` for the other variables, with
  `lc_x(q)(α) ≠ 0` and `q(x, α)` squarefree. Of a few such points, keep the one with
  fewest univariate factors. A point with one factor proves `q` irreducible.
- **Shift.** `y → y + α` moves the point to 0, so Taylor coefficients are plain
  coefficients and evaluating means dropping terms.
- **Leading coefficients.**
  - **Wang's distribution.** Factor `L = lc_x(q)` into irreducibles `tⱼ`. When each
    `tⱼ(α)` has a prime divisor shared with no other and not with the contents (Wang's
    non-divisor test), the `tⱼ` dividing each univariate factor's leading coefficient
    can be read off, and the lift runs on `q` itself. Points passing the test are
    preferred, but the values stay small, since a larger range densified the shifted
    polynomial and cost far more.
  - **Imposition.** Otherwise `L` is imposed on every factor of `L^(r-1)·q`, and
    primitive parts are taken at the end.
- **Lifting.** Multivariate Hensel lifting, one variable at a time, with multivariate
  Diophantine equations (Geddes–Czapor–Labahn 6.2–6.4), over ℚ. Working over ℚ needs no
  p-adic coefficient bound, and an exact product check decides success.
- **Spurious factors.** If `α` split a true factor, the full lift fails. Then two-factor
  lifts of subsets of the univariate factors, smallest first, find each irreducible
  factor.
- **Fallback.** Kronecker substitution (`x_i → t^(w_i)`, factoring the univariate image,
  decoding subset products) is the fallback if packed exponents would overflow, and the
  reference in the tests. Its image has degree about `Π (deg_i + 1)` and is capped at
  4000.
- **Speed.** Against Kronecker:

  | polynomial | Kronecker | EEZ |
  |:--|--:|--:|
  | `(wy + x + 1)(xz − y² + 3w)(x²y + z³ + w)` | 14.2 s | 2.5 ms |
  | `x⁶ − y⁶z⁶` | 3.6 s | 1.7 ms |
  | `(x⁵y³ + z⁴w² − 7)(x⁴ − y⁵z + w³ + 2)(xyzw + 1)` | beyond its cap | 38 ms |

**Timings.**

| polynomial | time |
|:--|--:|
| 31 linear factors | 3.3 ms |
| `x^60 - 1` | 5.6 ms |
| `x^105 - 1` | 45 ms |
| the degree-16 Swinnerton-Dyer polynomial, shown irreducible | 2.2 ms |
| small multivariate inputs | under 2 ms |

**Expressions.** `factor(ex)` treats atoms the way `cancel` does, and factors the
numerator and denominator of rational expressions.

### Partial fractions

`apart(ex, x)` decomposes a rational function of `x` with exact numeric coefficients
into a polynomial plus terms `c/pʲ`, with `p` irreducible and `deg c < deg p`.

1. Division gives the polynomial part and a proper remainder `r`.
2. `factor` gives `d = c·Π pᵢ^mᵢ`.
3. For each `i`, `Aᵢ = (r/c)·(d/pᵢ^mᵢ)⁻¹ mod pᵢ^mᵢ`, by extended Euclid over ℚ.
4. Expanding `Aᵢ` in powers of `pᵢ` splits `Aᵢ/pᵢ^mᵢ` into the terms.

**Symbolic coefficients.** When other symbols or parameters appear, the decomposition
is over the field `ℚ(a, b, …)`. For example,
`apart(1/((x - a)(x - b)), x) = 1/((b - a)(x - b)) - 1/((b - a)(x - a))`.

- **Factoring.** By Gauss's lemma, the irreducible factors of the denominator over
  `ℚ(params)[x]` are its multivariate factors over ℤ that involve `x`; the others are a
  unit there.
- **Arithmetic.** The steps above run in `Poly{RatFun}`. `RatFun` (`src/ratfun.jl`) is
  an exact rational function in the parameters, kept in a canonical form so that
  equality and zero tests are exact:
  - a constant denominator becomes 1, with rational numerator coefficients;
  - otherwise, a coprime integer numerator and denominator (the modular gcd removes the
    common factor, including integer content) with a positive leading coefficient.
- **Generic code.** `Poly`'s division, `rem`, powers and modular inverses are generic in
  the coefficient type, so they work over `RatFun` unchanged.
- **Display.** Each numerator `c(x)` is written over a factored common denominator.

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

## Parameters

`parameter(:a)` or `@parameters a b` makes a symbolic parameter: a named constant of
unknown value, as opposed to a variable. It is a node kind, `PARAMETER`, beside
`SYMBOL`, not a new Julia type. Expressions stay one concrete type, so containers
stay concretely typed and dispatch can't distinguish expressions that contain
parameters; `isparameter` and `free_parameters` answer that instead. The parameter
`a` and the variable `a` hash differently and are distinct.

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
  - Faster LLL for van Hoeij: floating-point Gram–Schmidt, and adding traces one at a
    time. Exact rationals are most of S₆'s 6.1 s.
  - Graded monomial orders.
  - Wider monomials (`UInt128` or exponent vectors) for many generators or high degree.
  - Evaluation by Horner's scheme.
  - Threads for large products.
- **Hash-consing.** Interning nodes would make equality pointer comparison and save
  memory, at the cost of a global, thread-safe table.
- **Further features:** series, `simplify`, `solve`,
  and symbolic matrices (noncommutative `MUL`).
