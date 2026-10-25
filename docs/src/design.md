# Design notes

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
  - *Not yet:* the bit and integer-representation functions.
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

`Poly{C,M}` is a sparse multivariate polynomial. Expressions convert with
`Poly(ex, x, y, …)` and convert back with `Sym(p)`. `Poly{C}` (any word `M`) is what
signatures use.

```julia
struct Poly{C,M<:Integer} <: AbstractSymbolic
    gens::Vector{Sym}       # generators, most significant first
    layout::_Layout         # field width per generator (with shifts and masks)
    monos::Vector{M}        # packed exponent vectors, strictly decreasing
    coefs::Vector{C}        # nonzero coefficients
end
```

- **Packed monomials.** Each exponent vector is packed into one integer word, with the
  first generator in the most significant field. Multiplying monomials is integer
  addition (degree checks rule out carries), and lex order is integer order. So terms are
  sorted parallel vectors, addition is a linear merge, `diff` is a subtraction, and
  changing generators means repacking. The all-ones field is reserved, so each exponent
  is at most `2^width − 2`.
- **Layouts.** The usual layout is uniform: `64 ÷ ngens` bits per field in a `UInt64`.
  4 generators allow degree 65 534 each, and 16 allow 14. All such polynomials share
  one precomputed layout, so mixing them never repacks.
  - **Uneven fields.** When the uniform layout can't hold the degrees, fields are sized
    to each generator's degree and the spare bits shared out as headroom. One generator
    of high degree among many linear ones then still fits one `UInt64`: `x^(2^40)` with
    15 other generators, or `x^70` with 10.
  - **Wider words.** Past 64 bits the polynomial is packed into a `UInt128`, then a
    `BigInt` (any width). Wider words are slower, but no degree or number of generators
    is refused.
  - **Construction.** `Poly(ex, …)` tries the uniform layout first, the common case at
    no extra cost. If a field overflows, it computes degree bounds for `ex` and builds
    once in a layout sized to them.
  - **Arithmetic keeps its word.** A product, power or division needing a wider field
    gets a new layout in the same word, which is type-stable. It never changes the
    word: when results could change type, multivariate `factor` ran 70% slower. Past
    the word's capacity it throws `_PackOverflow` (shown as an `OverflowError`).
    `BigInt` words just grow their fields.
  - **Entry points restart.** `gcd`, `lcm`, `factor`, `factor_squarefree`,
    `squarefree_part` on polynomials, and `cancel`, `factor`, `apart`, `gcd`, … on
    expressions, catch `_PackOverflow`. They rerun the computation with a raised minimum
    field width (`_with_widening`; task-local, so nested calls compose). Internally
    they call non-restarting versions (`_gcd_poly`, `_factor_poly`, …), so callers stay
    type-stable.
  - **Direct arithmetic.** `p * q`, `p^k` or `divrem` on polynomials built by hand raise
    rather than widen. Build the polynomial from an expression, or `convert` it to
    `Poly{C,UInt128}`.
  - **Limits elsewhere.**
    - The modular gcd (and so Zippel) works on `UInt64` and `UInt128` words. `BigInt`
      ones use the PRS gcd.
    - `RatFun` keeps its numerator and denominator in `UInt64` words (concrete fields
      keep parameter arithmetic fast), so parameters of very high degree are refused
      with an error.
  - **Bug fixed.** `cancel` used to return its input unchanged, silently, when the
    exponents didn't fit or there were more than 64 symbols.
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
  exponents) goes through the generic tree expansion. The
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
- **Overflow.** Lex division can raise the degrees of the later generators (`x^k` by
  `x + y^5` leaves `±y^(5k)`). Each quotient term is checked against the exponent limit:
  `BigInt` words widen and retry; fixed-size words throw `_PackOverflow` for an entry
  point to restart.

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

**Zippel's sparse interpolation** is chosen per call, with no global setting:
- **The keywords.** `gcd` and `lcm`, on polynomials and expressions, take
  `method = :auto | :dense | :zippel`. The functions built on gcds take the same values
  as `gcd_method`: `cancel`, `factor`, `factor_squarefree`, `squarefree_part`, `apart`.
  A different name, because `factor(ex; method = :zippel)` would suggest a factoring
  method.
- **The default, `:auto`,** estimates the dense method's work as about
  `Π (deg_i + 1)` univariate gcds over all but the main variable, and compares it with
  the inputs' term counts, a measure of sparsity. Zippel is chosen when the estimate
  exceeds both 256 and 3 times the larger term count.
  - Measured: dense inputs have ratios near 1, and there Brown is up to 30x faster.
    Sparse or many-variable inputs have ratios of 4 and up, and there Zippel is
    3–100x faster.
  - This removed a cliff in Brown's method, which is exponential in the number of
    variables. A 14-variable low-degree gcd took 1.1 s and now takes 0.05 s; 20
    variables took more than a minute and now take 0.18 s.
- **Internal gcds.** Gcds inside `RatFun` arithmetic (linear systems, `apart` with
  parameters) and inside EEZ use `:auto`.

It changes
how the images inside Brown's loop over `x_k` are found; everything else is the same.
- **The skeleton.** The first image comes from the full recursion. Its terms are taken
  as the gcd's support.
- **Later images.** The other variables are evaluated at powers `β^j` of a random
  point. Univariate gcds in the main variable (the first active one) then give equations
  for the skeleton's coefficients.
- **Scales.** Each univariate gcd is monic, so it equals the image only up to an unknown
  scale `m_j`. That scale is an extra unknown per evaluation (the LINZIP form), so the
  gcd's leading coefficient in the main variable needn't be a single term. One extra
  evaluation checks the solution.
- **Fallback.** A degree mismatch, or a singular or inconsistent system, means the
  support is wrong or the point is bad. That image is then computed densely, and the
  skeleton refreshed from it. This also covers a content in the other variables, which
  the scales would absorb. The final division over ℤ still guarantees the result.

It loses on dense gcds, hence the automatic choice. It needs about as many terms' worth
of evaluations as the gcd has, plus a linear solve.

**Performance.** `bench/gcd.jl` compares against Nemo, which wraps FLINT. The first six
inputs are `gcd(g·a, g·b)` with `g = (x+y+z+w+1)^k`; the last two have sparse gcds of
degree 10 in four variables:

| case | auto (default) | dense | Zippel | PRS | Nemo (FLINT) |
|:--|--:|--:|--:|--:|--:|
| `k = 10`, linear cofactors | 1.64 ms | 1.64 ms | 54 ms | 4.6 ms | 0.55 ms |
| `k = 4`, degree-4 cofactors | 0.23 ms | 0.23 ms | 0.34 ms | 2.1 ms | 0.31 ms |
| `k = 4`, degree-6 cofactors | 0.36 ms | 0.35 ms | 0.43 ms | 16 ms | 0.44 ms |
| `k = 6`, degree-6 cofactors | 0.87 ms | 0.87 ms | 2.8 ms | 80 ms | 0.79 ms |
| `k = 8`, degree-6 cofactors | 1.85 ms | 1.85 ms | 15 ms | 720 ms | 1.35 ms |
| big coefficients (`3^20`, `5^15`, `7^10`), cubed | 0.60 ms | 0.59 ms | 0.45 ms | 36 ms | 0.69 ms |
| sparse gcd, 5 terms | 0.086 ms | 1.31 ms | 0.085 ms | 0.28 ms | 0.30 ms |
| sparse gcd, non-monomial leading coefficient | 0.068 ms | 0.48 ms | 0.068 ms | 1.8 ms | 0.34 ms |

- **The automatic choice** picks the faster method in every case except big
  coefficients, where the two are close.
- **Dense gcds:** faster than FLINT on four of the first six cases and within 1.4x on
  `k = 8`. The weak case is a high-degree gcd with tiny cofactors (first row), where
  Brown's method needs about `(deg+1)^(nvars−1)` univariate gcds. There about 18% of the
  time is the final division check.
- **Sparse gcds:** Zippel is 6–16x faster than dense and 3.5–5x faster than FLINT.

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
    down so that true factors give short vectors. LLL is exact and all-integer (Cohen,
    Algorithm 2.6.7): it keeps the integers d_i and λ_ij in place of the rational
    Gram–Schmidt data, so every division is exact and there are no gcd normalizations.
    Exactness matters because a single short vector is taken as proof of irreducibility.
    Compared with exact rationals this is 13x faster on S₅ (0.98 s → 74 ms) and 8x on
    S₆ (6.1 s → 0.76 s). FLINT is still about 60x faster on S₅; closing that gap needs
    floating-point LLL and traces truncated to their leading bits (van Hoeij–Novocin).
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
  decoding subset products) is the fallback when the EEZ cannot be used, and the
  reference in the tests. (An EEZ whose intermediate exponents overflow the packing is
  rerun with wider words by `factor`.) Its image has degree about `Π (deg_i + 1)` and is capped at
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
  - **Except square roots of rationals.** When the reduced denominator involves roots
    `sᵢ = sqrt(cᵢ)` (up to three), both sides are multiplied by its other conjugates
    (each `sᵢ → ±sᵢ`), which makes it rational. The gcd of that denominator and every
    coefficient of the numerator, grouped by its `s`-monomial, is a common factor. When
    only the numerator has roots, the same with the roles swapped:
    `(x + sqrt(2)y)/(x^2 − 2y^2) → 1/(x − sqrt(2)y)`. The result is kept only if the
    denominator's degree in the other atoms goes down, so `1/(x − sqrt(2))` is left
    alone. These roots are ordered last among the generators, so the sign normalization
    (positive leading coefficient) is decided by the variables.
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

## Solving

`solve(eq, x)` extends `CommonSolve.solve`, so code written against CommonSolve's
interface works with SymJ equations.

- **Input.** The equation must be a `Relation` such as `x^2 ~ 4`. A bare expression is
  never taken to mean "= 0", on pedagogical grounds, and gives an error suggesting
  `ex ~ 0`.
- **Method.**
  - `lhs − rhs` is brought to lowest terms with `cancel`, and the numerator's roots are
    the solutions. Numerator and denominator share no factor, so no root is a pole:
    `(x² − 1)/(x − 1) ~ 0` gives `[-1]`.
  - The numerator is factored; other symbols and parameters may appear in the
    coefficients.
- **Each irreducible factor in `x`:**
  - degree 1 and 2: the usual formulas;
  - `a·xⁿ + b`: all `n` roots exactly, `(−b/a)^(1/n)·(−1)^(2k/n)`;
  - degree 3: Cardano. With numeric coefficients and three real roots it uses the
    trigonometric form instead, which is exact and real (Cardano's radicals would need
    complex numbers, the *casus irreducibilis*): `x³ − 3x + 1 → 2cos(2π/9), …`;
  - degree 4: Ferrari, with biquadratics solved as quadratics in `x²`;
  - degree ≥ 5 with numeric coefficients: numeric roots by Aberth's iteration, with no
    dependency; with symbolic coefficients, an error.
- **Output.** A vector of distinct roots. Numeric roots are sorted, real ones first.
- **Equations in one function of `x`.** When `x` appears only through one function
  application `g` (`sin(2x)`, `exp(x)`, `sqrt(x + 1)`, `2^x`), `u` replaces `g`, the
  polynomial in `u` is solved, and `g(t) = r` is inverted into equations for its
  argument, solved recursively.
  - Inversions: `sin(t) = r` gives `t = asin(r)` and `π − asin(r)` (principal solutions
    in one period); `cos` gives `±acos(r)`; `bᵗ = r` gives `log(r)/log(b)`;
    `t^(p/q) = r` gives `r^(q/p)`.
  - Candidates are checked by substitution, which drops extraneous roots
    (`sqrt(x) ~ −1` has none).
  - Several powers of one base count as one function: `exp(2x)` and `exp(x)` are
    powers of `u = exp(x)`, and `4^x` and `2^x` are powers of `2^x`.
  - **Several logarithms.** When `x` appears only inside logarithms and the equation is
    linear in them with rational coefficients, `Σ cᵢ log(fᵢ) + R = 0` (constant logs
    such as `log(6)` included, `R` free of logs) is exponentiated to
    `Π fᵢ^kᵢ = exp(−L·R)` with integers `kᵢ = L·cᵢ`. Exponentiating forgets the branch
    (`log(−2) + log(−3) = log(6) + 2πi`), so candidates are checked in the original
    equation: `log(x) + log(x − 1) ~ log(6)` gives `[3]`, not `[−2, 3]`.
  - A square root mixed with `x` is reduced to `A + B·sqrt(t) = 0` using `R² = t`,
    squared to the polynomial `A² = B²t`, and extraneous roots dropped:
    `sqrt(x + 1) ~ x − 1` gives `[3]`.
- **Lambert W.** When `x` appears both directly and inside an exponential
  `b^(γx + δ)`, with the equation linear in that exponential, it reduces to
  `(αx + β)·e^(γ′x + δ′) = K`. Its solution is
  `x = W(K·γ′/α·e^(γ′β/α − δ′))/γ′ − β/α`, using `lambertw` from LambertW.jl (a
  dependency, re-exported). The branch `W₋₁` is added when the argument is real in
  `(−1/ℯ, 0)`, where there are two real solutions. Examples: `x·eˣ ~ 1 → W(1)`,
  `2^x ~ 3x`, `eˣ ~ x + 2`.
- **Inequalities.** `solve(x^2 ≫ 1, x)` returns a `SymSet`, a union of disjoint
  intervals with exact endpoints: `(-∞, -1) ∪ (1, ∞)`.
  - **Types.** `SymSet` and `SymInterval` are deliberately not exported. An interval
    has `Sym` endpoints, with `nothing` for ±∞, and open/closed flags; a point is
    `[a, a]`, printed `{a}`.
  - **Method.** The critical points are the real roots of the numerator and
    denominator of `lhs − rhs`. Between consecutive ones the sign is constant, so one
    sample decides each open piece. Each critical point is decided exactly, and poles
    are excluded. Adjacent included pieces merge.
  - **Set operations.** `∪`, `∩` and `setdiff` reuse that construction over the
    operands' endpoints. `in` compares against the exact endpoints.
  - **Limits.** Numeric coefficients only: with parameters, the order of the endpoints
    depends on their values, which needs assumptions.
- **Systems.** `solve(eqs, vars)` returns a vector of solutions, each a `Dict` from
  variables to values, for linear and nonlinear systems alike. Every unknown is a key
  of every solution; one left free maps to itself. One rule means generic
  code can loop over the result; SymPy's return type depends on the kind of system.
- **Linear systems** have at most one solution in that vector.
  - **Underdetermined systems** give the pivot variables in terms of the free ones,
    which map to themselves: `solve([x + y + z ~ 1, x - z ~ 0], [x, y, z]) ==
    [Dict(x => z, y => 1 − 2z, z => z)]`.
  - **Inconsistent systems** give an empty vector.
  - **Method.** Gauss–Jordan elimination over an exact field: `SymNumber` when the
    coefficients are numeric, otherwise `RatFun` in the other symbols and parameters.
    Its canonical form makes the pivot zero tests exact.
  - **Parameters.** The result is the generic solution, assuming nonzero pivots (as
    SymPy does): `[a·x + y ~ 1, x − y ~ b] → x = (b + 1)/(a + 1)`.
  - **Linearity.** An equation must have no unknown in its denominator, and be of
    degree ≤ 1 in the unknowns.
- **Systems of inequalities.**
  - **One unknown:** each relation is solved and the sets intersected; equations
    contribute their real roots as points: `[x^2 ~ 4, x ≫ 0] → {2}`.
  - **Several unknowns, linear:** the equations are solved and substituted first. Then
    exact Fourier–Motzkin elimination removes the unknowns in the given order,
    recording each one's bounds in terms of the unknowns still to come. Constraints
    are normalized and deduplicated (the tighter one kept), and bounds beaten by a
    parallel tighter one are dropped. An infeasible system gives `∅`.
  - **Result:** a cylindrical description, `CylindricalSet` (not exported). The last
    unknown ranges over an interval; each earlier one lies between bounds in the later
    ones, using `max`/`min` for several bounds; unknowns fixed by equations come last.
    For example `[x + y ≪ 2, y ≪ 3] → y ∈ (-∞, 3), x ∈ (-∞, -y + 2)`.
  - **Calling a set.** A set substitutes when called, like expressions and relations.
    A `CylindricalSet` keeps the constraints it was solved from:
    - with every unknown given, `C(x => 1, y => 4)` is a `Bool`;
    - with some given, `C(y => 1)` re-solves exactly in the rest: a `SymSet` for one
      unknown, a `CylindricalSet` for several. That's simpler and more reliable than
      slicing the nested bounds;
    - `pt in C`, for a point given as a `Dict`, is the same test.

    A `SymSet` from `solve(ineq, x)` records `x`, so both `S(2)` and `S(x => 2)` test
    membership.
  - **Not yet:** nonlinear systems in several unknowns (cylindrical algebraic
    decomposition), `≠` in several unknowns (not convex), and symbolic coefficients.
- **`<` and friends.** `<`, `≤`, `>`, `≥` compare constant expressions as numbers,
  returning a `Bool`. With a symbol, they raise an error naming the symbolic operator
  to use instead (`≪`, `≦`, `≫`, `≧`).
- **Not yet:**
  - nonlinear systems;
  - several unrelated functions of `x` (`sin(x) + cos(x) ~ 1`);
  - exponents that aren't rational multiples of one another (`9^x ~ 3^(x+2)`);
  - Lambert-type equations beyond the linear form (`x·e^(x²) ~ 1`).

Two canonical-form additions came with it:

- **Positive constants come out of radicals.** `(c·m)^e = c^e·m^e` for `c > 0` on the
  principal branch, so `sqrt(4a) = 2sqrt(a)` and `sqrt(2a + 2b) = sqrt(2)·sqrt(a + b)`,
  where the content of a sum is the gcd of its coefficients.
- **Inverse trig functions have exact values** at the standard points (±1/2, ±√2/2,
  ±√3/2, and √3, √3/3 for `atan`).

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

## Simplification

- **`refine`** applies the identities that need facts (`sqrt(x^2) = abs(x)` for real `x`,
  `log(exp(x)) = x`, `sin(x + n*π)` for integer `n`), each where the assumptions decide
  its condition, bottom up.
- **`simplify`** refines, simplifies function arguments, then searches: specialized
  routines (`simplify_rational`, `simplify_powers`, `simplify_trig`, expansion,
  factoring) each give a
  candidate, and the one with the fewest nodes is kept, the input if none is smaller.
  As in SymPy and Mathematica, a front end over specialized routines; their names start
  with `simplify_`, and only `simplify` and `refine` are exported.
- **Rules are Julia functions** on the canonical form, not a general pattern matcher:
  the canonical form already orders sums and products, which is most of what a matcher's
  associative-commutative matching is for, and the rule sets are small.
- **`simplify_trig`.** Arguments are grouped into base angles (multiples of one
  another share a base), and sines and cosines of multiples written as polynomials in
  `S = sin g`, `C = cos g` (Chebyshev). Reducing modulo `S² + C² - 1` is a normal form;
  with one relation per base, whose leading term is `S²` or `C²`, the reduction is the
  substitution `S² → 1 - C²` (or the reverse), both tried. Linearizing through
  `S = (z - 1/z)/2i`, `C = (z + 1/z)/2`, `z = exp(i g)`, turns products into sines and
  cosines of sums of multiples of the bases (Fu's product-to-sum rule), which gives
  multiple and sum angles. The smallest candidate is kept, `sin/cos` shown as `tan`.
- **`simplify_hyperbolic`** is the same method in another mode: `sinh`, `cosh` with
  `C² - S² = 1`, linearized through `z = exp(g)`. In both modes exponentials join in
  (`exp(i k g) = cos(kg) + i sin(kg)`, `exp(kg) = cosh(kg) + sinh(kg)`), and the
  exponential form, `z` substituted back, is one more candidate.
- **Presentation.** Results write `sin/cos` as `tan` (a negative power as `cot`, and a
  partial match as `tan·sec`), and a remaining negative power of a cosine or sine as a
  power of `sec`, `csc`, `sech` or `csch`. The size counts a negative exponent as one
  more node (a division, as SymPy's `count_ops` counts it), so `csch(x)^2` is smaller
  than `1/sinh(x)^2`.
- **`simplify_log`** combines logarithms of positive arguments into one, or expands a
  logarithm over its positive factors (a positive factor does not change the argument
  of a product, so the rest may be anything), whichever is smaller.
- **`simplify_radicals`** denests square roots of `a + b sqrt(c)` when `a² - b²c` is a
  rational square, rationalizes denominators `a + b sqrt(c)` by the conjugate, and
  factors radicands so that `refine` can take a power out of the root.

## Summation

- **`SymJ.summation(f, k, a, b)`**, a chain of methods as in SymPy: linearity; direct
  addition for a short numeric range; an antidifference `F` (`F(k+1) - F(k) = f(k)`)
  giving `F(b+1) - F(a)`, by Bernoulli polynomials for polynomials and Gosper's
  algorithm for hypergeometric terms (the term ratio from polynomials, `r^k`,
  factorials, binomials, `gamma`; the Gosper form through `resultant` and `gcd`; the
  polynomial `x(k)` from its degree bound and a linear system); partial fractions for
  rational functions without one (harmonic numbers, polygamma, exact `ζ(2n)` by
  Bernoulli numbers); known sums (the binomial theorem; `₀F₀`, `₁F₀`, `₀F₁`
  hypergeometric series: `exp`, `cos`, `sin`, `cosh`, `sinh`).
- **Infinite sums** are limits of the antidifference (`limit`). A polynomial times
  `r^k` converges exactly when `-1 < r < 1`, decided by cases on a parameter.
- **Otherwise** the sum stays unevaluated, `Sum(f, k, a, b)`, which adds up once its
  bounds become numbers. The name is `summation`, not a method of `Base.sum`, since the
  call is different (an index and bounds, not an iterable).

## Linear algebra

- **Exact elimination, first-nonzero pivoting.** Julia's generic `det`, `inv` and `lu`
  pivot by comparing magnitudes, which has no meaning for expressions. SymJ eliminates
  over a field with an exact zero test: exact rationals for numbers, otherwise
  rational functions (`RatFun`) in the entries' atoms, whose canonical form makes a
  zero pivot certain.
- **`det`** is fraction-free (Bareiss) for polynomial entries, every division exact;
  **`inv`** is the adjugate (Bareiss cofactors over the determinant) for polynomial
  entries up to 6×6, Gauss–Jordan over `RatFun` otherwise. **`lu`** returns a
  `LinearAlgebra.LU`, so `F.L`, `F.U`, `F.p`, `det(F)` and `F \ b` work as usual.
- **Names.** Methods extend LinearAlgebra's functions (`det`, `inv`, `lu`, `qr`, `rank`,
  `nullspace`, `eigvals`, `eigvecs`, `eigen`) and Base's (`exp`, `sin`, `cos`, `^`);
  functions LinearAlgebra does not have stay unexported: `SymJ.rref`,
  `SymJ.jordan_form`, `SymJ.charpoly`.
- **`SymJ.rref`, `rank`, `nullspace`** come from the same elimination; symbolic pivots are
  taken to be nonzero (the generic case), as for linear systems. **`eigvals`** are the
  roots of the factored characteristic polynomial, repeated by multiplicity;
  **`eigvecs`** the null spaces of `A - v*I`; **`eigen`** a `LinearAlgebra.Eigen` for a
  diagonalizable matrix.
- **`qr`** is Gram–Schmidt with exact arithmetic (norms are square roots, simplified),
  symbols taken as real in the inner product; it returns a small `SymQR` (Julia's `QR`
  stores Householder reflectors). **`SymJ.jordan_form`** builds Jordan chains from the
  kernels of `(A - λI)^k`. **Matrix functions** (`exp`, `sin`, `cos`, `A^n` for a
  symbolic `n`) go through the Jordan form: on a block for `λ`, `f(J)` has
  `f^(j)(λ)/j!` on its `j`-th superdiagonal.
- **Algebraic constants.** As rational functions, `sqrt(2)` or `im` would be
  independent atoms, and `sqrt(2)^2 - 2` would not test zero. A matrix with such
  constants (radicals of numbers or of expressions, `im`) is eliminated over `Sym`
  values normalized after each operation (expanded, radicals combined, denominators
  rationalized), which the canonical form makes exact. A quadratic factor of a
  characteristic polynomial whose discriminant is `-u²` gets the roots `(-b ± i u)/2a`
  rather than ones with `sqrt(-u²)`.
- **`I` and scalars.** `A ± I` works (`Sym ± UniformScaling`); `UniformScaling` needs a
  `Number` scale, so a symbolic multiple is `λ*I(n)`. An expression times an array is
  elementwise, as for numbers.
- **LinearAlgebra** is a dependency (a standard library, loaded at startup anyway).

## Assumptions

Assumptions are passed per call or for a block, never stored on symbols. This keeps
canonical form context-free: `x` is always the same node with the same hash, and
construction never consults assumptions. The trade-off is that identities such as
`sqrt(x^2) = x` apply only inside the operations that take assumptions, never during
construction.

- **Facts.**
  - *Relations:* `x ≫ 0`, `a ≦ b + 1`, `a ≶ 0` (≠).
  - *Membership:* `n ∈ ℤ`, with `ℂ ⊃ ℝ ⊃ ℚ ⊃ ℤ`. On a number, `∈` is still a `Bool`.
  - *Combining:* a list of facts is their conjunction, and `!` negates a fact.
  - *Intervals:* `0 ≦ x ≪ 1` cannot be written as a chain, because Julia lowers chains
    through `&&`; it is two facts, or a chain inside `@assume`.
- **Domains.** Symbols are complex by default, so every identity that holds over ℂ is
  safe. Any order fact about a symbol makes it real.
- **Deciding power.**
  - Facts about one symbol, each a domain and a `SymSet`. Within one symbol, And, Not
    and Or all reduce to one set, so they cost nothing extra.
  - Propagation of sign, reality and integrality through expressions.
  - Linear relations among symbols, treating nonlinear monomials as atoms (`a*b > 0`
    is linear in the atom `a*b`), decided by Fourier–Motzkin elimination.
  - No Or across symbols, which would multiply cases, and no CAD.
  - Queries are three-valued: true, false or unknown.
- **`Assumptions`.** `Assumptions(facts...)` normalizes the facts once: per-symbol
  domains and sets, plus linear relations. It checks consistency, throwing on a
  contradiction, and caches sign queries, so reusing it is cheap.
- **Scope.** As in Mathematica:
  - per call: `solve(eq, x; assume = A)`, where a tuple of facts also works;
  - for a block: `assuming(facts...) do … end`, or `@assuming facts... begin … end`,
    through a `ScopedValue`, so it nests and is thread-safe;
  - a keyword's facts inside a block are added to the block's.
- **Cases.** As in GruntzLimit's `limit`: an operation returns a plain answer when the
  assumptions decide it, and otherwise a `Cases`. That is a vector of
  `conditions => answer` pairs; calling it with parameter values gives the answer for
  those values, and indexing with them gives their case.
- **First users.** `solve` (filtering roots, and splitting on the signs of
  parameters) and `limit`, which hands the assumptions to GruntzLimit as signs and
  returns SymJ's answers; `refine` and `simplify`.
- **`@assume`**, named after the keyword, writes facts with Julia's operators:
  `@assume(x > 0 && 1 < y <= 2, n in Integers)`. Chained comparisons are split, and
  `&&`, `||`, `!` combine. The number sets have ASCII names, `Integers`, `Rationals`,
  `Reals` and `Complexes`.

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
