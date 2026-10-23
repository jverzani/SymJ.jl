# Implementation

This page discusses the main design choices in SymJ and the reasons for them. The
[Design notes](@ref) record each decision in more detail, with measurements.

SymJ set out to be a symbolic kernel in pure Julia with the speed of SymEngine, exact
numbers, and expressions that behave like functions. Most of the choices below follow
from wanting speed without giving up exactness, and type stability without giving up
generality.

## One concrete node type

Every expression is a `Sym`, a single concrete type with a `kind` tag (number, symbol,
constant, sum, product, power, function application). This is SymEngine's
approach, and differs from a Julia type per kind of node:

```julia
mutable struct Sym <: AbstractSymbolic
    const kind::Kind              # NUMBER SYMBOL CONSTANT ADD MUL POW FN, …
    const hash::UInt              # structural hash, computed once
    const name::Symbol            # symbol, constant or function name
    const num::SymNumber          # value, constant term or coefficient
    const args::Vector{Sym}       # terms, bases, [base, exponent] or arguments
    const exps::Vector{Sym}       # MUL: the exponent of each base
    const coefs::Vector{SymNumber}# ADD: the coefficient of each term
end
```

- **Concretely typed children.** A sum's terms are a `Vector{Sym}`, so walking a tree
  never boxes a value, and the compiler infers everything.
- **Branching instead of dispatch.** The hot paths, `add`, `mul` and `pow`, branch on a
  byte; with a type per kind they would dispatch dynamically on the children.
- **A stored hash.** Each node computes its structural hash once. Equality compares
  hashes first, and dictionaries keyed by expressions never rehash a subtree.
- **Identity.** The struct is `mutable` with all fields `const`, so `===` is pointer
  identity, a free first test for equality, and storing a node copies a pointer, not a
  wide struct.

The cost is that Julia's dispatch cannot tell sums from products: functions inspect the
kind. Function applications are not separate types either. `sin(x)` is a node naming
`sin`, and its rules (exact values, symmetries, derivative) live in a registry
(`register_function!`), which package extensions add to for SpecialFunctions,
HypergeometricFunctions and PolyLog.

## Canonical form as you go

Sums collect like terms and products collect like bases when they are built:
`x + 2x` is `3x`, `x*y*x` is `x^2*y`, `exp(x)*exp(y)` is `exp(x + y)`. There is no
separate simplification pass for these identities, so `a - a == 0` holds for any `a`.

A sum is stored as a constant plus sorted parallel vectors of terms and coefficients, a
product as a coefficient plus sorted vectors of bases and exponents. This is SymEngine's
dictionary representation, laid out as sorted vectors instead of hash tables. Adding two
sums is a merge of sorted vectors, and equality is elementwise comparison.

**The storage order is by hash.** The canonical order compares the stored hashes first,
and the structure only when two hashes collide, which practically never happens. Sorting
a sum's terms is then a sort of integers. When the order compared structure, it took over
40% of the time to convert a large polynomial to an expression; sorting by hash made that
conversion 3.4x faster.

The order has no mathematical meaning, so printing computes its own (by degree, then by
name). It is deterministic within a Julia version, but `Symbol` hashes differ between
versions. So anything that depends on the order of the variables can differ between
versions, such as which term leads and therefore the sign of a gcd. The tests allow for
this.

## Exact numbers

Numbers are `SymNumber`, which holds a small rational in two `Int`s, a
`Rational{BigInt}`, or a `Float64`.

- **The fast path.** Arithmetic on small values uses `Int` with overflow checks. On
  overflow it moves to `BigInt`, and back when a result fits again. Most coefficient
  arithmetic never allocates.
- **Exactness.** No answer is ever wrong because of overflow. (In the CAS comparison,
  another Julia system's `expand((x + y)^100)` silently overflowed `Int`.)
- **Radicals.** A number to a rational power is reduced: `sqrt(8)` is `2*sqrt(2)`,
  `4^(1/4)` is `sqrt(2)`. Standard values of functions are exact: `sin(π/3)` is
  `sqrt(3)/2`.
- **Floats are contagious.** A float in an operation makes the result a float, so an
  exact answer is never silently replaced by an approximation and back.

## Expressions as functions

`Sym` is a subtype of `Function`. The goal is for expressions to be drop-in replacements
for functions in numeric code:

- **Calling with a number evaluates** with Julia's arithmetic on that number: if
  `u = x^2`, then `u(3) === 9`, `u(3.0) === 9.0` and `u(1//2) === 1//4`. Because this is
  generic, dual numbers (automatic differentiation), intervals and `BigFloat` work.
- **Calling with pairs substitutes.** A single value is accepted only when there is one
  free variable; with several, an implicit order would be confusing.
- **Calling with keywords names the values**: `ex(; x = 2, y = 3)`. It evaluates with
  Julia arithmetic when every variable gets a number, and substitutes
  otherwise. Names that don't appear are ignored, so one set of keywords can serve
  several expressions. The keywords are matched to the expression's variables by
  name, which is possible because a variable is identified by its name.
- **Calling with no arguments** gives the value of a constant expression, and returns the
  expression unchanged if it has a free variable. It does not throw, deliberately: a
  Julia exception costs about 40 μs (the backtrace), and generic code that reads numeric
  leaves this way (GruntzLimit does) would pay that for every symbol.
- **`lambdify` compiles.** Evaluation by calling interprets the tree (about 2.4 μs for
  `sin(x)^2 + 3x^5 - 2x + exp(-x/2)`). `lambdify` generates Julia code, with shared
  subexpressions computed once, and compiles it with RuntimeGeneratedFunctions, so it can
  be called immediately. The result runs as fast as a hand-written function (9.3 ns
  against 9.2 ns for that expression).


## Polynomials

Polynomial algorithms (expansion, division, gcd, factoring, partial fractions, solving)
run on a separate type, `Poly{C,M}`, a sparse multivariate polynomial: a vector of
generators, sorted vectors of monomials and coefficients, and a description of how
exponents are packed.

### Packed monomials

A monomial's exponents are packed into one machine word, with the first generator in the
most significant bits. For `x^3*y^2` in two generators with 32-bit fields, the word is
`3 << 32 + 2`. Then:

- multiplying monomials is adding words (degree checks rule out carries);
- lexicographic order is integer order, so sorting terms sorts integers;
- divisibility and differentiation are a few shifts and subtractions.

Products accumulate in a purpose-built hash table keyed by the packed monomial (about 3x
faster than a `Dict`), then are sorted once. When all coefficients are machine integers
and a bound shows the sums cannot overflow, they accumulate in `Int128`, which removes
nearly all of the coefficient cost. Powers use the multinomial theorem when the base is
sparse: `expand((a + x)^100)` takes 26 μs.

`expand` converts to `Poly` whenever an expression is a polynomial with numeric
coefficients, which is why SymJ expands large products far faster than tree-based
systems: `expand(e*(e + w))` with `e = expand((x + y + z + w)^15)` takes 5 ms, against
415 ms for SymEngine.jl.

### When the exponents don't fit

The usual packing gives each of `n` generators `64 ÷ n` bits: 4 generators allow degree
65 534 each, and 16 allow 14. Exponents that don't fit are handled in two stages:

1. **Uneven fields.** The fields are sized to each generator's degree, with the spare
   bits shared out. One generator of high degree among many linear ones still fits one
   64-bit word: `x^(2^40)` with 15 other generators, or `x^70` with 10.
2. **Wider words.** Past 64 bits the polynomial is packed into a `UInt128`, then a
   `BigInt`, which is slower but has no limit.

The word type is a type parameter, `Poly{C,M}`, and arithmetic never changes it. This is
a trade-off between generality and speed. A first version widened inside each
operation, so the result of `p * q` could be a different type from `p`. Julia could no
longer infer types through polynomial code, and code doing many small polynomial
operations slowed down: multivariate factoring by 70%. Now an operation that needs a
wider word signals it, and the public entry points (`Poly(ex, …)`, `gcd`, `factor`,
`cancel`, `apart`, …) catch the signal and rerun the computation with wider words. The
common case keeps its speed, and only direct arithmetic on hand-built polynomials can
still raise an error.

### Division

Multivariate division keeps the dividend in the same hash table used for products, with
a heap of its monomials. Each step takes the largest monomial; everything added to the
table is smaller, so no monomial is revisited and the quotient and remainder come out
already sorted, without forming intermediate polynomials.

## Greatest common divisors

Polynomial gcds are at the heart of rational function arithmetic (`cancel`), squarefree
and full factorization, partial fractions and solving, so their speed sets the speed of
much else.

### Why a modular algorithm

The classical method is a sequence of pseudo-remainders, made primitive at each step.
Its coefficients grow quickly, and it is very slow beyond small inputs: 720 ms for a
4-variable gcd that the modular method does in 1.9 ms. SymJ keeps it only as a reference
implementation for the tests, and for polynomials packed in `BigInt` words.

The default is **Brown's dense modular algorithm** (in the form of Geddes, Czapor and
Labahn):

- **Images mod primes.** The gcd is computed modulo primes below 2^31, so products fit in
  an `Int`, and the images are combined by the Chinese remainder theorem.
- **Variables by interpolation.** Modulo a prime, the last variable is evaluated at
  points, the gcds of the evaluated polynomials are computed recursively, and the
  variable is restored by Newton interpolation. The base case is Euclid's algorithm on
  dense vectors.
- **Unlucky primes and points.** A prime or point can give an image of too high degree.
  These are detected by comparing leading monomials and discarded.
- **Certainty.** The candidate is accepted only if it divides both inputs exactly. The
  answer is therefore certain, not probabilistic. The division is tried only once the
  reconstruction stops changing, since a wrong candidate has coefficients spread over
  the whole modulus.

The packed monomial layout pays off again: with the eliminated variable in the least
significant field, evaluation keeps terms sorted, and terms that share their other
exponents are contiguous runs, which are exactly the coefficients of a polynomial in that
variable. Evaluation and interpolation are linear scans.

Against FLINT (through Nemo), on gcds of `(x+y+z+w+1)^k` times two cofactors, this is
faster on four of six cases and within 1.4x on another. The weak case is a high-degree
gcd with small cofactors, where the method needs about `(deg + 1)^(nvars − 1)` small
gcds.

### Sparse gcds: Zippel's method

Brown's method interpolates every possible term, so its cost is exponential in the
number of variables, whatever the size of the answer. A gcd of 14 low-degree variables
took 1.1 s, and 18 variables more than a minute. Zippel's method learns the gcd's terms
from one image and then solves only for their coefficients:

- **The skeleton.** The first image comes from the full recursion. Its terms are taken
  as the gcd's support.
- **Later images.** The other variables are evaluated at powers of a random point.
  Univariate gcds then give a linear system for the skeleton's coefficients, with an
  unknown scale for each evaluation (the LINZIP form), so the leading coefficient can
  have several terms.
- **Fallback.** Any inconsistency (a degree mismatch, a singular system) means the
  support was wrong or the point bad. That image is computed densely instead, and the
  skeleton refreshed. The final division still guarantees the result.

Zippel's method loses on dense gcds, by up to 33x, and wins on sparse ones, by 6–16x.
So the default, `method = :auto`, chooses per call. It estimates Brown's work as the
product of `(degree + 1)` over the variables and compares it with the number of terms of
the inputs, a measure of sparsity. On the benchmarks, dense inputs have ratios near 1
and sparse ones 4 and up, so Zippel is chosen above 3 (and above an absolute minimum).
The 14-variable gcd now takes 0.05 s and the 20-variable one 0.18 s, while dense gcds
are unchanged. `method = :dense` and `method = :zippel` force a choice.

## Resultants

A resultant is the determinant of the Sylvester matrix, which is never formed. There are
two methods, and `resultant` chooses between them:

- **Modular, for dense inputs in few variables.** The resultant is computed modulo
  word-size primes, up to a bound on its coefficients. For each prime, the other
  variables are evaluated on a grid of points, each point takes Euclid's algorithm in
  `Z_p[x]`, and the values are interpolated back.
  - Keeping the degrees formal makes every prime and point usable. Where a leading
    coefficient vanishes, the Euclidean recurrence still gives the determinant of the
    full Sylvester matrix, so there are no unlucky cases to detect.
  - The coefficient bound is Hadamard's, with the 1-norms of the coefficient
    polynomials: they bound the polynomials on the unit torus, and the coefficients of
    the resultant are averages over it.
- **Subresultant PRS, for sparse inputs in many variables.** With symbolic coefficients
  the grid is far larger than the answer: two cubics `a x³ + …` and `e x³ + …` give
  4⁸ points and 34 terms. Here the subresultant algorithm on sparse polynomials is
  several hundred times faster.

The choice compares the grid with the product of the input sizes.

## Gröbner bases

`groebner` runs Buchberger's algorithm, with the Gebauer–Möller criteria to discard
pairs and the sugar strategy to choose them.

- **Monomials** are packed into one word, like `Poly`'s, with a guard bit above each
  exponent field. Divisibility is one subtraction tested against the guard bits, and an
  exponent that outgrows its field sets a guard bit, which reruns the computation with
  wider words.
- **The monomial order** is a second packed word, linear in the exponents, so
  multiplying by a monomial adds to both words.
  - grlex packs the degree above the exponents.
  - grevlex packs the degree above the prefix sums `e₁ + … + eₖ`. At equal degree, a
    smaller last exponent is a larger prefix sum.
- **Reduction is fraction-free** over the integers, taking primitive parts after each
  reduction. With parameters, the coefficients are integer polynomials in them.
- **Lex bases go through FGLM when there are finitely many solutions.** A lex basis
  computed directly suffers coefficient growth: Katsura-4 took over five minutes. A
  grevlex basis is cheap, and FGLM converts it by linear algebra on normal forms,
  computed with multiplication matrices. Katsura-4 then takes 3 ms, and Katsura-5
  0.12 s.

## Factorization

Factoring over the rationals starts with the **squarefree decomposition**, by Yun's
algorithm: a few gcds with a derivative separate the factors by multiplicity, and each
squarefree part is then factored into irreducibles.

**Univariate polynomials** use the Zassenhaus approach: factor modulo a small prime
chosen to keep the number of factors small, lift the factors by Hensel lifting past a
bound on the true factors' coefficients, then recombine the lifted factors into true
ones, checking each candidate by exact division.

Recombination by trying subsets is exponential in the number of modular factors, and
some irreducible polynomials split into many factors modulo every prime (the
Swinnerton-Dyer polynomials). From nine factors on, SymJ uses **van Hoeij's method**:
traces of the lifted factors form a lattice in which the true factors give short
vectors, found by lattice reduction (LLL). The reduction is exact, in integers (Cohen's
integral LLL), because a single short vector is taken as proof that the polynomial is
irreducible, and floating point could not guarantee that. Keeping it in integers instead
of rationals made it 8–13x faster: the degree-64 Swinnerton-Dyer polynomial is shown
irreducible in 0.76 s.

**Multivariate polynomials** use Wang's extended Zassenhaus method:

- evaluate all variables but one at a small integer point, chosen so that the
  univariate image keeps its degree and stays squarefree;
- factor the univariate image;
- lift the factors back to all variables by multivariate Hensel lifting.

The leading coefficients are predicted before lifting, by Wang's method of matching the
factors of the leading coefficient to the univariate factors through prime divisors.
Kronecker's substitution, which maps the problem to one univariate factorization,
remains as a fallback and as a reference for the tests; it is correct but far slower
(14 s against 2.5 ms on a 4-variable example).

## Rational functions and partial fractions

`cancel` brings an expression to the form `p/q` in lowest terms. Its polynomial
variables are the expression's atoms: symbols, but also function applications and
non-integer powers. Treating `sin(x)` as an independent variable is always sound, since
a polynomial identity in it remains true after substituting back. So
`cancel((sin(x)^2 - 1)/(sin(x) + 1))` is `sin(x) - 1`. The price is that relations
among atoms are not used: `(x - 1)/(sqrt(x) - 1)` is not reduced.

Square roots of rational numbers are the exception, since `sqrt(2)^2 = 2` is known.
When the denominator involves them, both sides are multiplied by the denominator's other
conjugates (`sqrt(2) → −sqrt(2)`), which makes the denominator rational; a common factor
of that denominator and the numerator's coefficients is then removed. So
`(x^2 + 2sqrt(2)x + 2)/(x + sqrt(2))` becomes `x + sqrt(2)`. The result is kept only when
the denominator gets simpler, so `1/(x − sqrt(2))` is not rewritten as
`(x + sqrt(2))/(x^2 − 2)`.

`apart` divides out the polynomial part, factors the denominator, and finds each
fraction's numerator by modular inverses. When other symbols appear, the
decomposition is over the field of rational functions in them, with arithmetic in an
exact rational function type kept in canonical form, so that its zero tests are exact.

## Solving

`solve` follows a few rules:

- **An equation is required.** It must be written `lhs ~ rhs`. A bare expression is not
  taken to mean "= 0", on pedagogical grounds.
- **Polynomials are factored first.** The left side minus the right is reduced with
  `cancel`, so no root is a pole, and the numerator is factored. Each irreducible factor
  is solved by formula up to degree four, numerically beyond.
  - A cubic with three real roots is solved in trigonometric form, which is exact and
    real, where Cardano's formula would need complex numbers.
- **One function of the unknown.** When the unknown appears through one function, such
  as `exp(x)` or `sqrt(x + 1)`, that function is solved for first, then inverted.
  Candidates are checked by substitution, which drops extraneous roots.
- **Several logarithms.** An equation linear in logarithms of the unknown, such as
  `log(x) + log(x − 1) ~ log(6)`, is exponentiated into `x(x − 1) = 6`. That can add
  roots where the logarithms are complex (`x = −2` here), and the check by substitution
  removes them.
- **Lambert W.** Equations mixing the unknown with an exponential of it, such as
  `x*exp(x) ~ 1`, are solved with the Lambert W function.
- **Inequalities** are solved exactly through the real roots of the numerator and
  denominator, between which the sign is constant. The result is a union of intervals.
- **Systems.** Linear systems are solved by Gauss–Jordan elimination over an exact
  field. Linear inequalities in several unknowns are solved by Fourier–Motzkin
  elimination, giving each unknown's bounds in terms of the unknowns after it.

`<` and friends are not used for symbolic inequalities: they keep returning `true` or
`false`, which generic code relies on, and raise an error naming the symbolic operator
(`≪`, `≦`, …) when given a symbol.

## Series

`series` substitutes `x = x0 + t` (or `x = 1/t` at infinity) and expands by truncated
power series arithmetic in `t`, rather than by repeated differentiation, which breaks
down at a pole or a branch point:

- **Exponents are rationals**, so Laurent and Puiseux series need nothing special.
- **`log(t)` is kept in the coefficients** as a placeholder symbol, which gives the
  logarithmic terms of `x^x = exp(x log x)`.
- **Each series carries the order to which it is correct.** A product, a quotient or a
  cancellation (as in `(sin(x) - x)/x^3`) lowers it, and `series` reruns with more terms
  until the result is correct to the order asked for.
- **Functions** are expanded by their Taylor coefficients at the constant term of their
  argument, from the derivative rules; `tan`, `cot`, `sec` and the like go through `sin`
  and `cos`, so poles come out of a division. An essential singularity, such as `exp(1/x)`
  at 0, is an error.

## Working with other packages

SymJ implements TermInterface, so generic code that walks terms with `operation` and
`arguments` and rebuilds them with `maketerm` works on its expressions. GruntzLimit and
RischIntegrate are written this way, for several symbolic systems.

Julia compiles code the first time it runs, so the first limit or integral would pay
for compiling those packages' code for SymJ's types: about 3 s for the first limit and
4 s for the first integral. SymJ has package extensions for GruntzLimit and
RischIntegrate that run a precompilation workload of limits and integrals on SymJ
expressions when those packages are installed. The compiled code then ships with the
extension, and the first limit takes about 1 ms.

## How it compares

The `bench/cas` directory compares SymJ (with GruntzLimit and RischIntegrate), Symbolics
(with the same two), SymPy and Giac on 311 problems in nine areas, checking every answer
numerically. SymJ verifies 309, Giac 305, and Symbolics and SymPy 296 each. The
[comparison page](https://rawcdn.githack.com/jverzani/SymJ.jl/d9f9703f64bef58d1aa11a55ba041a7611429fa2/bench/cas/results/cas_comparison.html)
has the details, the times, and each system's answer to each problem.
