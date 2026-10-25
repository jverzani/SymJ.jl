# Design note: pattern matching and rewrite rules

A proposal, written 2026-10-25, for pattern matching in SymJ: patterns with wildcards,
matching against expressions, and rules that rewrite with the matches. Nothing here is
implemented yet. The decisions marked **Open** are for the author to settle before the
prototype.

## Decision: inside SymJ, with the SymbolicUtils syntax

The matcher lives in SymJ and works on `Sym` nodes directly, not through TermInterface.
The surface syntax is the one AssociativeCommutativePatternMatching (ACPM) borrowed from
SymbolicUtils: expressions with `~x` slots make patterns, and macros compile them.

The reason is cost. Through TermInterface a matcher sees only `operation` and
`arguments`, and must assume nothing about them: sums may be nested, unsorted and
uncombined, so associative-commutative (AC) matching searches orderings of arguments and
splits of them into groups. That search is what makes ACPM's `_eachmatch` slow. In
ACPM's 62-case benchmark (`benchmark/matching`, run 2026-10-18), the first match takes:

| Matcher | Median | Slowest |
|:--|--:|--:|
| `ac`: `_eachmatch`, AC (MatchPy) | 13.1 μs | 628 μs (nested) |
| `c`: `_ceachmatch`, commutative | 1.8 μs | 15.1 μs |
| `cc`: commutative, pattern compiled once | 0.45 μs | 9.2 μs |
| `su`: SymbolicUtils `@rule` | 0.93 μs | 6.6 μs |

SymJ's canonical form removes most of the search before it starts:

- **Flat.** `ADD` and `MUL` are never nested in themselves, so associativity costs
  nothing for the subject.
- **Sorted and combined.** Arguments are in a canonical order and like terms are merged
  (`a + a` is `2a`, `a·a` is `a^2`). A pattern element with no slots is found by hash
  comparison, not by trying orderings.
- **Coefficients apart.** An `ADD` stores `(num, coefs, args)`: `2xy + 3` is
  `num = 3, args = [x*y], coefs = [2]`. A `MUL` stores `(num, args, exps)`: `x·e^(2x)` is
  `args = [x, ℯ], exps = [1, 2x]`. Coefficients, constants and exponents are read from
  fields, not found by splitting products.
- **Hashed.** Every node carries its structural hash, so a repeated slot (`sin(~x)cos(~x)`)
  is checked in constant time.
- **No allocation.** The matcher reads `Sym` fields in place. Through TermInterface,
  `arguments` on SymJ rebuilds a vector at every step.

A TermInterface-level matcher could only share this speed if every term type supplied
canonical-form hooks (is this operator AC, sorted arguments, the coefficient split).
That interface would be large and most term types could not honor it. If sharing is
wanted later, a narrow interface between the matcher and `Sym` can be added, with ACPM
as the generic fallback.

## Surface syntax

Slots, as in ACPM and SymbolicUtils:

| Syntax | With a predicate | Matches |
|:--|:--|:--|
| `~x` | `~x::pred` | exactly one argument (one term of a sum, one factor of a product) |
| `~!x` | `~!x::pred` | as `~x`, or a default when absent: 0 in a sum, 1 in a product or exponent |
| `~~x` | `~~x::pred` | zero or more arguments, bound to their sum or product (0 or 1 when empty) |
| `~~~x` | `~~~x::pred` | one or more arguments, bound to their sum or product |

A name may appear more than once (a non-linear pattern) but with one flavor. A predicate
gets the bound `Sym` and must return `true`, or the match is rejected; an error counts as
`false`.

The functions are unexported and called qualified, or they extend existing Base functions:

```julia
p = SymJ.@pattern ~a * sin(~x) + ~~r        # compiled where the macro is expanded
p = SymJ.pattern(:(~a * sin(~x) + ~~r))     # the same, built at run time from an Expr

match(p, ex)        # Dict{Symbol,Sym} of the first match, or nothing
eachmatch(p, ex)    # an iterator of all the matches
occursin(p, ex)     # a match anywhere in ex

r = SymJ.@rule sin(~x)^2 + cos(~x)^2 => 1
r = SymJ.@rule log(~x::ispos) + log(~y::ispos) => log(~x * ~y)
r = SymJ.@rule ~x^~n => ~n * ~x^(~n - 1)    # the right side is Julia, with ~x the bound values
r(ex)                                       # the rewritten expression, or nothing
SymJ.rewrite(ex, rules...; order = :post, fixpoint = true)
replace(ex, rules...)                       # every match anywhere in ex, once
```

As in SymbolicUtils, a right side that returns `nothing` rejects the match. In a quoted
pattern, predicate names are looked up in `Main` and `Base`, as in ACPM; the macros
resolve them, and the right side, in the calling module.

`Base.:~` with one argument is not defined for `Sym`, so `~x` outside a pattern is an
error rather than a silent bitwise `not`.

## Patterns are `Sym` values

A slot is a private function node, `_slot(x)` (like gammasimp's `_Γ`), with its flavor
recorded in the node. Predicates are kept in the pattern object, keyed by slot name. A
pattern is built with SymJ's ordinary arithmetic, so it is in the same canonical form as
the expressions it is matched against. That is the point of the design, but it changes
patterns, and the documented meaning of a pattern is its canonical form:

| Written | Canonical pattern | Consequence |
|:--|:--|:--|
| `~x + ~x` | `2*~x` | matches `2a`, with `x => a` |
| `~x * ~x` | `~x^2` | matches `a^2` |
| `(~x * ~y)^2` | `~x^2 * ~y^2` | matches `a^2 b^2` (also the subject's form of `(ab)^2`) |
| `exp(~x)` | `ℯ^~x` | matches the factor `ℯ^(2x)` inside `x·e^(2x)` |
| `sin(-~x)` | `-sin(~x)` | matches `sin(-a)` with `x => a`, not `x => -a` |
| `exp(log(~x))` | `~x` | matches everything: almost surely not what was meant |
| `~x - ~x` | `0` | the slot is gone |

The constructors flag the last two kinds. A pattern whose slot disappears is an error.
A pattern whose top-level head changed under canonicalization (from `exp` to a bare slot,
say) gives a warning.

**Open:** the warning, or an error, for a changed head.

## What a match means

Matching is structural on canonical forms. Two expressions that are equal but stored
differently (`(a + b)^2` and `a^2 + 2ab + b^2`) do not match each other's patterns;
rewriting to a normal form first is the caller's choice.

Each kind is matched as follows:

- **Leaves and literals.** A pattern part without slots matches an equal node, by hash
  and `_eq`.
- **Function applications.** The name and the number of arguments must agree, and the
  arguments are matched in order. A sequence slot may take a run of arguments:
  `f(~x, ~~rest)`.
- **Powers.** The base and exponent are matched, with `~!n` defaulting to 1, so
  `~x^~!n` also matches a plain `x`.
- **Sums and products.** An `ADD` is matched as the multiset of its terms `cᵢ·tᵢ` plus
  its constant, when that is not 0. A `MUL` is matched as the multiset of its factors
  `bᵢ^eᵢ` plus its coefficient, when that is not 1.
  - A slot `~x` takes one element, sequence slots take what is left, and order does not
    matter.
  - A subject that is not a sum is a one-element sum for this purpose, and likewise for
    products. So `~a*~x + ~!b` matches `sin(t)` with `a => 1`, `x => sin(t)`, `b => 0`.

Associativity comes from sequence slots, not from `~x`. `~x + ~y` against `a + b + c`
has no match; `~x + ~~y` has three. ACPM's `_eachmatch` instead lets every `~x` take a
group, which is where its cost goes. This is the semantics of SymbolicUtils and of ACPM's
commutative matcher. **Open:** an `associative = true` keyword that treats every `~x` in
a sum or product as `~~~x`, for the cases that want MatchPy's meaning.

**Numbers.** A number in a pattern matches only an equal number, as in Mathematica, with
one exception. A pattern term `c·P`, where `P` holds a factor slot `~x` directly, matches
a subject term `d·t` by binding `~x` to `(d/c)` times its share of `t`. Without this,
`~x - ~y` would not match `a - 2b`, because the term `-~y` has coefficient −1 and `-2b`
has −2. With it, `~x - ~y` matches any two-term sum, with `y => -b` against `a + b`.
`-sin(~x)` is unaffected, since its slot is inside `sin`: it matches only terms with
coefficient −1. **Open:** this exception (SymPy's `Wild` behaves this way, Mathematica
does not).

**Matches.** The bindings are a `Dict{Symbol,Sym}`. Matches come in a deterministic
order, following the canonical order of the subject's elements.

## Algorithm

A pattern is compiled once into a tree of matcher nodes, one per pattern node:

- a literal, compared by hash;
- a slot, which binds or checks its name;
- a function application, which checks the name and arity, then matches the arguments
  in order;
- a power, which matches the base and exponent;
- a sum or product, matched by the AC procedure below.

Matching is backtracking in continuation-passing style. The bindings live on a stack
that is cut back on failure, so a failed branch allocates nothing. `match` stops at the
first success and `eachmatch` continues.

**The AC procedure.** For a sum or product, the pattern's elements are matched in this
order, re-sorted as slots get bound:

1. Elements with no slots: found in the subject by hash, with no choice to make. A
   missing one fails the whole match at once.
2. Elements whose slots are all bound: these are literals now, handled as in step 1.
3. Elements with a fixed head (`sin(~x)`, `ℯ^~u`, `~x^2`): the candidates are only the
   subject elements with that head, tried in turn.
4. Single slots with predicates, then single slots without.
5. One sequence slot takes the remainder, with no choice to make.
6. Several sequence slots: the remaining elements are split among them. This is the
   only place the search is exponential in the subject's size, and patterns rarely need
   it.

The used subject elements are a bit mask (a `UInt64` up to 64 elements, a `BitVector`
beyond). An element is an index into the subject's vectors; it becomes a `Sym` only when
a slot binds it.

**Compilation.** The first stage builds the matcher tree at run time (the analogue of
ACPM's `cc`). `@pattern` and `@rule` build it when the module loads. Generating
specialized Julia code from the macros is a later stage, worth doing only if the
benchmark shows the tree's dispatch matters.

**Rule sets.** `rewrite` with many rules indexes them by the head of the left side (kind
and function name), so a node is tried only against rules that could match. A full
discrimination net waits until a large rule set, such as an integration table, needs
it.

## Rewriting

A rule's right side is evaluated with the bindings, and the result goes back through
the canonical constructors (`_rebuild_node`), so a rewritten expression is canonical
again.

A rule whose left side is a sum or product rewrites part of a larger one: `sin(~x)^2 +
cos(~x)^2 => 1` applied to `sin(t)^2 + cos(t)^2 + a` gives `1 + a`. The rule gets an
implicit trailing `~~rest`, and the result is the right side plus (or times) the rest.
**Open:** whether this is on by default. It is what trig identities want, but it means
`@rule ~x + ~y => …` matches every sum of two or more terms.

`SymJ.rewrite` walks bottom-up (`order = :post`) or top-down (`:pre`). With `fixpoint =
true` it repeats until nothing changes, up to a step limit, and gives a warning when the
limit is hit rather than looping.

## Relation to the rest of SymJ

- **Existing code.** SymJ's ad-hoc matchers (`_monomial_power`, `_integrate_linear_fn`,
  `_ode_linear_coefs`) stay as they are. They can be restated as rules later if that is
  as fast and clearer.
- **First users.** The first real users are SymJ's own rewrite tables: trig identities
  inside `simplify`'s candidates, and a rule-based integration layer if one comes.
- **TermInterface.** SymJ keeps its TermInterface methods, so external rewriters
  (SymbolicUtils, ACPM) keep working on SymJ terms.

## Testing and benchmarks

- **Tests.** Port the cases from ACPM's `matchpy-test.jl` and `commutative-test.jl` whose
  meaning agrees with the semantics above. The canonical-form cases in the table above
  get their own tests.
- **Benchmark.** Add a SymJ column to ACPM's 62-case benchmark (`benchmark/matching`),
  matching against SymJ terms.
- **Targets.** A first-match median at or below `cc`'s 0.45 μs, and the nested cases
  (628 μs for `ac`) within 10 μs. Cases that need MatchPy's associative meaning are
  compared only under `associative = true`.

## Stages

1. **Matcher.** Slots of all four flavors, predicates, defaults, the AC procedure,
   `SymJ.pattern`, `match`, `eachmatch` and `occursin`, plus the tests and the benchmark
   column.
2. **Rules.** `@pattern`, `@rule`, `SymJ.rewrite` and `replace`, with sub-sum and
   sub-product rewriting.
3. **Speed.** Indexing rule sets, and code generation from the macros if the benchmark
   calls for it.
4. **Use.** Trig rules among `simplify`'s candidates, measured against the current
   `simplify` for speed and result size.

## Open questions

1. **Changed heads:** a warning or an error when canonicalization changes a pattern's
   head (`exp(log(~x))`).
2. **Coefficients:** whether a numeric coefficient divides into a direct factor slot,
   so `~x - ~y` matches `a - 2b`.
3. **Associativity:** whether to offer `associative = true` (MatchPy's meaning for `~x`).
4. **Partial rewriting:** whether rules match part of a sum or product by default.
5. **Head slots:** whether to match function heads (`(~f)(~x)`, binding a name), as ACPM
   does. It is not needed for stage 1.
