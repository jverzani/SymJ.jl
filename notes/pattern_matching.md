# Design note: pattern matching and rewrite rules

A proposal, written 2026-10-25, for pattern matching in SymJ: patterns with wildcards,
matching against expressions, and rules that rewrite with the matches. Nothing here is
implemented yet. The questions the first draft left open were settled the same day;
the answers are folded in, and listed at the end.

## Decision: inside SymJ, with the SymbolicUtils syntax

The matcher lives in SymJ and works on `Sym` nodes directly, not through TermInterface.
The surface syntax is the one AssociativeCommutativePatternMatching (ACPM) borrowed from
SymbolicUtils: expressions with `~x` slots make patterns, and macros compile them.

The reason is cost. Through TermInterface a matcher sees only `operation` and
`arguments`, and must assume nothing about them: sums may be nested, unsorted and
uncombined, so associative-commutative (AC) matching searches orderings of arguments and
splits of them into groups. That search is what makes ACPM's `_eachmatch` slow. The
meaning of a match stays MatchPy's (see below); what SymJ removes is the work spent on
a subject that is not in canonical form. In ACPM's 62-case benchmark (`benchmark/matching`, run 2026-10-18), the first match takes:

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
| `~x` | `~x::pred` | one argument of a function; one or more terms of a sum or factors of a product |
| `~!x` | `~!x::pred` | as `~x`, or a default when absent: 0 in a sum, 1 in a product or exponent |
| `~~x` | `~~x::pred` | zero or more arguments |
| `~~~x` | `~~~x::pred` | one or more arguments |
| `(~f)(…)` | `(~f::pred)(…)` | the head of a function application |

`+` and `*` are associative, so inside a sum or product `~x` is a sequence slot: it takes
one or more terms (factors) and is bound to their sum (product). There `~x` and `~~~x`
mean the same, and `~~x` and `~!x` also allow none. In a function's arguments, `~x` takes
exactly one, and the sequence slots bind a vector of arguments.

A name may appear more than once (a non-linear pattern) but with one flavor. A predicate
gets the bound value and must return `true`, or the match is rejected; an error counts as
`false`.

The functions are unexported and called qualified, or they extend existing Base functions:

```julia
p = SymJ.@pattern ~a * sin(~x) + ~~r        # compiled where the macro is expanded
p = SymJ.pattern(:(~a * sin(~x) + ~~r))     # the same, built at run time from an Expr

match(p, ex)        # Dict{Symbol,Any} of the first match, or nothing
eachmatch(p, ex)    # an iterator of all the matches
occursin(p, ex)     # a match anywhere in ex

r = SymJ.@rule sin(~x)^2 + cos(~x)^2 + ~!r => 1 + ~r
r = SymJ.@rule log(~x::ispos) + log(~y::ispos) => log(~x * ~y)
r = SymJ.@rule (~f)(-~x) => (~f)(~x)       # a head slot, bound to the function
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
| `exp(log(~x))` | `~x` | matches everything, as `~x` does |
| `~x - ~x` | `0` | the slot is gone |

A pattern is taken as it canonicalizes, with no warning: `exp(log(~x))` is the pattern
`~x`. The one error is a slot that disappears (`~x - ~x`), since a rule's right side
could then refer to a name the match never binds.

## What a match means

Matching is structural on canonical forms. Two expressions that are equal but stored
differently (`(a + b)^2` and `a^2 + 2ab + b^2`) do not match each other's patterns;
rewriting to a normal form first is the caller's choice.

Each kind is matched as follows:

- **Leaves and literals.** A pattern part without slots matches an equal node, by hash
  and `_eq`.
- **Function applications.** The name and the number of arguments must agree, and the
  arguments are matched in order. A sequence slot may take a run of arguments:
  `f(~x, ~~rest)`. A head slot `(~f)(…)` matches any function application, so `(~f)(~x)`
  matches `sin(t)` with `f => sin`, `x => t`. It is bound to the function (a SymJ
  function such as `sin`, or a `SymFunction`), so a right side can call it. Sums,
  products and powers have their own rules below and are not matched by head slots.
- **Powers.** The base and exponent are matched, with `~!n` defaulting to 1, so
  `~x^~!n` also matches a plain `x`.
- **Sums and products.** An `ADD` is matched as the multiset of its terms `cᵢ·tᵢ` plus
  its constant, when that is not 0. A `MUL` is matched as the multiset of its factors
  `bᵢ^eᵢ` plus its coefficient, when that is not 1.
  - Each pattern element takes a group of subject elements, as a whole: a slot takes one
    or more (or none, for `~~x` and `~!x`), anything else takes one. Order does not
    matter.
  - A subject that is not a sum is a one-element sum for this purpose, and likewise for
    products. So `~a*~x + ~!b` matches `sin(t)` with `a => 1`, `x => sin(t)`, `b => 0`.

This is MatchPy's meaning, and so ACPM's `_eachmatch`: `~x + ~y` against `a + b + c` has
six matches, one for each way to split the three terms into two non-empty groups, among
them `x => a + b, y => c`. A pattern that is a sum therefore also matches part of a
larger sum without any special rule: `sin(~x)^2 + cos(~x)^2 + ~!r` matches
`sin(t)^2 + cos(t)^2 + a` with `r => a`, and the identity alone with `r => 0`.

**Numbers.** A number in a pattern matches an equal number, and a numeric coefficient
also divides into a slot, as SymPy's `Wild` does. A pattern term `c·P`, where `P` holds a
slot directly as a factor, matches a subject term `d·t` with that slot taking `d/c`
along with its share of `t`. So `~x - ~y` matches `a - 2b` with `y => 2b`, and `a + b`
with `y => -b` (or `-a`). `-sin(~x)` has its slot inside `sin`, so it matches only terms
with coefficient −1. In a product, `2*~x` against `6a` binds `x => 3a`.

**Matches.** The bindings are a `Dict{Symbol,Any}`: a `Sym` for an ordinary slot, a
vector of `Sym`s for a sequence slot among a function's arguments, a function for a head
slot. Matches come in a deterministic order, following the canonical order of the
subject's elements.

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
4. Slots with predicates, then slots without.
5. The last slot takes the remainder, with no choice to make.

Because `~x` takes groups in a sum or product, step 4 splits the remaining elements among
the free slots, and the number of splits grows exponentially with the subject's size.
That is the cost MatchPy pays, and it can't be removed, but here it is confined:

- Steps 1–3 run first and take no groups. Most patterns have only one or two free slots
  left by step 4.
- With one free slot there is nothing to split. With several, the first split tried
  (each slot one element, the last slot the rest) is a match unless a predicate or a
  repeated slot rejects it. So `match` stops at once in the common case, and the search
  only grows when matches are rejected.
- `eachmatch` costs in proportion to the matches it returns, which are themselves
  exponential in number (`~x + ~y` against `n` terms has `2ⁿ - 2`).
- A repeated slot that is already bound needs a group equal to its value. Since the
  value is canonical, that group is found from its elements by hash, not searched for.

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

There is no implicit rest. A rule that should rewrite part of a larger sum says so with
a slot, as in `sin(~x)^2 + cos(~x)^2 + ~!r => 1 + ~r`; associativity does the rest.
Without the slot, the rule matches only the whole identity.

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

- **Tests.** Port the cases from ACPM's `matchpy-test.jl`, whose meaning is the one
  above, and the cases from `commutative-test.jl` that don't depend on `~x` taking one
  term. The canonical-form cases in the table above, and the coefficient division, get
  their own tests. A match differs from ACPM's only where SymJ's canonical form stores
  the subject differently.
- **Benchmark.** Add a SymJ column to ACPM's 62-case benchmark (`benchmark/matching`),
  matching against SymJ terms. `ac` is the like-for-like comparison, since it has the
  same meaning.
- **Targets.** A first-match median at or below `cc`'s 0.45 μs, with the AC meaning that
  `cc` lacks. The nested cases (628 μs for `ac`) within 10 μs. Every match found by
  `eachmatch` costing well under `ac`'s time per match.

## Stages

1. **Matcher.** Slots of all four flavors and head slots, predicates, defaults,
   coefficient division, the AC procedure, `SymJ.pattern`, `match`, `eachmatch` and
   `occursin`, plus the tests and the benchmark column.
2. **Rules.** `@pattern`, `@rule`, `SymJ.rewrite` and `replace`.
3. **Speed.** Indexing rule sets, and code generation from the macros if the benchmark
   calls for it.
4. **Use.** Trig rules among `simplify`'s candidates, measured against the current
   `simplify` for speed and result size.

## Decisions

The first draft's open questions, as settled:

1. **Changed patterns:** a pattern is what it canonicalizes to. `exp(log(~x))` is `~x`,
   with no warning.
2. **Coefficients:** a numeric coefficient divides into a slot, as SymPy does, so
   `~x - ~y` matches `a - 2b`.
3. **Associativity:** `+` and `*` are associative, and in them `~x` takes one or more
   terms or factors, as a sequence slot does (MatchPy's meaning).
4. **Partial rewriting:** nothing implicit. Associativity already lets `~x + ~y` match
   with `x => sin(t)^2 + cos(t)^2, y => a`, and a rule for part of a sum carries its own
   slot for the rest.
5. **Head slots:** patterns match function heads, `(~f)(~x)`, binding the function.
