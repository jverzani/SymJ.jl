# Design note: pattern matching and rewrite rules

A proposal, written 2026-10-25, for pattern matching in SymJ: patterns with wildcards,
matching against expressions, and rules that rewrite with the matches. The questions
the first draft left open were settled the same day; the answers are folded in, and
listed at the end. Stages 1–3, the matcher, rules and speed, are implemented in
`src/patterns.jl`, and stage 4, rules in `simplify`, in `src/simprules.jl`; see the
"as built" sections at the end.

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
| `~!x` | `~!x::pred` | one argument; any number of terms or factors, none giving the default 0 (1); in an exponent, 1 when there is no power |
| `~~x` | `~~x::pred` | zero or more arguments |
| `~~~x` | `~~~x::pred` | one or more arguments |
| `(~f)(…)` | `(~f::pred)(…)` | the head of a function application |

`+` and `*` are associative, so inside a sum or product `~x` is a sequence slot: it takes
one or more terms (factors) and is bound to their sum (product). There `~x` and `~~~x`
mean the same. `~!x` is extended the same way: it takes none or one or more, and is
bound to their sum (product), which is the default 0 (1) when it takes none. There `~!x`
and `~~x` mean the same. In a function's arguments, `~x` and `~!x` take exactly one (a
function has no default to give), and the sequence slots bind a vector of arguments. In
an exponent, `~!n` is 1 when there is no power.

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

r = SymJ.@rule sin(~x)^2 + cos(~x)^2 + ~~r => 1 + ~r
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
  - Each pattern element takes a group of subject elements, as a whole: `~x` and `~~~x`
    take one or more, `~~x` and `~!x` take any number (none gives the default 0 or 1),
    and anything else takes one. Order does not matter.
  - A subject that is not a sum is a one-element sum for this purpose, and likewise for
    products. So `~!a*~x + ~!b` matches `sin(t)` with `a => 1`, `x => sin(t)`, `b => 0`.

This is MatchPy's meaning, and so ACPM's `_eachmatch`: `~x + ~y` against `a + b + c` has
six matches, one for each way to split the three terms into two non-empty groups, among
them `x => a + b, y => c`. A pattern that is a sum therefore also matches part of a
larger sum without any special rule: `sin(~x)^2 + cos(~x)^2 + ~~r` matches
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

**Which match comes first.** In a sum or product, `~!x` and `~~x` mean the same, so a
pattern can match in several ways. `~!a*~x + ~!b` against `3y + 5` gives `a => 3,
x => y, b => 5`, but also `a => 1, x => 3y, b => 5` and `a => 1, x => 3y + 5, b => 0`.
`eachmatch` gives them all. `match` gives the one found first, and the search is ordered
so that it is the expected one:

1. **Numbers go to default slots first.** A product's coefficient, or a sum's constant,
   is offered to a `~!` slot of that product or sum before any other slot. So
   `~!a*~x + ~!b` against `3y + 5` gives `a => 3`, `x => y`, `b => 5`, and
   `~!a*~x^~!n` against `2y^3` gives `a => 2`, `x => y`, `n => 3`.
2. **Otherwise, defaults take nothing first.** Apart from numbers, a `~!` slot first
   takes no elements, so that it has its default, and takes elements only when that
   fails. So `~x + ~!r` against `a + b` gives `x => a + b`, `r => 0`.

These rules only order the search. They don't change which matches exist, and with
them the common rule shapes match without backtracking.

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
a slot, as in `sin(~x)^2 + cos(~x)^2 + ~~r => 1 + ~r`; associativity does the rest.
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
6. **First match:** numbers go to `~!` slots first, and otherwise a `~!` slot takes
   nothing first (see "Which match comes first").

## Stage 1 as built

`src/patterns.jl`, 2026-10-25: `SymJ.pattern`, `SymJ.Pattern`, and methods of `match`,
`eachmatch` and `occursin`, with 43 tests. Beyond the design above:

- **Order of slots.** Slots are tried in the order they first appear in the written
  pattern (after the defaults, scaled and predicated slots), so the number in `3x` goes
  to `~b` in `~b*~x` whatever the slots' names. The other parts of a sum or product are
  tried in the same written order.
- **Ambiguous splits** follow SymJ's stored order of terms, which is not the printed
  one: `~x + ~y` against `a + b` gives `x => b` first. Sorting into printed order cost
  3–17 μs per match (it prints each term) and was dropped.
- **Defaults in powers.** `(~x)^(~!m)` against `x^3` matches twice: `x => x, m => 3`
  first, then `x => x^3, m => 1`. ACPM gives only the first.
- **Continuations** are a linked list of concrete goal records, not closures, which
  halved the search time.
- **Predicates** get `Sym`s. SymJ defines Base's `isinteger`, `iseven` and `isodd` for
  `Sym`, true for numbers only (a symbol `n` is not known to be an integer), so
  `~n::isinteger` takes integer values.

On ACPM's 62 cases (run in a scratch script; ACPM's checkout is unchanged), the first
match takes a median 0.83 μs and at most 6.1 μs, against `ac`'s 13.1 and 628 μs, `cc`'s
0.45 and 9.2 μs, and SymbolicUtils' 0.93 μs median. The match counts agree with `ac`'s
except in six cases, all from the decisions above: `~x + ~x` matches `2a` (SymbolicUtils
stores `a + a` as `2a` too, but `ac` doesn't match it) and `a + b` (`x => a/2 + b/2`);
`~!a` may take no factor (`~!a * ~x` against `3x` also gives `a => 1, x => 3x`); and
default exponents also match a power as a whole.

## Stage 2 as built

`src/patterns.jl`, 2026-10-25: `SymJ.@pattern`, `SymJ.@rule`, `SymJ.rule`, `SymJ.Rule`,
`SymJ.rewrite` and a method of `replace`, with 21 tests.

- **`@pattern` and `@rule`** quote the left side with each predicate interpolated, so
  predicates (and any `$v` the user writes) are evaluated where the macro is called.
  Names in a pattern are symbols: `a` is `Sym(:a)`, and a value goes in as `$a`. The
  right side of `@rule` is ordinary Julia code in the caller's scope, `~x` becoming the
  bound value. The pattern is built each time the macro's code runs, so a rule in a
  function body is better made once, as a `const`.
- **Rejected matches.** A right side that gives `nothing` sends the search on to the next
  match, inside the same backtracking search, so `@rule ~x * ~y => (isinteger(~x) ?
  nothing : ~x)` against `2a` gives `a`.
- **`SymJ.rule(lhs => rhs)`** makes a rule at run time. `rhs` is a function of the
  bindings, or a quoted template read like a pattern with the slots filled in: a
  sequence slot's arguments go in place (`f(~~r, ~x)`), and a head slot calls the bound
  function. Either way a right side naming a slot the left side lacks is an error when
  the rule is made.
- **`rewrite`** walks the terms of sums, factors of products and arguments, `order =
  :post` (default) or `:pre`; at each part the first rule with a result wins, and a
  result equal to the part counts as no change. `fixpoint = true` repeats whole passes,
  at most `maxsteps` (100), with a warning at the limit. `replace(ex, rules...)` is one
  `:post` pass.

## Stage 3 as built

`src/patterns.jl`, 2026-10-26. Times are best of five, on a machine with a load average
near 18.

- **Rule sets are indexed.** `rewrite` sorts its rules once into those whose left side is
  a function application (keyed by the name), head slots (any function node) and the
  rest; at a node it tries only the rules that could match there, in their original
  order. The merged list per function name is built with the set, so a set can be
  shared, as `simplify`'s is (stage 4).
- **A quick reject** runs before any search state is allocated: the head, a function's
  arity (unless a sequence slot is among its arguments), whether a power is a power, and
  whether a sum or product has as many terms or factors as the pattern needs.
- **Smaller savings:** no sorting of the slots when none is bound; the free elements are
  listed with the number first without a sort; a slot that takes all of a sum or product
  is bound to the subject itself; a rule's right side reads the bindings in place, so no
  `Dict` is built per match (a user function given to `SymJ.rule` still gets one).
- **Effect.** A 66-rule set (64 rules keyed by function names that never match, the
  Pythagorean identity and the log-of-a-product rule) rewriting a 30-term expression:
  4.3 ms before, 0.43 ms after. On ACPM's 62 cases the first match is a median 0.78 μs
  and at most 3.6 μs (stage 1: 0.83 and 6.1 μs).
- **Not done: generating code from the macros.** The matcher is now within 2x of ACPM's
  commutative-only `cc` (0.45 μs median), with associativity, which `cc` lacks. A
  generated matcher would remove the dispatch on node kinds and the goal records, but
  the remaining cost is mostly building the subject's element lists and the bound
  values, which generated code would still do. Not worth it until a real rule set shows
  the matcher as the bottleneck.

## Stage 4 as built

`src/simprules.jl` and `_simplify` in `src/simplify.jl`, 2026-10-26.

- **Gaps first.** A probe of 44 trigonometric and hyperbolic expressions found the
  existing `simplify_trig` (a normal form modulo `sin^2 + cos^2 = 1`, with products
  linearized) already handles most identities. What it misses: power reduction read
  backwards (`(1 - cos(2x))/2`), half angles (`sin(x)/(1 + cos(x))`), forms that need
  expanding first (`sin(x + y) - cos(y)*sin(x)`), and arguments of `exp` and of roots,
  which are powers and were not simplified as function arguments are.
- **The rules** are seven: four for power reduction (trig and hyperbolic) and three for
  half angles, as one `_RuleSet` built when the package loads. They are a candidate
  together with `expand_trig` and `simplify_rational ∘ expand_trig`, all only for
  expressions with a trig or hyperbolic function; a candidate is kept only when smaller.
  `tan(x) + cot(x) = 2csc(2x)` was tried and dropped: the right side is larger.
- **Lessons for writing rules.**
  - Coefficient division (decision 2) makes `cos(2 * ~u)` match `cos(x)` with
    `u => x/2`. Bottom up, that rewrote `1 + cos(x)` inside `sin(x)/(1 + cos(x))`
    before the half-angle rule could see it, and a rule `cos(2 * ~u) => …` loops. The
    rules use `cos(~v::_halvable)`, a predicate that the factor 2 is written.
  - A one-or-more slot takes no implicit coefficient 1: `~c - ~c * cos(…)` misses
    `1 - cos(2x)`, as `~c * cos(…)` is matched (with `c` free) before `~c`. The rules
    use `~!c`, whose default in a product is 1.
  - A constant is not split: `~!c - ~!c * cos(2x) + ~!r` misses `2 + y - cos(2x)`.
- **A matcher fix.** A bound slot whose value is 1 in a product now takes no factors
  (`h(~c, ~c * a)` matches `h(1, a)`). The sum analogue, a bound 0, is still turned away
  by the quick reject, which does not see bindings.
- **Effect.** Seven of the probe's expressions simplify further; none changes otherwise.
  The probe takes 12.6 ms instead of 11.1 ms (best of five, loaded machine): about
  11 μs per expression for the rules and 7 μs for the two expansions. Expressions
  without trig functions skip all three; a 17-case set of them is unchanged in time.
