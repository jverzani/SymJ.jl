# SymJ in RischIntegrate's Apostol benchmark

`worker_symj.jl` runs RischIntegrate on SymJ terms in the harness of
`RischIntegrate/benchmark/apostol` (175 problems from `test/test_apostol.jl`). It is
kept here until RischIntegrate is on GitHub; then it can move there.

## Running it

From a copy of `RischIntegrate/benchmark/apostol` (or the checkout itself):

1. Copy `worker_symj.jl` into `benchmark/apostol/`.
2. Make an environment with RischIntegrate and SymJ developed from their checkouts,
   and SpecialFunctions and PolyLog (as in `env_risch`; they load SymJ's extensions):

   ```sh
   julia --project=env_symj -e 'using Pkg
       Pkg.develop([PackageSpec(path = "path/to/RischIntegrate"), PackageSpec(path = "path/to/SymJ")])
       Pkg.add(["SpecialFunctions", "PolyLog"])'
   ```

3. Run it with the same per-problem limit as the report, then summarize:

   ```sh
   julia supervise.jl symj env_symj 120      # results/results_symj.tsv
   ```

   To include it in `results/apostol_comparison.html`, add
   `("symj", "RischIntegrate", "Risch", "SymJ")` to `METHODS` in `build.jl` and run
   `julia build.jl`.

## Result (2026-10-22)

RischIntegrate at 6b5dd5c (with an uncommitted `constant_fields.jl` change), SymJ at
the commit adding this folder, Julia 1.13, 120 s limit: 172 verified, 1 wrong, 2
errors; integration 47.5 s, wall 50.9 s, load 0.6 s.

- #19 and #172 fail on every term type (RischIntegrate: not implemented / not
  elementary).
- #55, `sqrt(1 + 3cos(x)^2) * sin(2x)`, is wrong on SymJ only: the u-substitution
  builds the constant `-(2/3) u^2` as a term, SymJ simplifies it to
  `-2cos(x)^2 - 2/3` (numbers distribute over sums), and `to_x` then no longer finds
  `u` by its printed form. Fix in RischIntegrate: map `U^k` to the variable inside
  `cterm` instead of recognizing `u^k` afterwards.
- Answers with constants such as `sqrt(pi)` come back as floats on SymJ:
  RischIntegrate's `term_output.jl` builds them with the `u(value, :)` placeholder
  convention, which SymJ does not support; `maketerm` (as in GruntzLimit) would keep
  them exact.
