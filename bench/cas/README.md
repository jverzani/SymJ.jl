# Computer algebra comparison

SymJ, with GruntzLimit for limits and RischIntegrate for integrals, against SymPy
(through SymPyPythonCall) and Giac (through Giac.jl), on the everyday tasks of a computer
algebra system. The result is a page, `results/cas_comparison.html`, with a summary
(counts and times) for each system, the same by area, and a grid of the problems of each
area; selecting a problem shows each system's answer.

| Area | Problems | Source |
|---|---|---|
| Expansion | 12 | `problems_basic.jl` |
| Factorization | 13 | `problems_basic.jl` |
| GCD | 8 | `problems_basic.jl` |
| Rational functions (`cancel`) | 8 | `problems_basic.jl` |
| Partial fractions | 6 | `problems_basic.jl` |
| Differentiation | 9 | `problems_basic.jl` |
| Equations (one unknown, linear systems) | 14 | `problems_basic.jl` |
| Limits | 66 | `problems_limits.jl`, from GruntzLimit's `benchmark/limits` |
| Integrals | 175 | `problems_integrals.jl`, from RischIntegrate's `benchmark/apostol` |

## Running it

SymJ, GruntzLimit and RischIntegrate come from the checkouts next to this one
(`../GruntzLimit`, `../RischIntegrate` beside the SymJ folder); everything else from the
registry. SymPy is installed by CondaPkg on the first `instantiate`.

```sh
cd bench/cas
julia --project=env_symj -e 'using Pkg; Pkg.instantiate()'    # once, or after changing a package
julia --project=env_sympy -e 'using Pkg; Pkg.instantiate()'
julia --project=env_giac -e 'using Pkg; Pkg.instantiate()'
julia run_all.jl            # or `julia run_all.jl 120` for a 120 s limit per problem
```

`run_all.jl` runs the systems one after another, so they do not compete for the CPU, and
then builds the page. Open `results/cas_comparison.html` in a browser; it is
self-contained. It was made for Julia 1.13.

To run or rerun one system, then rebuild the page:

```sh
julia supervise.jl symj env_symj 60
julia build.jl
```

To try a worker on some problems only: `julia --project=env_symj worker_symj.jl START OUTFILE STOP`.

## How it works

| File | Role |
|---|---|
| `problems_basic.jl`, `problems_limits.jl`, `problems_integrals.jl` | the problems, as quoted expressions; `problems.jl` puts them in one list of `(area, op, args, expect, source)` |
| `worker_symj.jl`, `worker_sympy.jl`, `worker_giac.jl` | one per system: loads the packages, makes the variables, and defines the operations (`EXPAND`, `FACTOR`, `GCD`, `CANCEL`, `APART`, `DIFF`, `SOLVE`, `LINSOLVE`, `LIMIT`, `INTEGRATE`) and what the checks need (`EVAL`, `NTERMS`, `DENDEG`, `UNSOLVED`, `CLEARCACHE`) |
| `core.jl` | shared by the workers: builds each problem with the worker's variables, runs and times it, checks the answer, writes a line per problem |
| `supervise.jl` | runs a worker; if one problem takes longer than the limit, kills the process, records a timeout and restarts at the next problem (from RischIntegrate's benchmark) |
| `build.jl`, `template.html` | turn the results into the page |
| `env_symj`, `env_sympy`, `env_giac` | the environments |

Each problem is built from its expression with the worker's own variables. Numeric
constants in the problems (`sqrt(2)`, `log(6)`) are made exact numbers of the system
first, so that Julia does not compute them in floating point.

What is recorded, for each problem and system:

- **Verified**: the answer's value, computed by the system, matches the problem evaluated
  in plain Julia in BigFloat (see `check` in `core.jl` for each operation), and the form
  is what was asked (the number of terms of an expansion or partial fraction
  decomposition, the number of irreducible factors, the degree of a reduced
  denominator, the number of real solutions).
- **Correct, incomplete**: the value is right but the form is not.
- **Answer, not checked**: the answer could not be evaluated numerically.
- **Wrong**: the value does not match.
- **Not done**: the problem came back unevaluated (`Integral(…)`, `limit(…)`), or no
  limit (`zoo`, `undef`).
- **Error**, **Timed out or crashed**.

Times: the first call in a fresh process (compilation included, for Julia packages), and a
repeated call (after clearing SymPy's cache) when the first took under 15 s.

The results file of each system, `results/results_NAME.tsv`, has a line per event:
`L` load time, `V` package versions, `S` problem started, `R` problem, status,
verification, first and repeat times, answer; `D` done, `T` wall time and limit.
