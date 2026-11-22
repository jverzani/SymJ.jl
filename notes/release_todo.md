# TODO before release v0.1.0

## Docs

- [x] Rework the intro and calculus sections (tutorial split into four pages, #206)
- [x] Attribution (Claude) on the design notes (design index, release-prep PR)
- [ ] Review the feature audit (`notes/feature_audit.md`, last updated 2026-11-19) so it
  is up to date, and organize it as a roadmap for future development: done, deferred by
  decision (Gruntz port, `==` on floats, parameters hidden, symbolic exponents step 2),
  open
- [ ] Organize the SymPy issue audit notes (`notes/sympy_review/`): a summary per round,
  linked from the roadmap
- [ ] README: status, installation, docs badge and links current
- [ ] Docs deployment: GitHub Pages enabled and the `DOCUMENTER_KEY` secret set, so
  `deploydocs` publishes (a local build skips deployment)

## Interface

- [ ] Review the exported functions (`src/SymJ.jl`)
  - [x] internal-looking names unexported: `AbstractSymbolic`, `SymNumber`,
    `SymFunction`, `Assumptions`, `evaluate`, `lambdify_expr`, `register_function!`,
    `numer_denom`, `factor_squarefree`, `squarefree_part`, `divexact`, `nterms`,
    `GruntzError`
  - [x] `SymJ.numer_denom`, `SymJ.factor_squarefree` in the tutorial pages
  - [ ] names that clash with companion packages (`Poly`, `degree`, `coefficients`,
    `gens`, `divides`, `Sym`, `@syms`, `subs`, `Eq`, `N`, `solve`), the one-letter `N`,
    and the operators `⩵`, `≷`

## Dependencies

- [x] Julia ≥ 1.12 (RischIntegrate main needs 1.12; CI's 1.10 job becomes 1.12)
- [x] Remove the `VERSION >= v"1.12"` test guards and the RischIntegrate 1.0 feature
  checks; RischIntegrate compat "1.5"
- [x] Bump RischIntegrate to 1.5.0 (with #184–#186: answers valid for every real
  parameter value); SymJ's #203 gives a wrong case with 1.4.0 (√(a+bx)√(a-bx), a < 0)
- [x] With that release, add the waiting tests: ∫₀¹ (√x - 1)/log x = log(3/2)
  (RischIntegrate #177) and SymPy #4231 (RischIntegrate #157)

## Release

- [x] Merge or close the open PRs (#203, #204, #205, #206)
- [ ] Merge the release-prep PR
- [ ] Register v0.1.0 in General (`Project.toml` is at 0.1.0; TagBot tags the release)
