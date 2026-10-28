# SymJ

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jverzani.github.io/SymJ.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jverzani.github.io/SymJ.jl/dev/)
[![Build Status](https://github.com/jverzani/SymJ.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jverzani/SymJ.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jverzani/SymJ.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jverzani/SymJ.jl)
[![CAS comparison, verified of 311: SymJ 309, Giac 305, Symbolics 296, SymPy 296](https://img.shields.io/badge/CAS%20comparison%2C%20verified%20of%20311-SymJ%20309%20%C2%B7%20Giac%20305%20%C2%B7%20Symbolics%20296%20%C2%B7%20SymPy%20296-2f8a4c)](https://rawcdn.githack.com/jverzani/SymJ.jl/d9f9703f64bef58d1aa11a55ba041a7611429fa2/bench/cas/results/cas_comparison.html)

A symbolic manipulation kernel for Julia: one concrete, hash-carrying node type, exact
numbers, canonical sums/products/powers built as you go, callable expressions, and the
TermInterface v2 contract. See [the design notes](docs/src/design.md).

The CAS badge compares SymJ (with GruntzLimit for limits and RischIntegrate for integrals),
Symbolics (with the same two), SymPy and Giac on 311 problems: expansion, factoring, gcds,
rational functions, partial fractions, derivatives, equations, limits and integrals. It
counts the answers that check out numerically. The [comparison](https://rawcdn.githack.com/jverzani/SymJ.jl/d9f9703f64bef58d1aa11a55ba041a7611429fa2/bench/cas/results/cas_comparison.html) has the numbers by
area, the times, and each system's answer to each problem; see `bench/cas`.

```julia
using SymJ
@syms x y f()

(x + y)^2 - (x + y)^2            # 0
sqrt(Sym(8))                     # 2*sqrt(2)
sin(Sym(π)/3)                    # sqrt(3)/2
expand((x + 1)^3)                # x^3 + 3*x^2 + 3*x + 1
diff(sin(x)^2, x)                # 2*cos(x)*sin(x)
(x^2 + y)(x => 2)                # y + 4
u = x^2 - 2
u(3), u(1.5), u(big(2.0))        # (7, 0.25, 2.0) -- callable like a function
f = lambdify(sin(x) * y, x, y)   # compiled: as fast as hand-written Julia
f(1.0, 2.0)
f(x, log(Sym(2), x))             # f(x, log(2, x))
cancel((x^2 - 1)/(x^2 + 2x + 1))  # (x - 1)/(x + 1)
solve(x^2 - 2 ~ 0, x)            # [-sqrt(2), sqrt(2)]
```

In a notebook that renders LaTeX (Jupyter, Pluto), `using Latexify` as well shows the
results as typeset mathematics.
