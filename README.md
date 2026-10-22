# SymJ

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jverzani.github.io/SymJ.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jverzani.github.io/SymJ.jl/dev/)
[![Build Status](https://github.com/jverzani/SymJ.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jverzani/SymJ.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jverzani/SymJ.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jverzani/SymJ.jl)

A symbolic manipulation kernel for Julia: one concrete, hash-carrying node type, exact
numbers, canonical sums/products/powers built as you go, callable expressions, and the
TermInterface v2 contract. See [the design notes](docs/src/design.md).

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
