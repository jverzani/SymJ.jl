# SymJ

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jverzani.github.io/SymJ.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jverzani.github.io/SymJ.jl/dev/)
[![Build Status](https://github.com/jverzani/SymJ.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jverzani/SymJ.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jverzani/SymJ.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jverzani/SymJ.jl)
[![CAS comparison, verified of 311: SymJ 311, Giac 305, Symbolics 297, SymPy 296](https://img.shields.io/badge/CAS%20comparison%2C%20verified%20of%20311-SymJ%20311%20%C2%B7%20Giac%20305%20%C2%B7%20Symbolics%20297%20%C2%B7%20SymPy%20296-2f8a4c)](https://rawcdn.githack.com/jverzani/SymJ.jl/f484f5df4e286ed814598c634bc8189dac523d52/bench/cas/results/cas_comparison.html)

SymJ is a Claude-programmed computer algebra system for Julia, along the lines of [SymPy](https://www.sympy.org/en/index.html) and [Giac](https://www-fourier.univ-grenoble-alpes.fr/~parisse/giac.html), both of which have Julia interfaces. Unlike SymPy and Giac, which Julia reaches through wrappers, SymJ is written entirely in Julia.

A few ideas run through the package:

* Variables are made with the `@syms` macro or related functions, and generic Julia code (arithmetic operations and mathematical functions) builds symbolic expressions from them.

* Expressions are kept in a canonical form as they are built, so like terms are combined.

* Numbers stay exact unless a floating point number is involved.

* Expressions are callable: an expression in one variable evaluates at a number to a Julia number (`u(3)`); pairs substitute (`(x^2 + y)(x => 2)` is `y + 4`).


```julia
using SymJ
@syms x y f()

(x + y)^2 - (x + y)^2             # 0
sqrt(Sym(8))                      # 2*sqrt(2)
sin(Sym(π)/3)                     # sqrt(3)/2
expand((x + 1)^3)                 # x^3 + 3*x^2 + 3*x + 1
diff(sin(x)^2, x)                 # 2*cos(x)*sin(x)
(x^2 + y)(x => 2)                 # y + 4
u = x^2 - 2
u(3), u(1.5), u(big(2.0))         # (7, 0.25, 2.0) -- callable like a function
f = lambdify(sin(x) * y, x, y)    # compiled: as fast as hand-written Julia
f(1.0, 2.0)
f(x, log(Sym(2), x))              # f(x, log(2, x))
cancel((x^2 - 1)/(x^2 + 2x + 1))  # (x - 1)/(x + 1)
solve(x^2 - 2 ~ 0, x)             # [-sqrt(2), sqrt(2)]
```

The CAS badge compares SymJ (with GruntzLimit for limits and
RischIntegrate for integrals), Symbolics (with the same two), SymPy
and Giac on 311 problems: expansion, factoring, gcds, rational
functions, partial fractions, derivatives, equations, limits and
integrals. It counts the answers that check out numerically. The
[comparison](https://rawcdn.githack.com/jverzani/SymJ.jl/f484f5df4e286ed814598c634bc8189dac523d52/bench/cas/results/cas_comparison.html)
has the numbers by area, the times, and each system's answer to each
problem; see `bench/cas`.

The development was greatly aided by a review of the large issue tracker of the SymPy project.
