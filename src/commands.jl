# The functions that are part of SymJ's interface but not exported: called qualified
# (`SymJ.dsolve`), or brought into scope together with `using SymJ.AllCommands`. One list
# feeds both the `public` declarations and the submodule's exports. Parameters
# (`SymJ.parameter`, …) are left out while their shared-name question is open, as are
# internals that happen to have docstrings (`add`, `kind`, the node kinds, `ONE`).
const _COMMANDS = (
    # solving
    :solveset, :fsolve, :rsolve, :dsolve, :infinitesimals, :linear_eq_to_matrix,
    :real_roots, :root_intervals,
    # calculus
    :summation, :harmonic, :residue, :Integral, :laplace_transform,
    :inverse_laplace_transform, :piecewise, :piecewise_fold, :heaviside, :dirac,
    :fresnels, :fresnelc, :gradient, :jacobian, :hessian, :divergence, :curl, :laplacian,
    # numbers, simplification and rewriting
    :evalf, :nsimplify, :unwrap_const, :rewrite, :expand_trig, :expand_log, :together,
    :horner, :gather, :cse, :decompose, :radsimp, :simplify_trig, :simplify_hyperbolic,
    :simplify_log, :simplify_powers, :simplify_radicals, :simplify_rational, :simplify_gamma,
    # expressions and polynomials
    :isconstant, :coeff, :numer_denom, :lambdify_expr, :minimal_polynomial,
    :factor_squarefree, :squarefree_part, :nterms, :divexact,
    # linear algebra
    :rref, :charpoly, :jordan_form, :permanent, :rotation_matrix, :hermite_normal_form,
    :hermite_normal_decomp, :smith_normal_form, :smith_normal_decomp,
    # patterns and rules
    Symbol("@rule"), Symbol("@pattern"), :rule, :pattern, :register_function!,
    # constants and the cases of an answer
    :PI, :E, :IM, :Cases,
)

eval(Expr(:public, _COMMANDS...))

"""
    SymJ.AllCommands

A module that exports SymJ's unexported functions, for those who would rather not
qualify them: after `using SymJ.AllCommands`, `dsolve`, `summation`, `rref`,
`simplify_trig` and the rest are in scope, as are the constants `PI`, `E` and `IM`. Plain
`using SymJ` keeps the smaller set of names; the others are `public` and can always be
written `SymJ.dsolve` and so on.

A name exported by another package loaded the same way (`gradient` in ForwardDiff, say)
is an error only where it is used unqualified; qualify it there.

```julia
using SymJ, SymJ.AllCommands
@syms x y()
dsolve(y'(x) ~ y(x), y(x))
```
"""
module AllCommands
import ..SymJ
eval(Expr(:using, Expr(:(:), Expr(:., :., :., :SymJ), (Expr(:., n) for n in SymJ._COMMANDS)...)))
eval(Expr(:export, SymJ._COMMANDS...))
end
