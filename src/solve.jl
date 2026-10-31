# Solving equations: CommonSolve's `solve` for a symbolic equation in one variable.
#
# The equation must be a Relation (`x^2 ~ 4`); a bare expression is not taken to mean
# "= 0". lhs - rhs is brought to lowest terms with `cancel`; the numerator's roots are
# the solutions (numerator and denominator then share no factor, so no root is a
# pole). The numerator is factored (other symbols and parameters may appear in the
# coefficients) and each irreducible factor solved:
#   degree 1 and 2: the usual formulas; a xⁿ + b: (-b/a)^(1/n) (-1)^(2k/n);
#   degree 3: Cardano; degree 4: Ferrari (biquadratic ones as quadratics in x²);
#   degree ≥ 5 with numeric coefficients: numeric roots (Aberth's iteration).

import CommonSolve: solve

"""
    solve(eq::Relation, x; assume = nothing, drop_empty = false)

The solutions of the equation `eq` (written with `~`, e.g. `x^2 ~ 4`) for `x`, as a
vector of distinct expressions. The equation must be polynomial or rational in `x`;
its coefficients may involve other symbols. Solutions are exact through
degree 4 and for equations that factor into such pieces; irreducible factors of higher
degree with numeric coefficients give numeric (floating point) roots. Numeric
solutions are sorted, real ones first. An equation with no solution gives an empty
vector. A solution must be in the equation's domain: `(x^2 - 1)/(x - 1) ~ 2` has none,
though cancelling gives `x + 1 = 2`. For a system, `solve([eqs...], [vars...])`, giving a vector of `Dict`s;
mixing the forms (`solve(eq, [x])`, `solve([eq], x)`) is an error. An unknown may also
be a function call or a derivative (`u(x)`, `diff(u(x), x)`), as in implicit
differentiation; other derivatives of it are held as they are.

An equation through `sin` or `cos` gives the solutions for its argument in `[0, 2π]`,
each once: `cos(x) ~ 0` gives `[π/2, 3π/2]`, `sin(x) ~ -1//2` gives `[7π/6, 11π/6]`, and
`cos(x) ~ 1` gives `[0, 2π]` (with the interval closed, 0 brings 2π). `tan(u) = r` gives
`atan(r)`.

```julia
solve(x^2 - 2 ~ 0, x)          # [-sqrt(2), sqrt(2)]
solve(x^2 ~ -1, x)             # [-im, im]
solve(a*x + b ~ 0, x)          # [-b/a]
solve(diff(x^2 + u(x)^2 - 1, x) ~ 0, diff(u(x), x))   # [-x/u(x)]
```

With assumptions (`assume = (x ≫ 0,)`, an [`Assumptions`](@ref), or inside
[`assuming`](@ref)), only the solutions they allow for `x` are kept: those in its
number set (`x ∈ ℝ`), its range (`x ≫ 0`), and consistent with its relations to other
symbols. Facts about the other symbols are used to decide that. When the answer
depends on something the assumptions leave open, the result is a `Cases` value: a
list of `conditions => solutions`, which can be called with values (`cs(a => 4)`).

```julia
solve(x^2 ~ 4, x; assume = (x ≫ 0,))                # [2]
solve(x^2 ~ a, x; assume = (x ∈ ℝ, a ∈ ℝ))         # a < 0 => [], a == 0 => [0], a > 0 => [-√a, √a]
solve(x^2 ~ a, x; assume = (x ≫ 0,))                # a > 0 => [√a], !(a > 0) => []
```

`drop_empty = true` drops the cases with no solutions, and when one case is left gives
its solutions alone: the conditions under which they hold are then not shown.

```julia
solve(x^2 ~ a, x; assume = (x ≫ 0,), drop_empty = true)   # [√a] (for a > 0)
```
"""
function solve(eq::Relation, x::Sym; assume = nothing, drop_empty::Bool = false)
    r = if !_isnamed(x)
        _solve_for_calls((es, vs, B) -> solve(only(es), only(vs); assume = B), [eq], [x];
                         assume = _assumptions(assume))
    else
        A = _assumptions(assume)
        A === nothing ? _solve_eq(eq, x) : _solve_assumed(eq, x, A)
    end
    drop_empty ? _drop_empty(r) : r
end

# the cases that have solutions; their solutions alone when one case is left
function _drop_empty(r)
    r isa Cases || return r
    cs = [c for c in r.cases if !(c.second isa AbstractVector && isempty(c.second))]
    isempty(cs) && return Sym[]
    length(cs) == 1 && return cs[1].second
    typeof(r)(cs)
end

## unknowns that are function calls: u(x), Derivative(u(x), x)
#
# Each such unknown is replaced by a fresh symbol (in the assumptions too), the equations
# solved for those, and the answer written back. Derivatives, integrals and sums that contain an unknown
# without being one are first replaced by symbols of their own, so that solving for u(x)
# leaves Derivative(u(x), x) as it is (rather than making it the derivative of a symbol).

const _SOLVE_HELD = (:Derivative, :Integral, :Sum)

function _solve_for_calls(f, eqs::Vector{<:Relation}, vars::Vector{Sym}; assume = nothing)
    for v in vars
        (_isnamed(v) || v.kind == FN) ||
            throw(ArgumentError("solve: $v is not a symbol or a function call (such as u(x) or a derivative)"))
    end
    n = Ref(0)
    fresh() = Sym(Symbol("_solve_u", n[] += 1))
    calls = [v for v in vars if !_isnamed(v)]
    # the operators holding an unknown, as symbols, outermost first
    held = Sym[]
    for eq in eqs, side in (eq.lhs, eq.rhs)
        _held_nodes!(held, side, calls)
    end
    hold = Pair{Sym,Sym}[h => fresh() for h in held]
    back = Dict{Sym,Sym}(t => h for (h, t) in hold)
    put = Pair{Sym,Sym}[]
    newvars = Sym[]
    for v in vars
        if _isnamed(v)
            push!(newvars, v)
        else
            t = fresh()
            push!(put, v => t)
            push!(newvars, t)
            back[t] = v
        end
    end
    # held operators first (so u(x) inside them stays), then the unknowns, larger first
    sort!(put; by = p -> -_nodes(p.first))
    sub(e) = subs(subs(e, hold...), put...)
    neweqs = [Relation(eq.op, sub(eq.lhs), sub(eq.rhs)) for eq in eqs]
    # facts about the unknowns (u(x) > 0) as facts about their symbols
    B = assume === nothing ? nothing :
        Assumptions((g isa Membership ? Membership(sub(g.x), g.set, g.negated) : _subsfact(g, sub) for g in assume.facts)...)
    _solve_back(f(neweqs, newvars, B), back)
end

function _held_nodes!(out, e::Sym, calls)
    if e.kind == FN && e.name in _SOLVE_HELD && !any(c -> _eq(c, e), calls) && any(c -> has(e, c), calls)
        any(h -> _eq(h, e), out) || push!(out, e)
        return out
    end
    foreach(a -> _held_nodes!(out, a, calls), e.args)
    foreach(a -> _held_nodes!(out, a, calls), e.exps)
    out
end

_solve_back(r::Sym, back) = subs(r, back)
_solve_back(r::AbstractVector, back) = [_solve_back(x, back) for x in r]
_solve_back(r::AbstractDict, back) = Dict(_solve_back(k, back) => _solve_back(v, back) for (k, v) in r)
_solve_back(r::Relation, back) = Relation(r.op, _solve_back(r.lhs, back), _solve_back(r.rhs, back))
_solve_back(r, back) = r isa Cases ? _map_cases(v -> _solve_back(v, back), r) : r   # (Cases: defined later)

# the solutions, without assumptions (what internal callers use)
function _solve_eq(eq::Relation, x::Sym)
    _isnamed(x) || throw(ArgumentError("solve: $x is not a symbol"))
    eq.op === :eq || return _solve_inequality(eq, x)
    _isnamed(x) || throw(ArgumentError("solve: $x is not a symbol"))
    n, _ = numer_denom(cancel(eq.lhs - eq.rhs))
    if !has(n, x)
        iszero(n) && throw(ArgumentError("solve: every value of $x satisfies $eq"))
        return Sym[]
    end
    atoms = Set{Sym}()
    exact = _poly_atoms!(atoms, n)
    bad = [a for a in atoms if !_eq(a, x) && has(a, x)]
    roots = if !isempty(bad)
        fs = _product_factors(eq.lhs - eq.rhs, n, x)
        fs === nothing ? _solve_through(eq, n, x, bad) : _solve_product(eq, fs, x)
    else
        exact ? _solve_exact(n, x, atoms) : _solve_inexact(n, x)
    end
    _sort_roots!(unique(_in_domain(roots, eq.lhs - eq.rhs, x)))
end

# the roots at which e is defined: cancel can remove a factor of a denominator, so
# (x^2 - 1)/(x - 1) = 2 has no solution, not x = 1. Checked through the expressions
# that must not vanish (denominators; log's argument; cos for tan, …), as tan(π/2)
# stays unevaluated with a large finite value.
function _in_domain(roots, e::Sym, x::Sym)
    isempty(roots) && return roots
    zs = unique(_zero_sets!(Sym[], e, x))
    isempty(zs) && return roots
    filter(r -> _defined_at(e, x, r) && !any(b -> _vanishes_at(b, x, r), zs), roots)
end

function _vanishes_at(b::Sym, x::Sym, r::Sym)
    v = try subs(b, x => r) catch err; err isa DivideError && return true; rethrow() end
    iszero(v) || iszero(cancel(v)) || (!_hassymbol(v) && _iszero_value(v))
end
# one equation takes one unknown and a system a vector of unknowns: the form of the call
# fixes the form of the result (values or Dicts)
const _SOLVE_FORMS = "solve(eq, x) for one equation (a vector of values), solve([eqs...], [vars...]) for a system (a vector of Dicts)"
solve(eq::Relation, x::Sym, rest...) = throw(ArgumentError("solve: give the unknowns of a system as a vector; use $_SOLVE_FORMS"))
solve(eq::Relation, vars::Union{AbstractVector,Tuple}) = throw(ArgumentError("solve: a single equation takes a single unknown; use $_SOLVE_FORMS"))
solve(eqs::AbstractVector{<:Relation}, x::Sym) = throw(ArgumentError("solve: a system takes a vector of unknowns; use $_SOLVE_FORMS"))
solve(ex::Sym, x::Sym) = throw(ArgumentError(
    "solve needs an equation, e.g. solve($ex ~ 0, $x); an expression is not taken to mean \"= 0\""))

# a product is zero where one of its factors is: solved factor by factor, a power b^p
# through its base (exp(x) has none). The candidates are checked in the whole equation
# and must leave every factor defined (x*log(x) = 0 has no root at 0, where log is not).
function _solve_product(eq::Relation, fs::Vector{Tuple{Sym,Sym}}, x::Sym)
    out = Sym[]
    for (b, p) in fs
        p.kind == NUMBER && p.num < N0 && continue                 # a denominator
        (has(b, x) || has(p, x)) && append!(out, _solve_eq(b ~ 0, x))
    end
    filter(r -> _satisfies(eq, x, r) && all(((b, p),) -> _defined_at(b^p, x, r), fs), out)
end

function _defined_at(f::Sym, x::Sym, r::Sym)
    v = try
        subs(f, x => r)
    catch err
        err isa DivideError && return false
        rethrow()
    end
    _hassymbol(v) && return true
    w = try N(v) catch; return true end
    isfinite(w)
end

# the factors (base, exponent) of e as a product, through which it is solved: e's own
# when it is one, else those of its numerator n, or the factors with x common to the
# terms of n with what remains; nothing when there are none
function _product_factors(e::Sym, n::Sym, x::Sym)
    e.kind == MUL && return _factor_pairs(e)
    n.kind == MUL && return _factor_pairs(n)
    n.kind == ADD && iszero(n.num) || return nothing
    terms = [Sym(c) * a for (c, a) in zip(n.coefs, n.args)]
    common = Dict{Sym,Sym}(b => p for (b, p) in _factor_pairs(terms[1]) if has(b, x) || has(p, x))
    for t in terms[2:end]
        ft = Dict{Sym,Sym}(_factor_pairs(t))
        for (b, p) in collect(common)
            q = get(ft, b, nothing)
            if q === nothing
                delete!(common, b)
            elseif !_eq(p, q)
                d = p - q
                d.kind == NUMBER ? (common[b] = d.num < N0 ? p : q) : delete!(common, b)
            end
        end
    end
    isempty(common) && return nothing
    C = prod((b^p for (b, p) in common); init = ONE)
    [Tuple{Sym,Sym}[(b, p) for (b, p) in common]; (_add_many(t / C for t in terms), ONE)]
end

_factor_pairs(t::Sym) = t.kind == MUL ? Tuple{Sym,Sym}[(b, p) for (b, p) in zip(t.args, t.exps)] :
                        t.kind == POW ? Tuple{Sym,Sym}[(t.args[1], t.args[2])] : Tuple{Sym,Sym}[(t, ONE)]

## equations in one function of x: f(g(x)) with f polynomial
#
# u stands for g(x); the polynomial in u is solved, then g(x) = r is inverted into an
# equation for g's argument and solved recursively. Inverting a power or a periodic
# function can introduce extraneous roots (sqrt(x) = -1), so candidates are checked
# by substitution. sin and cos give the solutions for their argument in [0, 2π].

function _solve_through(eq::Relation, n::Sym, x::Sym, bad::Vector{Sym})
    common = _common_power(bad)                     # exp(2x), exp(x): powers of one exp
    if common === nothing
        if length(bad) > 1
            r = _solve_logs(eq, n, x)
            r === nothing || return r
            # several radicals: one squared away at a time (the rest in the recursion)
            if all(_issqrt, bad) && length(bad) <= 4
                R = first(bad)
                u = Sym(Symbol("_u", objectid(R)))
                return _solve_radical(eq, subs(n, R => u), u, R, x)
            end
            throw(ArgumentError(
                "solve: $eq involves $x through several functions ($(join(bad, ", "))); not supported yet"))
        end
        g = only(bad)
        u = Sym(Symbol("_u", objectid(g)))
        nu = subs(n, g => u)
    else
        g, ks = common
        u = Sym(Symbol("_u", objectid(g)))
        nu = subs(n, (a => u^k for (a, k) in zip(bad, ks))...)
    end
    if has(nu, x)
        _issqrt(g) && return _solve_radical(eq, nu, u, g, x)
        r = _solve_lambert(eq, nu, u, g, x)
        r === nothing || return r
        throw(ArgumentError("solve: $eq involves both $x and $g; not supported (it may need, e.g., Lambert W)"))
    end
    out = Sym[]
    for r in _solve_eq(nu ~ 0, u), (arg, val) in _invert(g, r, x)
        append!(out, _solve_eq(arg ~ val, x))
    end
    filter(r -> _satisfies(eq, x, r), out)
end

# Σ cᵢ log(fᵢ) + R = 0, linear in the logarithms with rational cᵢ (constant logs such as
# log(6) included) and R free of logs: exponentiated to Π fᵢ^kᵢ = exp(-L R), kᵢ = L cᵢ
# integers. Exponentiating loses the branch information (log(-2) + log(-3) ≠ log(6)), so
# candidates are checked in the original equation. nothing if n is not of this form.
function _solve_logs(eq::Relation, n::Sym, x::Sym)
    atoms = Set{Sym}()
    _poly_atoms!(atoms, n)
    logs = [a for a in atoms if a.kind == FN && a.name === :log && length(a.args) == 1]
    any(a -> !_eq(a, x) && has(a, x) && !any(l -> _eq(l, a), logs), atoms) && return nothing
    P = Poly{Sym}(n, logs)
    cs, R = Sym[], ZERO
    ks = zeros(Int, length(logs))
    for (ev, c) in zip(exponent_vectors(P), coefficients(P))
        d = sum(ev)
        if d == 0
            has(c, x) && return nothing
            R = c
        elseif d == 1 && c.kind == NUMBER && isexact(c.num)
            push!(cs, c)
        else
            return nothing                              # log products or symbolic coefficients
        end
    end
    L = reduce(lcm, (big(_denominator(c.num)) for c in cs); init = big(1))
    for (ev, c) in zip(exponent_vectors(P), coefficients(P))
        sum(ev) == 1 && (ks[findfirst(==(1), ev)] = Int(_numerator(c.num) * (L ÷ _denominator(c.num))))
    end
    lhs = prod((logs[i].args[1]^ks[i] for i in eachindex(logs) if ks[i] != 0); init = ONE)
    cands = _solve_eq(lhs ~ exp(-Sym(SymNumber(L)) * R), x)
    filter(r -> _satisfies(eq, x, r), cands)
end

# atoms b^(e_i) with one base b free of x and exponents rational multiples of one
# another: (b^g, k_i) with b^(e_i) = (b^g)^(k_i), k_i integers; else nothing
function _common_power(atoms::Vector{Sym})
    length(atoms) > 1 || return nothing
    all(a -> a.kind == POW, atoms) || return nothing
    b = atoms[1].args[1]
    es = [a.args[2] for a in atoms]
    if !all(a -> _eq(a.args[1], b), atoms)
        # integer bases that are powers of one integer: 4^x = 2^(2x)
        bases = [a.args[1] for a in atoms]
        all(c -> c.kind == NUMBER && _isint(c.num) && c.num > SymNumber(1), bases) || return nothing
        r = _common_integer_base([big(_numerator(c.num)) for c in bases])
        r === nothing && return nothing
        b = Sym(SymNumber(r))
        es = [e * Sym(SymNumber(_ilog(big(_numerator(c.num)), r))) for (e, c) in zip(es, bases)]
    end
    e1 = es[1]
    ratios = [cancel(e / e1) for e in es]
    all(r -> r.kind == NUMBER && isexact(r.num), ratios) || return nothing
    q = reduce(lcm, (big(_denominator(r.num)) for r in ratios); init = big(1))
    g = e1 / Sym(SymNumber(q))
    ks = [Int(_numerator(r.num) * (q ÷ _denominator(r.num))) for r in ratios]
    (pow(b, g), ks)
end

# the least integer r > 1 of which every base is a power, or nothing
function _common_integer_base(bs::Vector{BigInt})
    m = minimum(bs)
    for k in 64:-1:1
        r = _perfect_root(m, k)
        r === nothing && continue
        all(c -> _ilog(c, r) !== nothing, bs) && return r
    end
    nothing
end
function _ilog(c::BigInt, r::BigInt)                 # k with r^k == c, or nothing
    k = 0
    while c > 1 && c % r == 0
        c ÷= r
        k += 1
    end
    c == 1 ? k : nothing
end

_issqrt(g::Sym) = g.kind == POW && g.args[2] == HALF

# an equation in x and R = sqrt(t): reduced with R^2 = t to A + B R = 0, squared to
# A^2 = B^2 t (polynomial), solved, and the extraneous roots of squaring dropped
function _solve_radical(eq::Relation, nu::Sym, u::Sym, R::Sym, x::Sym)
    t = R.args[1]
    atoms = Set{Sym}()
    _poly_atoms!(atoms, nu)
    P = Poly{Sym}(nu, [u])                               # coefficients in x
    A, B = ZERO, ZERO
    for (ev, c) in zip(exponent_vectors(P), coefficients(P))
        k = ev[1]
        iseven(k) ? (A += c * t^(k ÷ 2)) : (B += c * t^((k - 1) ÷ 2))
    end
    cands = _solve_eq(expand(A^2 - B^2 * t) ~ 0, x)
    filter(r -> _satisfies(eq, x, r), cands)
end

# x both directly and in an exponential g = b^(γx + δ), the equation linear in g:
# c1(x) g + c0(x) = 0 with one of c1, c0 linear and the other constant. As
# (αx + β) e^(γ'x + δ') = K (γ' = γ log b): y = γ'(x + β/α) satisfies y e^y =
# K γ'/α e^(γ'β/α - δ'), so x = W(K γ'/α e^(γ'β/α - δ'))/γ' - β/α. The other branch
# W₋₁ is added when the argument is real in (-1/e, 0) (two real solutions).
function _solve_lambert(eq::Relation, nu::Sym, u::Sym, g::Sym, x::Sym)
    g.kind == POW || return nothing
    b, e = g.args
    has(b, x) && return nothing
    γ = cancel(diff(e, x))
    has(γ, x) && return nothing                               # exponent linear in x
    δ = cancel(e - γ * x)
    lb = b === E ? ONE : log(b)
    γ, δ = γ * lb, δ * lb
    P = Poly{Sym}(nu, [u])
    degree(P, u) == 1 || return nothing
    cs = coefficients(P)
    c1, c0 = cs[1], length(cs) == 2 ? cs[2] : ZERO
    lin(c) = (Q = Poly{Sym}(c, [x]); degree(Q, x) <= 1 ? Q : nothing)
    L1, L0 = lin(c1), lin(c0)
    (L1 === nothing || L0 === nothing) && return nothing
    if degree(L1, x) == 1 && degree(L0, x) <= 0              # (αx + β) e^(γx+δ) = -c0
        α, β, K = _lin_coeffs(c1, x)..., -c0
    elseif degree(L1, x) <= 0 && degree(L0, x) == 1          # c1 e^(γx+δ) = -(αx + β)
        (α, β), K = _lin_coeffs(c0, x), -c1
        γ, δ = -γ, -δ
    else
        return nothing
    end
    iszero(K) && return filter(r -> _satisfies(eq, x, r), [-β / α])
    arg = K * γ / α * exp(γ * β / α - δ)
    ws = [lambertw(arg)]
    if !_hassymbol(arg)
        v = try N(arg) catch; nothing end
        v isa Real && -exp(-1) < v < 0 && push!(ws, lambertw(arg, -1))
    end
    filter(r -> _satisfies(eq, x, r), [w / γ - β / α for w in ws])
end
_lin_coeffs(c::Sym, x::Sym) = (α = cancel(diff(c, x)); (α, cancel(c - α * x)))

# g(x) = r as equations arg ~ value for an argument of g
function _invert(g::Sym, r::Sym, x::Sym)
    k = g.kind
    if k == FN && length(g.args) == 1
        t = g.args[1]
        name = g.name
        name === :sin && return [(t, a) for a in _one_period(asin(r), PI - asin(r))]
        name === :cos && return [(t, a) for a in _one_period(acos(r), 2PI - acos(r))]
        name === :tan && return [(t, atan(r))]
        name === :log && return [(t, exp(r))]
        name === :asin && return [(t, sin(r))]
        name === :acos && return [(t, cos(r))]
        name === :atan && return [(t, tan(r))]
        name === :sinh && return [(t, asinh(r))]
        name === :cosh && return [(t, acosh(r)), (t, -acosh(r))]
        name === :tanh && return [(t, atanh(r))]
    elseif k == POW
        b, e = g.args
        if !has(b, x)                                       # b^t = r: t = log(r)/log(b)
            iszero(r) && return Tuple{Sym,Sym}[]
            return [(e, b === E ? log(r) : log(r) / log(b))]
        elseif !has(e, x) && e.kind == NUMBER && isexact(e.num)   # t^(p/q) = r
            return [(b, r^inv(e))]
        end
    end
    throw(ArgumentError("solve: cannot invert $g for $x"))
end

# the angles of sin or cos in [0, 2π], each once (and 0 with 2π: cos(x) = 1 gives both)
function _one_period(a::Sym, b::Sym)
    (_hassymbol(a) || _hassymbol(b)) && return [a, b]
    out = Sym[]
    for c in (a, b)
        v = try N(c) catch; return [a, b] end
        v isa Real || return [a, b]
        c = c - 2PI * floor(Int, v / (2Float64(π)) + 1e-12)   # into [0, 2π)
        any(o -> isapprox(N(o), N(c); atol = 1e-12), out) || push!(out, c)
    end
    any(iszero, out) && push!(out, 2PI)
    sort!(out; by = N)
end

# does r satisfy eq? (exactly, or numerically for constants; kept when undecidable)
function _satisfies(eq::Relation, x::Sym, r::Sym)
    d = try
        subs(eq.lhs - eq.rhs, x => r)
    catch err
        err isa DivideError && return false
        rethrow()
    end
    (iszero(d) || iszero(cancel(d))) && return true
    _hassymbol(d) && return true
    v = try N(d) catch; return true end
    abs(v) <= 1e-9 * max(1, abs(N(r)))
end

function _solve_exact(n::Sym, x::Sym, atoms)
    gens = [x; sort!([a for a in atoms if !_eq(a, x)]; lt = _lt)]
    _, facs = factor(Poly{SymNumber}(n, gens))
    roots = Sym[]
    for (f, _) in facs
        degree(f, x) == 0 && continue
        append!(roots, _solve_irreducible(_xcoeffs(f), f))
    end
    roots
end

# coefficients c[k+1] of x^k (x the first generator), as expressions
_xcoeffs(f::Poly{SymNumber}) = [Sym(_coeff_of(f, 1, k)) for k in 0:_deg_in_all(f, 1)]

# numeric coefficients only (floats): numeric roots of the whole numerator
function _solve_inexact(n::Sym, x::Sym)
    c = try
        _xcoeffs(Poly{SymNumber}(n, [x]))
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("solve: inexact coefficients that are not numbers are not supported"))
    end
    _numeric_roots(c)
end

## an irreducible factor, by degree

function _solve_irreducible(c::Vector{Sym}, f)
    d = length(c) - 1
    d == 1 && return [-c[1] / c[2]]
    d == 2 && return _quadratic(c...)
    if all(iszero, c[2:(end-1)])                              # a xⁿ + b
        r = (-c[1] / c[end])^(1 // d)
        return [r * pow(MINUS_ONE, Sym(SymNumber(2k // d))) for k in 0:(d - 1)]
    end
    d == 3 && return _cardano(c...)
    d == 4 && return _ferrari(c...)
    r = _solve_decomposed(c)
    r === nothing || return r
    all(ci -> !_hassymbol(ci), c) || throw(ArgumentError(
        "solve: an irreducible factor of degree $d with symbolic coefficients"))
    _numeric_roots(c)
end

# f = g₁ ∘ ⋯ ∘ gₖ with every gᵢ of degree at most 4: the roots of g₁, then for each value
# v the roots of g₂ = v, and so on; nothing when f has no such decomposition
function _solve_decomposed(c::Vector{Sym})
    chain = _decompose(c)
    (length(chain) > 1 && all(g -> length(g) <= 5, chain)) || return nothing
    vals = _solve_irreducible(chain[1], nothing)
    for g in chain[2:end]
        vals = reduce(vcat, (_solve_irreducible([g[1] - v; g[2:end]], nothing) for v in vals))
    end
    vals
end

function _quadratic(c0, c1, c2)
    s = sqrt(c1^2 - 4 * c2 * c0)
    [(-c1 - s) / (2c2), (-c1 + s) / (2c2)]
end

# is a constant expression zero? (exactly, or numerically when it has no symbols)
function _iszero_value(e::Sym)
    iszero(e) && return true
    _hassymbol(e) && return false
    v = try N(e) catch; return false end
    abs(v) < 1e-12
end

const _OMEGA = (-1 + sqrt(Sym(-3))) / 2                      # a primitive cube root of 1

function _cardano(d, c, b, a)                                # a x³ + b x² + c x + d
    # three distinct real roots and numeric coefficients: the trigonometric form, exact
    # and real (Cardano's radicals would need complex numbers: casus irreducibilis)
    if !any(_hassymbol, (a, b, c, d))
        p = (3a * c - b^2) / (3a^2)
        q = (2b^3 - 9a * b * c + 27a^2 * d) / (27a^3)
        disc = -(4p^3 + 27q^2)
        v = N(disc)
        if v isa Real && v > 0
            θ = acos(3q / (2p) * sqrt(-3 / p)) / 3
            return [2sqrt(-p / 3) * cos(θ - 2PI * k / 3) - b / (3a) for k in 0:2]
        end
    end
    Δ0 = b^2 - 3a * c
    Δ1 = 2b^3 - 9a * b * c + 27a^2 * d
    s = sqrt(Δ1^2 - 4Δ0^3)
    C = ((Δ1 + s) / 2)^(1 // 3)
    _iszero_value(C) && (C = ((Δ1 - s) / 2)^(1 // 3))
    [-(b + _OMEGA^k * C + Δ0 / (_OMEGA^k * C)) / (3a) for k in 0:2]
end

function _ferrari(e, d, c, b, a)                             # a x⁴ + b x³ + c x² + d x + e
    p = (8a * c - 3b^2) / (8a^2)
    q = (b^3 - 4a * b * c + 8a^2 * d) / (8a^3)
    shift = -b / (4a)
    if _iszero_value(q)                                      # biquadratic in y = x - shift
        r = (-3b^4 + 256a^3 * e - 64a^2 * b * d + 16a * b^2 * c) / (256a^4)
        out = Sym[]
        for y2 in _quadratic(r, p, ONE), sgn in (-1, 1)
            push!(out, shift + sgn * sqrt(y2))
        end
        return out
    end
    Δ0 = c^2 - 3b * d + 12a * e
    Δ1 = 2c^3 - 9b * c * d + 27b^2 * e + 27a * d^2 - 72a * c * e
    Q = ((Δ1 + sqrt(Δ1^2 - 4Δ0^3)) / 2)^(1 // 3)
    _iszero_value(Q) && (Q = ((Δ1 - sqrt(Δ1^2 - 4Δ0^3)) / 2)^(1 // 3))
    S = sqrt(-2p / 3 + (Q + Δ0 / Q) / (3a)) / 2
    for _ in 1:2                                             # S ≠ 0 for some cube root of Q
        _iszero_value(S) || break
        Q *= _OMEGA
        S = sqrt(-2p / 3 + (Q + Δ0 / Q) / (3a)) / 2
    end
    u, v = sqrt(-4S^2 - 2p + q / S) / 2, sqrt(-4S^2 - 2p - q / S) / 2
    [shift - S - u, shift - S + u, shift + S - v, shift + S + v]
end

## numeric roots: Aberth–Ehrlich iteration

function _numeric_roots(c::Vector{Sym})
    a = ComplexF64[_evalc(ci) for ci in c]
    while length(a) > 1 && iszero(a[end])
        pop!(a)
    end
    n = length(a) - 1
    n == 0 && return Sym[]
    a ./= a[end]                                             # monic
    R = 1 + maximum(abs, a[1:(end-1)])                       # Cauchy's root bound
    z = [R * cis(2π * (k + 0.25) / n) for k in 0:(n - 1)]
    p(t) = evalpoly(t, a)
    dp(t) = evalpoly(t, [k * a[k+1] for k in 1:n])
    for _ in 1:500
        δmax = 0.0
        for k in 1:n
            w = p(z[k]) / dp(z[k])
            s = sum(1 / (z[k] - z[j]) for j in 1:n if j != k; init = 0.0im)
            δ = w / (1 - w * s)
            z[k] -= δ
            δmax = max(δmax, abs(δ) / max(1.0, abs(z[k])))
        end
        δmax < 1e-15 && break
    end
    map(z) do t
        abs(imag(t)) <= 1e-12 * max(1.0, abs(t)) ? Sym(real(t)) : Sym(real(t)) + Sym(imag(t)) * IM
    end
end

# numeric (constant) roots sorted: real ones first, by real then imaginary part
function _sort_roots!(roots::Vector{Sym})
    any(_hassymbol, roots) && return roots
    vals = map(r -> (v = try _evalc(r) catch; return nothing end; v), roots)
    any(isnothing, vals) && return roots
    key(v) = (abs(imag(v)) > 1e-10 * max(1, abs(v)), real(v), imag(v))
    roots[sortperm([key(v) for v in vals])]
end
