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
differentiation; other derivatives of it are held as they are. In a system, `abs`,
`real`, `imag` and `conj` of the unknowns are solved over the complex numbers, through
each unknown's real and imaginary parts (`[I*R ~ 220, abs(I) ~ 2]`: a family written
with a free part such as `R_im`); with the unknowns assumed real (`assume = (I ∈ ℝ,)`),
`abs` is `±` its argument and only real solutions are kept (`I = ±2`).

An equation through `sin` or `cos` gives the solutions for its argument in `[0, 2π]`,
each once: `cos(x) ~ 0` gives `[π/2, 3π/2]`, `sin(x) ~ -1//2` gives `[7π/6, 11π/6]`, and
`cos(x) ~ 1` gives `[0, 2π]` (with the interval closed, 0 brings 2π). `tan(u) = r` gives
`atan(r)`. When the assumptions keep `x` in a bounded interval and the argument is
linear in `x`, every solution in that interval is given (`cos(2x) ~ 1//2` with
`0 ≤ x ≤ 2π`: four). An equation `c*floor(u) + d ~ 0` (or `ceil`) holds on intervals: its solution
is a set, as for an inequality (`floor(x) ~ 2` gives `[2, 3)`).

```julia
solve(x^2 - 2 ~ 0, x)          # [-sqrt(2), sqrt(2)]
solve(x^2 ~ -1, x)             # [-im, im]
solve(diff(x^2 + u(x)^2 - 1, x) ~ 0, diff(u(x), x))   # [-x/u(x)]
```

With parameters, the solutions are correct for every value of them: where a
coefficient could vanish, the result is a `Cases` value, a list of
`conditions => solutions` (called with values, `cs(a => 4)`, it gives that case's
solutions). The cases come from the leading coefficients of a polynomial (the next one
when it vanishes), from denominators free of `x` (where they vanish the equation is not
defined: no solutions), from denominators of the solutions, and from a symbolic
exponent (`x^a = b` for `a = 0`; for `a ≠ 0` the answer is the principal root `b^(1/a)`,
not the other roots `b^(1/a) e^(2πik/a)`). A case in which every value of `x` is a
solution holds an error. Assumptions that settle a coefficient give the plain answer.

```julia
solve(a*x + b ~ 0, x)          # a == 0, b == 0 => (every x), a == 0, b != 0 => [], a != 0 => [-b/a]
solve(a*x + b ~ 0, x; assume = (a ≫ 0,))   # [-b/a]
solve(exp(a*x) ~ 2, x)         # a == 0 => [], a != 0 => [log(2)/a]
solve(x^2 ~ a, x)              # [-sqrt(a), sqrt(a)]
solve(x^a ~ 2, x)              # a == 0 => [], a != 0 => [2^(1/a)]
```

With assumptions (`assume = (x ≫ 0,)`, an [`Assumptions`](@ref), or inside
[`assuming`](@ref)), only the solutions they allow for `x` are kept: those in its
number set (`x ∈ ℝ`), its range (`x ≫ 0`), and consistent with its relations to other
symbols. Facts about the other symbols are used to decide that. For a real `x`, a
floating point root whose imaginary part is tiny beside its size (a double root split
by rounding, `264 ± 6e-6im`) is taken as real, with a warning. When the answer
depends on something the assumptions leave open, the result is again a `Cases` value.

```julia
solve(x^2 ~ 4, x; assume = (x ≫ 0,))                # [2]
solve(x^2 ~ a, x; assume = (x ∈ ℝ, a ∈ ℝ))         # a < 0 => [], a == 0 => [0], a > 0 => [-√a, √a]
solve(x^2 ~ a, x; assume = (x ≫ 0,))                # a > 0 => [√a], !(a > 0) => []
```

`drop_empty = true` drops the cases with no solutions (and those holding an error, such
as every value of `x` solving the equation), and when one case is left gives its
solutions alone: the conditions under which they hold are then not shown.

```julia
solve(x^2 ~ a, x; assume = (x ≫ 0,), drop_empty = true)   # [√a] (for a > 0)
```
"""
function solve(eq::Relation, x::Sym; assume = nothing, drop_empty::Bool = false)
    # c floor(u) + d = 0 (or ceil) holds on intervals of u: its solution set, as for an
    # inequality (SymPy #25001)
    eq.op === :eq && _isnamed(x) && (S = _solve_floor(eq, x, assume); S !== nothing) && return S
    A = something(_assumptions(assume), Assumptions())
    r = try
        if !_isnamed(x)
            with(_NOT_PARAMS => _named(x)) do
                _solve_for_calls((es, vs, B) -> solve(only(es), only(vs); assume = B), [eq], [x]; assume = A)
            end
        else
            _solve_assumed(eq, x, A)
        end
    catch err
        (err isa ArgumentError && occursin("inexact coefficients", err.msg)) || rethrow()
        r = _solve_floats_exactly(eq, x, assume)
        r === nothing ? rethrow() : r
    end
    drop_empty ? _drop_empty(r) : r
end

# floats where the exact methods need exact numbers (t^1.1, a float exponent): each as
# the shortest rational within its rounding, solved exactly, the solutions numbers again
# (SymPy #23284, #22584); nothing when that fails too
function _solve_floats_exactly(eq::Relation, x::Sym, assume)
    (_isexactex(eq.lhs) && _isexactex(eq.rhs)) && return nothing
    e2 = Relation(eq.op, _exactfloats(eq.lhs; short = true), _exactfloats(eq.rhs; short = true))
    r = try solve(e2, x; assume) catch err; err isa ArgumentError || rethrow(); return nothing end
    num(v) = v isa Sym ? (_hassymbol(v) ? _floatify(v) : (w = N(v); w isa Number ? Sym(w) : v)) : v
    r isa Cases && return _map_cases(vs -> vs isa AbstractVector ? map(num, vs) : vs, r)
    r isa AbstractVector ? map(num, r) : nothing
end

# floor(u) = n ⇔ n ≤ u < n + 1 and ceil(u) = n ⇔ n - 1 < u ≤ n, n an integer (no
# solutions otherwise), for an equation c floor(u) + d = 0 with c, d constants
function _solve_floor(eq::Relation, x::Sym, assume)
    _has_fn(eq.lhs - eq.rhs, (:floor, :ceil)) || return nothing     # (no expanding first)
    g = expand(eq.lhs - eq.rhs)
    F = nothing
    _replace_walk(g) do t
        t.kind == FN && t.name in (:floor, :ceil) && length(t.args) == 1 && has(t.args[1], x) &&
            (F === nothing || _eq(F, t)) && (F = t)
        nothing
    end
    F === nothing && return nothing
    u = Sym(Symbol("_floor_u"))
    h = expand(subs(g, F => u))
    has(h, x) && return nothing
    P = try Poly{Sym}(h, [u]) catch; return nothing end
    degree(P, u) == 1 || return nothing
    c, d = coeff(h, u, 1), coeff(h, u, 0)
    (_hassymbol(c) || _hassymbol(d)) && return nothing
    n = -d / c
    (n.kind == NUMBER && isexact(n.num)) || return nothing
    _isint(n.num) || return SymSet(SymInterval[])
    w = F.args[1]
    lo, hi = F.name === :floor ? (Relation(:ge, w, n), Relation(:lt, w, n + 1)) :
                                  (Relation(:gt, w, n - 1), Relation(:le, w, n))
    S1, S2 = solve(lo, x; assume), solve(hi, x; assume)
    (S1 isa SymSet && S2 isa SymSet) || return nothing
    S1 ∩ S2
end

# the generic solutions, coefficients with parameters taken to be nonzero (what internal
# callers use)
_solve_generic(eq::Relation, x::Sym) = _isnamed(x) ? _solve_eq(eq, x) :
    _solve_for_calls((es, vs, _) -> _solve_generic(only(es), only(vs)), [eq], [x])

# the cases that have solutions; their solutions alone when one case is left. Cases that
# hold an error (every value a solution, a case not supported) are dropped too.
function _drop_empty(r)
    r isa Cases || return r
    cs = [c for c in r.cases if !(c.second isa AbstractVector && isempty(c.second)) && !(c.second isa Exception)]
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
    # an equation left free of the unknowns but not of their variables (cos(2a) = 0 for
    # cos(a)): its other functions of them in terms of the unknowns (cos(2a) = 2cos(a)² - 1),
    # refused when that leaves the variables still; taken as a constant it gave no solutions (SymPy #12129)
    argsyms = Set{Sym}(s for c in calls for s in free_symbols(c))
    lost(q) = !any(t -> has(q.lhs, t) || has(q.rhs, t), newvars) && any(s -> has(q.lhs, s) || has(q.rhs, s), argsyms)
    if any(lost, neweqs)
        neweqs = [Relation(eq.op, sub(expand_trig(eq.lhs)), sub(expand_trig(eq.rhs))) for eq in eqs]
        # (all of them: 2cos(a)sin(a) = 0 for sin(a) has cos(a) = 0 too)
        stays(q) = any(s -> has(q.lhs, s) || has(q.rhs, s), argsyms)
        any(stays, neweqs) && throw(ArgumentError(
            "solve: the equations involve the variables of $(join(calls, ", ")) other than through them; not supported"))
    end
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
# (a Cases is a vector too: its values mapped, not its pairs left as they are; SymPy #23459)
_solve_back(r::AbstractVector, back) = r isa Cases ? _map_cases(v -> _solve_back(v, back), r) : [_solve_back(x, back) for x in r]
_solve_back(r::AbstractDict, back) = Dict(_solve_back(k, back) => _solve_back(v, back) for (k, v) in r)
_solve_back(r::Relation, back) = Relation(r.op, _solve_back(r.lhs, back), _solve_back(r.rhs, back))
_solve_back(r, back) = r isa Cases ? _map_cases(v -> _solve_back(v, back), r) : r   # (Cases: defined later)

# the solutions, without assumptions (what internal callers use)
function _solve_eq(eq::Relation, x::Sym)
    _isnamed(x) || throw(ArgumentError("solve: $x is not a symbol"))
    eq.op === :eq || return _solve_inequality(eq, x)
    _isnamed(x) || throw(ArgumentError("solve: $x is not a symbol"))
    r = _solve_piecewise(eq, x)
    r === nothing || return r
    n, d = numer_denom(cancel(eq.lhs - eq.rhs))
    # tan, cot, sec, csc of x in a denominator: their poles are zeros of the expression
    # (1/tan x = cos x/sin x vanishes at π/2), which a numerator taken over them loses;
    # solved in sin and cos instead (SymPy #12527)
    _has_recip_trig(d, x) && return _solve_eq(_recip_trig_sincos(eq.lhs - eq.rhs, x) ~ 0, x)
    if !has(n, x)
        iszero(n) && throw(_identity_error(x, eq))
        return Sym[]
    end
    atoms = Set{Sym}()
    exact = _poly_atoms!(atoms, n)
    bad = [a for a in atoms if !_eq(a, x) && has(a, x)]
    # floating point coefficients that are not a polynomial's numbers alone (those get
    # numeric roots): solved with the rationals they stand for, the solutions' numbers
    # made floats again (y = c + 0.5b, 175 = 1.5 log(x); SymPy #27528, #27712)
    if !_isexactex(n) && !_float_exponent(n) && (!isempty(bad) || any(a -> !_eq(a, x), atoms) || !exact)
        q = _exactfloats(eq.lhs - eq.rhs)
        if !_isexactex(q)
            # (an infinite or NaN coefficient is no complex one: say so; SymPy #11015)
            _has_nonfinite(q) &&
                throw(ArgumentError("solve: $eq has an infinite or NaN coefficient; no finite solutions are found from it"))
            throw(ArgumentError("solve: complex floating point coefficients are not supported"))
        end
        rs = _solve_eq(q ~ 0, x)
        return rs isa AbstractVector ? Sym[_floatify(r) for r in rs] : rs
    end
    roots = if !isempty(bad)
        fs = _product_factors(eq.lhs - eq.rhs, n, x)
        fs === nothing ? _solve_through(eq, n, x, bad) : _solve_product(eq, fs, x)
    else
        exact ? _solve_exact(n, x, atoms) : _solve_inexact(n, x)
    end
    # (each side's points: lhs - rhs may have cancelled them, 1/(x² - 4) + x = 1/(x² - 4) + 2
    # gave the pole x = 2; SymPy #9522)
    _sort_roots!(unique(_in_domain(_in_domain(roots, eq.lhs, x), eq.rhs, x)))
end

_recip_trig_sincos(ex::Sym, x::Sym) = _replace_walk(ex) do e
    (e.kind == FN && e.name in (:tan, :cot, :sec, :csc) && length(e.args) == 1 && has(e.args[1], x)) || return nothing
    s, c = sin(e.args[1]), cos(e.args[1])
    e.name === :tan ? s / c : e.name === :cot ? c / s : e.name === :sec ? 1 / c : 1 / s
end
_has_recip_trig(e::Sym, x::Sym) =
    (e.kind == FN && e.name in (:tan, :cot, :sec, :csc) && length(e.args) == 1 && has(e.args[1], x)) ||
    any(a -> _has_recip_trig(a, x), e.args) || any(a -> _has_recip_trig(a, x), e.exps)

# a floating point exponent somewhere in e (x^1.11: its exact rational has a huge
# denominator, so those keep their own path)
_float_exponent(e::Sym) = (e.kind == POW && e.args[2].kind == NUMBER && isfloat(e.args[2].num)) ||
    (e.kind == MUL && any(p -> p.kind == NUMBER && isfloat(p.num), e.exps)) ||
    any(_float_exponent, e.args) || any(_float_exponent, e.exps)

# the numbers of e as floats, except exponents and coefficients ±1 (b/2 + y as
# 0.5*b + y), and its constant parts as their values (exp(350/3), π/6)
function _floatify(e::Sym)
    f(c::SymNumber) = isexact(c) ? SymNumber(Float64(c)) : c
    if !_hassymbol(e) && e.kind != NUMBER       # a constant: its value (π/6, exp(350/3))
        v = try N(e) catch; nothing end
        # (a real root in complex radicals: its imaginary part of rounding size dropped,
        # not left as 1e-77im; SymPy #9691)
        if v isa Complex && isfinite(v) && abs(imag(v)) <= 1e-12 * abs(v) && _isreal_const(e) === true
            v = real(v)
        end
        v isa Number && isfinite(v) && return Sym(v)
    end
    k = e.kind
    k == NUMBER && return Sym(f(e.num))
    if k == ADD
        return Sym(f(e.num)) + _add_many((_isN1(abs(c)) ? Sym(c) : Sym(f(c))) * _floatify(a) for (c, a) in zip(e.coefs, e.args))
    elseif k == MUL
        c = _isN1(abs(e.num)) ? Sym(e.num) : Sym(f(e.num))
        # (exp(-600t) stays exp, not 2.718281828459045^(-600t); im = (-1)^(1/2) stays im, not
        # sqrt(-1.0), SymPy #14371)
        fb(a) = a === E || _eq(a, MINUS_ONE) ? a : _floatify(a)
        return c * prod((fb(a)^x for (a, x) in zip(e.args, e.exps)); init = ONE)
    elseif k == POW
        # (e^u: u's numbers floats too, exp(18014398509481984t/5404319552844595) was
        # exp(3.33t) written exactly; SymPy #11542)
        e.args[1] === E && return E^_floatify(e.args[2])
        return (_eq(e.args[1], MINUS_ONE) ? e.args[1] : _floatify(e.args[1]))^e.args[2]
    elseif k == FN
        # (an unknown function's arguments are its own: f(2) stays f(2), not f(2.0), which
        # is another expression; SymPy #22524)
        haskey(FUNCTION_RULES, e.name) || return e
        return apply_function(e.name, Sym[_floatify(a) for a in e.args])
    end
    e
end

# an infinite or NaN float among e's numbers (coefficients included)
_has_nonfinite(e::Sym) = (bad(c) = isfloat(c) && !isfinite(Float64(c));
    (e.kind == NUMBER && bad(e.num)) || (e.kind in (ADD, MUL) && bad(e.num)) || (e.kind == ADD && any(bad, e.coefs)) ||
    any(_has_nonfinite, e.args) || any(_has_nonfinite, e.exps))

# the error for an equation that every value satisfies, naming the points where it is
# not defined (x/(x - 1) - 1/(x - 1) = 1 holds for every x except 1; SymPy #23254)
function _identity_error(x::Sym, eq)
    e = eq isa Relation ? eq.lhs - eq.rhs : Sym(eq)
    pts = Sym[]
    for z in unique(_zero_sets!(Sym[], e, x))
        rs = try _solve_eq(z ~ 0, x) catch err
            err isa Union{ArgumentError,DivisionByZero} || rethrow()
            continue
        end
        (rs isa AbstractVector && !(rs isa Cases)) || continue
        for r in rs
            r isa Sym && !has(r, x) && !any(q -> _eq(q, r), pts) && push!(pts, r)
        end
    end
    isempty(pts) && return ArgumentError("solve: every value of $x satisfies $eq")
    ArgumentError("solve: every value of $x satisfies $eq, except " * join(("$x = $r" for r in pts), ", ") *
                  ", where it is not defined")
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
    v = try subs(b, x => r) catch err; _isdivzero(err) && return true; rethrow() end
    iszero(v) && return true
    # with parameters: nonzero at a sample of them is not identically zero, without a
    # cancel over a cubic's radicals (minutes for Welch's equation; SymPy #25920)
    _hassymbol(v) && _nonzero_somewhere(v) && return false
    iszero(cancel(v)) || (!_hassymbol(v) && _iszero_value(v))
end

# is v (with parameters) clearly nonzero at some exact sample of their values?
function _nonzero_somewhere(v::Sym)
    syms = collect(free_symbols(v))
    for k in 1:3
        vals = [s => Sym(_SAMPLE_VALUES[mod1(k + 3j, length(_SAMPLE_VALUES))]) for (j, s) in enumerate(syms)]
        w = try N(subs(v, vals...)) catch; continue end
        w isa Number && isfinite(w) && abs(w) > 1e-6 && return true
    end
    false
end
# one equation takes one unknown and a system a vector of unknowns: the form of the call
# fixes the form of the result (values or Dicts)
const _SOLVE_FORMS = "solve(eq, x) for one equation (a vector of values), solve([eqs...], [vars...]) for a system (a vector of Dicts)"
# (keywords accepted, so a call with assume = … gets the message too; SymPy #25783)
solve(eq::Relation, x::Sym, rest...; kws...) = throw(ArgumentError("solve: give the unknowns of a system as a vector; use $_SOLVE_FORMS"))
solve(eq::Relation, vars::Union{AbstractVector,Tuple}; kws...) = throw(ArgumentError("solve: a single equation takes a single unknown; use $_SOLVE_FORMS"))
solve(eqs::AbstractVector{<:Relation}, x::Sym; kws...) = throw(ArgumentError("solve: a system takes a vector of unknowns; use $_SOLVE_FORMS"))
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
        _isdivzero(err) && return false
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

const _WEIERSTRASS_FNS = (:sin, :cos, :tan, :cot, :sec, :csc)

# the one argument of trigonometric functions, all of them among sin, cos, tan, cot, sec,
# csc; nothing otherwise
function _multiple_angle_form(n::Sym, bad::Vector{Sym}, x::Sym)
    all(b -> b.kind == FN && b.name in _WEIERSTRASS_FNS && length(b.args) == 1, bad) || return nothing
    cs = SymNumber[]
    for b in bad
        c = cancel(b.args[1] / x)
        (c.kind == NUMBER && isexact(c.num) && !iszero(c.num)) || return nothing
        push!(cs, c.num)
    end
    L = reduce(lcm, (_denominator(c) for c in cs); init = big(1))
    ns = [BigInt(_numerator(c * SymNumber(L))) for c in cs]
    G = reduce(gcd, ns)
    ks = [k ÷ G for k in ns]
    maximum(abs, ks) <= 12 || return nothing
    W = Sym(Symbol("_mw", objectid(x)))
    w = x * Sym(SymNumber(G) / SymNumber(L))
    n2 = n
    for (b, k) in zip(bad, ks)
        n2 = subs(n2, b => expand_trig(apply_function(b.name, Sym[Sym(SymNumber(k)) * W])))
    end
    (n2, w, W) |> t -> (subs(t[1], t[3] => t[2]), t[2])
end

function _common_trig_arg(bad::Vector{Sym})
    all(b -> b.kind == FN && b.name in _WEIERSTRASS_FNS && length(b.args) == 1, bad) || return nothing
    w = bad[1].args[1]
    all(b -> _eq(b.args[1], w), bad) ? w : nothing
end

# n(sin w, cos w, tan w, …) = 0 through t = tan(w/2): sin w = 2t/(1 + t²), cos w =
# (1 - t²)/(1 + t²), a rational equation in t; w = 2 atan(t) in [0, 2π), and w = π
# (t = ∞) when it is a solution. Candidates are checked in the original equation.
function _solve_weierstrass(eq::Relation, n::Sym, x::Sym, w::Sym)
    t = Sym(Symbol("_t", objectid(w)))
    S, C = 2t / (1 + t^2), (1 - t^2) / (1 + t^2)
    sub = Dict(:sin => S, :cos => C, :tan => S / C, :cot => C / S, :sec => 1 / C, :csc => 1 / S)
    nt = _replace_walk(n) do e
        (e.kind == FN && e.name in _WEIERSTRASS_FNS && length(e.args) == 1 && _eq(e.args[1], w)) || return nothing
        sub[e.name]
    end
    num, _ = numer_denom(cancel(nt))
    num = expand(num)
    # x outside the trigonometric functions too (2x - 30 - 2(3 - sin x)cos x): no equation
    # in t alone, and its "roots" in t are no roots at all (SymPy #26253)
    has(num, x) && throw(ArgumentError(
        "solve: $eq involves $x both inside and outside trigonometric functions; not supported"))
    # in a solve by cases, a leading coefficient that may vanish is decided: where it does,
    # the degree drops and t = ∞ (w = π) is a root (a sin x + b cos x = c, b + c = 0)
    at_pi = false
    A = _SOLVE_ASSUME[]
    if A !== nothing && _ispoly_in(num, t)
        while has(num, t)
            k = degree(Poly{Sym}(num, [t]), t)
            lc = coeff(num, t, k)
            all(f -> _decide(Relation(:ne, f, ZERO), A), _param_factors(lc, (x, t))) && break
            num = expand(num - lc * t^k)
            at_pi = true
        end
        if !has(num, t)                           # a constant: zero for every w, or never
            (iszero(num) || !all(f -> _decide(Relation(:ne, f, ZERO), A), _param_factors(num, (x, t)))) &&
                throw(_identity_error(x, eq))
            return at_pi ? filter(r -> _defined_at(n, x, r), _eq(w, x) ? Sym[PI] : _solve_generic(w ~ PI, x)) : Sym[]
        end
    end
    # with algebraic numbers in the coefficients (√3), through the polynomial over Q, which
    # factors: the roots come out in radicals of low degree, not by Cardano
    lifted = _solve_lifted(Sym[num], Sym[t], true)
    ts = lifted === nothing ? _solve_generic(num ~ 0, t) : Sym[d[t] for d in lifted if haskey(d, t)]
    ts isa AbstractVector || throw(ArgumentError("solve: $eq: the equation in tan($w/2) has cases"))
    ws = Sym[PI]
    for r in ts
        # t = ±i: no finite w (tan(w/2) = ±i only as Im w → ±∞); atan raised a pole
        # (cos x + i sin x = 1 has t = 0 and t = i; SymPy #10809)
        (iszero(expand(r - IM)) || iszero(expand(r + IM))) && continue
        u = 2atan(r)
        # 2atan(r) as atan(2r/(1 - r²)) (± π past |r| = 1) when smaller: 2atan(√10/3 - 1/3)
        # is atan(3) (SymPy #19102)
        rv = try N(r) catch; nothing end
        if rv isa Real && !isapprox(abs(rv), 1) && !_hassymbol(r)
            # (in the field of r's radicals, denominators rationalized: √10 terms cancel to 3)
            q = try Sym(_AlgF(_algnorm(2r)) * inv(_AlgF(_algnorm(1 - r^2)))) catch err
                err isa InterruptException && rethrow(); nothing
            end
            if q isa Sym                                  # (not w: that is the argument)
                c = atan(q) + (abs(rv) < 1 ? ZERO : rv > 0 ? PI : -PI)
                _nodes(c) < _nodes(u) && (u = c)
            end
        end
        v = try N(u) catch; nothing end
        v isa Real && v < 0 && (u += 2PI)
        any(s -> _eq(s, u), ws) || push!(ws, u)
    end
    out = Sym[]
    for (i, u) in enumerate(ws)
        for r in (_eq(w, x) ? Sym[u] : _solve_generic(w ~ u, x))
            any(s -> _eq(s, r), out) && continue
            # w = π (t = ∞) is no root of the polynomial in t: kept only where the equation
            # is known to hold, not when parameters are left (a sin x + b cos x = c holds at
            # π only for b + c = 0)
            (i == 1 ? (at_pi || _holds_at_point(eq, x, r)) : _satisfies(eq, x, r)) || continue
            _defined_at(n, x, r) && push!(out, r)
        end
    end
    out
end

# eq holds at x = r: exactly, or numerically when no symbols are left
function _holds_at_point(eq::Relation, x::Sym, r::Sym)
    d = try subs(eq.lhs - eq.rhs, x => r) catch err; _isdivzero(err) && return false; rethrow() end
    (iszero(d) || iszero(cancel(d))) && return true
    _hassymbol(d) && return iszero(simplify(d))
    v = try N(d) catch; return false end
    abs(v) <= 1e-9 * max(1, abs(N(r)))
end

function _solve_through(eq::Relation, n::Sym, x::Sym, bad::Vector{Sym})
    # several hyperbolic functions (and exponentials): all as exponentials, sinh(8x) +
    # cosh(8x) = e^(8x) (SymPy #19489)
    hyp = (:sinh, :cosh, :tanh, :coth, :sech, :csch)
    if length(bad) > 1 && any(g -> g.kind == FN && g.name in hyp, bad) &&
       all(g -> (g.kind == FN && g.name in hyp) || (g.kind == POW && g.args[1] === E), bad)
        n2 = numer_denom(cancel(expand(rewrite(n, (sinh, cosh, tanh, coth, sech, csch) => exp))))[1]
        _has_fn(n2, hyp) || return filter(r -> _satisfies(eq, x, r), _solve_generic(n2 ~ 0, x))
    end
    common = _common_power(bad)                     # exp(2x), exp(x): powers of one exp
    if common === nothing
        if length(bad) > 1
            r = _solve_logs(eq, n, x)
            r === nothing || return r
            r = _solve_exp_bases(eq, x)
            r === nothing || return r
            r = _solve_exp_sum(eq, n, x)
            r === nothing || return r
            # cube roots among them (or a square root beside a cube root): each radical an
            # unknown w with w^q = R, the polynomial system solved (SymPy #21167, #21026)
            if !all(_issqrt, bad) && length(bad) <= 3 && all(b -> _issqrt(b) || _iscuberoot(b), bad)
                sols = _solve_radical_system(Relation[eq], [x])
                sols === nothing || return _sort_roots!(unique([s[x] for s in sols if !_hassymbol(s[x]) || !has(s[x], x)]))
            end
            # several radicals: one squared away at a time (the rest in the recursion)
            if all(_issqrt, bad) && length(bad) <= 4
                R = first(bad)
                u = Sym(Symbol("_u", objectid(R)))
                return _solve_radical(eq, subs(n, R => u), u, R, x)
            end
            # sin, cos, tan, … of one argument: through t = tan(w/2)
            w = _common_trig_arg(bad)
            w === nothing || return _solve_weierstrass(eq, n, x, w)
            # rational multiples c x of x: each f(c x) = f(k w), w = x·g for the largest
            # g with every c/g an integer k, expanded in sin w, cos w (sin x + cos(x/2) = 0
            # in w = x/2; SymPy #18315)
            mw = _multiple_angle_form(n, bad, x)
            mw === nothing || return _solve_weierstrass(eq, mw[1], x, mw[2])
            throw(ArgumentError(
                "solve: $eq involves $x through several functions ($(join(bad, ", "))); not supported yet"))
        end
        g = only(bad)
        u = Sym(Symbol("_u", objectid(g)))
        nu = subs(n, g => u)
    else
        g, ks, cs = common
        u = Sym(Symbol("_u", objectid(g)))
        nu = subs(n, (a => c * u^k for (a, k, c) in zip(bad, ks, cs))...)
    end
    if has(nu, x)
        _issqrt(g) && return _solve_radical(eq, nu, u, g, x)
        # log(x) and x: x = e^y, an equation in y and e^y (1/x = 1/log(x)²: y² = e^y; #26461)
        if g.kind == FN && g.name === :log && length(g.args) == 1 && _eq(g.args[1], x)
            y = Sym(Symbol("_y", objectid(g)))
            ys = _solve_eq(subs(nu, u => y, x => exp(y)) ~ 0, y)
            return filter(r -> _satisfies(eq, x, r), Sym[exp(r) for r in ys])
        end
        r = _solve_lambert(eq, nu, u, g, x)
        r === nothing || return r
        r = _solve_generator_first(eq, nu, u, g, x)
        r === nothing || return r
        r = _solve_inner_generator(eq, n, x, g)
        r === nothing || return r
        throw(ArgumentError("solve: $eq involves both $x and $g; not supported (it may need, e.g., Lambert W)"))
    end
    out = Sym[]
    for r in _solve_eq(nu ~ 0, u), (arg, val) in _invert(g, r, x)
        # x = val needs no solving: solving would expand val (a Cardano root to the 10th
        # power), whose numeric value then loses the root to cancellation
        _eq(arg, x) ? push!(out, _float_value(val)) : append!(out, _solve_eq(arg ~ val, x))
    end
    out = filter(r -> _satisfies(eq, x, r), out)
    # (roots through a power of x, (…)^6/19683 from x^(1/6): simplified when that makes
    # them smaller, and then the same root once; SymPy #23028)
    common === nothing || (out = _tidy_roots(out))
    out
end

# nu, a polynomial of degree ≥ 2 in u = g with coefficients in x: solved for u first,
# then g = each root (a²/log(b/x + c)² = x²: log(b/x + c) = ±a/x, each through Lambert
# W; SymPy #23253). nothing when nu is not such a polynomial or no root gives solutions.
function _solve_generator_first(eq::Relation, nu::Sym, u::Sym, g::Sym, x::Sym)
    P = try Poly{Sym}(expand(nu), [u]) catch; return nothing end
    degree(P, u) >= 2 || return nothing
    us = try _solve_eq(nu ~ 0, u) catch err; err isa ArgumentError || rethrow(); return nothing end
    (us isa AbstractVector && !(us isa Cases)) || return nothing
    out = Sym[]
    found = false
    for r in us
        rs = try _solve_eq(g ~ r, x) catch err; err isa ArgumentError || rethrow(); continue end
        (rs isa AbstractVector && !(rs isa Cases)) || continue
        found = true
        append!(out, rs)
    end
    found || return nothing
    unique(filter(r -> _satisfies(eq, x, r), out))
end

_iscuberoot(b::Sym) = (b.kind == FN && b.name === :cbrt && length(b.args) == 1) ||
    (b.kind == POW && b.args[2].kind == NUMBER && isexact(b.args[2].num) && _denominator(b.args[2].num) == 3)

function _tidy_roots(rs::Vector{Sym})
    length(rs) <= 12 || return rs
    out = Sym[]
    for r in rs
        if !_hassymbol(r) && _nodes(r) > 8
            s = try simplify(r) catch err; err isa InterruptException && rethrow(); r end
            _nodes(s) < _nodes(r) && (r = s)
        end
        any(q -> _eq(q, r), out) || push!(out, r)
    end
    out
end

# a piecewise function of x in the equation: solved piece by piece, a piece's roots kept
# where its condition holds and no earlier one does (SymPy #27037); nested ones in the
# recursion. Conditions left with parameters are decided (cases). nothing if there is none.
function _solve_piecewise(eq::Relation, x::Sym)
    # a piecewise inside a piecewise's piece: folded into one first (SymPy #18359)
    e0 = eq.lhs - eq.rhs
    if _has_node(e0, t -> t.kind == FN && t.name === :piecewise && any(a -> _hasfn(a, :piecewise), t.args))
        f = try piecewise_fold(e0, x) catch err; err isa ArgumentError || rethrow(); nothing end
        if f isa Sym && !_eq(f, e0)
            return _solve_piecewise(Relation(:eq, f, ZERO), x)
        end
    end
    P = nothing
    _replace_walk(eq.lhs - eq.rhs) do t
        P === nothing && t.kind == FN && t.name === :piecewise && has(t, x) && (P = t)
        nothing
    end
    P === nothing && return nothing
    holds(c, r) = begin
        v = subs(c, x => r)
        v.kind == NUMBER && return !iszero(v)
        _iscond(v) || throw(ArgumentError("solve: cannot decide the condition $v"))
        A = _SOLVE_ASSUME[]
        A === nothing && throw(ArgumentError("solve: cannot decide the condition $v"))
        v.name === :And ? all(a -> holds(a, r), v.args) : _decide(_relation_of(v), A)
    end
    out = Sym[]
    ps = _pieces(P)
    zero = falses(length(ps))
    for (i, (c, v)) in enumerate(ps)
        e = Relation(:eq, subs(eq.lhs, P => v), subs(eq.rhs, P => v))
        # a piece that is the equation identically: all of its interval (SymPy #21079)
        if !has(e.lhs - e.rhs, x) && iszero(simplify(e.lhs - e.rhs))
            zero[i] = true
            continue
        end
        rs = try
            has(e.lhs - e.rhs, x) ? _solve_eq(e, x) : Sym[]
        catch err
            err isa ArgumentError && occursin("every value", sprint(showerror, err)) || rethrow()
            throw(ArgumentError("solve: every $x where $c holds satisfies the equation"))
        end
        rs isa AbstractVector || return nothing
        for r in rs
            holds(c, r) && !any(j -> holds(ps[j][1], r), 1:(i - 1)) || continue
            any(o -> _eq(o, r), out) || push!(out, r)
        end
    end
    any(zero) || return out
    # the intervals where an identically satisfied piece is the one that applies, with
    # the other roots as points: a set (the conditions' breakpoints numbers)
    bps = _breakpoints(P, x)
    (bps === nothing || any(_hassymbol, bps)) &&
        throw(ArgumentError("solve: every $x where $(ps[findfirst(zero)][1]) holds satisfies the equation"))
    pts = [r for r in out if !_hassymbol(r) && _isreal_const(r) === true]
    length(pts) == length(out) ||
        throw(ArgumentError("solve: every $x where $(ps[findfirst(zero)][1]) holds satisfies the equation, and there are non-real roots"))
    function active(v)
        w = v isa Sym ? v : Sym(v)
        for (i, (c, _)) in enumerate(ps)
            t = _truth(subs(c, x => w))
            t == ONE && return i
            t == ZERO || return nothing
        end
        nothing
    end
    S = _set_from(unique(vcat(bps, pts)), (v, atpoint) -> begin
        i = active(v)
        (i !== nothing && zero[i]) || (atpoint && any(o -> _eq(o, v), pts))
    end)
    SymSet(S.intervals, x)
end

# Π bᵢ^uᵢ = Π cⱼ^vⱼ with positive constant bases (and positive constant factors): their
# logarithms, Σ uᵢ log bᵢ = Σ vⱼ log cⱼ (2^x 3^(x²) = 6: x log 2 + x² log 3 = log 6;
# SymPy #27330). Candidates checked in the original equation; nothing if not of this form.
function _solve_exp_bases(eq::Relation, x::Sym)
    # a T₁ + b T₂ = 0 with a, b of opposite signs: a T₁ = -b T₂, each side positive
    # (-2^(3x/11) e^(x/11) + π^(x/11) = 0; SymPy #17604)
    g = eq.lhs - eq.rhs
    if g.kind == ADD && iszero(g.num) && length(g.args) == 2 && !(iszero(eq.rhs) && eq.lhs.kind != ADD)
        (c1, c2), (t1, t2) = g.coefs, g.args
        if isexact(c1) && isexact(c2) && sign(c1) != sign(c2)
            p1, p2 = c1 > N0 ? (Sym(c1) * t1, -Sym(c2) * t2) : (Sym(c2) * t2, -Sym(c1) * t1)
            _eq(p1, eq.lhs) && _eq(p2, eq.rhs) || return _solve_exp_bases(Relation(:eq, p1, p2), x)
        end
    end
    function logsum(e::Sym)
        if e.kind == NUMBER
            (isexact(e.num) || isfloat(e.num)) && e.num > N0 || return nothing
            return _log_primes(e)
        elseif e.kind == POW
            b, p = e.args
            (!_hassymbol(b) && _ispositive_const(b)) || return nothing
            return p * _log_primes(b)
        elseif e.kind == MUL
            e.num > N0 || return nothing
            acc = _log_primes(Sym(e.num))
            for (b, p) in zip(e.args, e.exps)
                (!_hassymbol(b) && _ispositive_const(b)) || return nothing
                acc += p * _log_primes(b)
            end
            return acc
        end
        nothing
    end
    L, R = logsum(eq.lhs), logsum(eq.rhs)
    (L === nothing || R === nothing) && return nothing
    (has(L, x) || has(R, x)) || return nothing
    cands = try _solve_eq(expand(L - R) ~ 0, x) catch err; err isa ArgumentError || rethrow(); return nothing end
    cands isa AbstractVector || return nothing
    filter(r -> _satisfies(eq, x, r), cands)
end

# Σ cᵢ Π b^(α x + β), b positive rationals, α, β rational: each term's x part is r^(kᵢ x)
# for one r when the bases' prime exponents lie on a line (4^x + 6^x = 9^x: with t =
# (3/2)^x, 1 + t = t²; SymPy #24138). A polynomial in t = r^x, its positive roots giving
# x = log(t)/log(r). nothing if not of this form.
function _solve_exp_sum(eq::Relation, n::Sym, x::Sym)
    n = expand(n)
    n.kind == ADD || return nothing
    ts = [Sym(c) * a for (c, a) in zip(n.coefs, n.args)]
    iszero(n.num) || push!(ts, Sym(n.num))
    terms = Tuple{Sym,Dict{BigInt,Rational{BigInt}}}[]
    for t in ts
        fs = t.kind == MUL ? collect(zip(t.args, t.exps)) : t.kind == POW ? [(t.args[1], t.args[2])] : [(t, ONE)]
        c = t.kind == MUL ? Sym(t.num) : ONE
        v = Dict{BigInt,Rational{BigInt}}()
        for (b, p) in fs
            if !has(b, x) && !has(p, x)
                c *= b^p
                continue
            end
            (b.kind == NUMBER && isexact(b.num) && b.num > N0 && !has(b, x)) || return nothing
            P = try Poly{Sym}(expand(p), [x]) catch; return nothing end
            degree(P, x) == 1 || return nothing
            α, β = coeff(expand(p), x, 1), coeff(expand(p), x, 0)
            (α.kind == NUMBER && isexact(α.num) && !_hassymbol(β)) || return nothing
            c *= b^β
            a = Rational{BigInt}(_numerator(α.num), _denominator(α.num))
            for (q, k) in _prime_exponents(b.num)
                v[q] = get(v, q, 0 // 1) + a * k
            end
        end
        _hassymbol(c) && has(c, x) && return nothing
        push!(terms, (c, v))
    end
    length(terms) >= 2 || return nothing
    ps = sort!(unique(reduce(vcat, [collect(keys(v)) for (_, v) in terms])))
    vec(v) = [get(v, q, 0 // 1) for q in ps]
    v1 = vec(terms[1][2])
    ds = [vec(v) - v1 for (_, v) in terms]
    j = findfirst(d -> any(!iszero, d), ds)
    j === nothing && return nothing
    w = ds[j]
    ks = Rational{BigInt}[]
    for d in ds
        i = findfirst(!iszero, w)
        k = d[i] / w[i]
        d == k .* w || return nothing                    # not on one line
        push!(ks, k)
    end
    L = reduce(lcm, (denominator(k) for k in ks); init = big(1))
    ki = [Int(numerator(k * L)) for k in ks]
    m = minimum(ki)
    r = prod((Sym(SymNumber(q))^Sym(SymNumber(wq / L)) for (q, wq) in zip(ps, w)); init = ONE)
    T = Sym(Symbol("_expsum_t"))
    poly = _add_many(c * T^(k - m) for ((c, _), k) in zip(terms, ki))
    roots = try _solve_eq(poly ~ 0, T) catch err; err isa ArgumentError || rethrow(); return nothing end
    roots isa AbstractVector || return nothing
    out = Sym[]
    for τ in roots
        _vanishes(τ) && continue
        v = try N(τ) catch; continue end
        (v isa Real && v > 0) || continue                # real x: t = r^x > 0
        cand = simplify(log(τ) / log(r))
        _satisfies(eq, x, cand) && !any(o -> _eq(o, cand), out) && push!(out, cand)
    end
    out
end

# (prime, exponent) pairs of a positive rational
function _prime_exponents(c::SymNumber)
    out = Pair{BigInt,Int}[]
    for (m, sgn) in ((big(_numerator(c)), 1), (big(_denominator(c)), -1))
        p = big(2)
        while p * p <= m && p <= 10^6                # (a large cofactor kept as one "prime")
            k = 0
            while iszero(m % p)
                m ÷= p
                k += 1
            end
            k > 0 && push!(out, p => sgn * k)
            p += p == 2 ? 1 : 2
        end
        m > 1 && push!(out, m => sgn)
    end
    out
end

# x only through h, the argument of g, elsewhere too (-5s² + 5log(s² + x²) - 5x² + 7: in
# v = s² + x², -5v + 5log(v) + 7): the equation in v, then h = v for each value (SymPy
# #24418). nothing when x does not go away with h = v.
function _solve_inner_generator(eq::Relation, n::Sym, x::Sym, g::Sym)
    (g.kind == FN && length(g.args) == 1 && has(g.args[1], x) && !_eq(g.args[1], x)) || return nothing
    h = g.args[1]
    v = Sym(Symbol("_v", objectid(h)))
    rs = try _solve_eq(h ~ v, x) catch err; err isa ArgumentError || rethrow(); return nothing end
    (rs isa AbstractVector && !isempty(rs)) || return nothing
    m = try simplify(subs(n, x => rs[1])) catch; return nothing end
    has(m, x) && return nothing
    has(m, v) || return nothing
    vs = try _solve_eq(m ~ 0, v) catch err; err isa ArgumentError || rethrow(); return nothing end
    vs isa AbstractVector || return nothing
    out = Sym[]
    for v0 in vs
        xs = try _solve_eq(h ~ v0, x) catch err; err isa ArgumentError || rethrow(); continue end
        for r in xs
            ok = try _satisfies(eq, x, r) catch err; err isa ArgumentError || rethrow(); true end
            ok && !any(o -> _eq(o, r), out) && push!(out, r)
        end
    end
    out
end

_ispositive_const(b::Sym) = (v = try N(b) catch; return false end; v isa Real && v > 0)

# log(c) of a positive rational as Σ k log(p) over its primes (log 6 = log 2 + log 3, so
# the equation's coefficients factor), by trial division up to 10^6; else log(c)
function _log_primes(c::Sym)
    (c.kind == NUMBER && isexact(c.num)) || return log(c)
    n, d = big(_numerator(c.num)), big(_denominator(c.num))
    out = ZERO
    for (m, sgn) in ((n, 1), (d, -1))
        p = big(2)
        while p * p <= m && p <= 10^6
            k = 0
            while iszero(m % p)
                m ÷= p
                k += 1
            end
            k > 0 && (out += sgn * k * log(Sym(p)))
            p += p == 2 ? 1 : 2
        end
        m > 1 && (out += sgn * log(Sym(m)))
    end
    out
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
    # the coefficients over their gcd g: 72/5 log(y - c) - 72/5 log(y) is log((y - c)/y)
    # times 72/5, not a polynomial of degree 72 (an ODE's implicit solution; SymPy #22462)
    rs = Rational{BigInt}[Rational{BigInt}(_numerator(c.num), _denominator(c.num)) for c in cs]
    g = reduce(gcd, rs)
    for (ev, c) in zip(exponent_vectors(P), coefficients(P))
        sum(ev) == 1 && (ks[findfirst(==(1), ev)] = Int(Rational{BigInt}(_numerator(c.num), _denominator(c.num)) / g))
    end
    lhs = prod((logs[i].args[1]^ks[i] for i in eachindex(logs) if ks[i] != 0); init = ONE)
    cands = _solve_eq(lhs ~ exp(-R / Sym(SymNumber(g))), x)
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
        bases = [a.args[1] for a in atoms]
        if all(c -> c.kind == NUMBER && _isint(c.num) && c.num > SymNumber(1), bases)
            # integer bases that are powers of one integer: 4^x = 2^(2x)
            r = _common_integer_base([big(_numerator(c.num)) for c in bases])
            r === nothing && return nothing
            b = Sym(SymNumber(r))
            es = [e * Sym(SymNumber(_ilog(big(_numerator(c.num)), r))) for (e, c) in zip(es, bases)]
        else
            # constant bases that are small rational powers of the first: √(2 - √3) is
            # 1/√(2 + √3) (SymPy #20709)
            any(_hassymbol, bases) && return nothing
            ks = SymNumber[]
            for c in bases
                k = findfirst(k -> iszero(simplify(c - b^Sym(k))), _BASE_POWERS)
                k === nothing && return nothing
                push!(ks, _BASE_POWERS[k])
            end
            es = [e * Sym(k) for (e, k) in zip(es, ks)]
        end
    end
    # exponents with constant parts: b^(x + 1) = b·b^x (SymPy #20199)
    # (only an exponent with symbols: x^(1/3) and √x keep theirs)
    xfree(e) = (e.kind == ADD && _hassymbol(e)) ?
               Sym(e.num) + _add_many(Sym(c) * t for (c, t) in zip(e.coefs, e.args) if !_hassymbol(t)) : ZERO
    ds = [xfree(expand(e)) for e in es]
    es = [expand(e - d) for (e, d) in zip(es, ds)]
    e1 = es[1]
    ratios = [cancel(e / e1) for e in es]
    all(r -> r.kind == NUMBER && isexact(r.num), ratios) || return nothing
    q = reduce(lcm, (big(_denominator(r.num)) for r in ratios); init = big(1))
    g = e1 / Sym(SymNumber(q))
    ks = [Int(_numerator(r.num) * (q ÷ _denominator(r.num))) for r in ratios]
    (pow(b, g), ks, Sym[b^d for d in ds])
end
const _BASE_POWERS = SymNumber[1, -1, 2, -2, 1 // 2, -1 // 2, 3, -3]

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
    if L1 === nothing || L0 === nothing
        # c1 g = -c0 with one side κ(x + β)^k, k ≥ 2: its k-th roots, x + β = ω^j (…)^(1/k)
        # g^(1/k), each a Lambert equation (y² = e^y is y = ±e^(y/2); SymPy #26461)
        return _solve_lambert_power(eq, c1, c0, g, x)
    end
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

function _solve_lambert_power(eq::Relation, c1::Sym, c0::Sym, g::Sym, x::Sym)
    # κ(x + β)^k for a polynomial c of degree k ≥ 2 in x: (κ, β, k), or nothing
    function lin_power(c)
        has(c, x) || return nothing
        P = try Poly{Sym}(c, [x]) catch err; err isa ArgumentError || rethrow(); return nothing end
        k = degree(P, x)
        k >= 2 || return nothing
        κ = Sym(coefficients(P)[1])
        β = cancel(coeff(expand(c), x, k - 1) / (k * κ))
        iszero(expand(c - κ * (x + β)^k)) || return nothing
        κ, β, k
    end
    b, e = g.args
    out = Sym[]
    if (t = lin_power(c0)) !== nothing && !has(c1, x)         # c1 g = -κ (x + β)^k
        κ, β, k = t
        r = (-c1 / κ)^(1 // k)
        for j in 0:(k - 1)
            ω = cos(2PI * j / k) + IM * sin(2PI * j / k)
            append!(out, _solve_eq(x + β - ω * r * b^(e / k) ~ 0, x))
        end
    elseif (t = lin_power(c1)) !== nothing && !has(c0, x)     # κ (x + β)^k g = -c0
        κ, β, k = t
        r = (-c0 / κ)^(1 // k)
        for j in 0:(k - 1)
            ω = cos(2PI * j / k) + IM * sin(2PI * j / k)
            append!(out, _solve_eq((x + β) * b^(e / k) - ω * r ~ 0, x))
        end
    else
        return nothing
    end
    filter(r -> _satisfies(eq, x, r), unique(out))
end

# the assumptions of the solve in progress (abs inverts only for a real argument)
const _SOLVE_ASSUME = ScopedValue{Any}(nothing)

# g(x) = r as equations arg ~ value for an argument of g
function _invert(g::Sym, r::Sym, x::Sym)
    k = g.kind
    if k == FN && length(g.args) == 1
        t = g.args[1]
        name = g.name
        if name === :abs                                    # |t| = r: t = ±r, for real t
            A = _SOLVE_ASSUME[]
            if A !== nothing && _isreal(t, A) === true
                _decide(Relation(:ge, r, ZERO), A) || return Tuple{Sym,Sym}[]
                return iszero(r) ? [(t, r)] : [(t, r), (t, -r)]
            end
            # |b^q| = |b|^q for b, q real: |√x| = r is |x| = r² (SymPy #25779)
            if A !== nothing && t.kind == POW && _isreal(t.args[1], A) === true && _isreal(t.args[2], A) === true
                b, q = t.args
                _decide(Relation(:ge, r, ZERO), A) || return Tuple{Sym,Sym}[]
                return [(abs(b), r^inv(q))]
            end
            # |t| = 0 only at t = 0, in ℂ too (round-11 review); |t| = r > 0 is a curve
            iszero(r) && return [(t, ZERO)]
            throw(ArgumentError("solve: cannot invert abs($t) for $x: in the complex numbers |$t| = $r is a curve; assume $x ∈ ℝ"))
        end
        # sign(t) = 0 exactly where t = 0; another value is a set of t (an inequality)
        if name === :sign
            iszero(r) && return [(t, ZERO)]
            r.kind == NUMBER && !(_isN1(abs(r.num))) && return Tuple{Sym,Sym}[]
            throw(ArgumentError("solve: sign($t) = $r holds on a set of $x (an inequality in $t), not at points"))
        end
        name === :sin && return [(t, a) for a in _periods_in_range(t, _one_period(asin(r), PI - asin(r)), 2PI, x)]
        name === :cos && return [(t, a) for a in _periods_in_range(t, _one_period(acos(r), 2PI - acos(r)), 2PI, x)]
        name === :tan && return [(t, a) for a in _periods_in_range(t, [atan(r)], PI, x)]
        # csc t = r and sec t = r: sin t = 1/r, cos t = 1/r; never 0 (SymPy #17536)
        if name === :csc || name === :sec
            iszero(r) && return Tuple{Sym,Sym}[]
            return _invert(apply_function(name === :csc ? :sin : :cos, Sym[t]), inv(r), x)
        end
        # cot t = r: t = acot r, π/2 for r = 0 (SymPy #19115)
        name === :cot && return [(t, a) for a in _periods_in_range(t, [iszero(r) ? PI / 2 : acot(r)], PI, x)]
        # heaviside(t) = 1/2 at t = 0; 0 or 1 hold on half-lines (SymPy #18939)
        if name === :heaviside
            r == Sym(1 // 2) && return [(t, ZERO)]
            (iszero(r) || isone(r)) && throw(ArgumentError("solve: heaviside($t) = $r holds on a set of $x ($t $(isone(r) ? ">" : "<") 0), not at points; solve the inequality $t $(isone(r) ? "≫" : "≪") 0"))
            r.kind == NUMBER && return Tuple{Sym,Sym}[]
        end
        name === :log && return [(t, exp(r))]
        # W(t) = r: t = r e^r, where r is a value of the principal branch (real ones are
        # ≥ -1, decided by cases; SymPy #12783)
        if name === :lambertw
            A = _SOLVE_ASSUME[]
            if r.kind == NUMBER && isreal(r.num) && r.num < SymNumber(-1)
                return Tuple{Sym,Sym}[]
            elseif A !== nothing && _hassymbol(r) && _isreal(r, A) === true
                _decide(Relation(:ge, r, Sym(-1)), A) || return Tuple{Sym,Sym}[]
            end
            return [(t, r * exp(r))]
        end
        name === :conj && return [(t, conj(r))]               # (SymPy #20305)
        name === :asin && return [(t, sin(r))]
        name === :acos && return [(t, cos(r))]
        # (atan's values lie strictly between -π/2 and π/2: atan(t) = π/2 has no root, as
        # the case of 1/(C + cos x) with a breakpoint there asked; SymPy #23613)
        if name === :atan
            v = try tan(r) catch err; _isdivzero(err) || rethrow(); nothing end
            return v === nothing ? Tuple{Sym,Sym}[] : [(t, v)]
        end
        name === :sinh && return [(t, asinh(r))]
        name === :cosh && return [(t, acosh(r)), (t, -acosh(r))]
        name === :tanh && return [(t, atanh(r))]
    elseif k == POW
        b, e = g.args
        if !has(b, x)                                       # b^t = r: t = log(r)/log(b)
            # a base with parameters: b = 1 and b = 0 are cases (1^t = r and 0^t = r hold
            # for every t or none, not at log(r)/log(b); round-11 review)
            A = _SOLVE_ASSUME[]
            if iszero(b)                                    # 0^t = r
                (iszero(r) || (A !== nothing && _decide(Relation(:eq, r, ZERO), A))) &&
                    throw(ArgumentError("solve: 0^$e = 0 holds wherever the real part of $e is positive, not at points"))
                # 0^0 = 1 (0^x = 1 at x = 0 only; it gave none; SymPy #11536)
                (isone(r) || (A !== nothing && _hassymbol(r) && _decide(Relation(:eq, r, ONE), A))) && return [(e, ZERO)]
                return Tuple{Sym,Sym}[]
            end
            if A !== nothing && _hassymbol(b)
                if _decide(Relation(:eq, b, ONE), A)
                    _decide(Relation(:eq, r, ONE), A) && throw(_identity_error(x, g ~ r))
                    return Tuple{Sym,Sym}[]
                end
                if _decide(Relation(:eq, b, ZERO), A)
                    _decide(Relation(:eq, r, ZERO), A) &&
                        throw(ArgumentError("solve: 0^$e = 0 holds wherever the real part of $e is positive, not at points"))
                    _decide(Relation(:eq, r, ONE), A) && return [(e, ZERO)]     # 0^0 = 1
                    return Tuple{Sym,Sym}[]
                end
            end
            iszero(r) && return Tuple{Sym,Sym}[]
            # a negative base and a real t: b^t is real only at integers, where |b|^t = |r|
            # and the parity gives r's sign ((-2)^x = -8 at x = 3; log(-8)/log(-2) is not
            # real; SymPy #19816)
            A = _SOLVE_ASSUME[]
            if A !== nothing && b.kind == NUMBER && isexact(b.num) && b.num < N0 && r.kind == NUMBER &&
               isexact(r.num) && _isreal(e, A) === true
                n = log(abs(b), abs(r))
                (n.kind == NUMBER && _isint(n.num)) || return Tuple{Sym,Sym}[]
                return b^n == r ? [(e, n)] : Tuple{Sym,Sym}[]
            end
            # (log(b, r): an exact integer when there is one, log(8)/log(2) = 3)
            b === E && return [(e, log(r))]
            n = log(b, r)
            n.kind == NUMBER && return [(e, n)]
            # r a small power of a constant base: √(2 + √3)^t = 2 - √3 at t = -2 (SymPy #20709)
            if !_hassymbol(b) && !_hassymbol(r) && !(b.kind == NUMBER)
                for k in (1, -1, 2, -2, 3, -3, 4, -4, 6, -6)
                    v = try simplify(b^k - r) catch err; err isa InterruptException && rethrow(); continue end
                    iszero(v) && return [(e, Sym(k))]
                end
            end
            return [(e, log(r) / log(b))]
        elseif iszero(r) && has(e, x)                       # b^e = 0: b = 0 (where e > 0,
            return [(b, ZERO)]                              # as the roots are checked; SymPy #20205)
        elseif isone(r) && has(e, x)
            # b^e = 1 with x in both: b = 1, e = 0, or b = -1 at an even e (the roots are
            # checked in the equation; (x² - 7x + 11)^(x² - 13x + 42) = 1; SymPy #15731)
            out = Tuple{Sym,Sym}[]
            add!(v) = any(o -> iszero(expand(o[2] - v)), out) || push!(out, (x, v))
            ok(f) = r -> try f(r) catch err; err isa InterruptException && rethrow(); false end
            for v in _solve_generic(b ~ 1, x)
                add!(v)
            end
            for v in _solve_generic(e ~ 0, x)                # (not where b = 0: 0^0)
                ok(r -> !iszero(simplify(subs(b, x => r))))(v) && add!(v)
            end
            for v in _solve_generic(b ~ -1, x)               # e an even integer there
                ok(r -> (w = simplify(subs(e, x => r)); w.kind == NUMBER && _isint(w.num) && iseven(_toint(w.num))))(v) && add!(v)
            end
            return out
        elseif !has(e, x) && e.kind == NUMBER && isexact(e.num)   # t^(p/q) = r
            return [(b, r^inv(e))]
        elseif !has(e, x)                                   # t^1.11 = r, t^a = r: the principal root
            # an integer exponent n has n roots r^(1/n) e^(2πik/n), not a finite list for a
            # symbolic n: refused rather than the one principal root (SymPy #21014)
            A = _SOLVE_ASSUME[]
            _hassymbol(e) && A !== nothing && _isinteger(e, A) === true && throw(ArgumentError(
                "solve: $g = $r has the roots $(isone(r) ? "" : "($r)^(1/$e) ")exp(2πik/$e) for k = 0, …, |$e| - 1, which depend on $e; substitute a number for $e"))
            # a vanishing exponent decided by cases: t^0 = 1 for every t (x^a = 1 gave [1]
            # alone; SymPy #12774)
            if _hassymbol(e) && A !== nothing && !_decide(Relation(:ne, e, ZERO), A)
                _decide(Relation(:eq, r, ONE), A) && throw(_identity_error(x, Relation(:eq, g, r)))
                return Tuple{Sym,Sym}[]
            end
            return [(b, r^(1 / e))]
        end
    end
    throw(ArgumentError("solve: cannot invert $g for $x"))
end

# the angles as, of period P, of an argument t = k x + c (numbers), when the assumptions
# keep x in a bounded interval: every angle a + nP that t takes there (cos(2x) = 1/2 on
# [0, 2π]: π/6, 5π/6, 7π/6, 11π/6; the one period alone lost two; SymPy #20798); as
# otherwise
function _periods_in_range(t::Sym, as::Vector{Sym}, P::Sym, x::Sym)
    A = _SOLVE_ASSUME[]
    A === nothing && return as
    S = get(A.sets, x, nothing)
    (S === nothing || length(S.intervals) != 1) && return as
    I = only(S.intervals)
    (I.lo === nothing || I.hi === nothing || _hassymbol(I.lo) || _hassymbol(I.hi)) && return as
    any(_hassymbol, as) && return as
    k = diff(t, x)
    if k.kind == NUMBER && !iszero(k.num)
        c = expand(t - k * x)
        _hassymbol(c) && return as
        ends = sort([Float64(N(k * I.lo + c)), Float64(N(k * I.hi + c))])
    else
        # t monotone in x there: the angles between its values at the ends; t unbounded
        # (log x as x → 0+): infinitely many roots, not those of one period (sin(log x) = 0
        # on (0, 1] gave [1] alone; SymPy #17579)
        tl = try limit(t, x => I.lo; dir = :+) catch; return as end
        th = try limit(t, x => I.hi; dir = :-) catch; return as end
        vl, vh = try (N(tl), N(th)) catch; return as end
        (vl isa Real && vh isa Real) || return as
        (isinf(vl) || isinf(vh)) && throw(ArgumentError(
            "solve: $t is unbounded for $x in the assumed interval: infinitely many solutions"))
        # its range: the values at the ends and at the critical points inside
        vals = Float64[vl, vh]
        m = _signs(k, A)
        if !(m !== nothing && (m & _NEG == 0 || m & _POS == 0))
            cs = try _solve_generic(k ~ 0, x) catch; return as end
            cs isa AbstractVector || return as
            lo, hi = Float64(N(I.lo)), Float64(N(I.hi))
            for cp in cs
                _hassymbol(cp) && return as
                w = try _evalc(cp) catch; return as end
                abs(imag(w)) > 1e-12 * max(1, abs(w)) && continue
                lo <= real(w) <= hi || continue
                tv = try N(subs(t, x => cp)) catch; return as end
                tv isa Real || return as
                push!(vals, Float64(tv))
            end
        end
        ends = [minimum(vals), maximum(vals)]
    end
    p = Float64(N(P))
    out = Sym[]
    for a in as
        av = Float64(N(a))
        n1, n2 = ceil(Int, (ends[1] - av) / p - 1e-9), floor(Int, (ends[2] - av) / p + 1e-9)
        n2 - n1 > 1000 && return as
        for n in n1:n2
            w = iszero(n) ? a : a + n * P
            any(o -> isapprox(N(o), N(w); atol = 1e-12), out) || push!(out, w)
        end
    end
    sort!(out; by = N)
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

# a constant with floats in it as its numeric value ((1.02 + 0.07im)^(100/3) as a
# complex float); others as they are
function _float_value(v::Sym)
    (_hassymbol(v) || _isexactex(v)) && return v
    w = try _evalc(v) catch; return v end
    abs(imag(w)) <= 1e-12 * max(1.0, abs(w)) ? Sym(real(w)) : Sym(real(w)) + Sym(imag(w)) * IM
end

# does r satisfy eq? (exactly, or numerically for constants; kept when undecidable)
function _satisfies(eq::Relation, x::Sym, r::Sym)
    # |u| of a real parameter in the root (√(y²(4y² + 4)) is 2|y|√(y² + 1)): decided by
    # cases on u's sign, so each case's check sees only its own parameters (both roots
    # ±|y|/√(y² + 1) of x/√(1 - x²) = y held somewhere, each on one side; SymPy #11534)
    A = _SOLVE_ASSUME[]
    if A !== nothing && _hassymbol(r)
        rr = try _refine(r, A) catch err; err isa ArgumentError || rethrow(); r end
        for u in _abs_args(rr)
            (_hassymbol(u) && !has(u, x) && _isreal(u, A) === true) || continue
            _decide(Relation(:ge, u, ZERO), A)
        end
    end
    d = try
        subs(eq.lhs - eq.rhs, x => r)
    catch err
        _isdivzero(err) && return false
        rethrow()
    end
    iszero(d) && return true
    if !_hassymbol(d)
        # a constant: its value decides, as cancel would (cancel of a quartic's nested
        # radical roots ran for minutes; SymPy #13625)
        v = try N(d) catch; return true end
        return abs(v) <= 1e-9 * max(1, abs(N(r)))
    end
    iszero(cancel(d)) && return true
    if _hassymbol(d)
        _holds_somewhere(d, r) && return true
        # a residual c·u^p vanishes exactly where u = 0: a case of the parameters, not a
        # spurious root (√(z - x) + √(y - x) = 0 has z = x when y = x)
        u = _vanishing_base(d)
        A = _SOLVE_ASSUME[]
        (u === nothing || A === nothing) && return false
        return _decide(Relation(:eq, u, ZERO), A)
    end
    v = try N(d) catch; return true end
    abs(v) <= 1e-9 * max(1, abs(N(r)))
end

# With parameters left, a candidate is dropped only when the equation fails at every one
# of several sample values of the parameters (complex arithmetic, principal branches):
# a root valid for some parameters only (x = a² of sqrt(x) = a, for a ≥ 0) is kept, a
# spurious one from squaring (SymPy #27381) is not. Kept when no sample can be evaluated.
# The samples are exact, so N evaluates d in the precision it needs: with floats,
# 1.38e-23 in an exponent's denominator (E2 - E1 times 7·10²²) cancels to noise (#25792).
const _SAMPLE_VALUES = (37 // 100, -61 // 100, 143 // 100, -189 // 100, 731 // 10000, 271 // 100,
                        -317 // 100, 83 // 100, 53 // 10, -29 // 100)

# the arguments of abs in e
function _abs_args(e::Sym)
    out = Sym[]
    _replace_walk(e) do t
        t.kind == FN && t.name === :abs && length(t.args) == 1 && !any(o -> _eq(o, t.args[1]), out) && push!(out, t.args[1])
        nothing
    end
    out
end

# u when d is c·u^p (c a nonzero number, p > 0) and u a polynomial in its symbols
function _vanishing_base(d::Sym)
    e = d
    for _ in 1:4
        if e.kind == MUL && length(e.args) == 1 && _positive_number(e.exps[1])
            e = e.args[1]
        elseif e.kind == POW && _positive_number(e.args[2])
            e = e.args[1]
        else
            break
        end
    end
    (_eq(e, d) || !_hassymbol(e) || !_ispoly_in_all(e)) ? nothing : e
end
_ispoly_in_all(e::Sym) = try
    Poly{SymNumber}(e, sort!(collect(free_symbols(e)); lt = _lt)); true
catch err
    err isa ArgumentError || rethrow()
    false
end

function _holds_somewhere(d::Sym, r::Sym)
    syms = collect(free_symbols(d))
    tried = false
    # (samples in the sets the assumptions give the symbols: y > 0 assumed, a root valid
    # only for y < 0 is spurious; SymPy #11534)
    A = _SOLVE_ASSUME[]
    pool(s) = (S = A === nothing ? nothing : get(A.sets, s, nothing);
               S === nothing ? collect(_SAMPLE_VALUES) : [v for v in _SAMPLE_VALUES if (try Sym(v) in S catch; true end)])
    pools = [pool(s) for s in syms]
    any(isempty, pools) && return true
    for k in 1:8
        vals = [s => Sym(P[mod1(k + 3j, length(P))]) for (j, (s, P)) in enumerate(zip(syms, pools))]
        v, w = try
            N(subs(d, vals...)), N(subs(r, vals...))
        catch
            continue
        end
        (v isa Number && w isa Number && isfinite(v) && isfinite(w)) || continue
        tried = true
        abs(v) <= 1e-9 * max(1, abs(w)) && return true
    end
    !tried
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
    # a cyclotomic polynomial Φ_N: the primitive N-th roots of unity as cos + i sin of
    # rational multiples of π, not radicals or floats (x^7 + 1; SymPy #25192)
    N = _cyclotomic_index(c)
    N === nothing || return Sym[cos(Sym(2k // N) * PI) + IM * sin(Sym(2k // N) * PI) for k in 1:N if gcd(k, N) == 1]
    # the minimal polynomial of cos(2π/N): its roots cos(2πk/N), not radicals or floats
    # (32x⁵ + 16x⁴ - 32x³ - 12x² + 6x + 1 for N = 11, from an 11×11 circulant; SymPy #21631)
    N = _cos_minpoly_index(c)
    N === nothing || return Sym[cos(Sym(2k // N) * PI) for k in 1:(N ÷ 2) if gcd(k, N) == 1]
    d == 3 && return _cardano(c...)
    d == 4 && return _ferrari(c...)
    r = _solve_decomposed(c)
    r === nothing || return r
    all(ci -> !_hassymbol(ci), c) || throw(ArgumentError(
        "solve: an irreducible factor of degree $d with symbolic coefficients"))
    _numeric_roots(c)
end

# N when the polynomial with coefficients c (low to high) is a multiple of Φ_N, else nothing
function _cyclotomic_index(c::Vector{Sym})
    all(ci -> ci.kind == NUMBER && isexact(ci.num), c) || return nothing
    d = length(c) - 1
    lc = c[end].num
    q = [Rational{BigInt}(_numerator(ci.num / lc), _denominator(ci.num / lc)) for ci in c]
    all(isinteger, q) || return nothing
    for N in 3:(2d^2 + 2)
        _totient(N) == d || continue
        _cyclotomic(N) == q && return N
    end
    nothing
end
_totient(n::Int) = count(k -> gcd(k, n) == 1, 1:n)
# N when the polynomial with coefficients c is a multiple of the minimal polynomial of
# cos(2π/N), N ≥ 7 (degree φ(N)/2 ≥ 3), else nothing
function _cos_minpoly_index(c::Vector{Sym})
    all(ci -> ci.kind == NUMBER && isexact(ci.num), c) || return nothing
    d = length(c) - 1
    d >= 3 || return nothing
    lc = c[end].num
    q = [Rational{BigInt}(_numerator(ci.num / lc), _denominator(ci.num / lc)) for ci in c]
    for N in 7:(4d^2 + 2)
        _totient(N) == 2d || continue
        _cos_minpoly(N) == q && return N
    end
    nothing
end
# the monic minimal polynomial of cos(2π/N), low to high: Φ_N(z) = z^d Ψ(z + 1/z),
# z^j + z^(-j) = D_j(y) (D₀ = 2, D₁ = y, D_j = y D_(j-1) - D_(j-2)), y = 2x
function _cos_minpoly(N::Int)
    a = _cyclotomic(N)
    d = (length(a) - 1) ÷ 2
    padd(p, q) = (r = zeros(Rational{BigInt}, max(length(p), length(q))); r[1:length(p)] .+= p; r[1:length(q)] .+= q; r)
    Ds = Vector{Vector{Rational{BigInt}}}(undef, d + 1)
    Ds[1] = Rational{BigInt}[2]
    d >= 1 && (Ds[2] = Rational{BigInt}[0, 1])
    for j in 2:d
        Ds[j+1] = padd([zero(Rational{BigInt}); Ds[j]], -Ds[j-1])
    end
    psi = Rational{BigInt}[a[d+1]]
    for j in 1:d
        psi = padd(psi, a[d+j+1] .* Ds[j+1])
    end
    px = [psi[i] * big(2)^(i - 1) for i in eachindex(psi)]   # y = 2x
    px ./ px[end]
end

const _CYCLOTOMIC = Dict{Int,Vector{Rational{BigInt}}}()
# (a global cache: filled under a lock, reentrant for the recursion, so concurrent solves
# do not corrupt it; SymPy #10449)
const _CYCLOTOMIC_LOCK = ReentrantLock()
# Φ_N's coefficients, low to high: (x^N - 1) / Π_{d | N, d < N} Φ_d
_cyclotomic(N::Int) = lock(() -> _cyclotomic_locked(N), _CYCLOTOMIC_LOCK)
function _cyclotomic_locked(N::Int)
    get!(_CYCLOTOMIC, N) do
        p = zeros(Rational{BigInt}, N + 1)
        p[1], p[end] = -1, 1
        for d in 1:(N - 1)
            N % d == 0 || continue
            p = _udiv_exact(p, d == 1 ? Rational{BigInt}[-1, 1] : _cyclotomic(d))
        end
        p
    end
end
function _udiv_exact(a::Vector{Rational{BigInt}}, b::Vector{Rational{BigInt}})
    a = copy(a)
    db, da = length(b) - 1, length(a) - 1
    q = zeros(Rational{BigInt}, da - db + 1)
    for i in (da - db):-1:0
        q[i + 1] = a[i + db + 1] / b[end]
        for j in 0:db
            a[i + j + 1] -= q[i + 1] * b[j + 1]
        end
    end
    q
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
    Δ0 = expand(b^2 - 3a * c)
    Δ1 = expand(2b^3 - 9a * b * c + 27a^2 * d)
    # a triple root (Δ0 = Δ1 = 0), and Δ0 = 0, where C = Δ1^(1/3): the general formula
    # would divide 0 by 0 ((x - 2im*y)^3 has Δ0 = Δ1 = 0 only once expanded)
    if iszero(Δ0)
        iszero(Δ1) && return [-b / (3a)]
        C = _real_cbrt(Δ1)
        return [-(b + _OMEGA^k * C) / (3a) for k in 0:2]
    end
    s = sqrt(Δ1^2 - 4Δ0^3)
    C = _real_cbrt((Δ1 + s) / 2)
    _iszero_value(C) && (C = _real_cbrt((Δ1 - s) / 2))
    [-(b + _OMEGA^k * C + Δ0 / (_OMEGA^k * C)) / (3a) for k in 0:2]
end

# u^(1/3), the real cube root -(-u)^(1/3) for a negative number u: the first root is then
# real in its form, not through ω (6√78 - 53 < 0 in 4x³ - 4x² + x - 2; SymPy #26069)
function _real_cbrt(u::Sym)
    _hassymbol(u) && return u^(1 // 3)
    v = try N(u) catch; return u^(1 // 3) end
    v isa Real && v < 0 ? -(expand(-u))^(1 // 3) : u^(1 // 3)
end

function _ferrari(e, d, c, b, a)                             # a x⁴ + b x³ + c x² + d x + e
    p = (8a * c - 3b^2) / (8a^2)
    q = (b^3 - 4a * b * c + 8a^2 * d) / (8a^3)
    shift = -b / (4a)
    # symbolic coefficients in a solve: the degenerate points (q = 0, Δ0 = 0, S = 0) are
    # cases, not skipped; the general formula is wrong there (½x⁴ + x² - 1 from the
    # generic quartic gave ±i twice; SymPy #17563)
    A = _SOLVE_ASSUME[]
    sym = A !== nothing && any(_hassymbol, (a, b, c, d, e))
    iszero_(t) = sym && _hassymbol(t) ? _decide(Relation(:eq, numerator(cancel(t)), ZERO), A) : _iszero_value(t)
    if iszero_(q)                                            # biquadratic in y = x - shift
        r = (-3b^4 + 256a^3 * e - 64a^2 * b * d + 16a * b^2 * c) / (256a^4)
        out = Sym[]
        for y2 in _quadratic(r, p, ONE), sgn in (-1, 1)
            push!(out, shift + sgn * sqrt(y2))
        end
        return out
    end
    Δ0 = c^2 - 3b * d + 12a * e
    Δ1 = 2c^3 - 9b * c * d + 27b^2 * e + 27a * d^2 - 72a * c * e
    if sym && iszero_(Δ0)                                     # Q = Δ1^(1/3), no Δ0/Q
        Q = Δ1^(1 // 3)
        T = Q
    else
        Q = ((Δ1 + sqrt(Δ1^2 - 4Δ0^3)) / 2)^(1 // 3)
        _iszero_value(Q) && (Q = ((Δ1 - sqrt(Δ1^2 - 4Δ0^3)) / 2)^(1 // 3))
        T = Q + Δ0 / Q
    end
    S2 = -2p / 3 + T / (3a)                                  # (2S)²
    for _ in 1:2                                             # S ≠ 0 for some cube root of Q
        iszero_(S2) || break
        Q *= _OMEGA
        T = sym && iszero_(Δ0) ? Q : Q + Δ0 / Q
        S2 = -2p / 3 + T / (3a)
    end
    S = sqrt(S2) / 2
    u, v = sqrt(-4S^2 - 2p + q / S) / 2, sqrt(-4S^2 - 2p - q / S) / 2
    [shift - S - u, shift - S + u, shift + S - v, shift + S + v]
end

## numeric roots: Aberth–Ehrlich iteration

function _numeric_roots(c::Vector{Sym})
    setprecision(BigFloat, 256) do
        C = Complex{BigFloat}
        B = C[(v = _bigvalue(ci, 256); v === nothing ? C(_evalc(ci)) : C(v)) for ci in c]
        while length(B) > 1 && iszero(B[end])
            pop!(B)
        end
        n = length(B) - 1
        n == 0 && return Sym[]
        realcoefs = all(b -> iszero(imag(b)), B)
        m = B ./ B[end]                                      # monic, exactly enough
        # coefficients beyond Float64's range (420-digit ones overflowed to Inf and every
        # root was NaN; SymPy #15355): the iteration in BigFloat; else in Float64
        nzm = [abs(b) for b in m if !iszero(b)]
        wide = any(b -> !(1e-290 < b < 1e290), nzm) || maximum(nzm) / minimum(nzm) > 1e80 || n >= 16
        z = wide ? _aberth(m, 1e-60) : C.(_aberth(ComplexF64.(m), 1e-15))
        # each root polished by Newton's method at 256 bits on the exact coefficients: the
        # roots of (1 + x)^30 - 3/2 expanded were off by 5e-5 (SymPy #15296)
        dm = [k * m[k + 1] for k in 1:n]
        for k in eachindex(z)
            for _ in 1:8                           # (nested: a break ends this root only)
                d = evalpoly(z[k], dm)
                iszero(d) && break
                δ = evalpoly(z[k], m) / d
                z[k] -= δ
                abs(δ) <= big(2.0)^-240 * max(1, abs(z[k])) && break
            end
        end
        # a real polynomial's complex roots in exact conjugate pairs (the last digits of
        # 0.35247154603172637 and …626 kept dsolve from pairing them)
        if realcoefs
            used = falses(n)
            for k in 1:n
                (used[k] || imag(z[k]) <= 0) && continue
                js = [j for j in 1:n if !used[j] && j != k && imag(z[j]) < 0]
                isempty(js) && continue
                j = js[argmin([abs(z[j] - conj(z[k])) for j in js])]
                abs(z[j] - conj(z[k])) <= big(10.0)^-30 * max(1, abs(z[k])) || continue
                z[j] = conj(z[k]); used[j] = used[k] = true
            end
        end
        map(z) do t
            # (exact coefficients, irreducible: simple roots, polished to ~1e-70; floats may
            # carry a double root, whose imaginary residue stays near √eps)
            tol = realcoefs && all(_isexactex, c) ? big(10.0)^-40 : big(10.0)^-12
            r, i = Float64(real(t)), Float64(imag(t))
            abs(imag(t)) <= tol * max(1, abs(t)) ? Sym(r) : Sym(r) + Sym(i) * IM
        end
    end
end

# Aberth–Ehrlich iteration for the roots of the monic polynomial with coefficients a
# (constant first), to relative steps below tol
function _aberth(a::Vector{T}, tol) where {T}
    n = length(a) - 1
    # Fujiwara's root bound, 2 max |aₙ₋ₖ|^(1/k): Cauchy's, 1 + max |aₖ|, is far larger for a
    # high degree, and its n-th power overflows (u^181 + u^74 - 50 gave NaN)
    R = 2 * maximum(abs(a[n + 1 - k])^(1 / k) for k in 1:n)
    R > 0 || (R = one(R))
    z = T[R * cis(2 * real(T)(π) * (k + 0.25) / n) for k in 0:(n - 1)]
    da = [k * a[k+1] for k in 1:n]
    for _ in 1:500
        δmax = zero(real(T))
        for k in 1:n
            w = evalpoly(z[k], a) / evalpoly(z[k], da)
            s = sum(1 / (z[k] - z[j]) for j in 1:n if j != k; init = zero(T))
            δ = w / (1 - w * s)
            z[k] -= δ
            δmax = max(δmax, abs(δ) / max(1, abs(z[k])))
        end
        δmax < tol && break
    end
    z
end

# numeric (constant) roots sorted: real ones first, by real then imaginary part
function _sort_roots!(roots::Vector{Sym})
    any(_hassymbol, roots) && return roots
    vals = map(r -> (v = try _evalc(r) catch; return nothing end; v), roots)
    any(isnothing, vals) && return roots
    key(v) = (abs(imag(v)) > 1e-10 * max(1, abs(v)), real(v), imag(v))
    roots[sortperm([key(v) for v in vals])]
end
