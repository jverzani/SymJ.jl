# Ordinary differential equations: SymJ.dsolve.
#
# The equation is read as F(x, y, y', …, y⁽ⁿ⁾) = 0 with placeholders Y₀, …, Yₙ for the
# derivatives, and the methods are tried in turn:
#
# - a missing y (and its lower derivatives): v = y⁽ʲ⁾ solved, then integrated j times;
# - linear (F affine in the Yₖ): first order by an integrating factor; constant
#   coefficients by the roots of the characteristic polynomial, a particular solution by
#   undetermined coefficients or variation of parameters; Cauchy–Euler through x = eᵗ;
# - first order, y' = f(x, y): separable, Bernoulli, exact (or with an integrating factor
#   in x or in y alone), homogeneous (f a function of y/x), Riccati (made linear, or
#   from a particular solution);
# - autonomous, y'' = f(y, y'): p = y' as a function of y, then ∫ dy/p = x + C₂.
#
# With `series = n`, a series solution instead ("series solutions", at the end); with
# method = :laplace, or heaviside or dirac on the right side, Laplace transforms
# (laplace.jl).
#
# Systems, linear with constant coefficients, are at the end ("systems").
#
# A solution is explicit, y = S(x, C₁, …), or implicit, φ(x, y) = C₁. It is checked by
# substitution: symbolically, else numerically at a few points.

"""
    SymJ.dsolve(eq, y(x), ics...; assume = nothing)

Solve the ordinary differential equation `eq` (an equation, written with `~`) for the
function `y(x)`. Derivatives are written `diff(y(x), x)` or `y'(x)`. The result is an
equation: `y(x) ~ S` with constants `C1, C2, …` (one for each order), or an implicit
solution `φ(x, y(x)) ~ C1` when `y` cannot be solved for. An implicit solution quadratic
in `y` gives a vector, one equation for each branch: `y'(x) ~ (1 - x)/y(x)` gives
`y(x) ~ ±sqrt(C1 + 2x - x^2)`.

Initial conditions are pairs, one for each order: `y(0) => 1, y'(0) => 0`. They fix the
constants.

With assumptions (`assume = (k ≫ 0, m ≫ 0)`, an [`Assumptions`](@ref), or inside
[`assuming`](@ref)), a characteristic root they show to be complex gives cosines and
sines: `y'' = -(k/m) y` is `C1 cos(√(km) x/m) + C2 sin(√(km) x/m)` with `k, m > 0`, and
`C1 exp(-√(-km) x/m) + C2 exp(√(-km) x/m)` without.

An autonomous equation `y'' = f(y, y')` whose `y'` has several branches (`y' =
±sqrt(…)`, as for `y'' = -sin(y)`) gives a vector of solutions, one for each branch;
with initial conditions, the branch through them (`∫ᵥʸ dy/y' = x - a`, a definite
`Integral(f, (y, v, y(x)))`, when the integral has no closed form).

With `method = :laplace`, a linear equation with constant coefficients is solved by
Laplace transforms ([`SymJ.laplace_transform`](@ref)); this is the method used, without
asking, when the right side has [`SymJ.heaviside`](@ref) or [`SymJ.dirac`](@ref). The
initial conditions (constants for those not given) are at one point, and the solution
holds from there on.

With `series = n`, the solution is a series about `x0` (the point of the initial
conditions, else `0`) with error `O((x - x0)^n)`, returned as the truncated sum. At an
ordinary point it is the Taylor series, for any equation that can be solved for its
highest derivative, nonlinear ones included; a linear one is grouped by constant,
`C1*(…) + C2*(…)`. At a regular singular point of a second-order linear homogeneous
equation it is the Frobenius series, with a `log` term when the indicial roots are equal
or an integer apart.

Methods: first order linear, separable, Bernoulli, exact (with an integrating factor in
`x` or `y` alone), homogeneous, Riccati; linear with constant coefficients, of any order
(undetermined coefficients, else variation of parameters, for the right side);
Cauchy–Euler; equations without `y` (reduced to a lower order); second order without
`x` (autonomous: `p = y'` as a function of `y`). Integrals are done by
the RischIntegrate extension when it is loaded (an unevaluated `Integral` otherwise, or
when there is no elementary antiderivative).

Not exported.

```julia
@syms x y()
SymJ.dsolve(y'(x) ~ y(x), y(x))                          # y(x) ~ C1*exp(x)
SymJ.dsolve(y''(x) + y(x) ~ 0, y(x), y(0) => 0, y'(0) => 1)   # y(x) ~ sin(x)
SymJ.dsolve(y'(x) ~ x * y(x)^2, y(x))                    # y(x) ~ -2/(x^2 + 2C1)
SymJ.dsolve(x^2 * y''(x) - 2y(x) ~ 0, y(x))              # y(x) ~ C1/x + C2*x^2
SymJ.dsolve(y''(x) - x * y(x) ~ 0, y(x); series = 6)      # C1*(1 + x^3/6) + C2*(x + x^4/12)
```

    SymJ.dsolve([eq₁, eq₂, …], [y₁(t), y₂(t), …], ics...)

A system of linear equations with constant coefficients, of any order, in normal form
(the highest derivatives can be solved for). The result is a vector of equations, with
as many constants as the unknowns' orders add up to; the initial conditions fix them, one
for each unknown and each order below its highest. A right side `F(t)` is integrated in
closed form when it is made of polynomials, exponentials, sines and cosines.

```julia
@syms t u() v()
SymJ.dsolve([u'(t) ~ v(t), v'(t) ~ -u(t)], [u(t), v(t)])
# [u(t) ~ C1*sin(t) + C2*cos(t), v(t) ~ C1*cos(t) - C2*sin(t)]
SymJ.dsolve([u'(t) ~ v(t), v'(t) ~ -u(t)], [u(t), v(t)], u(0) => 1, v(0) => 0)
# [u(t) ~ cos(t), v(t) ~ -sin(t)]
```
"""
function dsolve(eq::Relation, yx::Sym, ics::Pair...; series::Union{Nothing,Integer} = nothing, x0 = nothing,
                method::Symbol = :auto, assume = nothing)
    assume === nothing || return with(() -> dsolve(eq, yx, ics...; series, x0, method), _ASSUMING => _assumptions(assume))
    method in (:auto, :laplace) || throw(ArgumentError("dsolve: method is :auto or :laplace, not $method"))
    P = _ode_parse(eq, yx)
    Cs = _ode_constants(P.n, eq)
    series === nothing || return yx ~ _ode_series(P, Cs, series, x0, [_ode_ic(p, P) for p in ics])
    # Laplace transforms: asked for, or a right side with heaviside or dirac
    if method === :laplace || (_hasfn(P.F, :heaviside) || _hasfn(P.F, :dirac))
        conds = [_ode_ic(p, P) for p in ics]
        isempty(ics) || length(ics) == P.n || throw(ArgumentError("dsolve: give $(P.n) initial condition$(P.n == 1 ? "" : "s"), one for each order"))
        S = try _ode_laplace(P, Cs, conds) catch err; err isa ArgumentError || rethrow(); nothing end
        if S !== nothing
            _ode_verified(P, :explicit, S, Cs) || throw(ArgumentError("dsolve: the solution found does not satisfy $(eq)"))
            return yx ~ S
        end
        method === :laplace && throw(ArgumentError("dsolve: no Laplace transform solution for $(eq) (linear, constant coefficients, transformable right side)"))
    end
    sol = _ode_solve(P, Cs)
    if (sol === nothing || sol[1] === :several) && !isempty(ics)
        # an autonomous equation with several branches of y': the constants and the
        # branch from the initial conditions
        length(ics) == P.n || throw(ArgumentError("dsolve: give $(P.n) initial condition$(P.n == 1 ? "" : "s"), one for each order"))
        bs = _ode_auto_branches(P, Cs; conds = [_ode_ic(p, P) for p in ics])
        if bs !== nothing
            rs = [k === :explicit ? (yx ~ S_) : (subs(S_, P.ys[1] => yx) ~ c_) for (k, S_, c_) in bs if _ode_verified(P, k, S_, Cs)]
            isempty(rs) && throw(ArgumentError("dsolve: the solution found does not satisfy $(eq)"))
            return length(rs) == 1 ? rs[1] : rs
        end
    end
    sol === nothing && throw(ArgumentError("dsolve: no method for $(eq)"))
    if sol[1] === :several                       # one solution for each branch
        rs = Relation[]
        for (k, S_) in sol[2]
            k === :explicit && (S_ = _ode_absorb(S_, Cs, P.x))
            _ode_verified(P, k, S_, Cs) || throw(ArgumentError("dsolve: the solution found does not satisfy $(eq)"))
            push!(rs, k === :explicit ? (yx ~ S_) : (subs(S_, P.ys[1] => yx) ~ Cs[P.n]))
        end
        return rs
    end
    kind, S = sol
    kind === :explicit && (S = _ode_absorb(S, Cs, P.x))
    _ode_verified(P, kind, S, Cs) || throw(ArgumentError("dsolve: the solution found does not satisfy $(eq)"))
    c = Cs[P.n]                                  # an implicit solution is φ = c
    if !isempty(ics)
        length(ics) == P.n || throw(ArgumentError("dsolve: give $(P.n) initial condition$(P.n == 1 ? "" : "s"), one for each order"))
        conds = [_ode_ic(p, P) for p in ics]
        r = try
            _ode_ics(P, kind, S, Cs, conds)
        catch err
            err isa ArgumentError || rethrow()
            # an autonomous equation's explicit solution can give a nonlinear system for
            # the constants; its implicit form, ∫ dy/p(y) - x = C₂, a simpler one
            alt = kind === :explicit ? (try _ode_autonomous(P, Cs; implicit = true) catch; nothing end) : nothing
            alt === nothing || alt[1] !== :implicit || !_ode_verified(P, alt..., Cs) ? rethrow() :
                _ode_ics(P, alt..., Cs, conds)
        end
        kind, S, c = r
    end
    if kind === :implicit && isempty(ics) && P.n == 1
        bs = _ode_quadratic_branches(P, S, c, Cs)
        bs === nothing || return [yx ~ b for b in bs]
    end
    kind === :explicit ? (yx ~ S) : (subs(S, P.ys[1] => yx) ~ c)
end

# an implicit solution φ = C quadratic in y: its two branches y = …, each a solution
function _ode_quadratic_branches(P, S::Sym, c::Sym, Cs)
    Y0 = P.ys[1]
    d2 = diff(S, Y0, 2)
    (iszero(d2) || has(d2, Y0)) && return nothing
    ys = try solve(S ~ c, Y0) catch; return nothing end
    (ys isa AbstractVector && length(ys) == 2) || return nothing
    bs = [_ode_absorb(_ode_tidy(y), Cs, P.x) for y in ys]
    all(b -> _ode_verified(P, :explicit, b, Cs), bs) ? bs : nothing
end
dsolve(eq, yx::Sym, ics::Pair...) = throw(ArgumentError("dsolve: give an equation, `lhs ~ rhs`"))

struct _ODE
    F::Sym              # the equation, F = 0, in the placeholders
    ys::Vector{Sym}     # Y₀ = y(x), Y₁ = y'(x), …
    x::Sym
    yx::Sym
    n::Int              # the order
end

function _ode_parse(eq::Relation, yx::Sym)
    (yx.kind == FN && length(yx.args) == 1 && yx.args[1].kind == SYMBOL) ||
        throw(ArgumentError("dsolve: the unknown is a function of one symbol, y(x); got $yx"))
    x = yx.args[1]
    F = lhs(eq) - rhs(eq)
    ders = Dict{Sym,Int}()
    order(e) = _eq(e, yx) ? 0 : (e.kind == FN && e.name === :Derivative && length(e.args) == 2 && _eq(e.args[2], x) &&
                                 (o = order(e.args[1]); o === nothing ? nothing : o + 1))
    function walk(e)
        o = e.kind == FN ? order(e) : nothing
        if o isa Int
            ders[e] = o
            return
        end
        foreach(walk, e.args)
        foreach(walk, e.exps)
    end
    walk(F)
    isempty(ders) && throw(ArgumentError("dsolve: $yx does not appear in $eq"))
    n = maximum(values(ders))
    n >= 1 || throw(ArgumentError("dsolve: no derivative of $yx in $eq"))
    ys = [Sym(Symbol("_ode_y", k)) for k in 0:n]
    Fs = subs(F, (d => ys[o + 1] for (d, o) in ders)...)
    _hasfn(Fs, yx.name) && throw(ArgumentError("dsolve: $(yx.name) appears other than as $yx and its derivatives"))
    _ODE(Fs, ys, x, yx, n)
end

# C1, …, Cn, skipping names in eq (or in the equations)
_ode_constants(n::Int, eq::Relation) = _ode_constants(n, [eq])
function _ode_constants(n::Int, eqs::AbstractVector)
    used = Set(s.name for eq in eqs for s in free_symbols(lhs(eq) - rhs(eq)))
    out = Sym[]
    i = 1
    while length(out) < n
        nm = Symbol("C", i)
        nm in used || push!(out, Sym(nm))
        i += 1
    end
    out
end

# the residual of a solution: explicit y = S, or implicit φ = C (y' = -φₓ/φ_y, and the
# higher derivatives by differentiating along the curve)
function _ode_residual(P::_ODE, kind::Symbol, S::Sym)
    if kind === :explicit
        return subs(P.F, (P.ys[k + 1] => diff(S, P.x, k) for k in 0:P.n)...)
    end
    subs(P.F, (P.ys[k + 1] => d for (k, d) in enumerate(_ode_implicit_derivs(P, S)))...)
end

# y', …, y⁽ⁿ⁾ on the curve φ(x, y) = C, as functions of x and y
function _ode_implicit_derivs(P::_ODE, S::Sym)
    x, Y0 = P.x, P.ys[1]
    yp = _ode_tidy(-diff(S, x) / diff(S, Y0))
    ds = [yp]
    for _ in 2:P.n
        push!(ds, _ode_tidy(diff(ds[end], x) + yp * diff(ds[end], Y0)))
    end
    ds
end

function _ode_verified(P::_ODE, kind::Symbol, S::Sym, Cs)
    r = _ode_residual(P, kind, S)
    z = try simplify(r) catch; r end
    iszero(z) && return true
    vars = free_symbols(r)
    ok = _ode_numzero(r, vars)
    ok === nothing ? true : ok          # no numeric value (an unevaluated integral): accepted
end

# r at a few rational points: true if ≈ 0 at all, false if clearly not, nothing if no
# point gives a value
function _ode_numzero(r::Sym, vars)
    seen = false
    for t in 1:3
        vals = [v => Sym((7i + 5t) // (3 + i + t)) for (i, v) in enumerate(vars)]
        v = try _evalc(isempty(vals) ? r : subs(r, vals...)) catch; continue end
        (v isa Number && isfinite(v)) || continue
        seen = true
        abs(v) < 1e-8 || return false
    end
    seen ? true : nothing
end

## solving

function _ode_solve(P::_ODE, Cs)
    for m in (_ode_reduce, _ode_linear, _ode_first_order, _ode_autonomous)
        r = try m(P, Cs) catch err; err isa InterruptException && rethrow(); nothing end
        r === nothing || return r
    end
    nothing
end

# exp(a + c log(u)) = u^c exp(a), and the solution simplified where smaller
function _ode_tidy(e::Sym)
    e = _ode_signs(_explog(e))
    s = try simplify(e) catch; e end
    _nodes(s) <= _nodes(e) ? s : e
end
# a product's sums with every coefficient negative, to an odd power, negated with it:
# -u/(-v - 1) is u/(v + 1)
function _ode_signs(e::Sym)
    e.kind == MUL || return e
    neg(b, p) = b.kind == ADD && b.num <= 0 && all(<(0), b.coefs) && p.kind == NUMBER && _isint(p.num) && isodd(_toint(p.num))
    count(i -> neg(e.args[i], e.exps[i]), eachindex(e.args)) > 0 || return e
    Sym(e.num) * prod((neg(b, p) ? (-b)^p * (-1) : b^p for (b, p) in zip(e.args, e.exps)); init = ONE)
end
function _explog(e::Sym)
    isempty(e.args) && return e
    e.kind == POW && _eq(e.args[1], E) && return _explog_exp(e.args[2])
    e.kind == MUL && return Sym(e.num) * prod((_eq(b, E) ? _explog_exp(p) : _explog(b)^_explog(p) for (b, p) in zip(e.args, e.exps)); init = ONE)
    _rebuild_node(e, Sym[_explog(a) for a in e.args])
end
function _explog_exp(p::Sym)
    begin
        u = _explog(p)
        terms = u.kind == ADD ? [(Sym(c), t) for (c, t) in zip(u.coefs, u.args)] : [(ONE, u)]
        num = u.kind == ADD ? Sym(u.num) : ZERO
        out, rest = ONE, num
        for (c, t) in terms
            l = t.kind == MUL && length(t.args) == 1 && _eq(t.exps[1], ONE) ? (Sym(t.num), t.args[1]) : (ONE, t)
            if l[2].kind == FN && l[2].name === :log && length(l[2].args) == 1
                out *= l[2].args[1]^(c * l[1])
            else
                rest += c * t
            end
        end
        return out * exp(rest)
    end
end

# a solution with its constants absorbed, where C appears once: k*C (k a number) is C,
# and e^C (in an exponent) is C. Given x, C in several exponents, each k*C + (free of C)
# with one k free of x, is C as well: e^(kC) is C
function _ode_absorb(S::Sym, Cs, x = nothing)
    for C in Cs
        x === nothing || _ode_count(S, C) == 1 || (S = _ode_absorb_exp(S, C, x))
        _ode_count(S, C) == 1 || continue
        k = _ode_ccoef(S, C)
        k === nothing || isone(k) || (S = subs(S, C => C / Sym(k)))
        if _ode_inexp(S, C)
            T = _explog(subs(S, C => log(C)))
            _ode_count(T, C) == 1 && !_ode_haslog(T, C) && (S = _ode_tidy(T))
        end
    end
    S
end
function _ode_absorb_exp(S::Sym, C::Sym, x::Sym)
    ks = Sym[]
    (_ode_expcoefs!(ks, S, C) && !isempty(ks)) || return S
    k = ks[1]
    (has(k, x) || iszero(k) || !all(kk -> iszero(_ode_tidy(kk - k)), ks)) && return S
    T = _ode_tidy(_explog(_ode_expandexp(subs(S, C => log(C) / k))))
    _ode_haslog(T, C) ? S : T
end
# the exponents of e expanded, so e^((log(C)/k + x) k) is C e^(k x)
function _ode_expandexp(e::Sym)
    isempty(e.args) && return e
    e.kind == POW && _eq(e.args[1], E) && return exp(expand(e.args[2]))
    e.kind == MUL && return Sym(e.num) * prod((_eq(b, E) ? exp(expand(p)) : _ode_expandexp(b)^p for (b, p) in zip(e.args, e.exps)); init = ONE)
    _rebuild_node(e, Sym[_ode_expandexp(a) for a in e.args])
end
# the coefficients k of C in the exponents of e that have it, each linear in C; false
# where C is elsewhere
function _ode_expcoefs!(ks, e::Sym, C::Sym)
    _eq(e, C) && return false
    function expo!(p)
        has(p, C) || return true
        k = diff(p, C)
        has(k, C) && return false
        push!(ks, k)
        true
    end
    if e.kind == POW && _eq(e.args[1], E)
        return expo!(e.args[2])
    elseif e.kind == MUL
        for (b, p) in zip(e.args, e.exps)
            _eq(b, E) ? (expo!(p) || return false) : (_ode_expcoefs!(ks, b, C) && _ode_expcoefs!(ks, p, C) || return false)
        end
        return true
    end
    all(a -> _ode_expcoefs!(ks, a, C), e.args) && all(a -> _ode_expcoefs!(ks, a, C), e.exps)
end
_ode_count(e::Sym, C::Sym) = _eq(e, C) ? 1 : sum(a -> _ode_count(a, C), e.args; init = 0) + sum(a -> _ode_count(a, C), e.exps; init = 0)
function _ode_ccoef(e::Sym, C::Sym)
    if e.kind == ADD
        for (c, t) in zip(e.coefs, e.args)
            (_eq(t, C) || (t.kind == MUL && any(i -> _eq(t.args[i], C) && _eq(t.exps[i], ONE), eachindex(t.args)))) && return c
        end
    elseif e.kind == MUL && length(e.args) == 1 && _eq(e.args[1], C) && _eq(e.exps[1], ONE)
        return e.num
    end
    for a in Iterators.flatten((e.args, e.exps))
        k = _ode_ccoef(a, C)
        k === nothing || return k
    end
    nothing
end
_ode_inexp(e::Sym, C::Sym) = (e.kind == POW && _eq(e.args[1], E) && has(e.args[2], C)) ||
    (e.kind == MUL && any(i -> _eq(e.args[i], E) && has(e.exps[i], C), eachindex(e.args))) ||
    any(a -> _ode_inexp(a, C), e.args) || any(a -> _ode_inexp(a, C), e.exps)
_ode_haslog(e::Sym, C::Sym) = (e.kind == FN && e.name === :log && has(e.args[1], C)) ||
    any(a -> _ode_haslog(a, C), e.args) || any(a -> _ode_haslog(a, C), e.exps)

# y missing: v = y⁽ʲ⁾, the equation of order n - j in v solved, then v integrated j times
function _ode_reduce(P::_ODE, Cs)
    j = findfirst(k -> has(P.F, P.ys[k + 1]), 0:P.n) - 1
    j == 0 && return nothing
    j == P.n && return _ode_integrate_only(P, Cs)
    m = P.n - j
    Q = _ODE(subs(P.F, (P.ys[k + 1] => P.ys[k - j + 1] for k in j:P.n)...), P.ys[1:(m + 1)], P.x, P.yx, m)
    sol = _ode_solve(Q, Cs[1:m])
    (sol === nothing || sol[1] === :several) && return nothing
    kind, v = sol
    if kind !== :explicit
        vs = try solve(v ~ Cs[1], Q.ys[1]) catch; Sym[] end
        length(vs) == 1 || return nothing
        v = vs[1]
    end
    for i in 1:j
        v = _integrate(v, P.x) + Cs[m + i]
    end
    (:explicit, _ode_tidy(v))
end

# y⁽ⁿ⁾ = f(x): integrated n times
function _ode_integrate_only(P::_ODE, Cs)
    Yn = P.ys[end]
    a = diff(P.F, Yn)
    has(a, Yn) && return nothing
    f = -subs(P.F, Yn => ZERO) / a
    has(f, Yn) && return nothing
    v = _ode_tidy(f)
    for i in 1:P.n
        v = _integrate(v, P.x) + Cs[i]
    end
    (:explicit, _ode_tidy(v))
end

## linear equations

# the coefficients aₖ of a linear equation Σ aₖ y⁽ᵏ⁾ = g, and g; nothing if not linear
function _ode_linear_coefs(P::_ODE)
    a = Sym[diff(P.F, Y) for Y in P.ys]
    any(c -> any(Y -> has(c, Y), P.ys), a) && return nothing
    g = -subs(P.F, (Y => ZERO for Y in P.ys)...)
    a, g
end

function _ode_linear(P::_ODE, Cs)
    lc = _ode_linear_coefs(P)
    lc === nothing && return nothing
    a, g = lc
    x = P.x
    if P.n == 1                                    # y' + p y = q: integrating factor
        p, q = _ode_tidy(a[1] / a[2]), _ode_tidy(g / a[2])
        μ = _ode_tidy(exp(_integrate(p, x)))
        S = (_integrate(_ode_tidy(μ * q), x) + Cs[1]) / μ
        return (:explicit, _ode_tidy(S))
    end
    if all(c -> !has(c, x), a)                     # constant coefficients
        return (:explicit, _ode_constcoeff(a, g, x, Cs))
    end
    # Cauchy–Euler: aₖ = cₖ xᵏ (up to a common factor), through x = eᵗ
    n = P.n
    cs = Sym[_ode_tidy(a[k + 1] * x^(n - k) / a[n + 1]) for k in 0:n]
    any(c -> has(c, x), cs) && return nothing
    t = Sym(:_ode_t)
    # x^k y⁽ᵏ⁾ = θ(θ - 1)…(θ - k + 1) y, θ = d/dt: the coefficients of the polynomial in θ
    r = Sym(:_ode_r)
    q = expand(sum(cs[k + 1] * prod((r - i for i in 0:(k - 1)); init = ONE) for k in 0:n))
    b = Sym[_ode_coeff(q, r, k) for k in 0:n]
    G = _explog(subs(_ode_tidy(g * x^n / a[n + 1]), x => exp(t)))
    S = _ode_constcoeff(b, G, t, Cs)
    (:explicit, _ode_tidy(_explog(subs(S, t => log(x)))))
end

_ode_coeff(q::Sym, r::Sym, k::Int) = (d = diff(q, r, k); subs(d, r => ZERO) / factorial(big(k)))

# Σ aₖ y⁽ᵏ⁾ = g with constant aₖ: the homogeneous solutions from the characteristic roots
# (e^(αx) cos βx and e^(αx) sin βx for a conjugate pair), times x^j for a multiple root,
# and a particular solution
function _ode_constcoeff(a::Vector{Sym}, g::Sym, x::Sym, Cs)
    basis = _ode_basis(a, x)
    yh = sum(C * u for (C, u) in zip(Cs, basis))
    iszero(g) && return _ode_tidy(yh)
    yp = _ode_undetermined(a, g, x)
    yp === nothing && (yp = _ode_variation(a, g, x, basis))
    _ode_tidy(yh + yp)
end

# the characteristic roots with multiplicities
function _ode_roots(a::Vector{Sym})
    r = Sym(:_ode_r)
    p = sum(a[k + 1] * r^k for k in 0:(length(a) - 1))
    fs = _ode_factors(factor(p))
    out = Tuple{Sym,Int}[]
    for (b, e) in fs
        has(b, r) || continue
        m = e.kind == NUMBER && _isint(e.num) ? _toint(e.num) : return nothing
        rs = solve(b ~ 0, r)
        length(rs) == degree(Poly{Sym}(expand(b), [r]), r) || throw(ArgumentError("dsolve: roots of $b not found"))
        append!(out, ((ρ, m) for ρ in rs))
    end
    sort!(out; by = t -> string(t[1]))          # the same order on every Julia version
end

function _ode_basis(a::Vector{Sym}, x::Sym)
    roots = _ode_roots(a)
    A = _ASSUMING[]
    basis = Sym[]
    done = Set{Sym}()
    pairs = Tuple{Sym,Sym}[]
    for (ρ, m) in roots
        ρ in done && continue
        re, im_ = _ode_reim(ρ, A)
        if im_ !== nothing && !iszero(im_)
            # a conjugate pair, once
            β = _ode_tidy(abs(im_))
            A === nothing || (β = refine(β; assume = A))
            any(p -> _eq(p[1], re) && _eq(p[2], β), pairs) && continue
            push!(done, ρ, conj(ρ))
            push!(pairs, (re, β))
            for j in 0:(m - 1)
                push!(basis, x^j * exp(re * x) * cos(β * x), x^j * exp(re * x) * sin(β * x))
            end
        else
            push!(done, ρ)
            append!(basis, (x^j * exp(ρ * x) for j in 0:(m - 1)))
        end
    end
    length(basis) == length(a) - 1 || throw(ArgumentError("dsolve: the characteristic roots were not all found"))
    basis
end

# the real and imaginary parts of a root, when it has numbers only (decidable), or when
# the assumptions decide them: sqrt(u) with u < 0 is i sqrt(-u), and the parts are then
# real
function _ode_reim(ρ::Sym, A = nothing)
    isempty(free_symbols(ρ)) && return _ode_tidy(real(ρ)), _ode_tidy(imag(ρ))
    A === nothing && return ρ, nothing
    r = expand(_ode_isqrt(ρ, A))
    re = subs(r, IM => ZERO)
    im_ = expand((r - re) / IM)
    has(im_, IM) && return ρ, nothing
    _isreal(re, A) === true && _isreal(im_, A) === true || return ρ, nothing
    _ode_tidy(re), _ode_tidy(im_)
end
function _ode_isqrt(e::Sym, A)
    isempty(e.args) && return e
    if e.kind == POW && _eq(e.args[2], Sym(1 // 2)) && _signs(e.args[1], A) === _NEG
        return IM * sqrt(-e.args[1])
    elseif e.kind == MUL
        return Sym(e.num) * prod((_eq(p, Sym(1 // 2)) && _signs(b, A) === _NEG ? IM * sqrt(-b) : _ode_isqrt(b, A)^p
                                  for (b, p) in zip(e.args, e.exps)); init = ONE)
    end
    _rebuild_node(e, Sym[_ode_isqrt(a, A) for a in e.args])
end

# the multiplicity of ρ as a characteristic root
function _ode_multiplicity(a::Vector{Sym}, ρ::Sym)
    r = Sym(:_ode_r)
    p = sum(a[k + 1] * r^k for k in 0:(length(a) - 1))
    s = 0
    while s < length(a) && iszero(_ode_tidy(subs(p, r => ρ)))
        s += 1
        p = diff(p, r)
    end
    s
end

# undetermined coefficients: g a sum of c xᵐ e^(αx) cos(βx) (or sin); nothing otherwise
function _ode_undetermined(a::Vector{Sym}, g::Sym, x::Sym)
    g = expand(g)
    terms = g.kind == ADD ? Sym[Sym(c) * t for (c, t) in zip(g.coefs, g.args)] : Sym[g]
    g.kind == ADD && !iszero(g.num) && push!(terms, Sym(g.num))
    groups = Dict{Tuple{Sym,Sym},Tuple{Int,Vector{Sym}}}()      # (α, β) => (degree, terms)
    for t in terms
        f = _ode_term_form(t, x)
        f === nothing && return nothing
        m, α, β, u = f
        d, ts = get(groups, (α, β), (0, Sym[]))
        groups[(α, β)] = (max(d, m), push!(ts, u))
    end
    yp = ZERO
    for ((α, β), (m, ts)) in groups
        s = _ode_multiplicity(a, iszero(β) ? α : α + IM * β)
        A = [Sym(Symbol("_ode_A", j)) for j in 0:m]
        B = [Sym(Symbol("_ode_B", j)) for j in 0:m]
        P1 = sum(A[j + 1] * x^j for j in 0:m)
        P2 = sum(B[j + 1] * x^j for j in 0:m)
        trial = iszero(β) ? x^s * exp(α * x) * P1 : x^s * exp(α * x) * (P1 * cos(β * x) + P2 * sin(β * x))
        unknowns = iszero(β) ? A : vcat(A, B)
        R = sum(a[k + 1] * diff(trial, x, k) for k in 0:(length(a) - 1)) - sum(ts)
        R = expand(_explog(R * exp(-α * x)))
        cv, sv = Sym(:_ode_cv), Sym(:_ode_sv)
        iszero(β) || (R = expand(subs(R, cos(β * x) => cv, sin(β * x) => sv)))
        vars = iszero(β) ? [x] : [x, cv, sv]
        Q = try Poly{Sym}(R, vars) catch; return nothing end
        any(c -> has(c, x), coefficients(Q)) && return nothing
        eqs = [c ~ 0 for c in coefficients(Q)]
        sol = solve(eqs, unknowns)
        isempty(sol) && return nothing
        yp += subs(trial, (u => get(sol[1], u, ZERO) for u in unknowns)...)
    end
    yp
end

# t = c xᵐ e^(αx) cos(βx + φ) (or sin, or neither) as (m, α, β, t), φ free of x; the
# phase is split off (cos(βx + φ) = cos φ cos βx - sin φ sin βx); nothing otherwise
_ode_factors(t::Sym) = t.kind == MUL ? [(b, p) for (b, p) in zip(t.args, t.exps)] :
    t.kind == POW ? [(t.args[1], t.args[2])] : [(t, ONE)]
function _ode_term_form(t::Sym, x::Sym)
    fs = _ode_factors(t)
    m, α, β = 0, ZERO, ZERO
    trig = false
    u = t.kind == MUL ? Sym(t.num) : ONE
    for (b, p) in fs
        if !has(b, x) && !has(p, x)
            u *= b^p
        elseif _eq(b, x) && p.kind == NUMBER && _isint(p.num) && p.num > N0
            m += _toint(p.num)
            u *= b^p
        elseif _eq(b, E) && _ode_islinear(p, x)
            α += diff(p, x)
            u *= b^p
        elseif b.kind == FN && b.name in (:sin, :cos) && _eq(p, ONE) && !trig && _ode_islinear(b.args[1], x)
            trig = true
            β = diff(b.args[1], x)
            φ = expand(b.args[1] - β * x)
            if iszero(φ)
                u *= b
            else
                s, c = sin(β * x), cos(β * x)
                u *= b.name === :sin ? s * cos(φ) + c * sin(φ) : c * cos(φ) - s * sin(φ)
            end
        else
            return nothing
        end
    end
    m, expand(α), expand(β), u
end
_ode_islinear(p::Sym, x::Sym) = (d = diff(p, x); has(p, x) && !has(d, x))

# variation of parameters: u'ᵢ = Wᵢ/W (Cramer), yₚ = Σ uᵢ yᵢ
function _ode_variation(a::Vector{Sym}, g::Sym, x::Sym, basis::Vector{Sym})
    n = length(basis)
    W = [diff(basis[j], x, i) for i in 0:(n - 1), j in 1:n]
    w = _ode_tidy(LinearAlgebra.det(W))
    rhs = [i == n ? g / a[end] : ZERO for i in 1:n]
    yp = ZERO
    for j in 1:n
        Wj = copy(W)
        Wj[:, j] = rhs
        uj = _ode_tidy(LinearAlgebra.det(Wj) / w)
        yp += _integrate(uj, x) * basis[j]
    end
    yp
end

## first order: y' = f(x, y)

function _ode_first_order(P::_ODE, Cs)
    P.n == 1 || return nothing
    Y0, Y1 = P.ys
    x = P.x
    a = diff(P.F, Y1)
    f = has(a, Y1) ? (fs = solve(P.F ~ 0, Y1); length(fs) == 1 ? fs[1] : return nothing) : -subs(P.F, Y1 => ZERO) / a
    f = _ode_tidy(f)
    C = Cs[1]
    fallback = nothing             # a solution with an unevaluated integral, if no other
    for m in (_ode_separable, _ode_bernoulli, _ode_exact, _ode_homogeneous, _ode_riccati)
        r = try m(P, f, a, C) catch err; err isa InterruptException && rethrow(); nothing end
        r === nothing && continue
        r = _ode_explicit(r, Y0, C)
        _hasfn(r[2], :Integral) || return r
        fallback === nothing && (fallback = r)
    end
    fallback
end

# φ = C solved for y when that has one solution
function _ode_explicit(r, Y0::Sym, C::Sym)
    kind, S = r
    kind === :explicit && return r
    ys = try solve(S ~ C, Y0) catch; Sym[] end
    length(ys) == 1 && return (:explicit, _ode_absorb(_ode_tidy(ys[1]), [C]))
    r
end

# f = g(x) h(y): ∫ dy/h = ∫ g dx + C
function _ode_separable(P::_ODE, f::Sym, a::Sym, C::Sym)
    Y0, x = P.ys[1], P.x
    g, h = _ode_split(f, x, Y0)
    g === nothing && return nothing
    (:implicit, _ode_tidy(_integrate(_ode_tidy(1 / h), Y0) - _integrate(g, x)))
end

# f as g(x) h(y), or (nothing, nothing)
function _ode_split(f::Sym, x::Sym, Y::Sym)
    has(f, Y) || return f, ONE
    has(f, x) || return ONE, f
    fa = try factor(f) catch; f end
    fs = _ode_factors(fa)
    g, h = fa.kind == MUL ? Sym(fa.num) : ONE, ONE
    for (b, p) in fs
        if _eq(b, E)                         # e^(u(x) + v(y))
            u = expand(p)
            parts = u.kind == ADD ? [Sym(c) * t for (c, t) in zip(u.coefs, u.args)] : [u]
            u.kind == ADD && (g *= exp(Sym(u.num)))
            for q in parts
                has(q, x) && has(q, Y) && return nothing, nothing
                has(q, Y) ? (h *= exp(q)) : (g *= exp(q))
            end
            continue
        end
        hx, hy = has(b, x) || has(p, x), has(b, Y) || has(p, Y)
        hx && hy && return nothing, nothing
        hy ? (h *= b^p) : (g *= b^p)
    end
    g, h
end

# y' + p(x) y = q(x) yⁿ: v = y^(1 - n) is linear, v' + (1 - n) p v = (1 - n) q
function _ode_bernoulli(P::_ODE, f::Sym, a::Sym, C::Sym)
    Y0, x = P.ys[1], P.x
    e = expand(f)
    terms = e.kind == ADD ? [Sym(c) * t for (c, t) in zip(e.coefs, e.args)] : [e]
    e.kind == ADD && !iszero(e.num) && return nothing
    byexp = Dict{Sym,Sym}()
    for t in terms
        k = _ode_ypower(t, Y0)
        k === nothing && return nothing
        c, n = k
        byexp[n] = get(byexp, n, ZERO) + c
    end
    ks = collect(keys(byexp))
    (length(ks) == 2 && ONE in ks) || return nothing
    n = only(filter(k -> k != ONE, ks))
    iszero(n) && return nothing
    p, q = -byexp[ONE], byexp[n]
    any(c -> has(c, Y0), (p, q)) && return nothing
    Q = _ODE(P.ys[2] + (1 - n) * p * Y0 - (1 - n) * q, P.ys, x, P.yx, 1)
    sol = _ode_linear(Q, [C])
    sol === nothing && return nothing
    (:explicit, _ode_tidy(sol[2]^(1 / (1 - n))))
end

# t = c(x) yⁿ as (c, n), n free of x and y
function _ode_ypower(t::Sym, Y::Sym)
    has(t, Y) || return (t, ZERO)
    _eq(t, Y) && return (ONE, ONE)
    t.kind == POW && _eq(t.args[1], Y) && !has(t.args[2], Y) && return (ONE, t.args[2])
    t.kind == MUL || return nothing
    c, n = Sym(t.num), ZERO
    for (b, p) in zip(t.args, t.exps)
        if _eq(b, Y) && !has(p, Y)
            n += p
        elseif has(b, Y) || has(p, Y)
            return nothing
        else
            c *= b^p
        end
    end
    c, n
end

# M + N y' = 0 with M_y = N_x (after an integrating factor μ(x) or μ(y)): φ = ∫M dx + h(y)
function _ode_exact(P::_ODE, f::Sym, a::Sym, C::Sym)
    Y0, Y1, x = P.ys[1], P.ys[2], P.x
    N, M = a, subs(P.F, Y1 => ZERO)
    has(N, Y1) && return nothing
    if !iszero(_ode_tidy(diff(M, Y0) - diff(N, x)))
        d = diff(M, Y0) - diff(N, x)
        q = _ode_tidy(d / N)
        if !has(q, Y0)
            μ = _ode_tidy(exp(_integrate(q, x)))
        else
            q = _ode_tidy(-d / M)
            has(q, x) && return nothing
            μ = _ode_tidy(exp(_integrate(q, Y0)))
        end
        M, N = _ode_tidy(μ * M), _ode_tidy(μ * N)
        iszero(_ode_tidy(diff(M, Y0) - diff(N, x))) || return nothing
    end
    I = _integrate(M, x)
    hp = _ode_tidy(N - diff(I, Y0))
    has(hp, x) && return nothing
    (:implicit, _ode_tidy(I + _integrate(hp, Y0)))
end

# f(x, y) = F(y/x): y = v x, x v' + v = F(v) is separable, ∫ dv/(F(v) - v) = log(x) + C
function _ode_homogeneous(P::_ODE, f::Sym, a::Sym, C::Sym)
    Y0, x = P.ys[1], P.x
    v = Sym(:_ode_v)
    Fv = _ode_tidy(subs(f, Y0 => v * x))
    has(Fv, x) && return nothing
    d = _ode_tidy(Fv - v)
    iszero(d) && return nothing
    φ = _integrate(_ode_tidy(1 / d), v) - log(x)
    (:implicit, _ode_tidy(subs(φ, v => Y0 / x)))
end

## Riccati: y' = q₀ + q₁ y + q₂ y²

# y = -u'/(q₂ u) makes u'' - (q₁ + q₂'/q₂) u' + q₀ q₂ u = 0, linear, solved by the
# methods above (constant coefficients, Cauchy–Euler); u = C u₁ + u₂ gives y with one
# constant. Else a particular solution y₁ = c xᵏ or c eᵏˣ, and y = y₁ + 1/v with v linear:
# v' = -(q₁ + 2 q₂ y₁) v - q₂.
function _ode_riccati(P::_ODE, f::Sym, a::Sym, C::Sym)
    Y0, x = P.ys[1], P.x
    q2 = _ode_tidy(diff(f, Y0, 2) / 2)
    q1 = _ode_tidy(subs(diff(f, Y0), Y0 => ZERO))
    q0 = _ode_tidy(subs(f, Y0 => ZERO))
    (iszero(q2) || any(q -> has(q, Y0), (q0, q1, q2))) && return nothing
    iszero(_ode_tidy(f - q0 - q1 * Y0 - q2 * Y0^2)) || return nothing
    r = _ode_riccati_linear(P, q0, q1, q2, C)
    r === nothing || return r
    _ode_riccati_particular(P, q0, q1, q2, C)
end

function _ode_riccati_linear(P::_ODE, q0, q1, q2, C)
    x = P.x
    U = [Sym(Symbol("_ode_u", k)) for k in 0:2]
    G = U[3] - _ode_tidy(q1 + diff(q2, x) / q2) * U[2] + _ode_tidy(q0 * q2) * U[1]
    C2 = Sym(:_ode_C2)
    sol = _ode_linear(_ODE(G, U, x, P.yx, 2), [C, C2])
    sol === nothing && return nothing
    u = subs(sol[2], C2 => ONE)
    (has(u, C) && !_hasfn(u, :Integral)) || return nothing
    (:explicit, _ode_tidy(-diff(u, x) / (q2 * u)))
end

function _ode_riccati_particular(P::_ODE, q0, q1, q2, C)
    x, c = P.x, Sym(:_ode_c)
    for y1 in (c, c * x, c / x, c * x^2, c / x^2, c * exp(x), c * exp(-x), c * exp(2x), c * exp(-2x))
        res = _ode_tidy(diff(y1, x) - q0 - q1 * y1 - q2 * y1^2)
        n, _ = numer_denom(res)
        has(n, c) || continue
        cs = try solve(expand(n) ~ 0, c) catch; continue end
        for cv in cs
            has(cv, x) && continue                    # c must be a constant
            y1v = _ode_tidy(subs(y1, c => cv))
            iszero(y1v) && continue
            iszero(_ode_tidy(subs(res, c => cv))) || continue
            # v' + (q₁ + 2 q₂ y₁) v = -q₂
            p = _ode_tidy(q1 + 2q2 * y1v)
            μ = _ode_tidy(exp(_integrate(p, x)))
            v = _ode_tidy((_integrate(_ode_tidy(-μ * q2), x) + C) / μ)
            return (:explicit, _ode_tidy(y1v + 1 / v))
        end
    end
    nothing
end

## autonomous: y'' = f(y, y')

# p = y' as a function of s = y, y'' = p dp/ds: a first-order equation for p(s), then
# ∫ ds/p(s) = x + C₂, solved for y where that gives one solution; several branches of p
# (p = ±sqrt(…)) give a solution each (:several).
function _ode_autonomous(P::_ODE, Cs; implicit::Bool = false)
    bs = _ode_auto_branches(P, Cs; implicit)
    bs === nothing && return nothing
    length(bs) == 1 ? bs[1][1:2] : (:several, [b[1:2] for b in bs])
end

# the solutions, one for each branch p = gᵢ(y) of y': (kind, S, c), S explicit or φ = c.
# With conditions y(a) = v, y'(a) = w, C₁ is fixed first, from the relation between p
# and y at (v, w) (so ∫ dy/p is often elementary when the general one is not), the
# branches through (v, w) kept, and C₂ from y(a) = v
function _ode_auto_branches(P::_ODE, Cs; implicit::Bool = false, conds = nothing)
    P.n == 2 || return nothing
    Y0, Y1, Y2 = P.ys
    x = P.x
    (has(P.F, x) || !has(P.F, Y0)) && return nothing
    s, p0, p1 = Sym(:_ode_s), Sym(:_ode_p0), Sym(:_ode_p1)
    Q = _ODE(expand(subs(P.F, Y0 => s, Y1 => p0, Y2 => p0 * p1)), [p0, p1], s, P.yx, 1)
    sol = _ode_solve(Q, Cs[1:1])
    sol === nothing && return nothing
    kind, g = sol
    C1, C2 = Cs[1], Cs[2]
    if conds !== nothing
        i0, i1 = findfirst(c -> c[1] == 0, conds), findfirst(c -> c[1] == 1, conds)
        (i0 === nothing || i1 === nothing) && return nothing
        _, a, v = conds[i0]
        _, a1, w = conds[i1]
        _eq(a, a1) || return nothing
        c1s = kind === :explicit ? (try solve(subs(g, s => v) ~ w, C1) catch; Sym[] end) :
              Sym[_ode_tidy(subs(g, s => v, p0 => w))]
        length(c1s) == 1 || return nothing
        g = _ode_tidy(subs(g, C1 => c1s[1]))
        gs = kind === :explicit ? [g] : (try solve(g ~ c1s[1], p0) catch; Sym[] end)
        gs = [gi for gi in gs if iszero(_ode_tidy(subs(gi, s => v) - w))]
    else
        gs = kind === :explicit ? [g] : (try solve(g ~ C1, p0) catch; Sym[] end)
    end
    # y real: sqrt(exp(2y)) is exp(y)
    gs = [refine(_ode_tidy(gi); assume = Assumptions(s ∈ ℝ)) for gi in gs if !iszero(gi)]
    isempty(gs) && return nothing
    out = Tuple{Symbol,Sym,Sym}[]
    for gi in gs
        h = _ode_tidy(1 / gi)
        H = _ode_tidy(_integrate(h, s))
        if conds !== nothing && _hasfn(H, :Integral)
            # no antiderivative: ∫ᵥʸ dy/p = x - a, a definite integral
            u = Sym(P.yx.name)                       # the bound variable named as the function
            any(z -> _eq(z, u), free_symbols(P.F)) && (u = Sym(:_ode_u))
            push!(out, (:implicit, _mkfn_node(:Integral, Sym[subs(h, s => u), u, v, Y0]) - x, -a))
            continue
        end
        c2 = conds === nothing ? C2 : _ode_tidy(subs(H, s => v) - a)
        ys = implicit || _hasfn(H, :Integral) ? Sym[] : (try solve(H ~ x + c2, s) catch; Sym[] end)
        conds === nothing || (ys = [yy for yy in ys if iszero(_ode_tidy(subs(yy, x => a) - v))])
        if length(ys) == 1
            push!(out, (:explicit, _ode_tidy(ys[1]), ZERO))
        else
            push!(out, (:implicit, _ode_tidy(subs(H, s => Y0) - x), c2))
        end
    end
    out
end

## initial conditions

# a condition as (order, point, value)
function _ode_ic(p::Pair, P::_ODE)
    c = _ode_ic_parts(Sym(p.first))
    c !== nothing && c[1] === P.yx.name && return (c[2], c[3], Sym(p.second))
    throw(ArgumentError("dsolve: an initial condition is y(a) => v or y'(a) => v; got $(p)"))
end

# y(a) or y⁽ᵒ⁾(a) as (name, order, point), or nothing
function _ode_ic_parts(k::Sym)
    k.kind == FN && length(k.args) == 1 && k.name !== :Derivative && return (k.name, 0, k.args[1])
    _isderivat(k) || return nothing
    f, t = k.args[1], k.args[2]
    o = 1
    while f.kind == FN && f.name === :Derivative && length(f.args) == 2 && _eq(f.args[2], t)
        o += 1
        f = f.args[1]
    end
    f.kind == FN && length(f.args) == 1 && _eq(f.args[1], t) ? (f.name, o, k.args[3]) : nothing
end

function _ode_ics(P::_ODE, kind::Symbol, S::Sym, Cs, conds)
    x, Y0 = P.x, P.ys[1]
    if kind === :explicit
        eqs = [_ode_tidy(subs(diff(S, x, o), x => a)) ~ v for (o, a, v) in conds]
        used = [C for C in Cs if has(S, C)]
        # one constant: the equation's solve, which takes C in exp, tan, …; a system's
        # needs it polynomial
        sols = length(eqs) == 1 && length(used) == 1 ? [[used[1] => c] for c in solve(eqs[1], used[1])] :
               solve(eqs, used)
        isempty(sols) && throw(ArgumentError("dsolve: no solution satisfies the initial conditions"))
        T = [_ode_tidy(subs(S, sol...)) for sol in sols]
        length(unique(T)) == 1 || throw(ArgumentError("dsolve: several solutions satisfy the initial conditions: $(join(T, ", "))"))
        return (:explicit, T[1], ZERO)
    end
    # implicit φ(x, y) = C: the constants from y(a) = v and the derivatives at a (by
    # implicit differentiation); y solved for where one branch passes through (a, v)
    i0 = findfirst(c -> c[1] == 0, conds)
    i0 === nothing && throw(ArgumentError("dsolve: an implicit solution takes a condition y(a) => v"))
    _, a, v = conds[i0]
    all(c -> _eq(c[2], a), conds) || throw(ArgumentError("dsolve: the conditions for an implicit solution are at one point"))
    Ci = Cs[P.n]
    if P.n == 1
        c = _ode_tidy(subs(S, x => a, Y0 => v))
    else
        ds = _ode_implicit_derivs(P, S)
        eqs = [_ode_tidy(subs(S, x => a, Y0 => v)) ~ Ci;
               [_ode_tidy(subs(ds[o], x => a, Y0 => v)) ~ w for (o, _, w) in conds if o > 0]]
        sols = solve(eqs, [C for C in Cs if has(S, C) || _eq(C, Ci)])
        length(sols) == 1 || throw(ArgumentError("dsolve: the initial conditions don't fix the constants"))
        S, c = _ode_tidy(subs(S, sols[1]...)), _ode_tidy(subs(Ci, sols[1]...))
    end
    ys = try solve(S ~ c, Y0) catch; Sym[] end
    good = [y for y in ys if (z = _ode_tidy(subs(y, x => a) - v); iszero(z) || _ode_numzero(z, Sym[]) === true)]
    length(good) == 1 && return (:explicit, _ode_tidy(good[1]), ZERO)
    (:implicit, S, c)
end

## systems: linear, constant coefficients
#
# The equations are put in normal form (each unknown's highest derivative solved for)
# and made a first-order system Z' = M Z + F(t), Z the unknowns and their lower
# derivatives. With M = P J P⁻¹ (Jordan form), the columns of P exp(J t) are solutions;
# a conjugate pair of complex columns is replaced by its real and imaginary parts. A
# particular solution is P exp(J t) ∫ exp(-J t) P⁻¹ F dt, the integrand taken to
# complex exponentials, where t^k e^(a t) integrates in closed form.

function dsolve(eqs::AbstractVector, yxs::AbstractVector, ics::Pair...; method::Symbol = :auto)
    method in (:auto, :laplace) || throw(ArgumentError("dsolve: method is :auto or :laplace, not $method"))
    all(e -> e isa Relation, eqs) || throw(ArgumentError("dsolve: give equations, `lhs ~ rhs`"))
    length(eqs) == length(yxs) || throw(ArgumentError("dsolve: give as many equations as unknowns"))
    yxs = Sym[Sym(y) for y in yxs]
    t = _odesys_var(yxs)
    Fs, Ys = _odesys_parse(eqs, yxs, t)
    M, F, Z = _odesys_first_order(Fs, Ys, t)
    N = length(Z)
    Cs = _ode_constants(N, eqs)
    if method === :laplace || any(f -> _hasfn(f, :heaviside) || _hasfn(f, :dirac), F)
        isempty(ics) || length(ics) == N || throw(ArgumentError("dsolve: give $N initial conditions, one for each unknown and order below its highest"))
        S = _odesys_laplace(M, F, Z, Ys, yxs, Cs, t, ics)
        _odesys_verified(Fs, Ys, S, t) || throw(ArgumentError("dsolve: the solution found does not satisfy the equations"))
        return [y ~ v for (y, v) in zip(yxs, S)]
    end
    Φ, Ψ, ΨinvF = _odesys_fundamental(M, F, t)
    yp = Ψ === nothing ? zeros(Sym, N) : _odesys_particular(Φ, Ψ, ΨinvF, F, t)
    sol = [_ode_tidy(sum(Φ[i, j] * Cs[j] for j in 1:N) + yp[i]) for i in 1:N]
    S = [sol[findfirst(z -> _eq(z, Ys[j][1]), Z)] for j in eachindex(yxs)]
    _odesys_verified(Fs, Ys, S, t) || throw(ArgumentError("dsolve: the solution found does not satisfy the equations"))
    if !isempty(ics)
        length(ics) == N || throw(ArgumentError("dsolve: give $N initial conditions, one for each unknown and order below its highest"))
        S = _odesys_ics(S, yxs, Cs, t, ics)
    end
    [y ~ s for (y, s) in zip(yxs, S)]
end

# the independent variable, shared by the unknowns y(t), z(t), …
function _odesys_var(yxs)
    for y in yxs
        (y.kind == FN && length(y.args) == 1 && y.args[1].kind == SYMBOL) ||
            throw(ArgumentError("dsolve: each unknown is a function of one symbol, y(t); got $y"))
    end
    t = yxs[1].args[1]
    all(y -> _eq(y.args[1], t), yxs) || throw(ArgumentError("dsolve: the unknowns are functions of the same variable"))
    length(unique(y.name for y in yxs)) == length(yxs) || throw(ArgumentError("dsolve: the unknowns are different functions"))
    t
end

# the equations F = 0 in placeholders: Ys[j][k + 1] for the k-th derivative of unknown j
function _odesys_parse(eqs, yxs, t)
    function order(e, y)                # e = y⁽ᵒ⁾(t): o, else nothing
        _eq(e, y) && return 0
        (e.kind == FN && e.name === :Derivative && length(e.args) == 2 && _eq(e.args[2], t)) || return nothing
        o = order(e.args[1], y)
        o === nothing ? nothing : o + 1
    end
    ders = [Dict{Sym,Int}() for _ in yxs]
    function walk(e)
        if e.kind == FN
            for (j, y) in enumerate(yxs)
                o = order(e, y)
                o isa Int && (ders[j][e] = o; return)
            end
        end
        foreach(walk, e.args)
        foreach(walk, e.exps)
    end
    Fs0 = Sym[lhs(eq) - rhs(eq) for eq in eqs]
    foreach(walk, Fs0)
    ns = [isempty(d) ? 0 : maximum(values(d)) for d in ders]
    for (j, n) in enumerate(ns)
        n >= 1 || throw(ArgumentError("dsolve: no derivative of $(yxs[j]) in the equations"))
    end
    Ys = [[Sym(Symbol("_ode_", yxs[j].name, "_", k)) for k in 0:ns[j]] for j in eachindex(yxs)]
    sb = [d => Ys[j][o + 1] for j in eachindex(yxs) for (d, o) in ders[j]]
    Fs = Sym[expand(subs(F, sb...)) for F in Fs0]
    for F in Fs, y in yxs
        _hasfn(F, y.name) && throw(ArgumentError("dsolve: $(y.name) appears other than as $y and its derivatives"))
    end
    Fs, Ys
end

# Z' = M Z + F: the highest derivatives solved for (the equations linear, with constant
# coefficients, and in normal form)
function _odesys_first_order(Fs, Ys, t)
    allY = reduce(vcat, Ys)
    for F in Fs, Y in allY
        c = diff(F, Y)
        (any(Z -> has(c, Z), allY) || has(c, t)) &&
            throw(ArgumentError("dsolve: a system is solved when linear with constant coefficients"))
    end
    m = length(Fs)
    tops = [Ys[j][end] for j in 1:m]
    B = Sym[diff(Fs[i], tops[j]) for i in 1:m, j in 1:m]
    R = Sym[expand(Fs[i] - sum(B[i, j] * tops[j] for j in 1:m)) for i in 1:m]
    iszero(_ode_tidy(LinearAlgebra.det(B))) &&
        throw(ArgumentError("dsolve: the highest derivatives can't be solved for (not in normal form)"))
    top = [expand(v) for v in -(inv(B) * R)]
    Z = reduce(vcat, [Y[1:(end - 1)] for Y in Ys])
    N = length(Z)
    M = zeros(Sym, N, N)
    F = zeros(Sym, N)
    row = 0
    for (j, Y) in enumerate(Ys)
        n = length(Y) - 1
        for k in 1:n
            row += 1
            if k < n                             # (y⁽ᵏ⁻¹⁾)' = y⁽ᵏ⁾
                M[row, findfirst(z -> _eq(z, Y[k + 1]), Z)] = ONE
            else                                 # (y⁽ⁿ⁻¹⁾)' = the solved y⁽ⁿ⁾
                for (c, z) in enumerate(Z)
                    M[row, c] = diff(top[j], z)
                end
                F[row] = _ode_tidy(subs(top[j], (z => ZERO for z in Z)...))
            end
        end
    end
    M, F, Z
end

# the real fundamental matrix Φ, and for the particular solution P exp(J t) and
# exp(-J t) P⁻¹ F (nothing when F = 0)
function _odesys_fundamental(M, F, t)
    P, J = jordan_form(M)
    P = map(_ode_cartesian, P)
    N = size(M, 1)
    EJ, EJm = zeros(Sym, N, N), zeros(Sym, N, N)   # exp(J t), exp(-J t)
    i = 1
    while i <= N
        k = 1
        while i + k <= N && J[i + k - 1, i + k] == ONE
            k += 1
        end
        λ = J[i, i]
        for r in i:(i + k - 1), c in r:(i + k - 1)
            d = c - r
            EJ[r, c] = t^d / factorial(d) * exp(λ * t)
            EJm[r, c] = (-t)^d / factorial(d) * exp(-λ * t)
        end
        i += k
    end
    Ψ = P * EJ
    Φ = zeros(Sym, N, N)
    skip = Dict{Sym,Int}()
    col = 0
    for j in 1:N
        λ = J[j, j]
        if get(skip, λ, 0) > 0
            skip[λ] -= 1
            continue
        end
        v = [expand(_ode_cis(Ψ[r, j])) for r in 1:N]
        re = [subs(e, IM => ZERO) for e in v]
        if all(r -> _eq(re[r], v[r]), 1:N)
            col += 1
            Φ[:, col] = [_ode_tidy(e) for e in v]
        else                                     # a complex column: its real and imaginary parts
            im_ = [expand((v[r] - re[r]) / IM) for r in 1:N]
            Φ[:, col + 1] = [_ode_tidy(e) for e in re]
            Φ[:, col + 2] = [_ode_tidy(e) for e in im_]
            col += 2
            c = subs(λ, IM => -IM)               # parameters taken to be real
            skip[c] = get(skip, c, 0) + 1
        end
    end
    col == N || throw(ArgumentError("dsolve: no real fundamental system found"))
    all(iszero, F) && return Φ, nothing, nothing
    Φ, Ψ, EJm * (inv(P) * F)
end

# a particular solution: Ψ ∫ v dt made real, when v is a sum of t^k e^(a t) terms;
# otherwise Φ ∫ Φ⁻¹ F dt with the real fundamental matrix, by _integrate
function _odesys_particular(Φ, Ψ, v, F, t)
    w = Sym[]
    for e in v
        r = _odesys_int(_ode_toexp(e), t)
        r === nothing && break
        push!(w, r)
    end
    if length(w) == length(v)
        return [_ode_tidy(expand(_ode_cis(expand(sum(Ψ[i, j] * w[j] for j in eachindex(w)))))) for i in 1:size(Ψ, 1)]
    end
    u = map(_ode_tidy, inv(Φ) * F)
    w = [_integrate(e, t) for e in u]
    [_ode_tidy(sum(Φ[i, j] * w[j] for j in eachindex(w))) for i in 1:size(Φ, 1)]
end

# ∫ e dt for e a sum of cf t^k e^(a t) terms, or nothing
function _odesys_int(e::Sym, t::Sym)
    e = expand(e)
    terms = e.kind == ADD ? [(Sym(c), a) for (c, a) in zip(e.coefs, e.args)] : [(ONE, e)]
    out = e.kind == ADD ? Sym(e.num) * t : ZERO
    for (c, a) in terms
        r = _odesys_expoly_int(a, t)
        r === nothing && return nothing
        out += c * r
    end
    out
end

# ∫ cf t^k e^(a t) dt (cf, a free of t, k ≥ 0 an integer), or nothing
function _odesys_expoly_int(term::Sym, t::Sym)
    cf, k, a = ONE, 0, ZERO
    fs = term.kind == MUL ? collect(zip(term.args, term.exps)) : [(term, ONE)]
    term.kind == MUL && (cf = Sym(term.num))
    for (b, p) in fs
        if !has(b, t) && !has(p, t)
            cf *= b^p
        elseif _eq(b, t) && p.kind == NUMBER && _isint(p.num) && p.num >= N0
            k += _toint(p.num)
        elseif b === E && (d = diff(p, t); !has(d, t))
            iszero(subs(p, t => ZERO)) || (cf *= exp(subs(p, t => ZERO)))
            a += d
        elseif b.kind == POW && b.args[1] === E && _eq(p, ONE) && (d = diff(b.args[2], t); !has(d, t))
            a += d
        else
            return nothing
        end
    end
    iszero(a) && return cf * t^(k + 1) / (k + 1)
    # ∫ t^k e^(a t) = e^(a t) Σ_j (-1)^j k!/(k-j)! t^(k-j) / a^(j+1)
    cf * exp(a * t) * sum((-1)^j * Sym(factorial(big(k)) // factorial(big(k - j))) * t^(k - j) / a^(j + 1) for j in 0:k)
end

# sin, cos, sinh, cosh as complex exponentials
function _ode_toexp(e::Sym)
    isempty(e.args) && return e
    e = _rebuild_node(e, Sym[_ode_toexp(a) for a in e.args])
    e.kind == MUL && return Sym(e.num) * prod((_ode_toexp(b)^_ode_toexp(p) for (b, p) in zip(e.args, e.exps)); init = ONE)
    (e.kind == FN && length(e.args) == 1) || return e
    u = e.args[1]
    e.name === :sin && return (exp(IM * u) - exp(-IM * u)) / (2IM)
    e.name === :cos && return (exp(IM * u) + exp(-IM * u)) / 2
    e.name === :sinh && return (exp(u) - exp(-u)) / 2
    e.name === :cosh && return (exp(u) + exp(-u)) / 2
    e
end

# a number with i as a + i b (no i in a denominator, as in 2/(i sqrt(3) - 1))
function _ode_cartesian(p::Sym)
    _eq(subs(p, IM => ZERO), p) && return p         # (parameters taken to be real)
    n, d = numer_denom(p)
    d̄ = subs(d, IM => -IM)                   # i appears as a number: its conjugate
    n, d = expand(n * d̄), expand(d * d̄)
    _eq(subs(d, IM => ZERO), d) || return p
    re = subs(n, IM => ZERO)
    _ode_tidy(re / d) + IM * _ode_tidy(expand((n - re) / IM) / d)
end

# exp(a + i b) as exp(a) (cos b + i sin b)
function _ode_cis(e::Sym)
    isempty(e.args) && return e
    if e.kind == POW && e.args[1] === E
        u = expand(_ode_cis(e.args[2]))
        re = subs(u, IM => ZERO)
        _eq(re, u) && return exp(u)
        b = expand((u - re) / IM)
        _eq(subs(b, IM => ZERO), b) || return exp(u)
        return exp(re) * (cos(b) + IM * sin(b))
    elseif e.kind == MUL
        return Sym(e.num) * prod((b === E ? _ode_cis(b^p) : _ode_cis(b)^_ode_cis(p) for (b, p) in zip(e.args, e.exps)); init = ONE)
    end
    _rebuild_node(e, Sym[_ode_cis(a) for a in e.args])
end

function _odesys_verified(Fs, Ys, S, t)
    for F in Fs
        r = subs(F, (Ys[j][k + 1] => diff(S[j], t, k) for j in eachindex(Ys) for k in 0:(length(Ys[j]) - 1))...)
        z = try simplify(r) catch; r end
        iszero(z) && continue
        ok = _ode_numzero(r, free_symbols(r))
        ok === false && return false
    end
    true
end

function _odesys_ics(S, yxs, Cs, t, ics)
    eqs = Relation[]
    for p in ics
        c = _ode_ic_parts(Sym(p.first))
        j = c === nothing ? nothing : findfirst(y -> y.name === c[1], yxs)
        j === nothing && throw(ArgumentError("dsolve: an initial condition is y(a) => v or y'(a) => v for an unknown y; got $(p)"))
        push!(eqs, _ode_tidy(subs(diff(S[j], t, c[2]), t => c[3])) ~ Sym(p.second))
    end
    sols = solve(eqs, Cs)
    isempty(sols) && throw(ArgumentError("dsolve: no solution satisfies the initial conditions"))
    length(sols) == 1 || throw(ArgumentError("dsolve: the initial conditions don't fix the constants"))
    [_ode_tidy(subs(s, sols[1]...)) for s in S]
end

## series solutions

# The solution as a series in t = x - x₀, terms below tⁿ. At an ordinary point (the
# equation solved for y⁽ᵏ⁾ = f, f with a series at x₀ in the lower derivatives), the
# Taylor coefficients: y⁽ʲ⁾(x₀) = Cⱼ₊₁ (or an initial condition) for j < k, then cₘ from
# the coefficient of tᵐ⁻ᵏ in f(x₀ + t, y, …, y⁽ᵏ⁻¹⁾), which needs only c₀, …, cₘ₋₁. At a
# regular singular point of a second-order linear homogeneous equation, Frobenius.
function _ode_series(P::_ODE, Cs, N::Integer, x0, conds)
    N >= 1 || throw(ArgumentError("dsolve: series = n takes n ≥ 1 terms"))
    a = x0 === nothing ? (isempty(conds) ? ZERO : conds[1][2]) : Sym(x0)
    all(c -> _eq(c[2], a), conds) || throw(ArgumentError("dsolve: initial conditions for a series solution are at its point, $a"))
    isempty(conds) || length(conds) == P.n || throw(ArgumentError("dsolve: give $(P.n) initial condition$(P.n == 1 ? "" : "s"), one for each order"))
    t = Sym(:_ode_t)
    n, x = P.n, P.x
    Yn = P.ys[end]
    c = diff(P.F, Yn)
    f = if has(c, Yn)
        fs = solve(P.F ~ 0, Yn)
        length(fs) == 1 || throw(ArgumentError("dsolve: the highest derivative can't be solved for uniquely"))
        fs[1]
    else
        -subs(P.F, Yn => ZERO) / c
    end
    ft = subs(f, x => a + t)
    vals = Sym[Cs[k + 1] for k in 0:(n - 1)]
    for (o, _, v) in conds
        vals[o + 1] = v
    end
    cs = Sym[vals[k + 1] / factorial(big(k)) for k in 0:(min(n, N) - 1)]
    for m in n:(N - 1)
        S = sum(cs[j + 1] * t^j for j in 0:(m - 1))
        e = subs(ft, (P.ys[k + 1] => diff(S, t, k) for k in 0:(n - 1))...)
        ps = try Base.invokelatest(SymJ.series, e, t, 0, m - n + 1).ps catch; nothing end
        if ps === nothing || any(<(0), ps.e) || any(r -> denominator(r) != 1, ps.e)
            isempty(conds) || throw(ArgumentError("dsolve: $a is a singular point; no series solution with initial conditions"))
            return _ode_frobenius(P, Cs, N, a)
        end
        i = findfirst(==(m - n), ps.e)
        cm = i === nothing ? ZERO : ps.c[i]
        push!(cs, _ode_tidy(expand(cm * factorial(big(m - n)) / factorial(big(m)))))
    end
    St = expand(sum(cs[j + 1] * t^j for j in 0:(length(cs) - 1)))
    _ode_series_check(P, St, t, a, N - n) || throw(ArgumentError("dsolve: the series found does not satisfy the equation"))
    _ode_series_group(St, Cs, t, x - a)
end

# Σ C (its series) + the rest, when the series is linear in the constants
function _ode_series_group(St::Sym, Cs, t::Sym, back::Sym)
    used = [C for C in Cs if has(St, C)]
    lin = all(C -> iszero(diff(St, C, 2)) && all(D -> iszero(diff(diff(St, C), D)), used), used)
    lin || return subs(St, t => back)
    rest = subs(St, (C => ZERO for C in used)...)
    sum((C * subs(diff(St, C), t => back) for C in used); init = subs(rest, t => back))
end

# the residual of a truncated series is O(tᵐ) (m the terms checked)
function _ode_series_check(P::_ODE, St::Sym, t::Sym, a::Sym, m::Integer, shift = 0)
    m <= 0 && return true
    R = subs(P.F, P.x => a + t, (P.ys[k + 1] => diff(St, t, k) for k in 0:P.n)...)
    ps = try Base.invokelatest(SymJ.series, R, t, 0, shift + m).ps catch; return true end
    all(i -> ps.e[i] >= shift + m || iszero(simplify(ps.c[i])), eachindex(ps.e))
end

# Frobenius at a regular singular point of a₂ y'' + a₁ y' + a₀ y = 0: t p = t a₁/a₂ and
# t² q = t² a₀/a₂ with series p₀ + p₁ t + …, q₀ + …; the indicial equation
# F(r) = r(r - 1) + p₀ r + q₀ = 0, and aₘ(r) F(r + m) = -Σₖ ((r + m - k) pₖ + qₖ) aₘ₋ₖ(r).
# Roots r₁ ≥ r₂: tʳ¹ Σ aₘ(r₁) tᵐ, and for the second solution tʳ² Σ aₘ(r₂) tᵐ when
# r₁ - r₂ is not an integer, else ∂/∂r of tʳ Σ aₘ(r) tᵐ (equal roots) or of
# (r - r₂) tʳ Σ aₘ(r) tᵐ (roots an integer apart) at r₂, which has the log(t) term.
function _ode_frobenius(P::_ODE, Cs, N::Integer, a::Sym)
    nope() = throw(ArgumentError("dsolve: $a is a singular point; a series there is found for a second-order linear homogeneous equation with a regular singular point"))
    P.n == 2 || nope()
    lc = _ode_linear_coefs(P)
    lc === nothing && nope()
    (A, g) = lc
    iszero(g) || nope()
    t, x = Sym(:_ode_t), P.x
    coefs(e) = begin
        ps = try Base.invokelatest(SymJ.series, _ode_tidy(subs(e, x => a + t)), t, 0, N).ps catch; nope() end
        all(r -> r >= 0 && denominator(r) == 1, ps.e) || nope()        # an irregular singular point
        v = zeros(Sym, N)
        for (e_, c_) in zip(ps.e, ps.c)
            e_ < N && (v[Int(e_) + 1] = c_)
        end
        v
    end
    p = coefs(_ode_tidy(subs(t, t => x - a) * A[2] / A[3]))
    q = coefs(_ode_tidy((x - a)^2 * A[1] / A[3]))
    r = Sym(:_ode_r)
    F(ρ) = ρ * (ρ - 1) + p[1] * ρ + q[1]
    rs = solve(expand(F(r)) ~ 0, r)
    isempty(rs) && nope()
    all(ρ -> _eq(subs(ρ, IM => ZERO), ρ), rs) ||
        throw(ArgumentError("dsolve: the indicial roots at $a are complex; not done"))
    r1, r2 = length(rs) == 1 ? (rs[1], rs[1]) : (rs[1], rs[2])
    d = _ode_tidy(r1 - r2)
    if d.kind == NUMBER && d.num < N0
        r1, r2, d = r2, r1, -d
    end
    function coeffs(ρ)
        out = Sym[ONE]
        for m in 1:(N - 1)
            s = sum(((ρ + m - k) * p[k + 1] + q[k + 1]) * out[m - k + 1] for k in 1:m)
            push!(out, cancel(-s / F(ρ + m)))
        end
        out
    end
    T = x - a
    ser(ρ, cs) = sum(_ode_tidy(c_) * T^(ρ + m) for (m, c_) in zip(0:(N - 1), cs))
    y1 = ser(r1, [_ode_tidy(subs(c_, r => r1)) for c_ in coeffs(r)])
    if !(d.kind == NUMBER && _isint(d.num))
        y2 = ser(r2, [_ode_tidy(subs(c_, r => r2)) for c_ in coeffs(r)])
    elseif iszero(d)
        ar = coeffs(r)
        y2 = y1 * log(T) + ser(r1, [_ode_tidy(subs(cancel(diff(c_, r)), r => r1)) for c_ in ar])
    else
        br = [cancel((r - r2) * c_) for c_ in coeffs(r)]
        y2 = log(T) * ser(r2, [_ode_tidy(subs(c_, r => r2)) for c_ in br]) +
             ser(r2, [_ode_tidy(subs(cancel(diff(c_, r)), r => r2)) for c_ in br])
    end
    # each solution's residual, over a₂/t², is O(t^(r₂ + n)) (when r₂ is a number)
    if r2.kind == NUMBER
        for y in (y1, y2)
            R = _ode_residual(P, :explicit, y) * (x - a)^2 / A[3]
            ps = try Base.invokelatest(SymJ.series, subs(R, x => a + t), t, 0, _jlnum(r2.num) + N).ps catch; continue end
            all(i -> ps.e[i] >= r2.num + N || iszero(simplify(ps.c[i])), eachindex(ps.e)) ||
                throw(ArgumentError("dsolve: the series found does not satisfy the equation"))
        end
    end
    Cs[1] * y1 + Cs[2] * y2
end

# Z' = M Z + F by Laplace transforms: Z(s) = (sI - M)⁻¹ (Z(a) + L{F}), Z(a) the initial
# conditions (constants for those not given), t shifted to their point a
function _odesys_laplace(M, F, Z, Ys, yxs, Cs, t, ics)
    N = length(Z)
    Z0 = Sym[Cs[i] for i in 1:N]
    a = ZERO
    for (k, p) in enumerate(ics)
        c = _ode_ic_parts(Sym(p.first))
        j = c === nothing ? nothing : findfirst(y -> y.name === c[1], yxs)
        j === nothing && throw(ArgumentError("dsolve: an initial condition is y(a) => v or y'(a) => v for an unknown y; got $(p)"))
        k == 1 && (a = c[3])
        _eq(c[3], a) || throw(ArgumentError("dsolve: the initial conditions are at one point"))
        i = findfirst(z -> _eq(z, Ys[j][c[2] + 1]), Z)
        i === nothing && throw(ArgumentError("dsolve: $(p.first) is not below the highest derivative of $(yxs[j])"))
        Z0[i] = Sym(p.second)
    end
    s = Sym(:_ode_s)
    Fs = Sym[laplace_transform(subs(f, t => t + a), t, s) for f in F]
    R = map(cancel, inv(s * Matrix{Sym}(LinearAlgebra.I, N, N) - M))
    Zs = R * (Z0 + Fs)
    [_ode_tidy(subs(inverse_laplace_transform(Zs[findfirst(z -> _eq(z, Ys[j][1]), Z)], s, t), t => t - a)) for j in eachindex(yxs)]
end

