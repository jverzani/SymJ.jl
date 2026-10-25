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
#   in x or in y alone), homogeneous (f a function of y/x).
#
# A solution is explicit, y = S(x, C₁, …), or implicit, φ(x, y) = C₁. It is checked by
# substitution: symbolically, else numerically at a few points.

"""
    SymJ.dsolve(eq, y(x), ics...)

Solve the ordinary differential equation `eq` (an equation, written with `~`) for the
function `y(x)`. Derivatives are written `diff(y(x), x)` or `y'(x)`. The result is an
equation: `y(x) ~ S` with constants `C1, C2, …` (one for each order), or an implicit
solution `φ(x, y(x)) ~ C1` when `y` cannot be solved for.

Initial conditions are pairs, one for each order: `y(0) => 1, y'(0) => 0`. They fix the
constants.

Methods: first order linear, separable, Bernoulli, exact (with an integrating factor in
`x` or `y` alone), homogeneous; linear with constant coefficients, of any order
(undetermined coefficients, else variation of parameters, for the right side);
Cauchy–Euler; equations without `y` (reduced to a lower order). Integrals are done by
the RischIntegrate extension when it is loaded (an unevaluated `Integral` otherwise, or
when there is no elementary antiderivative).

Not exported.

```julia
@syms x y()
SymJ.dsolve(y'(x) ~ y(x), y(x))                          # y(x) ~ C1*exp(x)
SymJ.dsolve(y''(x) + y(x) ~ 0, y(x), y(0) => 0, y'(0) => 1)   # y(x) ~ sin(x)
SymJ.dsolve(y'(x) ~ x * y(x)^2, y(x))                    # y(x) ~ -2/(x^2 + 2C1)
SymJ.dsolve(x^2 * y''(x) - 2y(x) ~ 0, y(x))              # y(x) ~ C1/x + C2*x^2
```
"""
function dsolve(eq::Relation, yx::Sym, ics::Pair...)
    P = _ode_parse(eq, yx)
    Cs = _ode_constants(P.n, eq)
    sol = _ode_solve(P, Cs)
    sol === nothing && throw(ArgumentError("dsolve: no method for $(eq)"))
    kind, S = sol
    kind === :explicit && (S = _ode_absorb(S, Cs))
    _ode_verified(P, kind, S, Cs) || throw(ArgumentError("dsolve: the solution found does not satisfy $(eq)"))
    c = Cs[1]                                    # an implicit solution is φ = c
    if !isempty(ics)
        length(ics) == P.n || throw(ArgumentError("dsolve: give $(P.n) initial condition$(P.n == 1 ? "" : "s"), one for each order"))
        kind, S, c = _ode_ics(P, kind, S, Cs, [_ode_ic(p, P) for p in ics])
    end
    kind === :explicit ? (yx ~ S) : (subs(S, P.ys[1] => yx) ~ c)
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

# C1, …, Cn, skipping names in eq
function _ode_constants(n::Int, eq::Relation)
    used = Set(s.name for s in free_symbols(lhs(eq) - rhs(eq)))
    out = Sym[]
    i = 1
    while length(out) < n
        nm = Symbol("C", i)
        nm in used || push!(out, Sym(nm))
        i += 1
    end
    out
end

# the residual of a solution: explicit y = S, or implicit φ = C (first order: y' = -φₓ/φ_y)
function _ode_residual(P::_ODE, kind::Symbol, S::Sym)
    if kind === :explicit
        return subs(P.F, (P.ys[k + 1] => diff(S, P.x, k) for k in 0:P.n)...)
    end
    Y0 = P.ys[1]
    yp = -diff(S, P.x) / diff(S, Y0)
    subs(P.F, P.ys[2] => yp)
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
    for m in (_ode_reduce, _ode_linear, _ode_first_order)
        r = try m(P, Cs) catch err; err isa InterruptException && rethrow(); nothing end
        r === nothing || return r
    end
    nothing
end

# exp(a + c log(u)) = u^c exp(a), and the solution simplified where smaller
function _ode_tidy(e::Sym)
    e = _explog(e)
    s = try simplify(e) catch; e end
    _nodes(s) <= _nodes(e) ? s : e
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
# and e^C (in an exponent) is C
function _ode_absorb(S::Sym, Cs)
    for C in Cs
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
    sol === nothing && return nothing
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
    basis = Sym[]
    done = Set{Sym}()
    for (ρ, m) in roots
        ρ in done && continue
        re, im_ = _ode_reim(ρ)
        if im_ !== nothing && !iszero(im_)
            # a conjugate pair, once
            push!(done, ρ, conj(ρ))
            β = _ode_tidy(abs(im_))
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

# the real and imaginary parts of a root, when it has numbers only (decidable)
function _ode_reim(ρ::Sym)
    isempty(free_symbols(ρ)) || return ρ, nothing
    _ode_tidy(real(ρ)), _ode_tidy(imag(ρ))
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
    for m in (_ode_separable, _ode_bernoulli, _ode_exact, _ode_homogeneous)
        r = try m(P, f, a, C) catch err; err isa InterruptException && rethrow(); nothing end
        r === nothing && continue
        return _ode_explicit(r, Y0, C)
    end
    nothing
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

## initial conditions

# a condition as (order, point, value)
function _ode_ic(p::Pair, P::_ODE)
    k, v = Sym(p.first), Sym(p.second)
    if k.kind == FN && k.name === P.yx.name && length(k.args) == 1
        return (0, k.args[1], v)
    end
    if _isderivat(k)
        f, t = k.args[1], k.args[2]
        o = 1
        while f.kind == FN && f.name === :Derivative && length(f.args) == 2 && _eq(f.args[2], t)
            o += 1
            f = f.args[1]
        end
        f.kind == FN && f.name === P.yx.name && length(f.args) == 1 && _eq(f.args[1], t) && return (o, k.args[3], v)
    end
    throw(ArgumentError("dsolve: an initial condition is y(a) => v or y'(a) => v; got $(p)"))
end

function _ode_ics(P::_ODE, kind::Symbol, S::Sym, Cs, conds)
    x, Y0 = P.x, P.ys[1]
    if kind === :explicit
        eqs = [_ode_tidy(subs(diff(S, x, o), x => a)) ~ v for (o, a, v) in conds]
        used = [C for C in Cs if has(S, C)]
        sols = solve(eqs, used)
        isempty(sols) && throw(ArgumentError("dsolve: no solution satisfies the initial conditions"))
        T = [_ode_tidy(subs(S, sol...)) for sol in sols]
        length(unique(T)) == 1 || throw(ArgumentError("dsolve: several solutions satisfy the initial conditions: $(join(T, ", "))"))
        return (:explicit, T[1], ZERO)
    end
    # implicit φ(x, y) = C: C = φ(a, v); y solved for where one branch passes through (a, v)
    o, a, v = only(conds)
    o == 0 || throw(ArgumentError("dsolve: an implicit solution takes a condition y(a) => v"))
    c = _ode_tidy(subs(S, x => a, Y0 => v))
    ys = try solve(S ~ c, Y0) catch; Sym[] end
    good = [y for y in ys if (z = _ode_tidy(subs(y, x => a) - v); iszero(z) || _ode_numzero(z, Sym[]) === true)]
    length(good) == 1 && return (:explicit, _ode_tidy(good[1]), ZERO)
    (:implicit, S, c)
end
