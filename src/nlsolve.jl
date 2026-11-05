# Polynomial systems: a lex Gröbner basis is triangular, and is solved by
# back-substitution from the last variable. At each level the polynomials whose largest
# variable is x_k, with the values found so far substituted, are equations in x_k
# alone: the one of least degree is solved (by `solve` for one equation) and its roots
# kept where the others vanish too. A level with no such polynomial is a free
# variable (infinitely many solutions), left symbolic (`solve` maps it to itself).
#
# The last variable's polynomial is factored first, and the basis recomputed with each
# irreducible factor: then the other polynomials are reduced modulo that factor, so the
# values substituted into them are as simple as they can be.
#
# Rational equations are solved through their numerators, and solutions where a
# denominator vanishes dropped. Floating point coefficients are converted to the exact
# rationals they stand for and the solutions converted back to floating point.

function _solve_nonlinear(eqs::Vector{<:Relation}, vars::Vector{Sym})
    r = _solve_abs_system(eqs, vars)
    r === nothing || return r
    r = _solve_radical_system(eqs, vars)
    r === nothing || return r
    try
        return _solve_polynomial_system(eqs, vars)
    catch err
        (err isa ArgumentError && occursin("must be polynomial or rational", sprint(showerror, err))) || rethrow()
        sols = _solve_triangular(eqs, vars)
        sols === nothing ? rethrow() : return sols
    end
end

# |u| in the equations, u with unknowns: the unknowns taken to be real (|u| of complex
# unknowns gives no finite set of solutions), the system with |u| = u where u ≥ 0 and
# |u| = -u where u ≤ 0, the real solutions kept that fall in their case (|x| + |y| = 2,
# x - y = 1; SymPy #24189). nothing when there is no such |u|.
function _solve_abs_system(eqs::Vector{<:Relation}, vars::Vector{Sym})
    a = nothing
    for eq in eqs
        _replace_walk(eq.lhs - eq.rhs) do t
            a === nothing && t.kind == FN && t.name === :abs && length(t.args) == 1 &&
                any(v -> has(t.args[1], v), vars) && (a = t)
            nothing
        end
        a === nothing || break
    end
    a === nothing && return nothing
    u = a.args[1]
    out = Dict{Sym,Sym}[]
    for (w, op) in ((u, :ge), (-u, :le))
        es = Relation[Relation(:eq, subs(e.lhs, a => w), subs(e.rhs, a => w)) for e in eqs]
        sols = _solve_generic(es, vars)
        for sol in sols
            all(v -> !_hassymbol(sol[v]) || _eq(sol[v], v), vars) || continue
            any(v -> _eq(sol[v], v), vars) && continue          # a free unknown: not decided here
            vals = [N(sol[v]) for v in vars]
            all(v -> isreal(v), vals) || continue
            uv = N(subs(u, (v => sol[v] for v in vars)...))
            (isreal(uv) && (op === :ge ? real(uv) >= -1e-12 : real(uv) <= 1e-12)) || continue
            any(o -> all(v -> _eq(o[v], sol[v]), vars), out) || push!(out, sol)
        end
    end
    out
end

# square roots √R of polynomials R in the unknowns: each a new unknown w with w² = R, the
# polynomial system solved, and the solutions kept where w is the principal root
# (Lagrange multipliers with a distance √((x+1)² + …); SymPy #24868, #24269). nothing
# when there is no such root.
function _solve_radical_system(eqs::Vector{<:Relation}, vars::Vector{Sym})
    roots = Sym[]
    for eq in eqs
        _replace_walk(eq.lhs - eq.rhs) do t
            if t.kind == POW && t.args[2].kind == NUMBER && abs(t.args[2].num) == SymNumber(1 // 2) &&
               any(v -> has(t.args[1], v), vars)
                R = t.args[1]
                any(o -> _eq(o, R), roots) || push!(roots, R)
            end
            nothing
        end
    end
    isempty(roots) && return nothing
    length(roots) <= 3 || return nothing
    atoms = Set{Sym}()
    all(R -> _poly_atoms!(atoms, R), roots) || return nothing
    any(a -> !any(v -> _eq(a, v), vars) && any(v -> has(a, v), vars), atoms) && return nothing
    ws = [Sym(Symbol("_w", i)) for i in eachindex(roots)]
    # R^(k/2) as w^k
    sub(e) = _replace_walk(e) do t
        t.kind == POW && t.args[2].kind == NUMBER && _isint(2 * t.args[2].num) || return nothing
        i = findfirst(R -> _eq(R, t.args[1]), roots)
        i === nothing ? nothing : ws[i]^_toint(2 * t.args[2].num)
    end
    es = Relation[Relation(:eq, sub(e.lhs), sub(e.rhs)) for e in eqs]
    any(e -> any(R -> any(v -> has(e.lhs - e.rhs, v), vars) && has(e.lhs - e.rhs, sqrt(R)), roots), es) && return nothing
    append!(es, [Relation(:eq, w^2, R) for (R, w) in zip(roots, ws)])
    sols = try _solve_generic(es, [vars; ws]) catch err; err isa ArgumentError || rethrow(); return nothing end
    out = Dict{Sym,Sym}[]
    for sol in sols
        ok = all(zip(roots, ws)) do (R, w)
            wv = sol[w]
            _hassymbol(wv) && return true                   # undecided: kept
            pr = sqrt(subs(R, (v => sol[v] for v in vars)...))
            z = try N(wv - pr) catch; return true end
            abs(z) <= 1e-9 * max(1, abs(N(wv)))
        end
        ok || continue
        d = Dict{Sym,Sym}(v => sol[v] for v in vars)
        any(o -> all(v -> _eq(o[v], d[v]), vars), out) || push!(out, d)
    end
    out
end

# a system that is not polynomial, solved an equation in one unknown at a time: that
# unknown's values substituted into the rest (T = h/c, then p from log(p) = …; SymPy
# #28064). nothing when no equation has a single unknown.
function _solve_triangular(eqs::Vector{<:Relation}, vars::Vector{Sym})
    # equations without unknowns: true ones dropped, a false one (generically) no solution
    free = [e for e in eqs if !any(v -> has(e.lhs - e.rhs, v), vars)]
    if !isempty(free)
        all(e -> iszero(simplify(e.lhs - e.rhs)), free) || return Dict{Sym,Sym}[]
        eqs = [e for e in eqs if any(v -> has(e.lhs - e.rhs, v), vars)]
        isempty(eqs) && return nothing
    end
    for (i, eq) in enumerate(eqs)
        vs = [v for v in vars if has(eq.lhs - eq.rhs, v)]
        length(vs) == 1 || continue
        v = only(vs)
        rs = try _solve_generic(eq, v) catch err; err isa ArgumentError || rethrow(); continue end
        rs isa AbstractVector || continue
        others = [u for u in vars if !_eq(u, v)]
        out = Dict{Sym,Sym}[]
        for r in rs
            es = Relation[]
            ok = true
            for (j, e) in enumerate(eqs)
                j == i && continue
                f = Relation(:eq, subs(e.lhs, v => r), subs(e.rhs, v => r))
                if any(u -> has(f.lhs - f.rhs, u), others)
                    push!(es, f)
                elseif !iszero(simplify(f.lhs - f.rhs))
                    ok = false
                    break
                end
            end
            ok || continue
            sub = isempty(others) || isempty(es) ? [Dict{Sym,Sym}()] : _solve_generic(es, others)
            for d in sub
                d = copy(d)
                d[v] = r
                push!(out, d)
            end
        end
        return out
    end
    nothing
end

function _solve_polynomial_system(eqs::Vector{<:Relation}, vars::Vector{Sym})
    inexact = false
    nums, dens = Sym[], Sym[]
    for eq in eqs
        ex = eq.lhs - eq.rhs
        if !_isexactex(ex)
            inexact = true
            ex = _exactfloats(ex)
        end
        n, d = numer_denom(cancel(ex))
        any(v -> has(d, v), vars) && push!(dens, d)
        atoms = Set{Sym}()
        _poly_atoms!(atoms, n) ||
            throw(ArgumentError("solve: complex floating point coefficients in a nonlinear system are not supported"))
        for a in atoms
            !any(v -> _eq(a, v), vars) && any(v -> has(a, v), vars) &&
                throw(ArgumentError("solve: nonlinear systems must be polynomial or rational in $(join(vars, ", ")), not involve $a"))
        end
        iszero(n) || push!(nums, n)
    end
    sols = isempty(nums) ? [Dict{Sym,Sym}()] : _solve_algebraic(nums, vars)
    # drop solutions at a pole
    filter!(s -> all(d -> !_vanishes(subs(d, s...)), dens), sols)
    if inexact
        sols = [Dict{Sym,Sym}(k => (_hassymbol(v) ? v : Sym(N(v))) for (k, v) in s) for s in sols]
    end
    _uniquesols([_polish(s, nums, vars) for s in sols])
end

# a floating point solution refined by Newton's method on the system, in extended
# precision: back-substitution loses the relative accuracy of small components
# (H = 2.6e-13 for 9.5e-14). The original is kept when the iteration does not converge.
function _polish(s::Dict{Sym,Sym}, nums::Vector{Sym}, vars::Vector{Sym})
    (isempty(nums) || any(v -> !haskey(s, v) || _hassymbol(s[v]), vars)) && return s
    all(v -> _isexactex(s[v]), vars) && return s
    z0 = try ComplexF64[_evalc(s[v]) for v in vars] catch; return s end
    names = Tuple(v.name for v in vars)
    J = [diff(n, v) for n in nums, v in vars]
    ev(e, z) = _ev(e, _Named(NamedTuple{names}(Tuple(z))), BigFloat, Ref{Union{Nothing,Sym}}(nothing))
    z = setprecision(BigFloat, 256) do
        z = Complex{BigFloat}.(z0)
        tol = big(2.0)^-200
        for _ in 1:30
            f = Complex{BigFloat}[ev(n, z) for n in nums]
            Jz = Complex{BigFloat}[ev(d, z) for d in J]
            δ = try (Jz' * Jz) \ (Jz' * f) catch; return nothing end
            all(isfinite, δ) || return nothing
            z -= δ
            maximum(abs, δ) <= tol * max(1, maximum(abs, z)) && return z
        end
        nothing
    end
    z === nothing && return s
    # the values Float64 rounding found: not a different solution
    all(isapprox.(ComplexF64.(z), z0; rtol = 1e-6, atol = 1e-6 * max(1, maximum(abs, z0)))) || return s
    Dict{Sym,Sym}(k => (i = findfirst(v -> _eq(v, k), vars); i === nothing ? w :
        (w0 = z0[i]; c = ComplexF64(z[i]); iszero(imag(w0)) ? Sym(real(c)) : Sym(real(c)) + Sym(imag(c)) * IM))
        for (k, w) in s)
end

# Algebraic constants in the coefficients (sqrt(3)) are not parameters: a Gröbner basis
# over Q(p) can have a leading coefficient that vanishes at p = sqrt(3) (p^2 - 3). Each
# b^(1/q) becomes a new last variable r with r^q - b in the system, and the solutions
# kept are those with r = b^(1/q). Only then: a new variable per constant is costly
# (10^(-189/20) gives r^20 - 10^-189), so the basis over Q(p) is tried first.
function _solve_algebraic(nums::Vector{Sym}, vars::Vector{Sym})
    try
        return _solve_groebner(nums, vars)
    catch e
        e isa ArgumentError && e.msg == _GROEBNER_ALGEBRAIC_MSG || rethrow()
    end
    out = _solve_lifted(nums, vars, false)
    out === nothing && throw(ArgumentError(_GROEBNER_ALGEBRAIC_MSG))
    out
end

# the solutions with each algebraic constant b^(1/q) of the coefficients a new variable r
# (r^q - b in the system), kept where r = b^(1/q); nothing when there are no such
# constants. With first = true the new variables come first in the lex order, so they are
# eliminated: a univariate polynomial over Q(√3) becomes one over Q, which factors.
function _solve_lifted(nums::Vector{Sym}, vars::Vector{Sym}, first::Bool)
    atoms = Set{Sym}()
    foreach(n -> _poly_atoms!(atoms, n), nums)
    algs = sort!([a for a in atoms if _isalgconst(a)]; lt = _lt)
    isempty(algs) && return nothing
    rs = Sym[]
    for i in eachindex(algs)
        k = i
        r = Sym(Symbol("_alg", k))
        while any(n -> has(n, r), nums) || any(v -> _eq(v, r), vars)
            k += length(algs)
            r = Sym(Symbol("_alg", k))
        end
        push!(rs, r)
    end
    sub = [a => r for (a, r) in zip(algs, rs)]
    eqs = Sym[expand(subs(n, sub...)) for n in nums]
    for (a, r) in zip(algs, rs)
        b, q = a.args[1], _denominator(a.args[2].num)
        push!(eqs, r^q - b^_numerator(a.args[2].num))
    end
    sols = _solve_groebner(eqs, first ? [rs; vars] : [vars; rs])
    out = Dict{Sym,Sym}[]
    for s in sols
        all(((a, r),) -> haskey(s, r) && _vanishes(s[r] - a), sub) || continue
        back = [r => a for (a, r) in sub]
        push!(out, Dict{Sym,Sym}(k => simplify(subs(v, back...)) for (k, v) in s if !any(r -> _eq(k, r), rs)))
    end
    out
end

# b^(p/q), b and p/q rational: an algebraic number
_isalgconst(a::Sym) = a.kind == POW && a.args[1].kind == NUMBER && isexact(a.args[1].num) &&
    a.args[2].kind == NUMBER && isexact(a.args[2].num) && !_isint(a.args[2].num)

# Back-substitution divides by the leading coefficient of each level's polynomial in its
# variable. Where that coefficient involves free unknowns (x1 (x5 + 1) - x2: x5 free),
# the component on which it vanishes is lost (SymPy #27008): each such factor is added to
# the system and that system solved too. Components already solved (the same basis
# reached through other factors) are skipped; past `_MAX_COMPONENTS` an error, rather
# than an answer missing solutions.
const _MAX_COMPONENTS = 256
const _COMPONENTS = ScopedValue{Union{Nothing,Set{Vector{Sym}}}}(nothing)

function _solve_groebner(nums::Vector{Sym}, vars::Vector{Sym})
    _COMPONENTS[] === nothing && return with(() -> _solve_groebner(nums, vars), _COMPONENTS => Set{Vector{Sym}}())
    # an equation that is a product, x (x - 1) (…) = 0: the system for each factor, whose
    # bases are smaller (one basis of the product took minutes; SymPy #24500)
    for (i, n) in enumerate(nums)
        fs = _unknown_factors([n], vars)
        length(fs) >= 2 || continue
        rest = Sym[nums[j] for j in eachindex(nums) if j != i]
        out = Dict{Sym,Sym}[]
        for f in fs
            append!(out, _solve_groebner([rest; f], vars))
        end
        return _uniquesols(out)
    end
    G = groebner(nums, vars...)
    G == [ONE] && return Dict{Sym,Sym}[]
    seen = _COMPONENTS[]
    G in seen && return Dict{Sym,Sym}[]
    push!(seen, G)
    length(seen) > _MAX_COMPONENTS && throw(ArgumentError("solve: the system splits into more than " *
        "$_MAX_COMPONENTS components (where divisors of its solutions vanish); not supported"))
    degen = Sym[]
    out = _solve_groebner_basis(G, nums, vars, degen)
    for f in _unknown_factors(degen, vars)
        append!(out, _solve_groebner([nums; f], vars))
    end
    isempty(degen) ? out : _uniquesols(out)
end

# the irreducible factors (up to constants) of the expressions in D that involve vars:
# through products and powers (sqrt(-y^2) vanishes where y does)
function _unknown_factors(D::Vector{Sym}, vars::Vector{Sym})
    out = Sym[]
    function add!(b)
        if b.kind == MUL
            for (a, p) in zip(b.args, b.exps)
                p.kind == NUMBER && p.num > N0 && add!(a)
            end
            return
        elseif b.kind == POW
            p = b.args[2]
            p.kind == NUMBER && p.num > N0 && add!(b.args[1])
            return
        end
        any(v -> has(b, v), vars) || return
        atoms = Set{Sym}()
        _poly_atoms!(atoms, b)
        all(_isnamed, atoms) || return            # a polynomial (no radicals of values found)
        f = try factor(b) catch err; err isa ArgumentError || rethrow(); b end
        f.kind in (MUL, POW) && return add!(f)
        f = _content_free(f)
        any(o -> cancel(o / f).kind == NUMBER, out) || push!(out, f)
    end
    foreach(c -> add!(numer_denom(c)[1]), D)
    out
end

function _solve_groebner_basis(G::Vector{Sym}, nums::Vector{Sym}, vars::Vector{Sym}, degen::Vector{Sym})
    _DIVISORS[] === nothing || foreach(g -> _gb_divisors(g, vars), G)
    # split on the irreducible factors of the last variable's polynomial
    n = length(vars)
    last = [g for g in G if _mainvar(g, vars) == n]
    if length(last) == 1 && _onlyvar(last[1], vars[n])
        _, facs = factor(Poly{SymNumber}(last[1], [vars[n]]))
        if length(facs) > 1
            out = Dict{Sym,Sym}[]
            for (f, _) in facs
                degree(f) >= 1 || continue
                Gf = groebner([G; Sym(f)], vars...)
                Gf == [ONE] || append!(out, _backsub(Gf, vars, n, Dict{Sym,Sym}(), degen))
            end
            return out
        end
    end
    _backsub(G, vars, n, Dict{Sym,Sym}(), degen)
end

# a basis element's leading coefficient (lex in vars) and the denominators of its
# coefficients: the basis over Q(parameters) divided by them
function _gb_divisors(g::Sym, vars::Vector{Sym})
    P = try Poly{Sym}(g, vars) catch; return end
    evs, cs = collect(exponent_vectors(P)), collect(coefficients(P))
    isempty(cs) && return
    _divisor!(cs[argmax(evs)])
    foreach(c -> _divisor!(numer_denom(c)[2]), cs)
end

# p a polynomial in x alone, with rational coefficients
_onlyvar(p::Sym, x::Sym) = (atoms = Set{Sym}(); _poly_atoms!(atoms, p) && all(a -> _eq(a, x), atoms))

# the index of the first of vars that p involves (0 if none)
_mainvar(p::Sym, vars::Vector{Sym}) = something(findfirst(v -> has(p, v), vars), 0)

function _backsub(G::Vector{Sym}, vars::Vector{Sym}, k::Int, sol::Dict{Sym,Sym}, degen::Vector{Sym} = Sym[])
    k == 0 && return [sol]
    x = vars[k]
    pairs = [v => sol[v] for v in vars[(k+1):end] if haskey(sol, v)]
    qs = Sym[]
    for g in G
        _mainvar(g, vars) == k || continue
        q = isempty(pairs) ? g : expand(subs(g, pairs...))
        if !has(q, x)
            # a leading coefficient vanished entirely; a nonzero number: no solution here
            _vanishes(q) && continue
            _hassymbol(q) && continue
            return Dict{Sym,Sym}[]
        end
        push!(qs, q)
    end
    isempty(qs) && return _backsub(G, vars, k - 1, sol, degen)      # x is free
    sort!(qs; by = q -> _xdegree(q, x))
    # where a factor common to q's coefficients vanishes, x is free there
    ct = _xcontent(qs[1], x)
    ct !== nothing && any(v -> has(ct, v), vars) && push!(degen, ct)
    out = Dict{Sym,Sym}[]
    for r in _solve1(qs[1], x)
        all(q -> _satisfies(q ~ 0, x, r), qs[2:end]) || continue
        # where r's denominator vanishes, the formula does not give the solutions
        den = numer_denom(r)[2]
        any(v -> has(den, v), vars) && push!(degen, den)
        s = copy(sol)
        s[x] = r
        append!(out, _backsub(G, vars, k - 1, s, degen))
    end
    out
end

# the gcd of q's coefficients in x, c for q = c x^d (nothing when q is not a polynomial
# in x, or when they have no common factor)
function _xcontent(q::Sym, x::Sym)
    P = try Poly{Sym}(q, [x]) catch; return nothing end
    cs = Sym[Sym(c) for c in coefficients(P)]
    isempty(cs) && return nothing
    g = length(cs) == 1 ? cs[1] : try reduce(gcd, cs) catch err; err isa ArgumentError || rethrow(); return nothing end
    _hassymbol(g) ? g : nothing
end

# the roots of q in x: numeric when floating point values were substituted into it
# (complex ones too, which `solve` does not take as coefficients)
function _solve1(q::Sym, x::Sym)
    (_isexactex(q) || _hassymbol(subs(q, x => ZERO))) && return _solve_eq(q ~ 0, x)
    P = Poly{Sym}(q, [x])
    c = fill(ZERO, degree(P, x) + 1)
    for (ev, cf) in zip(exponent_vectors(P), coefficients(P))
        c[ev[1]+1] = cf
    end
    _numeric_roots(c)
end

_xdegree(q::Sym, x::Sym) = try degree(Poly{Sym}(q, [x]), x) catch; typemax(Int) end

# zero, exactly or numerically (to high precision) when it is a constant
function _vanishes(e::Sym)
    (iszero(e) || iszero(cancel(e))) && return true
    _hassymbol(e) && return false
    v = try N(e) catch; return false end
    abs(v) <= 1e-12 * max(1, abs(v))
end

function _uniquesols(sols::Vector{Dict{Sym,Sym}})
    out = Dict{Sym,Sym}[]
    for s in sols
        any(t -> keys(t) == keys(s) && all(k -> _vanishes(t[k] - s[k]), keys(s)), out) || push!(out, s)
    end
    out
end

# exact: no floating point numbers anywhere in the expression
_isexactex(e::Sym) = e.kind == NUMBER ? isexact(e.num) :
    (e.kind == ADD ? isexact(e.num) && all(isexact, e.coefs) : e.kind == MUL ? isexact(e.num) : true) &&
    all(_isexactex, e.args) && all(_isexactex, e.exps)

# floating point numbers replaced by the rationals they represent exactly
_exactnum(n::SymNumber) = isexact(n) ? Sym(n) :
    (v = N(Sym(n)); v isa Real ? Sym(SymNumber(Rational{BigInt}(v))) :
     throw(ArgumentError("solve: complex floating point coefficients in a nonlinear system are not supported")))
# (short = true: the shortest rational within rounding, 0.3 as 3/10, not 5404319552844595/2^54)
function _exactfloats(e::Sym; short::Bool = false)
    k = e.kind
    q(c) = short && !isexact(c) && (v = N(Sym(c)); v isa Real && isfinite(v)) ? Sym(SymNumber(Rational{BigInt}(rationalize(BigInt, Float64(v))))) :
        _exactnum(c)
    r(a) = _exactfloats(a; short)
    k == NUMBER && return q(e.num)
    k == ADD && return q(e.num) + _add_many(q(c) * r(a) for (c, a) in zip(e.coefs, e.args))
    k == MUL && return q(e.num) * prod((r(a)^r(x) for (a, x) in zip(e.args, e.exps)); init = ONE)
    k == POW && return r(e.args[1])^r(e.args[2])
    k == FN && return apply_function(e.name, Sym[r(a) for a in e.args])
    e
end
