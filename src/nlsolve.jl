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
    # factors e^u, never 0, dropped from an equation that is a product ((x² + 2x - y²)e^x
    # = 0; SymPy #21908)
    stripped = Relation[_drop_exp_factors(e, vars) for e in eqs]
    any(i -> stripped[i] !== eqs[i], eachindex(eqs)) && return _solve_generic(stripped, vars)
    # c₁ Π bᵢ^eᵢ = c₂ Π bⱼ^eⱼ in positive numbers: their logs, linear in the exponents
    # (32·2^m·2^n = 4^n, 27·3^m = 3^(-n): m = -4, n = 1; SymPy #19587)
    logged = Relation[something(_exp_products_logged(e, vars), e) for e in eqs]
    any(i -> logged[i] !== eqs[i], eachindex(eqs)) && return _solve_generic(logged, vars)
    les = [_log_eq(e, vars) for e in eqs]
    if any(!isnothing, les)
        es = Relation[l === nothing ? e : l[1] for (e, l) in zip(eqs, les)]
        As = Sym[l[2] for l in les if l !== nothing]
        sols = _solve_generic(es, vars)
        sols isa AbstractVector || return sols
        return filter(s -> all(A -> (a = try subs(A, s...) catch; return false end; !_vanishes(a)), As), sols)
    end
    r = _solve_abs_system(eqs, vars)
    r === nothing || return r
    r = _solve_radical_system(eqs, vars)
    r === nothing || return r
    try
        return _solve_polynomial_system(eqs, vars)
    catch err
        (err isa ArgumentError && occursin("must be polynomial or rational", sprint(showerror, err))) || rethrow()
        sols = _solve_triangular(eqs, vars)
        sols === nothing && (sols = _solve_by_elimination(eqs, vars))
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
    # the radicals: (R, q, real) for √R (principal, q = 2), R^(k/3) (principal, q = 3) and
    # cbrt(R) (the real cube root); cube roots too (cbrt(x - 1) + cbrt(x) + cbrt(x + 1) = 0;
    # SymPy #21167, #21026)
    roots = Tuple{Sym,Int,Bool}[]
    add!(R, q, rl) = any(o -> _eq(o[1], R) && o[2] == q && o[3] == rl, roots) || push!(roots, (R, q, rl))
    for eq in eqs
        _replace_walk(eq.lhs - eq.rhs) do t
            if t.kind == POW && t.args[2].kind == NUMBER && isexact(t.args[2].num) &&
               _denominator(t.args[2].num) in (2, 3) && any(v -> has(t.args[1], v), vars)
                add!(t.args[1], Int(_denominator(t.args[2].num)), false)
            elseif t.kind == MUL
                for (b, p) in zip(t.args, t.exps)
                    p.kind == NUMBER && isexact(p.num) && _denominator(p.num) in (2, 3) && any(v -> has(b, v), vars) &&
                        add!(b, Int(_denominator(p.num)), false)
                end
            elseif t.kind == FN && t.name === :cbrt && length(t.args) == 1 && any(v -> has(t.args[1], v), vars)
                add!(t.args[1], 3, true)
            end
            nothing
        end
    end
    isempty(roots) && return nothing
    length(roots) <= 3 || return nothing
    atoms = Set{Sym}()
    all(r -> _poly_atoms!(atoms, r[1]), roots) || return nothing
    any(a -> !any(v -> _eq(a, v), vars) && any(v -> has(a, v), vars), atoms) && return nothing
    ws = [Sym(Symbol("_w", i)) for i in eachindex(roots)]
    # R^(k/q) as w^k, cbrt(R) as w
    sub(e) = _replace_walk(e) do t
        if t.kind == FN && t.name === :cbrt && length(t.args) == 1
            i = findfirst(r -> _eq(r[1], t.args[1]) && r[3], roots)
            return i === nothing ? nothing : ws[i]
        end
        t.kind == POW && t.args[2].kind == NUMBER && isexact(t.args[2].num) || return nothing
        q = Int(_denominator(t.args[2].num))
        i = findfirst(r -> _eq(r[1], t.args[1]) && r[2] == q && !r[3], roots)
        i === nothing ? nothing : ws[i]^_toint(q * t.args[2].num)
    end
    es = Relation[Relation(:eq, sub(e.lhs), sub(e.rhs)) for e in eqs]
    # (a radical of the unknowns left after the substitution: not this method's)
    left(g) = any(t -> any(v -> has(t.args[1], v), vars), _radical_nodes(g)) ||
              (_hasfn(g, :cbrt) && any(v -> has(g, v), vars))
    any(e -> left(e.lhs - e.rhs), es) && return nothing
    append!(es, [Relation(:eq, w^q, R) for ((R, q, _), w) in zip(roots, ws)])
    sols = try
        _solve_generic(es, [vars; ws])
    catch err
        err isa ArgumentError || rethrow()
        # (a basis that did not finish: that reason, not "must be polynomial" after it)
        occursin("did not finish", sprint(showerror, err)) && rethrow()
        return nothing
    end
    out = Dict{Sym,Sym}[]
    for sol in sols
        ok = all(zip(roots, ws)) do ((R, q, rl), w)
            wv = sol[w]
            _hassymbol(wv) && return true                   # undecided: kept
            Rv = subs(R, (v => sol[v] for v in vars)...)
            z = try
                r = rl ? (v = N(Rv); v isa Real ? sign(v) * abs(v)^(1 / 3) : nothing) : N(Rv^(Sym(1) // q))
                r === nothing ? nothing : N(wv) - r
            catch
                return true
            end
            z === nothing && return false
            abs(z) <= 1e-9 * max(1, abs(N(wv)))
        end
        ok || continue
        d = Dict{Sym,Sym}(v => sol[v] for v in vars)
        any(o -> all(v -> _eq(o[v], d[v]), vars), out) || push!(out, d)
    end
    out
end

# the POW nodes of e with a non-integer exponent
function _radical_nodes(e::Sym)
    out = Sym[]
    _replace_walk(e) do t
        t.kind == POW && t.args[2].kind == NUMBER && !_isint(t.args[2].num) && push!(out, t)
        nothing
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
    # unknowns that an equation of their own determines last in the lex order, solved
    # first: d² = 8, p² = 3, x = d + p gives x = d + p, not x from x⁴ - 22x² + 25 with
    # d, p back from x (nested radicals; SymPy #23322). (Here, for the unknowns: the
    # variables for algebraic constants are placed on purpose.)
    own = [v for v in vars if any(n -> has(n, v) && all(u -> _eq(u, v) || !has(n, u), vars), nums)]
    ovars = vars
    if !isempty(own) && length(own) < length(vars) && !all(v -> v in own, vars[end-length(own)+1:end])
        ovars = [[v for v in vars if !(v in own)]; own]
    end
    sols = isempty(nums) ? [Dict{Sym,Sym}()] : _solve_algebraic(nums, ovars)
    # drop solutions at a pole
    filter!(s -> all(d -> !_vanishes(subs(d, s...)), dens), sols)
    if inexact
        # (with parameters too: their numbers floats again, not the 53-bit rationals the
        # floats were taken as; SymPy #24193)
        sols = [Dict{Sym,Sym}(k => (_hassymbol(v) ? _floatify(v) : Sym(N(v))) for (k, v) in s) for s in sols]
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
    # (the new variables last, else first: with im and parameters, a x + i a y = c, the
    # basis with r last met the same vanishing coefficient; SymPy #23889)
    out = try
        _solve_lifted(nums, vars, false)
    catch e
        e isa ArgumentError && e.msg == _GROEBNER_ALGEBRAIC_MSG || rethrow()
        _solve_lifted(nums, vars, true)
    end
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
    _parametric_degree_guard(nums, vars)
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

# With parameters, the system at sample values of them: when the last unknown's
# polynomial there has an irreducible factor of degree > 4 that does not decompose, so has
# the generic one (specializing only splits factors), whose roots solve cannot give; an
# error then, before the basis over the parameters' rational functions, which can run for
# hours (two quartics with 10 parameters; SymPy #23604)
function _parametric_degree_guard(nums::Vector{Sym}, vars::Vector{Sym})
    length(vars) >= 2 || return
    ps = Sym[]
    for n in nums, a in _named(n)
        any(v -> _eq(v, a), vars) || any(q -> _eq(q, a), ps) || push!(ps, a)
    end
    isempty(ps) && return
    atoms = Set{Sym}()
    all(n -> _poly_atoms!(atoms, n), nums) || return
    all(_isnamed, atoms) || return
    vals = [p => Sym(SymNumber((3i + 4) // (2i + 7))) for (i, p) in enumerate(ps)]
    G = try groebner([expand(subs(n, vals...)) for n in nums], vars...) catch err
        err isa ArgumentError || rethrow(); return
    end
    last = [g for g in G if _mainvar(g, vars) == length(vars)]
    length(last) == 1 && _onlyvar(last[1], vars[end]) || return
    _, facs = factor(Poly{SymNumber}(last[1], [vars[end]]))
    for (f, _) in facs
        d = degree(f)
        # a cubic or quartic in many parameters: Cardano's or Ferrari's formula over
        # coefficients in all of them, and the basis to get there, are too large to use
        if 3 <= d <= 4 && length(ps) >= _PARAMETRIC_RADICAL_MAX
            throw(ArgumentError("solve: the solutions are roots of an irreducible polynomial of degree $d in $(vars[end]) whose coefficients involve $(length(ps)) parameters; a closed form is too large to compute or use. Substitute numbers for the parameters"))
        end
        d > 4 || continue
        c = _xcoeffs(f)
        chain = _decompose(c)
        (length(chain) > 1 && all(g -> length(g) <= 5, chain)) && continue
        throw(ArgumentError("solve: with the parameters, the solutions are roots of an irreducible polynomial of degree $d in $(vars[end]), which have no closed form; substitute numbers for the parameters (the roots are then found numerically)"))
    end
end
const _PARAMETRIC_RADICAL_MAX = 6

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
    _drop_special(out)
end

# solutions that are special cases of a family among the others dropped: {x = 0, y = 0}
# lies in {x = 0} (y free, as an unknown missing from the Dict or y => y is), the
# family's free unknowns taking s's values (x y = x z = x v = 0: two families, not
# eight; SymPy #23110)
function _drop_special(sols::Vector{Dict{Sym,Sym}})
    2 <= length(sols) <= 200 || return sols
    function within(s, t)
        s === t && return false
        for (k, w) in t
            if haskey(s, k)
                d = try subs(w, s...) - s[k] catch; return false end
                _vanishes(d) || return false
            else
                _eq(w, k) || return false
            end
        end
        true
    end
    keep = [!any(t -> within(s, t), sols) for s in sols]
    # (two families each within the other, the same set: the first kept)
    for i in eachindex(sols), j in 1:(i - 1)
        keep[i] || keep[j] || !(within(sols[i], sols[j]) && within(sols[j], sols[i])) || (keep[j] = true)
    end
    sols[keep]
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


# log(A) = log(B) as A = B: the principal logarithm is one to one on the nonzero
# numbers, so the two are equivalent where A ≠ 0 (SymPy #23662). The equation, and the
# A that must not vanish, or nothing when eq is not of that form.
function _log_eq(eq::Relation, vars::Vector{Sym})
    d = expand(eq.lhs - eq.rhs)
    d.kind == ADD && iszero(d.num) && length(d.args) == 2 || return nothing
    (c1, t1), (c2, t2) = zip(d.coefs, d.args)
    c1 == -c2 || return nothing
    islog(t) = t.kind == FN && t.name === :log && length(t.args) == 1
    (islog(t1) && islog(t2)) || return nothing
    A, B = t1.args[1], t2.args[1]
    any(v -> has(A, v) || has(B, v), vars) || return nothing
    Relation(:eq, A, B), A
end

# a system that is not polynomial, by elimination: an unknown that an equation gives
# as a rational function of the others (x = 0.0015 z) substituted into the rest; or the
# equations that are polynomial solved first, their solutions (families too) put into
# the others (2x/y - 2xy = 0 gives y = ±1, then log(1/x) + 1 = 0; SymPy #22723,
# #23318). nothing when neither applies.
function _solve_by_elimination(eqs::Vector{<:Relation}, vars::Vector{Sym})
    length(vars) >= 2 || return nothing
    # floats: the system with the rationals they stand for, the solutions floats again
    if !all(e -> _isexactex(e.lhs) && _isexactex(e.rhs), eqs)
        ex = Relation[Relation(:eq, _exactfloats(e.lhs; short = true), _exactfloats(e.rhs; short = true)) for e in eqs]
        r = _solve_by_elimination(ex, vars)
        r === nothing && return nothing
        num(v) = _hassymbol(v) ? _floatify(v) : (w = N(v); w isa Number ? Sym(w) : v)
        return [Dict{Sym,Sym}(k => num(v) for (k, v) in s) for s in r]
    end
    ds = [e.lhs - e.rhs for e in eqs]
    function ispoly(e)
        atoms = Set{Sym}()
        n, d = numer_denom(cancel(_exactfloats(e)))         # (floats count as numbers)
        _poly_atoms!(atoms, n) && _poly_atoms!(atoms, d) &&
            !any(a -> !any(v -> _eq(a, v), vars) && any(v -> has(a, v), vars), atoms)
    end
    # an unknown of degree 1 in a polynomial equation, its coefficient free of the
    # unknowns: v = -c₀/c₁
    for (i, e) in enumerate(ds)
        ispoly(e) || continue
        n, d = numer_denom(cancel(e))
        for v in vars
            has(n, v) || continue
            P = try Poly{Sym}(n, [v]) catch; continue end
            degree(P, v) == 1 || continue
            c1, c0 = coeff(n, v, 1), coeff(n, v, 0)
            any(u -> has(c1, u), vars) && continue
            iszero(c1) && continue
            r = cancel(-c0 / c1)
            others = [u for u in vars if !_eq(u, v)]
            rest = Relation[Relation(:eq, subs(eqs[j].lhs, v => r), subs(eqs[j].rhs, v => r)) for j in eachindex(eqs) if j != i]
            rest = [f for f in rest if !iszero(simplify(f.lhs - f.rhs))]
            sub = isempty(rest) ? [Dict{Sym,Sym}()] :
                  try _solve_generic(rest, others) catch err; err isa ArgumentError || rethrow(); return nothing end
            sub isa AbstractVector || return nothing
            out = Dict{Sym,Sym}[]
            for s in sub
                s = copy(s)
                s[v] = cancel(subs(r, s...))
                _defined_in(eqs, s, vars) && push!(out, s)
            end
            return out
        end
    end
    # the polynomial equations first
    P = [i for i in eachindex(ds) if ispoly(ds[i])]
    (isempty(P) || length(P) == length(ds)) && return nothing
    pvars = [v for v in vars if any(i -> has(ds[i], v), P)]
    psols = try _solve_generic(eqs[P], pvars) catch err; err isa ArgumentError || rethrow(); return nothing end
    psols isa AbstractVector || return nothing
    out = Dict{Sym,Sym}[]
    for s in psols
        fixed = Dict{Sym,Sym}(k => v for (k, v) in s if !_eq(k, v))
        isempty(fixed) && return nothing
        rest = Relation[]
        ok = true
        for j in eachindex(eqs)
            j in P && continue
            f = try Relation(:eq, subs(eqs[j].lhs, fixed...), subs(eqs[j].rhs, fixed...)) catch err
                _isdivzero(err) || rethrow()
                ok = false; break
            end
            g = f.lhs - f.rhs
            if any(v -> has(g, v), vars)
                push!(rest, f)
            else
                # (undefined, log(1/0), or false: not a solution)
                _defined_at_all(g) && iszero(simplify(g)) || (ok = false; break)
            end
        end
        ok || continue
        left = [v for v in vars if !haskey(fixed, v)]
        sub = isempty(rest) ? [Dict{Sym,Sym}()] :
              try _solve_generic(rest, left) catch err; err isa ArgumentError || rethrow(); continue end
        sub isa AbstractVector || continue
        for t in sub
            d = Dict{Sym,Sym}(k => cancel(subs(v, t...)) for (k, v) in fixed)
            merge!(d, t)
            _defined_in(eqs, d, vars) && push!(out, d)
        end
    end
    out
end
_defined_at_all(g::Sym) = (v = try N(g) catch; return true end; !(v isa Number) || isfinite(v))

function _drop_exp_factors(eq::Relation, vars::Vector{Sym})
    g = eq.lhs - eq.rhs
    g.kind == MUL || return eq
    keep = [!(b === E && any(v -> has(p, v), vars)) for (b, p) in zip(g.args, g.exps)]
    all(keep) && return eq
    Relation(:eq, Sym(g.num) * prod((b^p for (b, p, k) in zip(g.args, g.exps, keep) if k); init = ONE), ZERO)
end

# lhs - rhs as a - b with a, b positive products of positive numbers to powers with the
# unknowns: log a = log b, expanded (linear when the exponents are); nothing otherwise
function _exp_products_logged(eq::Relation, vars::Vector{Sym})
    e = expand(eq.lhs - eq.rhs)
    (e.kind == ADD && iszero(e.num) && length(e.args) == 2) || return nothing
    (c1, t1), (c2, t2) = zip(e.coefs, e.args)
    (c1 > N0) == (c2 > N0) && return nothing
    # log of a positive rational with its prime powers apart: log 32 = 5 log 2
    function lg(r)
        r.kind == NUMBER && isexact(r.num) || return log(r)
        out = ZERO
        for (v, sg) in ((_numerator(r.num), 1), (_denominator(r.num), -1))
            v == 1 && continue
            p, k = _int_root_power(BigInt(v))
            out += sg * k * log(Sym(SymNumber(p)))
        end
        out
    end
    function logof(c, t)
        s = lg(Sym(abs(c)))
        fs = t.kind == MUL ? collect(zip(t.args, t.exps)) : t.kind == POW ? [(t.args[1], t.args[2])] : [(t, ONE)]
        t.kind == MUL && (s += lg(Sym(t.num)))
        for (b, q) in fs
            (b.kind == NUMBER && isexact(b.num) && b.num > N0) || return nothing
            s += q * lg(b)
        end
        s
    end
    any(v -> has(t1, v), vars) && any(v -> has(t2, v), vars) || return nothing
    l1, l2 = logof(c1, t1), logof(c2, t2)
    (l1 === nothing || l2 === nothing) && return nothing
    Relation(:eq, expand(l1 - l2), ZERO)
end

