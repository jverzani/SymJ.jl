# Linear systems: Gauss–Jordan elimination over an exact field, SymNumber for numeric
# coefficients and RatFun (rational functions in the other symbols and parameters)
# otherwise; RatFun's canonical form makes the pivot zero tests exact. With parameters
# the result is the generic solution (pivots such as `a` assumed nonzero), split by
# cases on those pivots by `solve`.

"""
    solve(eqs::AbstractVector{<:Relation}, vars)

Solve a system of equations (written with `~`) for the variables `vars`. The result is a
vector of solutions, each a `Dict` with every one of `vars` as a key; an inconsistent system
gives an empty vector. Both arguments must be vectors (or tuples): one equation in one
unknown is `solve(eq, x)`, which gives a vector of values instead.

- **Linear systems** are solved by Gauss–Jordan elimination. An underdetermined system
  has one solution, with the pivot variables given in terms of the free ones, which map
  to themselves.
- **Polynomial (and rational) systems** are solved through a lex Gröbner basis
  ([`groebner`](@ref)), which is triangular, by back-substitution. The values are exact
  where the equations in one variable along the way can be solved exactly, numeric
  otherwise. With infinitely many solutions, the variables left free map to
  themselves.

Coefficients may involve parameters. Where the solution divides by an expression in
them (a pivot of the elimination, a coefficient of the Gröbner basis, a denominator of
a solution) that could vanish, the result is a `Cases` value, the case where it
vanishes solved again with it substituted; a case that cannot be solved so holds an
error. `assume` for a system states facts about parameters and unknowns: `x ∈ ℝ` makes
`abs(x)`, `real(x)` real-valued (and keeps only real solutions), `x ≫ 0` keeps the
solutions with `x > 0`; without it the unknowns are complex.

```julia
solve([a*x + y ~ 1, x - y ~ 2], [x, y])
# a == 0 => [Dict(x => 3, y => 1)], a == -1 => [], a != 0, a != -1 => [Dict(x => 3/(a + 1), y => (1 - 2a)/(a + 1))]
```

```julia
solve([x + y ~ 3, x - y ~ 1], [x, y])          # [Dict(x => 2, y => 1)]
solve([x + y + z ~ 1, x - z ~ 0], [x, y, z])    # [Dict(x => z, y => 1 - 2z, z => z)]
solve([x^2 + y^2 ~ 1, x ~ y], [x, y])          # [Dict(x => -√2/2, y => -√2/2), Dict(x => √2/2, y => √2/2)]
solve([x*y ~ 1, x + y ~ 3], [x, y])
```
"""
function solve(eqs::AbstractVector{<:Relation}, vars::Union{AbstractVector,Tuple}; assume = nothing)
    # assumptions: about the unknowns, they choose the domain (x ∈ ℝ: abs(x) as ±x, real
    # solutions only) and filter the solutions; about parameters, they decide cases
    if assume !== nothing
        A = _assumptions(assume)
        A === nothing && (A = Assumptions())
        vs = Sym[Sym(v) for v in vars]
        r = with(() -> solve(eqs, vars), _SYSTEM_ASSUME => A)
        keep(sols) = sols isa AbstractVector ? filter(s -> !(s isa AbstractDict) || _system_sol_ok(s, vs, A), sols) : sols
        return r isa Cases ? _map_cases(keep, r) : keep(r)
    end
    vars = Sym[Sym(v) for v in vars]
    all(_isnamed, vars) || return with(_NOT_PARAMS => unique!(reduce(vcat, map(_named, vars)))) do
        _solve_for_calls((es, vs, _) -> solve(es, vs), collect(Relation, eqs), vars)
    end
    all(e -> e.op === :eq, eqs) || return _solve_ineq_system(eqs, vars)
    r = _solve_float_constants(eqs, vars)
    r === nothing || return r
    A0 = _SYSTEM_ASSUME[]
    _by_cases(A0 === nothing ? Assumptions() : A0) do B
        _solve_system_case(collect(Relation, eqs), vars, B)
    end
end

# a solution against the assumptions about the unknowns: a real one's value real, one in a
# set within it (numbers only; one with symbols kept)
function _system_sol_ok(s::AbstractDict, vars::Vector{Sym}, A)
    for v in vars
        w = get(s, v, v)
        (_eq(w, v) || _hassymbol(w)) && continue
        z = try _evalc(w) catch; continue end
        if _isreal(v, A) === true
            abs(imag(z)) > 1e-10 * max(1, abs(z)) && return false
        end
        S = get(A.sets, v, nothing)
        if S !== nothing && abs(imag(z)) <= 1e-10 * max(1, abs(z))
            inS = try w in S catch; true end
            inS || return false
        end
    end
    true
end

# floats among parameters: each float an opaque constant c, the system solved exactly in
# the c's, then the floats put back: conditions of the c's alone decided by their
# values, and the answer's numbers floats. (Exactly, a float such as 2^(1/6) became a
# 53-bit rational whose powers nearly cancelled, and cases split on that; SymPy #24193.)
# nothing when there are no such floats.
function _solve_float_constants(eqs, vars::Vector{Sym})
    exprs = Sym[e.lhs - e.rhs for e in eqs]
    any(e -> !_isexactex(e), exprs) || return nothing
    any(e -> any(p -> !any(v -> _eq(p, v), vars), _named(e)), exprs) || return nothing   # no parameters
    floats = SymNumber[]
    for e in exprs
        _replace_walk(e) do t
            # (by absolute value: -1.12 and 1.12 are one constant and its negative)
            cs0 = t.kind == NUMBER ? [t.num] : t.kind == ADD ? [t.num; t.coefs] : t.kind == MUL ? [t.num] : SymNumber[]
            for c in cs0
                isfloat(c) && !iszero(c) && !isone(abs(c)) && !any(f -> f == abs(c), floats) && push!(floats, abs(c))
            end
            nothing
        end
    end
    isempty(floats) && return nothing
    length(floats) <= 4 || return nothing
    cs = [Sym(Symbol("_fc", i)) for i in eachindex(floats)]
    sw(c::SymNumber) = (i = isfloat(c) ? findfirst(f -> f == abs(c), floats) : nothing;
                        i === nothing ? Sym(c) : (c < N0 ? -cs[i] : cs[i]))
    swap(e) = _replace_walk(e) do t
        if t.kind == NUMBER && isfloat(t.num)
            return sw(t.num)
        elseif t.kind == ADD && (any(isfloat, t.coefs) || isfloat(t.num))
            return sw(t.num) + _add_many(sw(c) * a for (c, a) in zip(t.coefs, t.args))
        elseif t.kind == MUL && isfloat(t.num)
            return sw(t.num) * (t / Sym(t.num))
        end
        nothing
    end
    es = Relation[Relation(:eq, swap(e.lhs), swap(e.rhs)) for e in eqs]
    any(e -> !_isexactex(e.lhs - e.rhs), es) && return nothing
    r = _by_cases(Assumptions()) do B
        _solve_system_case(es, vars, B)
    end
    back = [c => Sym(f) for (c, f) in zip(cs, floats)]
    putback(v) = v isa Sym ? _floatify(subs(v, back...)) :
                 v isa AbstractDict ? Dict{Sym,Sym}(k => putback(w) for (k, w) in v) :
                 v isa AbstractVector ? map(putback, v) : v
    r isa Cases || return putback(r)
    out = Pair{Vector{Any},Any}[]
    for (conds, v) in r.cases
        keep = true
        cs2 = Any[]
        for c in conds
            c isa Relation || (push!(cs2, c); continue)
            cr = Relation(c.op, subs(c.lhs, back...), subs(c.rhs, back...))
            d = expand(cr.lhs - cr.rhs)
            if !_hassymbol(d)
                holds = _REL_INFO[cr.op][2](_csign(d) === nothing ? sign(real(N(d))) : _csign(d), 0)
                holds || (keep = false; break)
            else
                push!(cs2, cr)
            end
        end
        keep && push!(out, cs2 => putback(v))
    end
    isempty(out) && return Dict{Sym,Sym}[]
    length(out) == 1 && isempty(out[1].first) && return out[1].second
    Cases{Any}(out)
end

# the divisions the generic solution made (pivots, Gröbner coefficients), collected
const _DIVISORS = ScopedValue{Union{Nothing,Vector{Sym}}}(nothing)
_divisor!(c) = (D = _DIVISORS[]; D === nothing || push!(D, Sym(c)); nothing)

# the system in the case B: no solutions where a denominator free of the unknowns
# vanishes; the generic solutions where none of the divisors does; else the system again
# with the vanishing one substituted
function _solve_system_case(eqs::Vector{Relation}, vars::Vector{Sym}, B)
    sub = [k => w for (k, w) in B.subst if !any(v -> _eq(v, k), vars)]
    es = try
        isempty(sub) ? eqs : Relation[Relation(:eq, subs(e.lhs, sub...), subs(e.rhs, sub...)) for e in eqs]
    catch err
        _isdivzero(err) || rethrow()
        return Dict{Sym,Sym}[]
    end
    for e in es, f in _param_factors(numer_denom(cancel(e.lhs - e.rhs))[2], vars)
        _decide(Relation(:ne, f, ZERO), B) || return Dict{Sym,Sym}[]
    end
    # an equation free of the unknowns is a condition on the parameters: no solutions
    # where it fails, and where it holds the parameter values it gives substituted, each
    # a case ((k + 1)² = 4(k + 4) beside a quadratic in x: k = -3 and k = 5; as a
    # coefficient field it was a nonzero constant, so no solutions; SymPy #19564)
    for e in es
        any(v -> has(e.lhs, v) || has(e.rhs, v), vars) && continue
        g = numer_denom(cancel(e.lhs - e.rhs))[1]
        iszero(g) && continue
        _hassymbol(g) || return Dict{Sym,Sym}[]
        # each factor decided: one linear in a parameter that vanishes is then substituted
        # by the case itself (k = -3), so the equation is gone when the system is solved
        # again; another's roots are the alternatives
        for f in _param_factors(g, vars)
            _decide(Relation(:ne, f, ZERO), B) && continue
            ps = _degenerate_subst([f], B, vars)
            ps === nothing && return Dict{Sym,Sym}[]
            rest = Relation[q for q in es if q !== e]
            return _alternatives(ps) do p
                _solve_system_case(Relation[Relation(:eq, subs(q.lhs, p), subs(q.rhs, p)) for q in rest], vars, B)
            end
        end
        return Dict{Sym,Sym}[]
    end
    D = Sym[]
    sols = with(() -> _solve_generic(es, vars), _DIVISORS => D)
    for s in sols, v in vars
        push!(D, numer_denom(s[v])[2])
    end
    ps = _degenerate_subst(D, B, vars)
    ps === nothing && return sols
    _alternatives(ps) do p
        _solve_system_case(Relation[Relation(:eq, subs(e.lhs, p), subs(e.rhs, p)) for e in es], vars, B)
    end
end

# the generic solutions, coefficients with parameters taken to be nonzero (what internal
# callers use)
function _solve_generic(eqs::AbstractVector{<:Relation}, vars::Union{AbstractVector,Tuple})
    vars = Sym[Sym(v) for v in vars]
    all(_isnamed, vars) ||
        return _solve_for_calls((es, vs, _) -> _solve_generic(es, vs), collect(Relation, eqs), vars)
    all(e -> e.op === :eq, eqs) || return _solve_ineq_system(eqs, vars)
    sols = if _islinear(eqs, vars)
        sol = _solve_linear(eqs, vars)
        sol === nothing ? Dict{Sym,Sym}[] : [sol]
    else
        _solve_nonlinear(collect(Relation, eqs), vars)
    end
    # every solution has every unknown as a key; a free one stands for itself
    for s in sols, v in vars
        haskey(s, v) || (s[v] = v)
    end
    filter!(s -> _defined_in(eqs, s, vars), sols)
end

# is every equation defined at the solution s? cancel can remove a factor of a
# denominator, so (x^2 - 1)/(x - 1) = 2 with y = x has no solution, not x = y = 1.
# (Only decided where the expressions that must not vanish become constants: with
# x free, x ≠ 1 is not recorded.)
function _defined_in(eqs, s::Dict{Sym,Sym}, vars::Vector{Sym})
    ps = [v => get(s, v, v) for v in vars]
    for eq in eqs
        e = eq.lhs - eq.rhs
        zs = Sym[]
        foreach(v -> _zero_sets!(zs, e, v), vars)
        for b in unique(zs)
            v = try subs(b, ps...) catch err; _isdivzero(err) && return false; rethrow() end
            (iszero(v) || (!_hassymbol(v) && _iszero_value(v))) && return false
        end
    end
    true
end

# every equation of degree at most 1 in the unknowns, with no unknown in a denominator
# or inside a function
function _islinear(eqs, vars::Vector{Sym})
    for eq in eqs
        n, d = numer_denom(cancel(eq.lhs - eq.rhs))
        any(v -> has(d, v), vars) && return false
        atoms = Set{Sym}()
        _poly_atoms!(atoms, _exactfloats(n)) || return false
        any(a -> !any(v -> _eq(a, v), vars) && any(v -> has(a, v), vars), atoms) && return false
        P = Poly{Sym}(n, vars)
        all(ev -> sum(ev) <= 1, exponent_vectors(P)) || return false
    end
    true
end

# (denominators smaller by trigonometric identities, cos² + sin² = 1, reduced; SymPy #23140)
function _solve_linear(eqs, vars::Vector{Sym})
    r = _solve_linear0(eqs, vars)
    (r isa AbstractDict && any(v -> _has_fn(v, _CIRCFNS), values(r))) || return r
    Dict{Sym,Sym}(k => _trig_den(v) for (k, v) in r)
end
function _solve_linear0(eqs, vars::Vector{Sym})
    A, b = _linear_system(eqs, vars)
    # square, with parameters: Cramer's rule, whose one divisor is det A (the pivots of
    # Gauss–Jordan are divisors too, each a case to decide, and its rational functions
    # took a gcd at every step: hours for a circuit in 13 parameters; SymPy #24135)
    if size(A, 1) == size(A, 2) <= _CRAMER_MAX && any(_hassymbol, A)
        r = _cramer_sym(hcat(A, b), size(A, 1))
        if r !== nothing
            X, D = r
            Dt = _trig_reduced(D)
            _divisor!(Dt)
            return Dict{Sym,Sym}(v => X[i, 1] for (i, v) in enumerate(vars))
        end
    end
    if all(e -> e.kind == NUMBER, Iterators.flatten((A, b)))
        _gauss_jordan(map(e -> e.num, A), map(e -> e.num, b), vars, Sym)
    elseif !any(_hassymbol, Iterators.flatten((A, b))) && any(e -> !_isexactex(e), Iterators.flatten((A, b)))
        # complex numbers with floats (1.5 + 2.0im): solved in complex floating point
        # (SymPy #21445)
        _gauss_jordan(map(e -> ComplexF64(N(e)), A), map(e -> ComplexF64(N(e)), b), vars, Sym)
    else
        atoms = Set{Sym}()
        all(e -> _poly_atoms!(atoms, numer_denom(e)[1]) && _poly_atoms!(atoms, numer_denom(e)[2]),
            Iterators.flatten((A, b))) ||
            throw(ArgumentError("solve: inexact coefficients together with symbols are not supported"))
        # radicals among the coefficients: over the field they generate, its elements
        # normalized (√6√11 and √2√33 are one number; as independent generators a zero
        # pivot looked nonzero and the solutions were 0/0; SymPy #22852)
        if any(_isalgebraic_atom, atoms)
            return _gauss_jordan(map(e -> _AlgF(_algnorm(e)), A), map(e -> _AlgF(_algnorm(e)), b), vars, Sym)
        end
        gens = sort!(collect(atoms); lt = _lt)
        function rf(e)
            n, d = numer_denom(e)
            _ratfun(Poly{SymNumber}(n, gens), Poly{SymNumber}(d, gens))
        end
        _gauss_jordan(map(rf, A), map(rf, b), vars, Sym)
    end
end
solve(eqs::Tuple{Vararg{Relation}}, vars; kws...) = solve(collect(Relation, eqs), vars; kws...)

# d simplified by trigonometric identities when that makes it smaller, else d itself
function _trig_reduced(d::Sym)
    _has_fn(d, _CIRCFNS) || return d
    t = try simplify_trig(d) catch err; err isa InterruptException && rethrow(); return d end
    _nodes(t) < _nodes(d) ? t : d
end

function _trig_den(v::Sym)
    n, d = numer_denom(v)
    dt = _trig_reduced(d)
    dt === d ? v : cancel(n / dt)
end

# A x = b from the equations: each lhs - rhs in lowest terms, no unknown in the
# denominator, of degree at most 1 in the unknowns
function _linear_system(eqs, vars::Vector{Sym})
    m, k = length(eqs), length(vars)
    A = fill(ZERO, m, k)
    b = fill(ZERO, m)
    for (i, eq) in enumerate(eqs)
        n, d = numer_denom(cancel(eq.lhs - eq.rhs))
        any(v -> has(d, v), vars) && throw(ArgumentError("solve: $eq is not linear in $(join(vars, ", "))"))
        P = Poly{Sym}(n, vars)
        for (ev, c) in zip(exponent_vectors(P), coefficients(P))
            s = sum(ev)
            s <= 1 || throw(ArgumentError("solve: $eq is not linear in $(join(vars, ", ")); nonlinear systems are not supported yet"))
            c = c / d
            s == 0 ? (b[i] = -c) : (A[i, findfirst(==(1), ev)] = c)
        end
    end
    A, b
end

# Gauss–Jordan over the field of the entries (SymNumber or RatFun)
function _gauss_jordan(A::Matrix{C}, b::Vector{C}, vars::Vector{Sym}, ::Type{Sym}) where {C}
    m, k = size(A)
    M = hcat(A, b)
    # floating point entries: a residue of rounding (5/3.0 - 1.6666666666666667) is zero,
    # not 0 = nonzero (SymPy #26223); the largest pivot is taken
    # (complex floats, C = ComplexF64, likewise; SymPy #21445)
    cx = C === ComplexF64
    fl = cx || (C === SymNumber && any(e -> isfloat(e), M))
    fabs(e) = cx ? abs(e) : abs(Float64(e))
    scale = fl ? maximum(fabs, M; init = 0.0) : 0.0
    # (an infinite or NaN entry: no tolerance means anything; every entry was zero below
    # it and the system "solved" with all unknowns free; SymPy #15113)
    fl && !isfinite(scale) && throw(ArgumentError("solve: a coefficient of the linear system is not finite"))
    nz(e) = !iszero(e) && !(fl && (cx || isfloat(e)) && fabs(e) <= 1e-12 * scale)
    r = 1
    pivots = Tuple{Int,Int}[]
    for col in 1:k
        r > m && break
        p = fl ? (is = [i for i in r:m if nz(M[i, col])];
                  isempty(is) ? nothing : is[argmax([fabs(M[i, col]) for i in is])] - r + 1) :
                 findfirst(i -> nz(M[i, col]), r:m)
        p === nothing && continue
        p += r - 1
        if p != r
            M[r, :], M[p, :] = M[p, :], M[r, :]
        end
        C === SymNumber || cx || _divisor!(M[r, col])
        piv = inv(M[r, col])
        M[r, :] = [e * piv for e in M[r, :]]
        for i in 1:m
            (i == r || !nz(M[i, col])) && continue
            f = M[i, col]
            M[i, :] = [M[i, j] - f * M[r, j] for j in 1:(k + 1)]
            fl && (M[i, col] = zero(M[i, col]))
        end
        fl && (M .= [nz(e) ? e : zero(e) for e in M])
        push!(pivots, (r, col))
        r += 1
    end
    if any(i -> nz(M[i, k+1]), r:m)                            # 0 = nonzero: inconsistent
        C === SymNumber || cx || foreach(i -> iszero(M[i, k+1]) || _divisor!(M[i, k+1]), r:m)
        return nothing
    end
    pivcols = Set(c for (_, c) in pivots)
    free = [j for j in 1:k if !(j in pivcols)]
    sol = Dict{Sym,Sym}()
    for (i, c) in pivots
        v = Sym(M[i, k+1])
        for j in free
            iszero(M[i, j]) || (v -= Sym(M[i, j]) * vars[j])
        end
        sol[vars[c]] = cancel(v)
    end
    sol
end
