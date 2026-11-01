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
    atoms = Set{Sym}()
    foreach(n -> _poly_atoms!(atoms, n), nums)
    algs = sort!([a for a in atoms if _isalgconst(a)]; lt = _lt)
    isempty(algs) && throw(ArgumentError(_GROEBNER_ALGEBRAIC_MSG))
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
    sols = _solve_groebner(eqs, [vars; rs])
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

function _solve_groebner(nums::Vector{Sym}, vars::Vector{Sym})
    G = groebner(nums, vars...)
    G == [ONE] && return Dict{Sym,Sym}[]
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
                Gf == [ONE] || append!(out, _backsub(Gf, vars, n, Dict{Sym,Sym}()))
            end
            return out
        end
    end
    _backsub(G, vars, n, Dict{Sym,Sym}())
end

# p a polynomial in x alone, with rational coefficients
_onlyvar(p::Sym, x::Sym) = (atoms = Set{Sym}(); _poly_atoms!(atoms, p) && all(a -> _eq(a, x), atoms))

# the index of the first of vars that p involves (0 if none)
_mainvar(p::Sym, vars::Vector{Sym}) = something(findfirst(v -> has(p, v), vars), 0)

function _backsub(G::Vector{Sym}, vars::Vector{Sym}, k::Int, sol::Dict{Sym,Sym})
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
    isempty(qs) && return _backsub(G, vars, k - 1, sol)            # x is free
    sort!(qs; by = q -> _xdegree(q, x))
    out = Dict{Sym,Sym}[]
    for r in _solve1(qs[1], x)
        all(q -> _satisfies(q ~ 0, x, r), qs[2:end]) || continue
        s = copy(sol)
        s[x] = r
        append!(out, _backsub(G, vars, k - 1, s))
    end
    out
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
function _exactfloats(e::Sym)
    k = e.kind
    k == NUMBER && return _exactnum(e.num)
    k == ADD && return _exactnum(e.num) + _add_many(_exactnum(c) * _exactfloats(a) for (c, a) in zip(e.coefs, e.args))
    k == MUL && return _exactnum(e.num) * prod((_exactfloats(a)^_exactfloats(x) for (a, x) in zip(e.args, e.exps)); init = ONE)
    k == POW && return _exactfloats(e.args[1])^_exactfloats(e.args[2])
    k == FN && return apply_function(e.name, Sym[_exactfloats(a) for a in e.args])
    e
end
