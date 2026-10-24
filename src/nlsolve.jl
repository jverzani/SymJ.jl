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
    sols = isempty(nums) ? [Dict{Sym,Sym}()] : _solve_groebner(nums, vars)
    # drop solutions at a pole
    filter!(s -> all(d -> !_vanishes(subs(d, s...)), dens), sols)
    if inexact
        sols = [Dict{Sym,Sym}(k => (_hassymbol(v) ? v : Sym(N(v))) for (k, v) in s) for s in sols]
    end
    _uniquesols(sols)
end

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
    (_isexactex(q) || _hassymbol(subs(q, x => ZERO))) && return solve(q ~ 0, x)
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
