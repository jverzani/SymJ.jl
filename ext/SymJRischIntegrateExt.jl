# RischIntegrate is generic over term types; this workload runs it on SymJ terms, so the
# compiled code ships with the extension instead of being compiled on the first integral
# (about 4 s otherwise).
module SymJRischIntegrateExt

using SymJ
using RischIntegrate: RischIntegrate, integrate
using PrecompileTools: @setup_workload, @compile_workload

# piecewise functions, heaviside and dirac: interval by interval (see SymJ's piecewise.jl);
# abs, sign, max and min, which the integrator does not take, as piecewise functions
# (for real x; conj(x), real(x) and imag(x) are x, x and 0, and a radical's real
# factors to even powers come out, √(x²(b² + 1)) = |x|√(b² + 1)). The assumptions
# (assume = …, or an @assuming block) decide the signs of parameter expressions the
# answer's form depends on, through the integrator's sign oracle, when the installed
# RischIntegrate has one; signs they leave open take its default.
# (floats in the integrand, which the integrator takes as the rationals they stand for:
# the answer's numbers floats again, not 27-digit fractions; found in the review)
function RischIntegrate.integrate(ex::Sym, x::Sym; assume = nothing, kw...)
    r = _unwrap(_indefinite(ex, x; assume, kw...))
    # exp(c log u) as u^c, the integrand's own form (xⁿ log x had exp(n log x) back)
    r isa Sym && (r = SymJ._exp_log_powers(r))
    r isa Sym && !SymJ._isexactex(ex) ? SymJ._floatify(r) : r
end

# an answer by cases (an inner definite integral's) integrated case by case, each with
# its conditions among the assumptions; an infinite value or an error stays (SymPy #23610)
function RischIntegrate.integrate(cs::SymJ.Cases, x::Sym; assume = nothing, kw...)
    A0 = something(SymJ._assumptions(assume), SymJ.Assumptions())
    SymJ.Cases{Any}(Any[c.first => _case_integral(c.second, c.first, x, A0; kw...) for c in cs.cases])
end
function _case_integral(v, conds, x::Sym, A0; kw...)
    v isa Exception && return v
    v = Sym(v)
    SymJ._isinfpt(v) && return v
    A = A0
    for c in conds
        B = try SymJ._with_fact(A, c) catch err; err isa ArgumentError || rethrow(); nothing end
        B === nothing || (A = B)
    end
    try
        RischIntegrate.integrate(v, x; assume = A, kw...)
    catch err
        err isa InterruptException && rethrow()
        err
    end
end

# (the definite integrals take the answer as it is: a ₂F₁ term holds only where its
# radicand is positive, which they check numerically)
function _indefinite(ex::Sym, x::Sym; assume = nothing, kw...)
    # a factor free of x outside (1/√(2π) beside √2 made a constant field the integrator
    # does not take; SymPy #22075)
    c, t = SymJ._split_kfree(ex, x)
    if !isone(c) && SymJ._hassymbol(t) && !(c.kind == SymJ.NUMBER)
        r = try _indefinite(t, x; assume, kw...) catch err; err isa InterruptException && rethrow(); nothing end
        r isa Sym && return c * r
        r isa SymJ.Cases && return SymJ._map_cases(v -> v isa Sym ? c * v : v, r)
    end
    A = SymJ._assumptions(assume)
    if A !== nothing && _has_oracle() && !haskey(kw, :sign)
        kw = (; sign = SymJ._sign_oracle(A), kw...)
    end
    # values of unknown functions free of x (f(1)) are constants: symbols in their place
    # while integrating (RischIntegrate takes no SymFunction; SymPy #14796)
    ucs = SymJ._unknown_consts!(Sym[], ex, x)
    if !isempty(ucs)
        ps = [c => Sym(Symbol("ω_f", i)) for (i, c) in enumerate(ucs)]
        r = _indefinite(subs(ex, ps...), x; assume, kw...)
        back = [p.second => p.first for p in ps]
        r isa SymJ.Cases && return SymJ._map_cases(v -> v isa Sym ? subs(v, back...) : v, r)
        return r isa Sym ? subs(r, back...) : r
    end
    # terms with an unknown function of x: held as ∫ f(x) dx, the rest integrated
    rest, held = SymJ._split_opaque(ex, x)
    if !isempty(held)
        # (∫ u'(x) dx = u(x): dsolve(f' ~ a g', f(x)) kept ∫ g'(x) dx; SymPy #13691)
        H = SymJ._add_many(c * (SymJ._isderiv_in(g, x) ? g.args[1] : SymJ._mkfn_node(:Integral, Sym[g, x])) for (c, g) in held)
        return iszero(rest) ? H : RischIntegrate.integrate(rest, x; assume, kw...) + H
    end
    integ = f -> RischIntegrate.risch_integrate(f, x; kw...)
    ex = SymJ._realsyms(ex, (x,))
    Ax = SymJ._with_fact(something(A, SymJ.Assumptions()), x ∈ SymJ.ℝ)
    Ax === nothing || (ex = SymJ._split_radicals(ex, Ax))
    # real and imaginary parts the assumptions resolve (|exp(w)|² = exp(2 re w) with
    # re w a rational function; SymPy #26102)
    if Ax !== nothing && SymJ._has_fn(ex, (:real, :imag, :conj, :abs))
        ex = SymJ._refine(ex, Ax)
    end
    SymJ._haspieces(ex, x) && return SymJ._integrate_pieces(ex, x, integ)
    try
        integ(ex)
    catch err
        err isa InterruptException && rethrow()
        pw = SymJ.piecewise(ex)
        if !SymJ._haspieces(pw, x)
            # |z| of a complex z with parameters, as for the pieces of |x - y|, for real
            # parameters: |exp(2πix) + exp(2πiy)| is √2 √(cos(2π(x - y)) + 1). (Only with
            # RischIntegrate's answers for √(1 + cos x) continuous across its zeros, from
            # its sign_jump_correction: before, a definite integral over part of a period
            # took the jump, ∫₀^½ at y = 0.8 came out -0.889, not 0.384.)
            isdefined(RischIntegrate, :sign_jump_correction) || rethrow()
            Ar = Ax
            for s in SymJ.free_symbols(ex)
                Ar === nothing && break
                Ar = SymJ._with_fact(Ar, s ∈ SymJ.ℝ)
            end
            (Ar === nothing || !SymJ._has_fn(ex, (:abs,))) && rethrow()
            er = SymJ._refine(ex, Ar)
            isequal(er, ex) && rethrow()
            return integ(er)
        end
        SymJ._integrate_pieces(pw, x, integ)
    end
end

# an answer with no root sums is its elementary part, an expression like any other
# (x*₂F₁(…) for 1/sqrt(1 + x^3), whose limit failed on the wrapper; SymPy #24879)
# A ₂F₁ term (HypTerm, c α^q x^(m+1)/(m+1) ₂F₁(-q, (m+1)/n; (m+1)/n + 1; -βxⁿ/α)) is
# such an expression too.
function _unwrap(r)
    r isa RischIntegrate.IntegralWithRootSums || return r
    isdefined(RischIntegrate, :HypTerm) || return r                      # (RischIntegrate 1.0)
    (r.elementary isa Union{Sym,Number} && all(t -> t isa RischIntegrate.HypTerm, r.rootsums)) || return r
    Sym(r.elementary) + SymJ._add_many(_hypterm(h) for h in r.rootsums)
end
function _hypterm(h)
    q(v) = Sym(SymJ.SymNumber(v))
    x = Sym(Symbol(h.var))
    a = q(h.m + 1) / h.n
    q(h.c) * q(h.α)^q(h.q) * x^(h.m + 1) / (h.m + 1) *
        SymJ.apply_function(:_₂F₁, Sym[-q(h.q), a, a + 1, -q(h.β) / q(h.α) * x^h.n])
end

# the ₂F₁ form of an answer on (lo, hi), where each term's radicand α + βxⁿ is positive
# (no real root inside, positive in the middle); nothing otherwise
function _symbolic_on(F, x::Sym, lo::Sym, hi::Sym)
    (F isa RischIntegrate.IntegralWithRootSums && isdefined(RischIntegrate, :HypTerm)) || return nothing
    (all(t -> t isa RischIntegrate.HypTerm, F.rootsums) && !isempty(F.rootsums)) || return nothing
    (SymJ._hassymbol(lo) || SymJ._hassymbol(hi)) && return nothing
    l, h = Float64(N(lo)), Float64(N(hi))
    for t in F.rootsums
        α, β, n = Float64(t.α), Float64(t.β), t.n
        r(v) = α + β * v^n
        # the real roots of α + β xⁿ
        q = -α / β
        roots = isodd(n) ? [sign(q) * abs(q)^(1 / n)] : q > 0 ? [q^(1 / n), -q^(1 / n)] : Float64[]
        any(z -> l < z < h, roots) && return nothing
        m = isfinite(l) && isfinite(h) ? (l + h) / 2 : isfinite(l) ? l + 1 : isfinite(h) ? h - 1 : 0.0
        r(m) > 0 || return nothing
    end
    _unwrap(F)
end

# does risch_integrate take a sign oracle (RischIntegrate with PR #70)?
const _ORACLE = Ref{Union{Nothing,Bool}}(nothing)
function _has_oracle()
    _ORACLE[] === nothing && (_ORACLE[] = any(m -> :sign in Base.kwarg_decl(m), methods(RischIntegrate.risch_integrate)))
    _ORACLE[]::Bool
end

# definite integrals: split at the integrand's singular points, one-sided limits at them
# and at ±∞, by cases on parameters; principal_value = true takes Cauchy principal values
# at the singular points inside and at ±∞ (see SymJ's definite.jl)
# (an inner integral by cases, as the integrand of an outer one: the cases the outer
# interval allows). drop_removable_singularities = true merges a case p = q into the
# general formula when its value is the formula's limit there (SymJ._drop_removable).
function RischIntegrate.integrate(ex::Union{Sym,SymJ.Cases}, xab::Tuple, more::Tuple...; assume = nothing,
                                  principal_value::Bool = false, drop_removable_singularities::Bool = false, kw...)
    antider = (f, x, A) -> _indefinite(f, x; assume = A, kw...)
    # each integral with the outer variables in their intervals (numeric ends): the inner
    # one's cases for z < -y - 2 drop out for y, z in (0, 1) (SymPy #14096)
    function within(ts)
        facts = Any[]
        for t in ts
            length(t) == 3 || continue
            y, a, b = Sym(t[1]), Sym(t[2]), Sym(t[3])
            (SymJ._hassymbol(a) || SymJ._hassymbol(b)) && continue
            lo, hi = SymJ._less(b, a, SymJ.Assumptions()) ? (b, a) : (a, b)
            SymJ._isinfpt(lo) || push!(facts, SymJ.Relation(:ge, y, lo))
            SymJ._isinfpt(hi) || push!(facts, SymJ.Relation(:le, y, hi))
        end
        isempty(facts) && return assume
        A = SymJ._assumptions(assume)
        A === nothing && (A = SymJ.Assumptions())
        B = SymJ._with_fact(A, SymJ.FactAnd(facts))
        B === nothing ? assume : B                      # (contradicting ones: as given)
    end
    v = SymJ._definite_integral(ex, xab, within(more), antider; principal_value)
    for (i, t) in enumerate(more)
        v = SymJ._definite_integral(v, t, within(more[(i + 1):end]), antider; principal_value)
    end
    v = drop_removable_singularities ? SymJ._drop_removable(v) : v
    # (floats in the integrand: the value's numbers floats, as for the antiderivative)
    if ex isa Sym && !SymJ._isexactex(ex)
        fl(w) = w isa Sym ? SymJ._floatify(w) : w
        v = v isa SymJ.Cases ? SymJ._map_cases(fl, v) : fl(v)
    end
    v
end

# RischIntegrate's special functions on expressions, so answers with them are expressions
# (not an IntegralWithRootSums of Julia Exprs): Clausen functions, Re/Im of polylogarithms
# and of ₂F₁, the hyperbolic sine and cosine integrals. Derivatives in the angle φ:
#   Cl′ₙ = Cl_{n-1} (n even), -Cl_{n-1} (n odd), Cl₁(φ) = -log(2 - 2cos φ)/2;
#   ∂φ Re Liₙ(ρe^{iφ}) = -Im Li_{n-1}, ∂φ Im Liₙ = Re Li_{n-1}, with
#   Li₁ = -log(1 - z): Re = -log(1 - 2ρ cos φ + ρ²)/2, Im = atan(ρ sin φ/(1 - ρ cos φ)).
# (RischIntegrate 1.2 on; older versions, as on Julia 1.10, have none of them)
const RI_FUNCTIONS = [:clausen2 => 1, :clausen => 2, :repolylog => 3, :impolylog => 3,
                      :rehyp2f1 => 4, :imhyp2f1 => 4, :sinhint => 1, :coshint => 1]
_hasri(f) = isdefined(RischIntegrate, f)
for (f, n) in RI_FUNCTIONS, sig in SymJ._symbolic_signatures(n)
    _hasri(f) || continue
    args = [Symbol(:a, i) for i in 1:n]
    @eval RischIntegrate.$f($([:($a::$T) for (a, T) in zip(args, sig)]...)) =
        SymJ.apply_function($(QuoteNode(f)), Sym[$(args...)])
end

_cl1(φ) = -log(2 - 2cos(φ)) / 2
_cl(n, φ) = n == 1 ? _cl1(φ) : n == 2 ? RischIntegrate.clausen2(φ) : RischIntegrate.clausen(n, φ)
_isposint(a::Sym) = a.kind == SymJ.NUMBER && SymJ.isexact(a.num) && SymJ._isint(a.num) && a.num > SymJ.SymNumber(0)
function _dclausen(args, i)
    n, φ = args
    i == 2 && _isposint(n) && n.num > SymJ.SymNumber(1) || return nothing
    k = SymJ._toint(n.num)
    iseven(k) ? _cl(k - 1, φ) : -_cl(k - 1, φ)
end
_reli(n, ρ, φ) = n == 1 ? -log(1 - 2ρ * cos(φ) + ρ^2) / 2 : RischIntegrate.repolylog(n, ρ, φ)
_imli(n, ρ, φ) = n == 1 ? atan(ρ * sin(φ) / (1 - ρ * cos(φ))) : RischIntegrate.impolylog(n, ρ, φ)
function _dpolylog(re::Bool)
    (args, i) -> begin
        n, ρ, φ = args
        i == 3 && _isposint(n) && n.num > SymJ.SymNumber(1) || return nothing
        k = SymJ._toint(n.num)
        re ? -_imli(k - 1, ρ, φ) : _reli(k - 1, ρ, φ)
    end
end

function __init__()
    SymJ._SYMBOLIC_ON[] = _symbolic_on
    reg(name; kw...) = _hasri(name) && SymJ.register_function!(name; julia = getfield(RischIntegrate, name), kw...)
    # Cl₂(kπ) = 0 for integers k (sin(jkπ) = 0 in its series; ∫₀^(π/2) log sin x came out
    # with clausen2(π) left in; SymPy #16084)
    reg(:clausen2; diff = (args, i) -> _cl1(args[1]),
        eval = args -> (r = SymJ.cancel(args[1] / SymJ.PI); r.kind == SymJ.NUMBER && SymJ.isexact(r.num) && SymJ._isint(r.num) ? SymJ.ZERO : nothing))
    reg(:clausen; diff = _dclausen)
    reg(:repolylog; diff = _dpolylog(true))
    reg(:impolylog; diff = _dpolylog(false))
    reg(:rehyp2f1)
    reg(:imhyp2f1)
    reg(:sinhint; diff = (args, i) -> sinh(args[1]) / args[1])
    reg(:coshint; diff = (args, i) -> cosh(args[1]) / args[1])
    # the integrator SymJ uses inside (dsolve)
    SymJ._INTEGRATOR[] = (f, x) -> integrate(f, x)
end

@setup_workload begin
    @compile_workload begin
        x, a, b, c, d = Sym(:x), Sym(:a), Sym(:b), Sym(:c), Sym(:d)
        for f in (x^3 - 2x + 1, 1 / (x^2 - 1), (2x + 3) / (x^2 + 2x + 5), 1 / (1 + x^3),
                  sqrt(1 + 3x), x / sqrt(4 - x^2), x^2 * (1 + x)^(1 // 3),    # algebraic
                  sqrt(5 - x^2), sqrt(7 - x^2), cos(x) / sqrt(sin(x)^3),
                  sin(x)^2 * cos(x), 1 / (2 + cos(x)), tan(x),                # trigonometric
                  1 / (3 + sin(x) + 2cos(x)), 1 / (7 + 3sin(x) - cos(x)), (sin(x) + 2cos(x)) / (sin(x) - cos(x)),
                  x * exp(2x), exp(x) * sin(x), log(x)^2, x * log(1 + x),     # exp and log
                  sin(x) * exp(cos(x)), exp(atan(x)) / (1 + x^2),
                  atan(x), x * asin(x),                                     # inverse trigonometric
                  1 / (a + x^2), 1 / (b + a * x^2), x^a * log(x),             # parameters
                  x^a * log(b * x), x^(a - 1) * cos(x^a), exp(a * x) * sin(b * x),
                  (a^2 - x^2)^(3 // 2), 1 / (1 + a * sin(x)),
                  (a * sin(x) + cos(x)) / (sin(x) + b * cos(x)), 1 / (a^2 * cos(x)^2 + sin(x)^2),
                  (a * sin(x) + b * cos(x)) / (c * sin(x) + d * cos(x)))
            # (an older RischIntegrate, the one Julia 1.10 gets, cannot do all of them)
            try integrate(f, x) catch end
        end
        P, H = SymJ.piecewise(x ≪ 0 => -x, true => x^2), SymJ.heaviside(x - 1)
        try integrate(P, x); integrate(x * H, x); integrate(cos(x) * SymJ.dirac(x - 1), x) catch end
        # answers with special functions (Ei, ₁F₁): not every term type has them
        for f in (exp(x) / x, exp(a * x) / x, 1 / log(x)^(2 + a))
            try integrate(f, x) catch end
        end
    end
end

end
