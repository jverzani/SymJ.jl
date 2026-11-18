# Numeric quadrature for the values of unevaluated definite integrals (N of
# Integral(f, x, a, b), the value of a definite integral proven to have no elementary
# antiderivative). Adaptive Gauss–Kronrod (7–15 points), from a dependency-free version by
# the package's author (its constants checked against QUADPACK's qk15). The nodes are inside
# each interval, so an integrable singularity at an end is never evaluated; infinite ends
# by a change of variable onto a finite interval.

const _GK_X = (0.0, 0.207784955007898467600689403773245, 0.405845151377397166906606412076961,
               0.586087235467691130294144845693013, 0.741531185599394439863864773280788,
               0.864864423359769072789712788640926, 0.949107912342758524526189684047851,
               0.991455371120812639206854697526329)
const _GK_WK = (0.209482141084727828012999174891714, 0.204432940075298892414161999234649,
                0.190350578064785409913256402421014, 0.169004726639267902826583426598550,
                0.140653259715525918745189590510238, 0.104790010322250183839876322541518,
                0.063092092629978553290700663189204, 0.022935322010529224963732008058970)
const _GK_WG = (0.417959183673469387755102040816327, 0.0, 0.381830050505118944950369775488975, 0.0,
                0.279705391489276667901467771423780, 0.0, 0.129484966168869693270611432679082, 0.0)

# the Kronrod value of ∫ₐᵇ f and its difference from the Gauss value (an error estimate)
function _gk15(f, a::Float64, b::Float64)
    m, Δ = (a + b) / 2, (b - a) / 2
    fc = f(m)
    ik, ig = fc * _GK_WK[1], fc * _GK_WG[1]
    for i in 2:8
        dx = Δ * _GK_X[i]
        s = f(m - dx) + f(m + dx)
        ik += s * _GK_WK[i]
        ig += s * _GK_WG[i]
    end
    ik * Δ, abs(ik - ig) * abs(Δ)
end

"""
    SymJ._quadgk(f, a, b; rtol = 1e-10, atol = 0.0, maxevals = 200_000)

∫ₐᵇ f(t) dt for a function `f` of a `Float64` (real or complex values), `a` and `b` real,
possibly infinite: `(value, error estimate)`. The interval with the largest error estimate
is halved until the sum of the estimates is within `max(atol, rtol |value|)`, or
`maxevals` evaluations are spent. Infinite ends by t = a + s/(1 - s), t = s/(1 - s²).
"""
function _quadgk(f, a::Real, b::Real; rtol::Real = 1e-10, atol::Real = 0.0, maxevals::Integer = 200_000)
    a, b = Float64(a), Float64(b)
    a == b && return (0.0, 0.0)
    a > b && return ((v, e) = _quadgk(f, b, a; rtol, atol, maxevals); (-v, e))
    if isinf(a) && isinf(b)
        g = s -> (d = 1 - s^2; f(s / d) * (1 + s^2) / d^2)
        return _quadgk_finite(g, -1.0, 1.0, rtol, atol, maxevals)
    elseif isinf(b)
        g = s -> (d = 1 - s; f(a + s / d) / d^2)
        return _quadgk_finite(g, 0.0, 1.0, rtol, atol, maxevals)
    elseif isinf(a)
        g = s -> (d = 1 - s; f(b - s / d) / d^2)
        return _quadgk_finite(g, 0.0, 1.0, rtol, atol, maxevals)
    end
    _quadgk_finite(f, a, b, rtol, atol, maxevals)
end

function _quadgk_finite(f, a::Float64, b::Float64, rtol, atol, maxevals)
    I0, E0 = _gk15(f, a, b)
    segs = [(a, b, I0, E0)]
    total, err = I0, E0
    evals = 15
    while err > max(atol, rtol * abs(total)) && evals < maxevals
        _, i = findmax(s -> s[4], segs)
        l, r, Ii, Ei = segs[i]
        m = (l + r) / 2
        (m == l || m == r) && break                   # no more room in Float64
        I1, E1 = _gk15(f, l, m)
        I2, E2 = _gk15(f, m, r)
        evals += 30
        segs[i] = (l, m, I1, E1)
        push!(segs, (m, r, I2, E2))
        total = sum(s -> s[3], segs)
        err = sum(s -> s[4], segs)
    end
    isfinite(total) || throw(ArgumentError("the integrand is not finite on the interval"))
    total, err
end

# Integral(f, v, a, b) by quadrature; `a` (_Named or a value) and `var` as in _ev, so
# the values of an enclosing integral's variables reach the integrand of an inner one
function _ev_integral(ex::Sym, a, var)
    f, v, lo, hi = ex.args
    if a isa _Named
        named = a.vals
    else
        fs = free_symbols(ex)
        if isempty(fs)
            named = NamedTuple()
        elseif length(fs) == 1
            s = only(fs)
            if var[] === nothing
                var[] = s
            elseif !_eq(var[]::Sym, s)
                throw(ArgumentError("expression has more than one free symbol; call it with pairs, e.g. ex(x => 1, y => 2)"))
            end
            named = NamedTuple{(s.name,)}((a,))
        else
            throw(ArgumentError("expression has more than one free symbol; call it with pairs, e.g. ex(x => 1, y => 2)"))
        end
    end
    la = _ev(lo, _Named(named), Float64, Ref{Union{Nothing,Sym}}(nothing))
    lb = _ev(hi, _Named(named), Float64, Ref{Union{Nothing,Sym}}(nothing))
    g(t) = _ev(f, _Named(merge(named, NamedTuple{(v.name,)}((t,)))), Float64, Ref{Union{Nothing,Sym}}(nothing))
    val, err = try
        if la isa Real && lb isa Real
            _quadgk(g, la, lb)
        else                                         # a segment of the complex plane
            za, zb = complex(la), complex(lb)
            (isfinite(za) && isfinite(zb)) || throw(ArgumentError("complex infinite limits"))
            _quadgk(s -> g(za + s * (zb - za)) * (zb - za), 0.0, 1.0)
        end
    catch err
        err isa ArgumentError || rethrow()
        throw(ArgumentError("no numeric value for $ex: $(err.msg)"))
    end
    err <= 1e-6 * max(1.0, abs(val)) ||
        throw(ArgumentError("no numeric value for $ex: the quadrature did not converge (estimated error $err)"))
    val
end
