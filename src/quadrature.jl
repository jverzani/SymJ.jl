# Numeric quadrature for the values of unevaluated definite integrals (N of
# Integral(f, x, a, b), the value of a definite integral proven to have no elementary
# antiderivative), by QuadGK's adaptive Gauss–Kronrod: complex values, infinite limits by a
# change of variable, nodes inside each interval (an integrable singularity at an end is
# never evaluated), and any float type: in BigFloat with a higher order, so evalf takes an
# integral to the precision asked.

import QuadGK

"""
    SymJ._quadgk(f, a, b; rtol = nothing, maxevals = 10^6)

∫ₐᵇ f(t) dt by QuadGK, `a` and `b` real (possibly infinite) of one float type `T`:
`(value, error estimate)`. The relative tolerance is `1e-12` in `Float64`, and about
`2^-(precision - 16)` in `BigFloat`, where the order of the rule rises with the precision.
"""
function _quadgk(f, a::Real, b::Real; rtol = nothing, maxevals::Integer = 10^6)
    T = promote_type(float(typeof(a)), float(typeof(b)))
    a, b = T(a), T(b)
    if T === BigFloat
        p = precision(BigFloat)
        tol = rtol === nothing ? BigFloat(2)^(16 - p) : rtol
        return QuadGK.quadgk(f, a, b; rtol = tol, maxevals, order = max(7, p ÷ 12))
    end
    QuadGK.quadgk(f, a, b; rtol = something(rtol, 1e-12), maxevals)
end

# Integral(f, v, a, b) by quadrature, in F's precision (Float64 when F is not a BigFloat
# type); `a` (_Named or a value) and `var` as in _ev, so the values of an enclosing
# integral's variables reach the integrand of an inner one
function _ev_integral(ex::Sym, a, var, F = Float64)
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
    R = F !== nothing && real(F) === BigFloat ? BigFloat : Float64
    nores() = Ref{Union{Nothing,Sym}}(nothing)
    la = _ev(lo, _Named(named), R, nores())
    lb = _ev(hi, _Named(named), R, nores())
    g(t) = _ev(f, _Named(merge(named, NamedTuple{(v.name,)}((t,)))), R, nores())
    val, err = try
        if la isa Real && lb isa Real
            _quadgk(g, la, lb)
        else                                         # a segment of the complex plane
            za, zb = complex(la), complex(lb)
            (isfinite(za) && isfinite(zb)) || throw(ArgumentError("complex infinite limits"))
            _quadgk(s -> g(za + s * (zb - za)) * (zb - za), zero(R), one(R))
        end
    catch err
        err isa Union{ArgumentError,DomainError} || rethrow()
        throw(ArgumentError("no numeric value for $ex: $(sprint(showerror, err))"))
    end
    tol = R === BigFloat ? BigFloat(2)^(24 - precision(BigFloat)) : 1e-8
    (isfinite(val) && err <= tol * max(1, abs(val))) ||
        throw(ArgumentError("no numeric value for $ex: the quadrature did not converge (estimated error $(Float64(err)))"))
    val
end
