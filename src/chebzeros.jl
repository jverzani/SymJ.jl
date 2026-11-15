# Numeric zeros of a function of one real variable by a Chebyshev-proxy rootfinder (Boyd
# 2002, 2014; the method of Chebfun's roots): `fsolve`'s method (fsolve.jl) for an
# equation in x alone that is not polynomial. Not exported; SymJ.ChebZeros.
#
# f is interpolated at Chebyshev points on [a, b], with n doubling to 128 until the
# Chebyshev coefficients decay (or level off at a noise plateau); a piece that doesn't
# resolve is halved, the piece nearest the middle first. The roots of each piece's
# interpolant are the eigenvalues of its colleague matrix. Real candidates are checked
# against f itself: bracketed and bisected where f changes sign (poles and jumps, where
# |f| doesn't shrink, rejected), else kept only as an even-order zero where |f| ≈ 0. An
# infinite interval is mapped onto [-1, 1] (x = c + s t/(1 - t²), or a + s(1 + t)/(1 - t))
# and f's values squashed by y/√(σ² + y²). The search stops at about 120 zeros (nearest
# the middle), so an f that oscillates forever gives those, with a warning.
#
# Complex zeros: a root λ of a piece's interpolant (ρ its Bernstein-ellipse parameter, n
# the degree, ε the relative size of the discarded coefficients) is trusted where ε ρⁿ is
# small, and then polished by secant steps on f (dropped if they don't converge). On an
# infinite interval they are looked for on a window around the real zeros.

module ChebZeros

using LinearAlgebra: eigvals

export cheb_zeros, cheb_all_zeros

# ---------------------------------------------------------------------------
# Chebyshev machinery
# ---------------------------------------------------------------------------

# Chebyshev coefficients c_0..c_n of the degree-n interpolant through values v[j+1] at
# x_j = cos(πj/n), j = 0..n. Direct O(n²) DCT-I: fine for n ≤ 256 and avoids FFTW.
function chebcoeffs(v::AbstractVector{T}) where {T}
    n = length(v) - 1
    c = zeros(T, n + 1)
    @inbounds for k in 0:n
        s = (v[1] + (iseven(k) ? v[end] : -v[end])) / 2
        for j in 1:n-1
            s += v[j+1] * cospi(T(mod(j * k, 2n)) / n)
        end
        c[k+1] = 2s / n
    end
    c[1] /= 2
    c[end] /= 2
    c
end

# Roots of sum_k c[k+1] T_k(x), with c[end] ≠ 0, as eigenvalues of the colleague matrix.
function colleague_roots(c::AbstractVector{T}) where {T}
    d = length(c) - 1
    d == 0 && return Complex{T}[]
    d == 1 && return [complex(-c[1] / c[2])]
    C = zeros(T, d, d)
    C[1, 2] = one(T)                      # x T_0 = T_1
    for k in 2:d-1                        # x T_k = (T_{k-1} + T_{k+1})/2
        C[k, k-1] = C[k, k+1] = T(1) / 2
    end
    C[d, d-1] = T(1) / 2
    for j in 1:d                          # eliminate T_d using p(x) = 0
        C[d, j] -= c[j] / (2c[end])
    end
    eigvals(C)
end

# ---------------------------------------------------------------------------
# Adaptive construction + recursive subdivision
# ---------------------------------------------------------------------------

struct Candidate{T}
    x::T        # approximate root
    vscale::T   # max |f| on the piece it came from (sets the scale for "|f| ≈ 0")
    w::T        # local length scale near x (b - a, or the stretched equivalent for a map)
end
Candidate(x, vscale) = Candidate(x, vscale, oftype(x, NaN))   # w filled in by the driver

# Resolve f on one piece [a, b]: push candidates; return the degree of the polynomial
# kept (0 if none) and whether the piece is done.
function _piece!(cands, f, a::T, b::T, opts, lastlevel::Bool) where {T}
    mid, rad = (a + b) / 2, (b - a) / 2
    pt(j, n) = mid + rad * cospi(T(j) / n)
    n = 16
    v = T[f(pt(j, n)) for j in 0:n]

    # f undefined (NaN, e.g. log of a negative number) at every sample: the piece is
    # outside f's domain, with no zeros
    all(isnan, v) && return 0, true
    while true
        all(isfinite, v) || break
        vscale = maximum(abs, v)
        if iszero(vscale)                         # f ≡ 0 on the piece
            @warn "f appears to vanish identically on [$a, $b]"
            return 0, true
        end
        c = chebcoeffs(v)
        K = findlast(ck -> abs(ck) > opts.tol * vscale, c)
        tail = (n + 1) - K
        converged = tail ≥ max(3, n ÷ 5)          # coefficients decayed below tol
        if !converged && n ≥ 32
            # Plateau test (a crude version of Chebfun's standardChop): if f is noisy, e.g.
            # because evaluating it is ill-conditioned, the coefficients level off above
            # tol. If they have levelled off at a small level and the second half of the
            # coefficients shows no further decay, the interpolant is as good as f allows.
            env = reverse(accumulate(max, reverse(abs.(c)))) ./ vscale   # tail envelope
            mid_level, end_level = env[(n+1)÷2], env[end-2]
            if mid_level ≤ opts.plateau && end_level ≥ mid_level / 10
                K = findfirst(≤(2mid_level), env)
                converged = true
            end
        end
        if converged
            λs = colleague_roots(c[1:K])
            if opts.cplx !== nothing
                # relative size of what was thrown away (or the noise plateau)
                ε = max(opts.tol, maximum(abs, view(c, K+1:n+1); init=zero(T)) / vscale)
                push!(opts.cplx, (; a, b, λs, ε, d=K - 1, vscale))
            end
            for λ in λs
                if abs(imag(λ)) ≤ opts.imtol && -1 - opts.imtol ≤ real(λ) ≤ 1 + opts.imtol
                    x = clamp(mid + rad * real(λ), a, b)
                    push!(cands, Candidate(x, vscale))
                end
            end
            return K - 1, true
        end
        2n > opts.maxdeg && break
        # double n, reusing samples: new even-index points are the old points
        w = Vector{T}(undef, 2n + 1)
        w[1:2:end] .= v
        for j in 1:2:2n
            w[j+1] = f(pt(j, 2n))
        end
        v, n = w, 2n
    end

    # Not resolved. On the last level, fall back to sign changes of the samples.
    lastlevel || return n, false
    # (at the edge of f's domain, where some samples are NaN, there is nothing to warn of)
    any(isnan, v) || push!(opts.unresolved, a, b)
    xs = [pt(j, n) for j in n:-1:0]
    fs = reverse(v)
    for i in 1:n
        if isfinite(fs[i]) && isfinite(fs[i+1]) && sign(fs[i]) * sign(fs[i+1]) ≤ 0
            push!(cands, Candidate((xs[i] + xs[i+1]) / 2, max(abs(fs[i]), abs(fs[i+1]))))
        end
    end
    0, true
end

# Refinement of [a, b], nearest `focus` first: the unresolved piece nearest `focus` is
# always the next one tried, and is halved (slightly off-center, as Chebfun does, to avoid
# symmetry traps) if it can't be resolved. Two limits stop the search early; the pieces
# left are then skipped (recorded in opts.skipped), so the zeros found are those nearest
# `focus`:
#   maxzeros    – the number of (real) candidate zeros found so far;
#   maxtotaldeg – the total degree of the polynomials kept (the size of the piecewise
#                 proxy), which bounds the work when f oscillates without crossing zero.
# Failed attempts on a piece don't count toward either; their cost is bounded by maxdepth.
function _proxy!(cands, f, a::T, b::T, opts, focus::T=(a + b) / 2) where {T}
    dist((l, u, _)) = l ≤ focus ≤ u ? zero(T) : min(abs(l - focus), abs(u - focus))
    queue = [(a, b, 0)]
    used = 0
    while !isempty(queue)
        k = argmin(dist.(queue))
        l, u, depth = queue[k]
        deleteat!(queue, k)
        if used ≥ opts.maxtotaldeg || length(cands) ≥ opts.maxzeros
            push!(opts.skipped, l, u)
            continue
        end
        tiny = (u - l) ≤ 1000 * eps(max(abs(l), abs(u)))
        deg, done = _piece!(cands, f, l, u, opts, depth ≥ opts.maxdepth || tiny)
        if done
            used += deg
        else
            m = (l + u) / 2 - (u - l) / 2 * T(0.004849834917525)
            push!(queue, (l, m, depth + 1), (m, u, depth + 1))
        end
    end
    cands
end

# ---------------------------------------------------------------------------
# Polishing and validation against f
# ---------------------------------------------------------------------------

# Bisection to machine precision on a sign-change bracket. (In Roots.jl this would
# just be `find_zero(f, (l, u), Bisection())`.)
function _bisect(f, l::T, u::T, fl, fu) where {T}
    while true
        m = l + (u - l) / 2
        (m == l || m == u) && break
        fm = f(m)
        iszero(fm) && return m, m, fm, fm
        if sign(fm) == sign(fl)
            l, fl = m, fm
        else
            u, fu = m, fm
        end
    end
    l, u, fl, fu
end

# Returns (root, bracketed::Bool) or nothing if the candidate is rejected.
function _polish(f, c::Candidate{T}, lo::T, hi::T, a::T, b::T, opts) where {T}
    w = c.w
    x = c.x
    fx = f(x)
    # An exact zero counts only if it is isolated: f ≡ 0 nearby (e.g. underflow in a
    # decaying tail) is not a root.
    δ = max(1e-6 * w, 64eps(abs(x)))
    isolated(p) = !iszero(f(max(p - δ, a))) || !iszero(f(min(p + δ, b)))
    iszero(fx) && return isolated(x) ? (x, true) : nothing
    # grow a bracket [x-h, x+h] ∩ [lo, hi] until f changes sign
    h = max(8eps(abs(x)), 1e-14 * w)
    hmax = max(x - lo, hi - x)
    while true
        l, u = max(x - h, lo), min(x + h, hi)
        fl, fu = f(l), f(u)
        iszero(fl) && return isolated(l) ? (l, true) : nothing
        iszero(fu) && return isolated(u) ? (u, true) : nothing
        if isfinite(fl) && isfinite(fu)
            for (p, q, fp, fq) in ((l, x, fl, fx), (x, u, fx, fu))
                if sign(fp) * sign(fq) < 0
                    l1, u1, fl1, fu1 = _bisect(f, p, q, fp, fq)
                    # Reject poles and jumps: at a true zero |f| at the final bracket is
                    # small relative to f's scale (at worst, evaluation noise), or at least
                    # much smaller than where bisection started. At a pole or a jump it
                    # stays comparable to the scale of f, or grows.
                    fend = min(abs(fl1), abs(fu1))
                    (fend ≤ opts.plateau * c.vscale || fend ≤ min(abs(fp), abs(fq)) / 2) ||
                        return nothing
                    return (abs(fl1) ≤ abs(fu1) ? l1 : u1, true)
                end
            end
        end
        h ≥ hmax && break
        h *= 10
    end
    # No sign change: accept only as an even-multiplicity zero, i.e. |f(x)| is tiny both
    # relative to the piece's scale and relative to |f| at the ends of the search region
    # (the second test matters when |f| varies over many orders of magnitude).
    # (look at least pairtol·w away, past the twin eigenvalue of a double root)
    H = max(x - lo, hi - x, opts.pairtol * w)
    fends = min(abs(f(max(x - H, a))), abs(f(min(x + H, b))))
    (abs(fx) ≤ opts.ztol * c.vscale && abs(fx) ≤ fends / 1000) ? (x, false) : nothing
end

# Squash f's values into (-1, 1) without moving zeros or changing signs:
# y ↦ (y/σ)/√(1+(y/σ)²) is smooth and odd, so simple zeros stay simple. σ is a typical
# size of |f| (a too-small σ would turn a steep zero crossing into a near-step). Used for infinite
# intervals, where f may grow without bound. DomainErrors (e.g. log of a negative
# number) become NaN, which makes the piece be subdivided.
function _safe(f)
    x -> begin
        y = try
            f(x)
        catch e
            e isa DomainError ? oftype(float(x), NaN) : rethrow()
        end
        # (a complex value at a real point: not a real zero there)
        y isa Complex ? (iszero(imag(y)) ? float(real(y)) : oftype(float(x), NaN)) : float(y)
    end
end
function _squash(f, σ)
    fs = _safe(f)
    x -> begin
        y = fs(x) / σ
        isnan(y) ? y : abs(y) > 1e150 ? sign(y) : y / sqrt(1 + y^2)
    end
end

function _defaults(T; maxdeg=128, tol=100eps(T), maxdepth=16, imtol=1e-5, ztol=1e3 * eps(T),
                   pairtol=1e-4, plateau=1e-8, scale=1, center=0, ctol=1e-6, cplx=nothing,
                   maxzeros=120, maxtotaldeg=5000)
    (; maxdeg, tol=T(tol), maxdepth, imtol=T(imtol), ztol=T(ztol), pairtol=T(pairtol),
       plateau=T(plateau), scale=T(scale), center=T(center), ctol=T(ctol), cplx,
       maxzeros, maxtotaldeg, unresolved=T[], skipped=T[])
end

"""
    cheb_zeros(f, a, b; kwargs...)
    cheb_zeros(f; kwargs...)             # whole real line, same as cheb_zeros(f, -Inf, Inf)

Find all zeros of a continuous function `f` on `[a, b]` with a Chebyshev-proxy rootfinder.
Returns a sorted vector.

`a` and/or `b` may be infinite. Then the interval is mapped onto [-1, 1]:

    (-∞, ∞):  x = center + scale·t/(1-t²)
    [a, ∞):   x = a + scale·(1+t)/(1-t)
    (-∞, b]:  x = b - scale·(1-t)/(1+t)

and f's values are squashed by y ↦ (y/σ)/√(1+(y/σ)²), σ a typical |f| from a few samples
(same zeros and signs, but bounded). This
works when f is "tame at infinity": it tends to a limit, or to ±∞ like a power or an
exponential. If f oscillates forever (sin(x)) there are infinitely many zeros; you get
those that could be resolved and a warning. `scale` (default 1) should be roughly the
size of x where the action is. A zero far out is still found (and then polished to full
accuracy in x), but at a cost in function evaluations.

Keywords:
* `maxdeg`   – largest interpolant degree on a piece before it is split (128)
* `tol`      – relative size (vs. max|f| on the piece) below which coefficients count as zero
* `maxdepth` – maximum number of subdivisions
* `maxzeros`  – stop once about this many zeros are found (120); with `maxtotaldeg`, a
               budget on the total degree of the piecewise polynomial proxy (5000) for an f
               that oscillates without crossing zero. When either runs out, the zeros
               nearest the middle (or the origin, or the finite end of a half-line) are
               returned, with a warning
* `imtol`    – how far from the real axis an eigenvalue may be and still be a candidate
* `ztol`     – for zeros without a sign change, accept if |f(x)| ≤ ztol·max|f|
* `plateau`  – coefficients that level off below plateau·max|f| count as noise (noisy f)
* `pairtol`  – unbracketed candidates closer than pairtol·w are one even-order zero
               (w = b - a, or its local equivalent under the map)
* `scale`, `center` – the maps above (infinite intervals only)
"""
cheb_zeros(f; kw...) = cheb_zeros(f, -Inf, Inf; kw...)

function cheb_zeros(f, a::Real, b::Real; kw...)
    a, b = float.(promote(a, b))
    a > b && ((a, b) = (b, a))
    first(_zeros(f, a, b, _defaults(typeof(a); kw...)))
end

# The real zeros, and the map x = X(t) (with derivative dX) from the variable the proxy
# pieces live in (identity for a finite interval) and the function F whose zeros they are
# (f, or f squashed for an infinite interval).
function _zeros(f, a::T, b::T, opts) where {T}

    if isfinite(a) && isfinite(b)
        f = _safe(f)                      # a DomainError (log of a negative) is a NaN
        cands = _proxy!(Candidate{T}[], f, a, b, opts)
        _warn_skipped(opts, (a + b) / 2, identity)
        cands = [Candidate(c.x, c.vscale, b - a) for c in cands]
        return _finish(f, cands, a, b, opts), identity, one, f
    end

    # --- infinite interval: change of variable x = X(t), t ∈ [-1, 1] -----------------
    s, c0 = opts.scale, opts.center
    if isinf(a) && isinf(b)
        X  = t -> c0 + s * t / (1 - t^2)
        dX = t -> s * (1 + t^2) / (1 - t^2)^2
    elseif isinf(b)
        X  = t -> a + s * (1 + t) / (1 - t)
        dX = t -> 2s / (1 - t)^2
    else
        X  = t -> b - s * (1 - t) / (1 + t)
        dX = t -> 2s / (1 + t)^2
    end
    tmax = prevfloat(one(T))                      # the endpoints t = ±1 are x = ±∞
    # σ = median |f| over a coarse sample of the mapped interval
    fs = _safe(f)
    ys = filter(y -> isfinite(y) && !iszero(y),
                [abs(fs(X(clamp(cospi(T(j) / 32), -tmax, tmax)))) for j in 0:32])
    σ = isempty(ys) ? one(T) : sort(ys)[(length(ys) + 1) ÷ 2]
    F = _squash(f, σ)
    G = t -> F(X(clamp(t, -tmax, tmax)))
    focus = isinf(a) && isinf(b) ? zero(T) : isinf(b) ? -one(T) : one(T)   # x = center, a or b
    tc = _proxy!(Candidate{T}[], G, -one(T), one(T), opts, focus)
    _warn_skipped(opts, focus, t -> X(clamp(t, -tmax, tmax)))
    if !isempty(opts.unresolved)
        xs = X.(clamp.(opts.unresolved, -tmax, tmax))
        R = 1e6 * max(s, abs(c0))
        if any(x -> abs(x) > R, xs)
            @warn "f is not resolved for |x| > $R (oscillating forever, or features far " *
                  "beyond `scale`?): zeros out there may be missing or incomplete"
        end
        near = filter(x -> abs(x) ≤ R, xs)
        if !isempty(near)
            @warn "f could not be resolved near x ∈ [$(minimum(near)), $(maximum(near))] " *
                  "(singularity, discontinuity, or edge of f's domain?)"
        end
        empty!(opts.unresolved)                   # already reported
    end
    cands = Candidate{T}[]
    for c in tc
        t = clamp(c.x, -tmax, tmax)
        push!(cands, Candidate(X(t), c.vscale, 2 * dX(t)))   # [-1,1] has width 2
    end
    _finish(F, cands, a, b, opts), X, dX, F
end

# The budget ran out: say where zeros were looked for (X maps a piece end to x).
function _warn_skipped(opts, focus, X)
    isempty(opts.skipped) && return
    ends = opts.skipped
    left = [ends[i+1] for i in 1:2:length(ends) if ends[i+1] ≤ focus]
    right = [ends[i] for i in 1:2:length(ends) if ends[i] ≥ focus]
    lo = isempty(left) ? "" : string(X(maximum(left)))
    hi = isempty(right) ? "" : string(X(minimum(right)))
    @warn "search stopped early (maxzeros = $(opts.maxzeros), maxtotaldeg = " *
          "$(opts.maxtotaldeg); f has very many zeros, or oscillates forever?): zeros " *
          "were searched for only between about " *
          "$(isempty(lo) ? "the left end" : lo) and $(isempty(hi) ? "the right end" : hi)"
    empty!(opts.skipped)
end

# Polish candidates against f, validate, merge duplicates.
function _finish(f, cands::Vector{Candidate{T}}, a::T, b::T, opts) where {T}
    if !isempty(opts.unresolved)
        u = opts.unresolved
        @warn "f could not be resolved on $(length(u)÷2) small piece(s) in " *
              "[$(minimum(u)), $(maximum(u))] (singularity or discontinuity?); " *
              "sign changes of samples were used there"
    end
    filter!(c -> isfinite(c.x), cands)
    sort!(cands; by=c -> c.x)

    # polish each candidate, bracket search limited to halfway to its neighbors
    # (and to a few thousand local length scales on an infinite side)
    found = Tuple{T,Bool,T}[]
    for (i, c) in enumerate(cands)
        lo = i == 1 ? (isfinite(a) ? a : c.x - 1000c.w) : (cands[i-1].x + c.x) / 2
        hi = i == length(cands) ? (isfinite(b) ? b : c.x + 1000c.w) : (c.x + cands[i+1].x) / 2
        lo, hi = min(lo, c.x), max(hi, c.x)
        r = _polish(f, c, lo, hi, a, b, opts)
        r === nothing || push!(found, (r..., c.w))
    end
    sort!(found; by=first)

    # Merge duplicates. Bracketed roots are at machine precision, so they merge only when
    # a few ulps apart. An even-multiplicity zero shows up as a pair of eigenvalues
    # x0 ± O(sqrt(tol)), so two unbracketed candidates that close together are one zero,
    # and their average is a much better estimate than either.
    out = T[]
    lastbr, lastw = false, zero(T)
    for (x, br, w) in found
        if !isempty(out)
            d = x - out[end]
            W = max(w, lastw)
            if br && lastbr
                d ≤ 4eps(max(abs(x), abs(out[end]))) && continue
            elseif !br && !lastbr
                if d ≤ opts.pairtol * W
                    out[end] = (out[end] + x) / 2
                    continue
                end
            elseif d ≤ sqrt(eps(T)) * W
                br && (out[end] = x; lastbr = true)    # keep the bracketed one
                continue
            end
        end
        push!(out, x)
        lastbr, lastw = br, w
    end
    out
end

# ---------------------------------------------------------------------------
# Complex zeros inside the trusted Bernstein ellipses
# ---------------------------------------------------------------------------

# Bernstein-ellipse parameter of z relative to [-1, 1]: the ellipse with foci ±1 through z
# has semi-axes summing to ρ(z) ≥ 1 (ρ = 1 is the interval itself).
function bernstein_rho(z)
    w = z + sqrt(z - 1) * sqrt(z + 1)
    max(abs(w), inv(abs(w)))
end

# Derivative-free complex polishing (secant). Returns (z, converged).
function _csecant(f, z0::Complex{T}, scale::T; maxit=60) where {T}
    z1 = z0 + max(1e-8 * scale, 1e-8 * abs(z0))
    f0, f1 = f(z0), f(z1)
    for _ in 1:maxit
        iszero(f1) && return z1, true
        Δf = f1 - f0
        iszero(Δf) && return z1, false
        z2 = z1 - f1 * (z1 - z0) / Δf
        isfinite(z2) || return z1, false
        abs(z2 - z1) ≤ 8eps(T) * max(abs(z2), scale) && return z2, true
        z0, f0, z1, f1 = z1, f1, z2, f(z2)
    end
    z1, false
end

"""
    cheb_all_zeros(f, a, b; ctol=1e-6, polish=true, kwargs...)
    cheb_all_zeros(f; kwargs...)        # whole real line

Real zeros of `f` on [a, b] (as `cheb_zeros`; a and b may be infinite), plus the complex
zeros of the Chebyshev proxy that can be trusted. Returns a named tuple `(real, complex)`.

`complex` holds the non-real zeros, each a named tuple `(z, rho, errest, polished)`:

* `rho` – Bernstein-ellipse parameter of the root for the piece it came from (1 = on the
  piece; larger = farther out)
* `errest` – estimated |p - f| / max|f| there for the proxy p: ε·ρⁿ, with ε the relative
  size of the discarded Chebyshev coefficients and n the degree. A root is kept only if
  `errest ≤ ctol`. The root's own error is roughly errest·max|f| / |f′(z)|.
* `polished` – true if complex secant steps on f converged. If f accepts complex input,
  roots that don't polish are dropped; if it doesn't, the proxy's roots are returned
  unpolished.

Each piece has its own ellipse, so the trusted region is a union of ellipses hugging the
real axis. Near-real conjugate pairs that are really a double real zero are dropped
(that zero is in `real`).

On an infinite interval the real zeros come from the mapped search, and complex zeros are
looked for (on the finite interval) in a window around them: |x - c| ≤ R, with
R = max(10·scale, 2·max|r - c|) over the real zeros r, and c the center (whole line) or
the finite end (half-line).
"""
cheb_all_zeros(f; kw...) = cheb_all_zeros(f, -Inf, Inf; kw...)

function cheb_all_zeros(f, a::Real, b::Real; ctol=1e-6, polish::Bool=true, kw...)
    a, b = float.(promote(a, b))
    a > b && ((a, b) = (b, a))
    T = typeof(a)
    pieces = NamedTuple[]
    opts = _defaults(T; ctol, cplx=pieces, kw...)
    realzeros, _ = _zeros(f, a, b, opts)
    if !(isfinite(a) && isfinite(b))
        # Infinite interval: complex zeros are looked for, by the finite-interval method,
        # in a window around the real zeros: |x - c| ≤ R, R = max(10·scale, 2·max|r - c|),
        # with c = center (whole line) or the finite end (half-line). (In the mapped
        # variable the pieces are too narrow for their ellipses to reach off the axis.)
        c = isinf(a) && isinf(b) ? opts.center : isinf(b) ? a : b
        R = max(10 * opts.scale, 2 * maximum(r -> abs(r - c), realzeros; init=zero(T)))
        lo, hi = isinf(a) ? c - R : a, isinf(b) ? c + R : b
        res = cheb_all_zeros(f, lo, hi; ctol, polish, kw...)
        return (; real=realzeros, complex=res.complex)
    end

    # can f be evaluated at complex points?
    canpolish = polish && try
        isfinite(f(complex((a + b) / 2, (b - a) / 7)))
    catch
        false
    end

    out = NamedTuple{(:z, :rho, :errest, :polished, :w),Tuple{Complex{T},T,T,Bool,T}}[]
    for pc in pieces
        mid, rad = (pc.a + pc.b) / 2, (pc.b - pc.a) / 2
        for λ in pc.λs
            abs(imag(λ)) > 0 || continue                  # real ones handled above
            ρ = bernstein_rho(λ)
            errest = pc.ε * ρ^pc.d
            errest ≤ opts.ctol || continue                # outside the trusted ellipse
            z, w = complex(mid + rad * real(λ), rad * imag(λ)), rad
            polished = false
            if canpolish
                zp, ok = _csecant(f, z, w)
                # accept the polished root if it converged and didn't wander off; if it
                # didn't, it is not a zero of f (where |f| is far below max|f| the
                # proxy's roots are noise: x e^(x²) - 1 on [-10, 10] near |x| = 9)
                (ok && abs(zp - z) ≤ 1e-2 * w) || continue
                z, polished = zp, true
            end
            push!(out, (; z, rho=ρ, errest, polished, w))
        end
    end

    # drop "complex" roots that are really (double) real zeros, then duplicates from
    # overlapping pieces, keeping the one with the smaller error estimate
    near_real(r) = abs(imag(r.z)) ≤ 1e3eps(T) * max(r.w, abs(r.z)) ||
                   any(x -> abs(r.z - x) ≤ 10sqrt(eps(T)) * r.w, realzeros)
    filter!(r -> !near_real(r), out)
    sort!(out; by=r -> r.errest)
    kept = similar(out, 0)
    for r in out
        any(k -> abs(k.z - r.z) ≤ 1e-8 * max(r.w, abs(r.z)), kept) || push!(kept, r)
    end
    sort!(kept; by=r -> (real(r.z), imag(r.z)))
    (; real=realzeros, complex=[(; z=r.z, rho=r.rho, errest=r.errest, polished=r.polished) for r in kept])
end

end # module
