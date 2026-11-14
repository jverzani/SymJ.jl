# Real and imaginary parts, conjugates and arguments: real, imag, conj, angle (and reim).
#
# Symbols are complex, and construction never consults assumptions, so only identities
# that hold for every complex value apply as these are built:
#   - constants are evaluated exactly: real(2 + 3im) = 2, conj(π) = π, angle(-2) = π;
#   - real and imag are linear over sums and real coefficients, and real(im*u) =
#     -imag(u);
#   - conj distributes over sums, products, integer powers and functions without a
#     branch cut (exp, sin, cos, sinh, cosh, tan, ...), but not over log or fractional
#     powers: conj(log(z)) is not log(conj(z)) on the negative real axis;
#   - real(x), imag(x), conj(x) of a symbol stay as they are.
# `refine` takes expressions apart into real and imaginary parts under assumptions
# (real(x) = x for real x, real(exp(im*x)) = cos(x)); see simplify.jl.

const _REALFNS = (:real, :imag, :abs, :angle)              # real-valued functions
const _CONJ_COMMUTES = (:sin, :cos, :tan, :cot, :sec, :csc, :sinh, :cosh, :tanh, :coth,
                        :sech, :csch, :abs, :real, :imag, :angle)

# the real and imaginary parts of a constant, exactly where its structure gives them
# (a + b im, products and exponentials of those); nothing otherwise
function _const_parts(e::Sym)
    _hassymbol(e) && return nothing
    k = e.kind
    (k == NUMBER || k == CONSTANT) && return (e, ZERO)
    _eq(e, IM) && return (ZERO, ONE)
    if k == ADD
        re, im = Sym(e.num), ZERO
        for (c, t) in zip(e.coefs, e.args)
            p = _const_parts(t)
            p === nothing && return nothing
            re += Sym(c) * p[1]
            im += Sym(c) * p[2]
        end
        return (re, im)
    elseif k == MUL
        re, im = Sym(e.num), ZERO
        for (b, q) in zip(e.args, e.exps)
            p = _const_parts(b^q)
            p === nothing && return nothing
            re, im = expand(re * p[1] - im * p[2]), expand(re * p[2] + im * p[1])
        end
        return (re, im)
    elseif k == POW
        b, q = e.args
        if _eq(b, E)                                        # exp(a + b im)
            p = _const_parts(q)
            p === nothing && return nothing
            return (exp(p[1]) * cos(p[2]), exp(p[1]) * sin(p[2]))
        end
        # b^q, q not an integer, b negative or complex or q complex: exp(q log b) with log b
        # in parts (i^i = e^(-π/2), √i = e^(iπ/4); SymPy #16348)
        if !(q.kind == NUMBER && isexact(q.num) && isone(_denominator(q.num))) &&
           (has(b, IM) || (b.kind == NUMBER && b.num < N0) || has(q, IM))
            lb = log(b)
            pl = _hasfn(lb, :log) && !_realconst(lb) ? nothing : _const_parts(lb)
            if pl !== nothing
                w = exp(expand(q * lb))
                _eq(w, e) || return _const_parts(w)
            end
        end
        if q.kind == NUMBER && isexact(q.num) && isone(_denominator(q.num)) && !iszero(q.num)
            p = _const_parts(b)                             # integer powers
            p === nothing && return nothing
            n = abs(_numerator(q.num))
            # (by squaring; past 8 only parts that are numbers, (1 + i)^100 is -2^50, as
            # surds' powers grow; SymPy #9890)
            n <= 8 || (p[1].kind == NUMBER && p[2].kind == NUMBER && n <= 4096) || return nothing
            mul((a, b), (c, d)) = (expand(a * c - b * d), expand(a * d + b * c))
            r, sq, m = (ONE, ZERO), p, Int(n)
            while m > 0
                isodd(m) && (r = mul(r, sq))
                m >>= 1
                m > 0 && (sq = mul(sq, sq))
            end
            q.num > N0 && return r
            # a negative power through the conjugate: 1/(a + b i) = (a - b i)/(a² + b²)
            # (real(1/(4 + 3i)) is 4/25; SymPy #10380)
            d = expand(r[1]^2 + r[2]^2)
            iszero(d) && return nothing
            return (expand(r[1] / d), expand(-r[2] / d))
        end
    end
    _isreal_const(e) === true ? (e, ZERO) : nothing
end

_realconst(e::Sym) = !_hassymbol(e) && (p = _const_parts(e); p !== nothing && iszero(expand(p[2])))
_isfn1(e::Sym, f::Symbol) = e.kind == FN && e.name === f && length(e.args) == 1

# real(u) or imag(u) (part = 1, 2): linear over sums and real coefficients
function _part(u::Sym, part::Int)
    p = _const_parts(u)
    p === nothing || return p[part]
    f = part == 1 ? :real : :imag
    if u.kind == ADD
        return (part == 1 ? Sym(u.num) : ZERO) + _add_many(Sym(c) * apply_function(f, Sym[t]) for (c, t) in zip(u.coefs, u.args))
    elseif u.kind == MUL
        rest = u / Sym(u.num)
        isone(u.num) || return Sym(u.num) * apply_function(f, Sym[rest])
        # im * w: real(im w) = -imag(w), imag(im w) = real(w)
        i = findfirst(j -> _eq(u.args[j]^u.exps[j], IM), eachindex(u.args))
        if i !== nothing
            w = expand(u / IM)
            return part == 1 ? -apply_function(:imag, Sym[w]) : apply_function(:real, Sym[w])
        end
        # a real constant factor comes out
        for (b, q) in zip(u.args, u.exps)
            t = b^q
            _realconst(t) && return t * apply_function(f, Sym[expand(u / t)])
        end
    elseif u.kind == FN && length(u.args) == 1 && u.name in _REALFNS
        return part == 1 ? u : ZERO                         # real(abs(z)) = abs(z)
    elseif u.kind == FN && u.name === :conj && length(u.args) == 1
        return part == 1 ? apply_function(:real, u.args) : -apply_function(:imag, u.args)
    elseif u.kind == FN && u.name === :log && length(u.args) == 1 && !_hassymbol(u.args[1])
        # log z of a complex constant: log|z| and arg z (principal; SymPy #17948)
        z = u.args[1]
        zp = _const_parts(z)
        if zp !== nothing && !iszero(expand(zp[2]))
            return part == 1 ? log(abs(z)) : angle(z)
        end
    end
    # real for every value of its symbols (abs(z)^2, conj(z) z): no assumptions
    _isreal(u, Assumptions()) === true && return part == 1 ? u : ZERO
    nothing
end

function _conj(u::Sym)
    p = _const_parts(u)
    p === nothing || return expand(p[1] - IM * p[2])
    k = u.kind
    k == ADD && return Sym(u.num) + _add_many(Sym(c) * conj(t) for (c, t) in zip(u.coefs, u.args))
    if k == MUL
        out = Sym(u.num)
        rest = Sym[]
        for (b, q) in zip(u.args, u.exps)
            t = b^q
            if q.kind == NUMBER && isexact(q.num) && isone(_denominator(q.num)) || _eq(b, E)
                out *= _eq(b, E) ? exp(conj(q)) : conj(b)^q
            else
                push!(rest, t)
            end
        end
        isempty(rest) && return out
        length(rest) == length(u.args) && return nothing
        return out * apply_function(:conj, Sym[prod(rest)])
    end
    if k == POW
        b, q = u.args
        _eq(b, E) && return exp(conj(q))
        q.kind == NUMBER && isexact(q.num) && isone(_denominator(q.num)) && return conj(b)^q
        return nothing
    end
    if k == FN && length(u.args) == 1
        u.name === :conj && return u.args[1]
        u.name in _REALFNS && return u
        u.name in _CONJ_COMMUTES && return apply_function(u.name, Sym[conj(u.args[1])])
    end
    nothing
end

function _angle(u::Sym)
    p = _const_parts(u)
    p === nothing && return nothing
    re, im = expand(p[1]), expand(p[2])
    sr, si = _csign(re), _csign(im)
    (sr === nothing || si === nothing) && return nothing
    si == 0 && return sr >= 0 ? ZERO : PI
    sr == 0 && return si > 0 ? PI / 2 : -PI / 2
    a = atan(im / re)
    sr > 0 ? a : si > 0 ? a + PI : a - PI
end

register_function!(:real; julia = real, eval = args -> length(args) == 1 ? _part(args[1], 1) : nothing)
register_function!(:imag; julia = imag, eval = args -> length(args) == 1 ? _part(args[1], 2) : nothing)
register_function!(:conj; julia = conj, eval = args -> length(args) == 1 ? _conj(args[1]) : nothing)
register_function!(:angle; julia = angle, eval = args -> length(args) == 1 ? _angle(args[1]) : nothing)

# abs of a constant from its parts: abs(3 + 4im) = 5 (the rule abs had, otherwise)
const _ABS_EVAL = FUNCTION_RULES[:abs].eval
register_function!(:abs; julia = abs, diff = FUNCTION_RULES[:abs].diff,
    eval = function (args)
        a = args[1]
        # |b^q| = |b|^re(q) e^(-im(q) arg b), b off the positive axis and q nonreal, constants
        # (SymPy #14274)
        if length(args) == 1 && a.kind == POW && !_hassymbol(a) && !_eq(a.args[1], E)
            b, q = a.args
            pb, pq, ab = _const_parts(b), _const_parts(q), _angle(b)
            if pb !== nothing && pq !== nothing && ab !== nothing && !iszero(ab) && !iszero(expand(pq[2]))
                return abs(b)^expand(pq[1]) * exp(-expand(pq[2]) * ab)
            end
        end
        # |b^n| = |b|^n, n an integer (abs((1 + i)^100) is 2^50; SymPy #9890)
        if length(args) == 1 && a.kind == POW && !_hassymbol(a) && (q = a.args[2]; q.kind == NUMBER &&
                isexact(q.num) && isone(_denominator(q.num)))
            ab = abs(a.args[1])
            _hasfn(ab, :abs) || return ab^q
        end
        if length(args) == 1 && a.kind != NUMBER && !_hassymbol(a)
            p = _const_parts(a)
            p !== nothing && !iszero(expand(p[2])) && return sqrt(expand(p[1]^2 + p[2]^2))
        end
        _ABS_EVAL(args)
    end)

# conj(w)^k w^k as abs(w)^2k throughout e (a candidate of simplify)
function _conj_pairs_to_abs(e::Sym)
    isempty(e.args) && return e
    if e.kind == MUL
        ts = [_conj_pairs_to_abs(b^p) for (b, p) in zip(e.args, e.exps)]
        rest = _unpair_conj(ts)
        length(rest) == length(ts) && return Sym(e.num) * prod(ts)
        pairs = ONE
        for t in ts
            any(r -> r === t, rest) && continue
            b, k = t.kind == POW ? (t.args[1], t.args[2]) : (t, ONE)
            b.kind == FN && b.name === :conj && (pairs *= abs(b.args[1])^(2k))
        end
        return Sym(e.num) * pairs * prod(rest; init = ONE)
    end
    _rebuild_node(e, Sym[_conj_pairs_to_abs(a) for a in e.args])
end

Base.real(x::Sym) = apply_function(:real, Sym[x])
Base.imag(x::Sym) = apply_function(:imag, Sym[x])
Base.conj(x::Sym) = apply_function(:conj, Sym[x])::Sym   # (inferred: adjoint of a Matrix{Sym} is one of Sym)
Base.angle(x::Sym) = apply_function(:angle, Sym[x])
Base.reim(x::Sym) = (real(x), imag(x))
# for arrays of expressions: A' and transpose(A)
Base.adjoint(x::Sym) = conj(x)::Sym
Base.transpose(x::Sym) = x

# the symbols vs taken to be real: conj(v), real(v), imag(v) are v, v, 0 (an integral's
# variable; dot(u, v) of symbol vectors conjugates its first argument)
function _realsyms(ex::Sym, vs)
    _has_fn(ex, (:conj, :real, :imag)) || return ex
    ps = Pair{Sym,Sym}[]
    for v in vs
        push!(ps, conj(v) => v, real(v) => v, imag(v) => ZERO)
    end
    subs(ex, ps...)
end
_realsyms(ex, vs) = ex
