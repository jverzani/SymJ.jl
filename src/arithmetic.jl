# Canonicalizing constructors for sums, products and powers.
#
# Invariants (what "canonical" means):
#   ADD: no NUMBER or ADD terms; terms are coefficient-free; coefficients are nonzero;
#        at least two terms, or one term and a nonzero constant.
#   MUL: no NUMBER (except as a base with non-integer exponent) or MUL bases; bases are
#        distinct; exponents are nonzero; not (coefficient 1 and a single factor);
#        never a numeric coefficient times a single sum (that is distributed).
#   POW: exponent is not 0 or 1.
# Terms of a sum and factors of a product are sorted by `_cmp`, so two operands can be
# combined with a linear merge.

## small vector helpers (never mutate inputs)

function _insert(v::Vector{T}, i::Int, x::T) where {T}
    r = Vector{T}(undef, length(v) + 1)
    copyto!(r, 1, v, 1, i - 1)
    r[i] = x
    copyto!(r, i + 1, v, i, length(v) - i + 1)
    r
end

function _deleteat(v::Vector{T}, i::Int) where {T}
    r = Vector{T}(undef, length(v) - 1)
    copyto!(r, 1, v, 1, i - 1)
    copyto!(r, i, v, i + 1, length(v) - i)
    r
end

# binary search in a sorted vector: (index, found)
function _search(v::Vector{Sym}, t::Sym)
    lo, hi = 1, length(v)
    @inbounds while lo <= hi
        mid = (lo + hi) >>> 1
        c = _cmp(v[mid], t)
        c == 0 && return (mid, true)
        c < 0 ? (lo = mid + 1) : (hi = mid - 1)
    end
    (lo, false)
end

## decompositions

# x (not NUMBER, not ADD) as coefficient * key, key coefficient-free
@inline function _split_coef(x::Sym)
    (x.kind != MUL || _isN1(x.num)) && return (N1, x)
    (x.num, _mkmul(N1, x.args, x.exps))
end

# x (not MUL) as base^exponent
@inline _split_pow(x::Sym) = x.kind == POW ? (@inbounds(x.args[1]), @inbounds(x.args[2])) : (x, ONE)

## assembly of canonical nodes from canonical parts

_mkpow(b::Sym, e::Sym) = _isone(e) ? b : _mkpow_node(b, e)

function _mkmul(c::SymNumber, bases::Vector{Sym}, exps::Vector{Sym})
    iszero(c) && return ZERO
    n = length(bases)
    n == 0 && return Sym(c)
    if n == 1
        _isN1(c) && return _mkpow(bases[1], exps[1])
        bases[1].kind == ADD && _isone(exps[1]) && return _scale_add(bases[1], c)
    end
    _mkmul_node(c, bases, exps)
end

# coefficient * key for a key from a sum
function _with_coef(c::SymNumber, key::Sym)
    _isN1(c) && return key
    key.kind == MUL && return _mkmul_node(c, key.args, key.exps)
    b, e = _split_pow(key)
    _mkmul_node(c, Sym[b], Sym[e])
end

function _mkadd(c::SymNumber, terms::Vector{Sym}, coefs::Vector{SymNumber})
    n = length(terms)
    n == 0 && return Sym(c)
    n == 1 && iszero(c) && return _with_coef(coefs[1], terms[1])
    _mkadd_node(c, terms, coefs)
end

## scaling by a number

_scale_add(A::Sym, c::SymNumber) = _mkadd(A.num * c, A.args, SymNumber[x * c for x in A.coefs])

function _scale(x::Sym, c::SymNumber)
    iszero(c) && return ZERO
    _isN1(c) && return x
    k = x.kind
    k == NUMBER && return Sym(x.num * c)
    k == ADD && return _scale_add(x, c)
    k == MUL && return _mkmul(x.num * c, x.args, x.exps)
    b, e = _split_pow(x)
    _mkmul_node(c, Sym[b], Sym[e])
end

## addition

"""
    add(a::Sym, b::Sym)

Canonical sum of two expressions (what `a + b` calls).
"""
function add(a::Sym, b::Sym)
    ka, kb = a.kind, b.kind
    if ka == NUMBER
        kb == NUMBER && return Sym(a.num + b.num)
        return _add_num(b, a.num)
    end
    kb == NUMBER && return _add_num(a, b.num)
    if ka == ADD
        kb == ADD && return _merge_add(a, b)
        return _add_term(a, b)
    end
    kb == ADD && return _add_term(b, a)
    ca, ta = _split_coef(a)
    cb, tb = _split_coef(b)
    c = _cmp(ta, tb)
    if c == 0
        s = ca + cb
        return iszero(s) ? ZERO : _with_coef(s, ta)
    end
    c < 0 ? _mkadd_node(N0, Sym[ta, tb], SymNumber[ca, cb]) :
            _mkadd_node(N0, Sym[tb, ta], SymNumber[cb, ca])
end

function _add_num(x::Sym, n::SymNumber)
    iszero(n) && return x
    x.kind == ADD && return _mkadd(x.num + n, x.args, x.coefs)
    c, t = _split_coef(x)
    _mkadd_node(n, Sym[t], SymNumber[c])
end

function _add_term(A::Sym, x::Sym)
    c, t = _split_coef(x)
    i, found = _search(A.args, t)
    if found
        s = A.coefs[i] + c
        iszero(s) && return _mkadd(A.num, _deleteat(A.args, i), _deleteat(A.coefs, i))
        coefs = copy(A.coefs)
        coefs[i] = s
        return _mkadd_node(A.num, A.args, coefs)
    end
    _mkadd(A.num, _insert(A.args, i, t), _insert(A.coefs, i, c))
end

function _merge_add(A::Sym, B::Sym)
    a, b, ac, bc = A.args, B.args, A.coefs, B.coefs
    na, nb = length(a), length(b)
    terms = Vector{Sym}(undef, na + nb)
    coefs = Vector{SymNumber}(undef, na + nb)
    i, j, k = 1, 1, 0
    @inbounds while i <= na && j <= nb
        c = _cmp(a[i], b[j])
        if c < 0
            k += 1; terms[k] = a[i]; coefs[k] = ac[i]; i += 1
        elseif c > 0
            k += 1; terms[k] = b[j]; coefs[k] = bc[j]; j += 1
        else
            s = ac[i] + bc[j]
            if !iszero(s)
                k += 1; terms[k] = a[i]; coefs[k] = s
            end
            i += 1; j += 1
        end
    end
    @inbounds while i <= na
        k += 1; terms[k] = a[i]; coefs[k] = ac[i]; i += 1
    end
    @inbounds while j <= nb
        k += 1; terms[k] = b[j]; coefs[k] = bc[j]; j += 1
    end
    resize!(terms, k)
    resize!(coefs, k)
    _mkadd(A.num + B.num, terms, coefs)
end

# n-ary sum through a hash table; used where many terms are combined at once
function _accum!(d::Dict{Sym,SymNumber}, t::Sym, c::SymNumber)
    d[t] = get(d, t, N0) + c
    d
end

function _accum_sym!(d::Dict{Sym,SymNumber}, c::SymNumber, x::Sym, s::SymNumber)
    k = x.kind
    if k == NUMBER
        c += x.num * s
    elseif k == ADD
        c += x.num * s
        for i in eachindex(x.args)
            _accum!(d, x.args[i], x.coefs[i] * s)
        end
    else
        cc, t = _split_coef(x)
        _accum!(d, t, cc * s)
    end
    c
end

function _add_from_dict(c::SymNumber, d::Dict{Sym,SymNumber})
    terms = Sym[]
    sizehint!(terms, length(d))
    for (t, v) in d
        iszero(v) || push!(terms, t)
    end
    permute!(terms, _canonical_perm(terms))
    _mkadd(c, terms, SymNumber[d[t] for t in terms])
end

function _add_many(xs)
    c = N0
    d = Dict{Sym,SymNumber}()
    for x in xs
        c = _accum_sym!(d, c, Sym(x), N1)
    end
    _add_from_dict(c, d)
end

## multiplication

"""
    mul(a::Sym, b::Sym)

Canonical product of two expressions (what `a * b` calls).
"""
function mul(a::Sym, b::Sym)
    ka, kb = a.kind, b.kind
    if ka == NUMBER
        kb == NUMBER && return Sym(a.num * b.num)
        return _scale(b, a.num)
    end
    kb == NUMBER && return _scale(a, b.num)
    if ka == MUL
        kb == MUL && return _merge_mul(a, b)
        return _mul_factor(a, b)
    end
    kb == MUL && return _mul_factor(b, a)
    ba, ea = _split_pow(a)
    bb, eb = _split_pow(b)
    c = _cmp(ba, bb)
    c == 0 && return pow(ba, add(ea, eb))
    c < 0 ? _mkmul(N1, Sym[ba, bb], Sym[ea, eb]) : _mkmul(N1, Sym[bb, ba], Sym[eb, ea])
end

# Record the canonical value `p` of `b^e` (just computed by `pow`) into a product under
# construction. Usually `p` is `b^e'` and becomes a factor; a numeric `p` folds into the
# coefficient; anything else (e.g. 2^(3/2) -> 2*sqrt(2)) is returned as a leftover to be
# multiplied in afterwards.
@inline function _push_power!(bases, exps, b::Sym, p::Sym, c::SymNumber, left)
    p.kind == NUMBER && return (c * p.num, left)
    pb, pe = _split_pow(p)
    if p.kind != MUL && _eq(pb, b)
        push!(bases, b)
        push!(exps, pe)
        return (c, left)
    end
    left === nothing && (left = Sym[])
    push!(left, p)
    (c, left)
end

function _finish_mul(c, bases, exps, left)
    r = _mkmul(c, bases, exps)
    if left !== nothing
        for p in left
            r = mul(r, p)
        end
    end
    r
end

function _mul_factor(M::Sym, x::Sym)
    b, e = _split_pow(x)
    i, found = _search(M.args, b)
    found || return _mkmul(M.num, _insert(M.args, i, b), _insert(M.exps, i, e))
    p = pow(b, add(M.exps[i], e))
    if p.kind == NUMBER
        return _mkmul(M.num * p.num, _deleteat(M.args, i), _deleteat(M.exps, i))
    end
    pb, pe = _split_pow(p)
    if p.kind != MUL && _eq(pb, b)
        exps = copy(M.exps)
        exps[i] = pe
        return _mkmul(M.num, M.args, exps)
    end
    mul(_mkmul(M.num, _deleteat(M.args, i), _deleteat(M.exps, i)), p)
end

function _merge_mul(A::Sym, B::Sym)
    a, b, ae, be = A.args, B.args, A.exps, B.exps
    na, nb = length(a), length(b)
    bases = Sym[]
    exps = Sym[]
    sizehint!(bases, na + nb)
    sizehint!(exps, na + nb)
    c = A.num * B.num
    left = nothing
    i, j = 1, 1
    @inbounds while i <= na && j <= nb
        k = _cmp(a[i], b[j])
        if k < 0
            push!(bases, a[i]); push!(exps, ae[i]); i += 1
        elseif k > 0
            push!(bases, b[j]); push!(exps, be[j]); j += 1
        else
            c, left = _push_power!(bases, exps, a[i], pow(a[i], add(ae[i], be[j])), c, left)
            i += 1; j += 1
        end
    end
    @inbounds while i <= na
        push!(bases, a[i]); push!(exps, ae[i]); i += 1
    end
    @inbounds while j <= nb
        push!(bases, b[j]); push!(exps, be[j]); j += 1
    end
    _finish_mul(c, bases, exps, left)
end

## powers

"""
    pow(b::Sym, e::Sym)

Canonical power (what `b^e` calls). Integer powers of products and powers are
distributed; numeric powers are evaluated exactly, extracting perfect powers from
radicals (`8^(1/2) == 2*2^(1/2)`).
"""
function pow(b::Sym, e::Sym)
    bk = b.kind
    if e.kind == NUMBER
        en = e.num
        iszero(en) && return ONE
        _isN1(en) && return b
        bk == NUMBER && return _numpow(b.num, en)
        isfloat(en) && (v = _floatpow(b, en)) !== nothing && return v
        if _isint(en)
            bk == MUL && return _mulpow_int(b, e)
            bk == POW && return pow(b.args[1], mul(b.args[2], e))
        elseif isexact(en)
            # (b^p)^q = b^(p q) for b > 0 (a number or e, π) and p real: (3^(1/3))^(1/4)
            # is 3^(1/12)
            if bk == POW && b.args[2].kind == NUMBER && isexact(b.args[2].num) &&
               ((b.args[1].kind == NUMBER && isexact(b.args[1].num) && b.args[1].num > N0) || b.args[1].kind == CONSTANT)
                return pow(b.args[1], mul(b.args[2], e))
            end
            # a positive constant factor comes out: (c m)^e = c^e m^e for c > 0 on the
            # principal branch (scaling by c > 0 leaves the argument unchanged), so
            # sqrt(4a) = 2 sqrt(a) and sqrt(2a + 2b) = sqrt(2) sqrt(a + b)
            c = _positive_content(b)
            c === nothing || return mul(_numpow(c, en), pow(_scale(b, inv(c)), e))
        end
        return _mkpow_node(b, e)
    end
    if bk == NUMBER
        _isN1(b.num) && return ONE
        # n^(u p/q) with n a perfect q-th power: 256^(x/2) = 16^x (SymPy #24090)
        if e.kind == MUL && isexact(e.num) && !_isint(e.num) && isexact(b.num) && _isint(b.num) && b.num > N1
            q = _denominator(e.num)
            r = _int_root(big(_numerator(b.num)), q)
            r === nothing || return pow(Sym(SymNumber(r)), mul(Sym(SymNumber(q)), e))
        end
    elseif b === E && e.kind == FN && e.name === :log && length(e.args) == 1
        return e.args[1]                        # exp(log(z)) == z
    elseif b === E
        # exp(iπ r), r rational: r reduced into (-1, 1]; exp(10iπ) = 1, exp(12iπ/5) =
        # exp(2iπ/5) (SymPy #24082)
        r = _ipi_coeff(e)
        if r !== nothing
            q = Rational{BigInt}(_numerator(r), _denominator(r))
            q2 = q - 2 * floor((q + 1) / 2)              # in [-1, 1)
            q2 == -1 && (q2 = 1 // 1)
            iszero(q2) && return ONE
            isone(q2) && return MINUS_ONE
            q2 == 1 // 2 && return IM
            q2 == -1 // 2 && return mul(MINUS_ONE, IM)
            q2 == q || return pow(E, mul(mul(Sym(SymNumber(q2)), IM), PI))
        end
    end
    _mkpow_node(b, e)
end

# the integer r with r^q == n (n > 1), or nothing
function _int_root(n::BigInt, q::Integer)
    q <= 1 && return nothing
    q > 64 * ndigits(n) && return nothing
    r = round(BigInt, BigFloat(n; precision = 4 * ndigits(n, base = 2) + 64)^(1 / BigFloat(q; precision = 4 * ndigits(n, base = 2) + 64)))
    for c in (r - 1, r, r + 1)
        c > 1 && c^q == n && return c
    end
    nothing
end

# r for e = r i π, r an exact rational; else nothing
function _ipi_coeff(e::Sym)
    e.kind == MUL && isexact(e.num) && length(e.args) == 2 || return nothing
    hasi = hasp = false
    for (a, p) in zip(e.args, e.exps)
        if a.kind == NUMBER && _isN1(-a.num) && p.kind == NUMBER && p.num == SymNumber(1 // 2)
            hasi = true
        elseif a === PI && p.kind == NUMBER && _isN1(p.num)
            hasp = true
        end
    end
    hasi && hasp ? e.num : nothing
end

# a constant to a floating point power, evaluated as _numpow does a number's (exp(0.5),
# π^1.5); nothing when b has symbols or is negative under a fractional power
function _floatpow(b::Sym, e::SymNumber)::Union{Sym,Nothing}
    (_hassymbol(b) || _hasparameter(b)) && return nothing
    # (behind an inference barrier: N reaches the evaluator and every function's rules, and
    # back into arithmetic; inferring *(Int64, Sym) ran past 2501 frames on Julia 1.13)
    bv = try Base.inferencebarrier(N)(b) catch; return nothing end
    bv isa Real || return nothing
    ef = Float64(e)
    bv < 0 && !isinteger(ef) && return nothing
    Sym(Float64(bv)^ef)
end

# the positive rational content of a product's coefficient or of a sum's coefficients
# (gcd of numerators over lcm of denominators), when it is not 1; else nothing
function _positive_content(b::Sym)
    k = b.kind
    if k == MUL
        isexact(b.num) || return nothing
        c = abs(b.num)
        return _isN1(c) ? nothing : c
    elseif k == ADD
        cs = iszero(b.num) ? b.coefs : vcat(b.coefs, b.num)
        all(isexact, cs) || return nothing
        num = reduce(gcd, (big(_numerator(c)) for c in cs); init = big(0))
        den = reduce(lcm, (big(_denominator(c)) for c in cs); init = big(1))
        c = SymNumber(num // den)
        return _isN1(c) ? nothing : c
    end
    nothing
end

function _mulpow_int(M::Sym, e::Sym)
    # a coefficient whose power is beyond any memory stays a power, as a number's does
    # ((10^9 x)^(10^9); SymPy #27278)
    if isexact(M.num) && !_isN1(M.num) && _toobig_pow(M.num, e.num)
        return mul(_mkpow_node(Sym(M.num), e), pow(_scale(M, inv(M.num)), e))
    end
    c = M.num^_toint(e.num)
    bases = Sym[]
    exps = Sym[]
    left = nothing
    for i in eachindex(M.args)
        b = M.args[i]
        c, left = _push_power!(bases, exps, b, pow(b, mul(M.exps[i], e)), c, left)
    end
    _finish_mul(c, bases, exps, left)
end

# (-1)^e with e reduced into (-1, 1], then a negative exponent written with the sign
# outside, (-1)^r = -(-1)^(r + 1): im^3 and 1/im are -im, as 1/sqrt(2) is sqrt(2)/2
function _minus_one_pow(e::SymNumber)
    r = _rat(e)
    r -= 2 * ceil(BigInt, (r - 1) / 2)
    isone(r) && return MINUS_ONE
    iszero(r) && return ONE
    r < 0 && return mul(MINUS_ONE, _mkpow_node(MINUS_ONE, Sym(_frombig(r + 1))))
    _mkpow_node(MINUS_ONE, Sym(_frombig(r)))
end

function _numpow(b::SymNumber, e::SymNumber)
    if isfloat(b) || isfloat(e)
        bf, ef = Float64(b), Float64(e)
        # an exact base to an infinite power: 0 stays exact (0^∞, (1/2)^∞, 2^-∞), as in limits
        !isfloat(b) && isinf(ef) && abs(bf) != 1 && (abs(bf) < 1) == (ef > 0) && return ZERO
        bf < 0 && !isinteger(ef) && return _mkpow_node(Sym(b), Sym(e))
        return Sym(bf^ef)
    end
    if _isint(e)
        iszero(b) && sign(e) < 0 && throw(DivisionByZero("division by zero"))
        iszero(e) && return ONE
        iszero(b) && return ZERO
        isone(b) && return ONE
        _isN1(-b) && return iseven(_numerator(e)) ? ONE : MINUS_ONE
        # a result beyond any memory (5^(5^3125)) stays a power
        _toobig_pow(b, e) && return _mkpow_node(Sym(b), Sym(e))
        return Sym(b^_toint(e))
    end
    if iszero(b)
        sign(e) > 0 && return ZERO
        throw(DivisionByZero("0 to the negative power $e"))
    end
    _isN1(b) && return ONE
    sign(b) < 0 && return mul(_minus_one_pow(e), _numpow(-b, e))
    _isint(b) || return mul(_numpow(SymNumber(_numerator(b)), e), _numpow(SymNumber(_denominator(b)), -e))
    _introot_pow(_numerator(b), e)
end

# b^e (e an integer) with more than about 10^7 bits in its numerator or denominator
const _MAX_POW_BITS = 10^7
function _toobig_pow(b::SymNumber, e::SymNumber)
    n = abs(_numerator(e))
    n <= 64 && return false
    bits = max(log2(abs(big(_numerator(b)))), log2(big(_denominator(b))))
    n * bits > _MAX_POW_BITS
end

# n^(p/q) for an integer n > 1 and q > 1
function _introot_pow(n::Integer, e::SymNumber)
    p, q = _numerator(e), _denominator(e)
    # a root index past Int (2^(1/10^50)): no perfect power to find, the power as it is
    # (Int(q) failed; SymPy #18526)
    q > typemax(Int) && return _mkpow_node(Sym(SymNumber(n)), Sym(e))
    k, r = fldmod(p, q)
    coef = SymNumber(n)^k
    a, c = _extract_power(n, Int(q))
    coef *= SymNumber(a)^r
    isone(c) && return Sym(coef)
    # c may still be a perfect j-th power with j | q, j at most c's bit length (not every
    # j up to a denominator of 5·10¹⁰; SymPy #21860)
    for j in min(Int(q), Base.top_set_bit(BigInt(c))):-1:2
        q % j == 0 || continue
        d = _perfect_root(c, j)
        d === nothing || return mul(Sym(coef), _numpow(SymNumber(d), SymNumber(r * j // q)))
    end
    _mkmul(coef, Sym[Sym(SymNumber(c))], Sym[Sym(SymNumber(r // q))])
end

## operators

Base.:+(a::Sym, b::Sym) = add(a, b)
Base.:*(a::Sym, b::Sym) = mul(a, b)
Base.:-(a::Sym) = _scale(a, NM1)
Base.:+(a::Sym) = a                  # a one-term sum or product, as for numbers
Base.:*(a::Sym) = a
Base.:-(a::Sym, b::Sym) = add(a, _scale(b, NM1))
Base.:/(a::Sym, b::Sym) = mul(a, pow(b, MINUS_ONE))
Base.://(a::Sym, b::Sym) = a / b
Base.:^(a::Sym, b::Sym) = pow(a, b)
Base.:^(a::Sym, b::Integer) = pow(a, Sym(b))
Base.inv(a::Sym) = pow(a, MINUS_ONE)

for op in (:+, :-, :*, :/, ://, :^)
    @eval Base.$op(a::Sym, b::Number) = $op(a, Sym(b))
    @eval Base.$op(a::Number, b::Sym) = $op(Sym(a), b)
end

Base.sum(xs::AbstractArray{Sym}) = _add_many(xs)

"""
    sqrt(x::Sym)

`x^(1//2)`, exact for numbers: `sqrt(Sym(12)) == 2*sqrt(Sym(3))`.
"""
Base.sqrt(x::Sym) = pow(x, HALF)
Base.exp(x::Sym) = pow(E, x)
