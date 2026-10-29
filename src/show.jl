# Printing. Output is valid Julia syntax where possible (`im`, `π`, `ℯ`, `sqrt`, `exp`).

const P_ADD, P_MUL, P_POW, P_ATOM = 1, 2, 3, 4

Base.show(io::IO, x::Sym) = _show(io, x, 0)
Base.show(io::IO, ::MIME"text/plain", x::Sym) = show(io, x)

_negexp(e::Sym) = (e.kind == NUMBER || e.kind == MUL) && sign(e.num) < 0
_ishalf(e::Sym) = e.kind == NUMBER && e.num == SymNumber(1 // 2)

# precedence of the printed form of x
function _prec(x::Sym)
    k = x.kind
    if k == NUMBER
        sign(x.num) < 0 && return P_ADD
        return _isint(x.num) || isfloat(x.num) ? P_ATOM : P_MUL
    end
    k == ADD && return P_ADD
    k == MUL && return sign(x.num) < 0 ? P_ADD : P_MUL
    if k == POW
        b, e = x.args
        _negexp(e) && return P_MUL
        (b === E || _ishalf(e)) && return P_ATOM
        return P_POW
    end
    P_ATOM
end

function _show(io::IO, x::Sym, p::Int)
    wrap = _prec(x) < p
    wrap && print(io, "(")
    k = x.kind
    if k == NUMBER
        _show_num(io, x.num)
    elseif k == SYMBOL || k == PARAMETER
        print(io, x.name)
    elseif k == CONSTANT
        print(io, x.name === :pi ? "π" : x.name === :E ? "ℯ" : x.name)
    elseif k == ADD
        _show_add(io, x)
    elseif k == MUL
        _show_mul(io, x.num, x.args, x.exps)
    elseif k == POW
        b, e = x.args
        _negexp(e) ? _show_mul(io, N1, x.args[1:1], x.args[2:2]) : _show_pow(io, b, e)
    elseif x.name === :piecewise
        _show_piecewise(io, x)
    elseif x.name === :DifferenceRoot && length(x.args) >= 8
        _show_droot(io, x)
    elseif _iscond(x)
        _show_cond(io, x)
    elseif x.name === :_slot
        print(io, ("~", "~!", "~~", "~~~")[_toint(x.args[2].num)], x.args[1].name)
    elseif x.name === :_hslot
        print(io, "(~", x.args[1].name, ")(")
        join(io, (sprint(_show, a, 0) for a in x.args[2:end]), ", ")
        print(io, ")")
    elseif (d = _deriv_at(x)) !== nothing
        print(io, d[1], "'"^d[2], "(", sprint(_show, d[3], 0), ")")
    elseif _isdefint(x)                         # as integrate takes it: Integral(f, (y, a, b))
        f, v, a, b = (sprint(_show, t, 0) for t in x.args)
        print(io, "Integral(", f, ", (", v, ", ", a, ", ", b, "))")
    else
        print(io, x.name === :atan2 ? :atan : x.name, "(")
        join(io, (sprint(_show, a, 0) for a in x.args), ", ")
        print(io, ")")
    end
    wrap && print(io, ")")
    nothing
end

# Derivative(f(v), v, a), Derivative(Derivative(f(v), v), v, a), …: (f, k, a), printed
# f'(a), f''(a) as the SymFunction derivative notation
function _deriv_at(x::Sym)
    (x.kind == FN && x.name === :Derivative && length(x.args) == 3) || return nothing
    f, v, a = x.args
    k = 1
    while f.kind == FN && f.name === :Derivative && length(f.args) == 2 && _eq(f.args[2], v)
        f, k = f.args[1], k + 1
    end
    (f.kind == FN && length(f.args) == 1 && _eq(f.args[1], v)) || return nothing
    f.name, k, a
end

function _show_num(io::IO, n::SymNumber)
    if isfloat(n)
        print(io, n.f)
    elseif _isint(n)
        print(io, _numerator(n))
    else
        print(io, _numerator(n), "/", _denominator(n))
    end
end

function _show_pow(io::IO, b::Sym, e::Sym)
    if b === E
        print(io, "exp(")
        _show(io, e, 0)
        print(io, ")")
    elseif _ishalf(e)
        if b == MINUS_ONE
            print(io, "im")
        else
            print(io, "sqrt(")
            _show(io, b, 0)
            print(io, ")")
        end
    else
        _show(io, b, P_ATOM)
        print(io, "^")
        _show(io, e, P_ATOM)
    end
end

_factor_str(b::Sym, e::Sym) = _isone(e) ? sprint(_show, b, P_MUL) : sprint(_show_pow, b, e)

function _show_mul(io::IO, c::SymNumber, bases::Vector{Sym}, exps::Vector{Sym})
    num, den = String[], String[]
    for (b, e) in _display_sort_pairs(bases, exps)
        _negexp(e) ? push!(den, _factor_str(b, -e)) : push!(num, _factor_str(b, e))
    end
    neg = sign(c) < 0
    c = abs(c)
    if isfloat(c)
        isone(c) || pushfirst!(num, string(c.f))
    else
        isone(_numerator(c)) || pushfirst!(num, string(_numerator(c)))
        isone(_denominator(c)) || pushfirst!(den, string(_denominator(c)))
    end
    neg && print(io, "-")
    print(io, isempty(num) ? "1" : join(num, "*"))
    if !isempty(den)
        print(io, "/")
        length(den) == 1 ? print(io, den[1]) : print(io, "(", join(den, "*"), ")")
    end
end

function _show_add(io::IO, x::Sym)
    terms = _display_sort(Sym[_with_coef(x.coefs[i], x.args[i]) for i in eachindex(x.args)])
    iszero(x.num) || push!(terms, Sym(x.num))
    for (i, t) in enumerate(terms)
        s = sprint(_show, t, P_ADD)
        if i == 1
            print(io, s)
        elseif startswith(s, "-")
            print(io, " - ", chop(s; head = 1, tail = 0))
        else
            print(io, " + ", s)
        end
    end
end

## display order: higher degree first, then by printed form

function _degree(x::Sym)::Float64
    k = x.kind
    k == SYMBOL && return 1.0
    k == PARAMETER && return 0.0      # order terms by the variables: a*x^2 + b*x + c
    (k == NUMBER || k == CONSTANT || k == FN) && return 0.0
    k == ADD && return maximum(_degree, x.args)
    if k == MUL
        return sum(_degree(_mkpow(x.args[i], x.exps[i])) for i in eachindex(x.args))
    end
    b, e = x.args
    _degree(b) * (e.kind == NUMBER ? Float64(e.num) : 1.0)
end

function _display_sort(xs::Vector{Sym})
    keys = [(-_degree(x), _sortstring(x)) for x in xs]
    xs[sortperm(keys)]
end

function _display_sort_pairs(bases::Vector{Sym}, exps::Vector{Sym})
    keys = [(bases[i].kind == NUMBER ? 0 : 1, sprint(_show, bases[i], 0)) for i in eachindex(bases)]
    p = sortperm(keys)
    zip(bases[p], exps[p])
end

# sort key that ignores a leading coefficient, so `2x` sorts with `x`
function _sortstring(x::Sym)
    x.kind == MUL || return sprint(_show, x, 0)
    sprint(_show_mul, N1, x.args, x.exps)
end
