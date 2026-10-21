# The expression node.
#
# Every expression is a `Sym`: one concrete, immutable-by-convention node type with a
# `kind` tag (in the spirit of SymEngine's TypeID). A single concrete type keeps child
# vectors concretely typed (`Vector{Sym}`) and lets the hot paths branch on a byte
# instead of dispatching dynamically. The structural hash is computed once at
# construction, so equality and dictionary lookups are cheap.

@enum Kind::UInt8 NUMBER SYMBOL CONSTANT ADD MUL POW FN

"""
    AbstractSymbolic <: Function

Supertype of SymJ's symbolic objects. Symbolic objects are callable: calling an
expression substitutes into it (see [`subs`](@ref)) and calling a [`SymFunction`](@ref)
builds an unevaluated function application.
"""
abstract type AbstractSymbolic <: Function end

"""
    Sym

A symbolic expression. Construct leaves with `Sym(:x)`, `Sym(2)`, `Sym(1//3)`,
`Sym(π)`, or the [`@syms`](@ref) macro; build compound expressions with the usual
arithmetic operators and functions. Expressions are kept in a canonical form as they are
built: like terms are collected in sums, like bases are collected in products, and
numeric powers such as `sqrt(Sym(8))` are reduced exactly (to `2*sqrt(2)`).

Layout by `kind`:

| kind       | `num`              | `args`               | `exps`         | `coefs`                |
|:-----------|:-------------------|:---------------------|:---------------|:-----------------------|
| `NUMBER`   | the value          |                      |                |                        |
| `SYMBOL`   |                    |                      |                |                        |
| `CONSTANT` |                    |                      |                |                        |
| `ADD`      | constant term      | terms (sorted)       |                | coefficient of each    |
| `MUL`      | coefficient        | bases (sorted)       | exponent of each |                      |
| `POW`      |                    | `[base, exponent]`   |                |                        |
| `FN`       |                    | arguments            |                |                        |

`name` holds the symbol, constant, or function name. The vectors must never be mutated
once a node is built; nodes share them freely.
"""
mutable struct Sym <: AbstractSymbolic
    const kind::Kind
    const hash::UInt
    const name::Symbol
    const num::SymNumber
    const args::Vector{Sym}
    const exps::Vector{Sym}
    const coefs::Vector{SymNumber}
end

const NOSYMS = Sym[]
const NONUMS = SymNumber[]
const NONAME = Symbol("")

const H_NUM   = 0x6a09e667f3bcc908 % UInt
const H_SYM   = 0xbb67ae8584caa73b % UInt
const H_CONST = 0x3c6ef372fe94f82b % UInt
const H_ADD   = 0xa54ff53a5f1d36f1 % UInt
const H_MUL   = 0x510e527fade682d1 % UInt
const H_POW   = 0x9b05688c2b3e6c1f % UInt
const H_FN    = 0x1f83d9abfb41bd6b % UInt

## Raw node constructors: no simplification, inputs must already be canonical.

_mknum(n::SymNumber) = Sym(NUMBER, hash(n, H_NUM + n.tag), NONAME, n, NOSYMS, NOSYMS, NONUMS)
_mksym(s::Symbol) = Sym(SYMBOL, hash(String(s), H_SYM), s, N0, NOSYMS, NOSYMS, NONUMS)
_mkconst(s::Symbol) = Sym(CONSTANT, hash(String(s), H_CONST), s, N0, NOSYMS, NOSYMS, NONUMS)

function _mkadd_node(c::SymNumber, terms::Vector{Sym}, coefs::Vector{SymNumber})
    h = hash(c, H_ADD)
    @inbounds for i in eachindex(terms)
        h = hash(coefs[i], h ⊻ terms[i].hash)
    end
    Sym(ADD, h, NONAME, c, terms, NOSYMS, coefs)
end

function _mkmul_node(c::SymNumber, bases::Vector{Sym}, exps::Vector{Sym})
    h = hash(c, H_MUL)
    @inbounds for i in eachindex(bases)
        h = hash(exps[i].hash, h ⊻ bases[i].hash)
    end
    Sym(MUL, h, NONAME, c, bases, exps, NONUMS)
end

_mkpow_node(b::Sym, e::Sym) = Sym(POW, hash(e.hash, b.hash ⊻ H_POW), NONAME, N0, Sym[b, e], NOSYMS, NONUMS)

function _mkfn_node(name::Symbol, args::Vector{Sym})
    h = hash(String(name), H_FN)
    for a in args
        h = hash(a.hash, h)
    end
    Sym(FN, h, name, N0, args, NOSYMS, NONUMS)
end

## Leaves

const _INTCACHE = [_mknum(SymNumber(i)) for i in -64:64]

function Sym(n::SymNumber)
    if n.tag == T_SMALL && n.d == 1 && -64 <= n.n <= 64
        return @inbounds _INTCACHE[n.n + 65]
    end
    _mknum(n)
end

Sym(x::Sym) = x
Sym(x::Union{Integer,Rational,AbstractFloat}) = Sym(SymNumber(x))
Sym(s::Symbol) = _mksym(s)
Sym(s::AbstractString) = Sym(Symbol(s))

const ZERO = Sym(N0)
const ONE = Sym(N1)
const MINUS_ONE = Sym(NM1)
const HALF = Sym(SymNumber(1 // 2))
const PI = _mkconst(:pi)
const E = _mkconst(:E)
const IM = _mkpow_node(MINUS_ONE, HALF)   # sqrt(-1)

Sym(::Irrational{:π}) = PI
Sym(::Irrational{:ℯ}) = E
Sym(z::Complex) = Sym(real(z)) + Sym(imag(z)) * IM

## Structural order, equality, hashing

Base.hash(x::Sym, h::UInt) = hash(x.hash, h)

function _cmpvec(a::Vector{Sym}, b::Vector{Sym})
    la, lb = length(a), length(b)
    la != lb && return la < lb ? -1 : 1
    @inbounds for i in 1:la
        c = _cmp(a[i], b[i])
        c != 0 && return c
    end
    0
end

"""
    _cmp(a, b)

Total order on expressions used for canonical storage (hash first, structure on ties).
It is deterministic across sessions but has no mathematical meaning; display order is
computed separately by the printer.
"""
function _cmp(a::Sym, b::Sym)
    a === b && return 0
    a.hash != b.hash && return a.hash < b.hash ? -1 : 1
    a.kind != b.kind && return a.kind < b.kind ? -1 : 1
    k = a.kind
    if k == NUMBER
        return _numcmp(a.num, b.num)
    elseif k == SYMBOL || k == CONSTANT
        return a.name === b.name ? 0 : (isless(a.name, b.name) ? -1 : 1)
    elseif k == ADD
        c = _numcmp(a.num, b.num)
        c != 0 && return c
        c = _cmpvec(a.args, b.args)
        c != 0 && return c
        for i in eachindex(a.coefs)
            c = _numcmp(a.coefs[i], b.coefs[i])
            c != 0 && return c
        end
        return 0
    elseif k == MUL
        c = _numcmp(a.num, b.num)
        c != 0 && return c
        c = _cmpvec(a.args, b.args)
        return c != 0 ? c : _cmpvec(a.exps, b.exps)
    elseif k == POW
        return _cmpvec(a.args, b.args)
    else
        a.name !== b.name && return isless(a.name, b.name) ? -1 : 1
        return _cmpvec(a.args, b.args)
    end
end

@inline _lt(a::Sym, b::Sym) = _cmp(a, b) < 0
@inline _eq(a::Sym, b::Sym) = a === b || (a.hash == b.hash && _cmp(a, b) == 0)

Base.:(==)(a::Sym, b::Sym) = _eq(a, b)
Base.:(==)(a::Sym, b::Number) = _eq(a, Sym(b))
Base.:(==)(a::Number, b::Sym) = _eq(Sym(a), b)
Base.isequal(a::Sym, b::Sym) = _eq(a, b)

## Introspection

"""
    kind(x::Sym)

The node kind: one of `NUMBER`, `SYMBOL`, `CONSTANT`, `ADD`, `MUL`, `POW`, `FN`
(accessible as `SymJ.ADD`, etc.).
"""
kind(x::Sym) = x.kind

isnumber(x::Sym) = x.kind == NUMBER
issymbol(x::Sym) = x.kind == SYMBOL

@inline _isone(x::Sym) = x.kind == NUMBER && _isN1(x.num)

Base.iszero(x::Sym) = x.kind == NUMBER && iszero(x.num)
Base.isone(x::Sym) = x.kind == NUMBER && isone(x.num)
Base.zero(::Type{Sym}) = ZERO
Base.zero(::Sym) = ZERO
Base.one(::Type{Sym}) = ONE
Base.one(::Sym) = ONE

Base.convert(::Type{Sym}, x::Sym) = x
Base.convert(::Type{Sym}, x::Number) = Sym(x)
Base.promote_rule(::Type{Sym}, ::Type{<:Number}) = Sym
