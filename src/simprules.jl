# Rules used by simplify: trigonometric and hyperbolic identities that the normal form
# of simplify_trig does not reach (it reduces modulo sin^2 + cos^2 = 1 and linearizes
# products, so it does not go back to half angles or squares). simplify offers the
# rewritten expression as one more candidate, kept only when it is smaller.

# whether v is 2u with the factor 2 showing: an even multiple, or a sum of them
_halvable(v::Sym) = v.kind == MUL ? _iseven_num(v.num) :
    v.kind == ADD && (iszero(v.num) || _iseven_num(v.num)) && all(_iseven_num, v.coefs)
_iseven_num(c) = _isint(c) && iseven(_numerator(c))

# whether e has an application of one of the named functions
_has_fn(e::Sym, names) = (e.kind == FN && e.name in names) || any(a -> _has_fn(a, names), e.args) ||
    (e.kind == MUL && any(p -> _has_fn(p, names), e.exps))
const _CIRCFNS = (_TRIGMODE.fns..., _HYPMODE.fns...)

const _SIMP_TRIG_RULES = _RuleSet(Rule[
    # power reduction, read backwards: for 2u with the 2 written (not cos(x) as
    # cos(2 * (x/2))), and ~!c so that the coefficient can be 1 (a one-or-more slot
    # takes no implicit 1 from cos(2u) alone)
    @rule(~!c - ~!c * cos(~v::_halvable) + ~!r => 2 * ~c * sin(~v / 2)^2 + ~r),
    @rule(~!c + ~!c * cos(~v::_halvable) + ~!r => 2 * ~c * cos(~v / 2)^2 + ~r),
    @rule(~!c * cosh(~v::_halvable) - ~!c + ~!r => 2 * ~c * sinh(~v / 2)^2 + ~r),
    @rule(~!c + ~!c * cosh(~v::_halvable) + ~!r => 2 * ~c * cosh(~v / 2)^2 + ~r),
    # half angles
    @rule(~!r * sin(~u) / (1 + cos(~u)) => ~r * tan(~u / 2)),
    @rule(~!r * (1 - cos(~u)) / sin(~u) => ~r * tan(~u / 2)),
    @rule(~!r * (1 - cos(~u)) * csc(~u) => ~r * tan(~u / 2)),
])

_simp_trig_rules(e::Sym) = _rewrite(_SIMP_TRIG_RULES, e)
