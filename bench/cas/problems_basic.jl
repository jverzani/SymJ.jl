# Problems for the everyday tasks of a computer algebra system, as quoted expressions so
# each worker builds them with its own variables. Each entry is
#   (area, op, args, expect, source)
# `op` and `args`:
#   :expand     (f,)            expand products and powers
#   :factor     (f,)            factor over the rationals; f is expanded first (untimed)
#   :gcd        (f, g)          polynomial gcd; f and g are expanded first (untimed)
#   :cancel     (f,)            reduce a rational function to lowest terms
#   :apart      (f, x)          partial fractions in x
#   :diff       (f, x, n)       the n-th derivative in x
#   :solve      (lhs, rhs, x)   the solutions of lhs = rhs in x
#   :linsolve   (eqs, vars)     a linear system, eqs a vector of (lhs, rhs)
# `expect`, checked after the value (see core.jl):
#   :expand, :apart   the number of terms of the result
#   :factor           the number of irreducible factors, counted with multiplicity
#   :gcd              the gcd, up to a constant factor
#   :cancel           the degree in x of the reduced denominator
#   :diff             nothing (the value is the check)
#   :solve            the number of distinct real solutions
#   :linsolve         nothing (the residuals are the check)
# Variables: x y z w t a b c d, and v1 … v12 for many variables.
[
    # ---- expansion
    ("Expansion", :expand, (:((x + y)^10),), 11, "binomial"),
    ("Expansion", :expand, (:((x + y)^100),), 101, "binomial, big coefficients"),
    ("Expansion", :expand, (:((a + b + c)^6),), 28, "trinomial"),
    ("Expansion", :expand, (:((x + y + z + 1)^8),), 165, "four terms"),
    ("Expansion", :expand, (:((x + y + z + w + 1)^15),), 3876, "five terms, large"),
    ("Expansion", :expand, (:((x + 1)^20 * (x - 1)^5),), 26, "product of powers"),
    ("Expansion", :expand, (:((x + y)^3 * (x - y)^3),), 4, "cancellation"),
    ("Expansion", :expand, (:((1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9)^10),), 91, "dense base"),
    ("Expansion", :expand, (:((a*x + b)^4 * (c*x + d)^2),), 15, "symbolic coefficients"),
    ("Expansion", :expand, (:((sqrt(2)*x + 1)^4),), 5, "radical coefficients"),
    ("Expansion", :expand, (:((x + 1/x)^6),), 7, "Laurent polynomial"),
    ("Expansion", :expand, (:((sin(x) + cos(x))^2),), 3, "function arguments"),

    # ---- factorization
    ("Factorization", :factor, (:(x^4 - 1),), 3, "difference of squares"),
    ("Factorization", :factor, (:(x^12 - 1),), 6, "cyclotomic"),
    ("Factorization", :factor, (:(x^15 - 1),), 4, "cyclotomic"),
    ("Factorization", :factor, (:(x^4 + 4),), 2, "Sophie Germain"),
    ("Factorization", :factor, (:(2x^4 + 4x^3 + 2x^2),), 4, "content and repeated factors"),
    ("Factorization", :factor, (:((x + 1)^3 * (x^2 + x + 1)^2 * (x^3 - 2)),), 6, "multiplicities"),
    ("Factorization", :factor, (:((123456789x + 987654321) * (x^2 - 1234567)),), 2, "large coefficients"),
    ("Factorization", :factor, (:(x^8 - 40x^6 + 352x^4 - 960x^2 + 576),), 1, "irreducible (Swinnerton-Dyer)"),
    ("Factorization", :factor, (:(x^2 * y^2 - 1),), 2, "bivariate"),
    ("Factorization", :factor, (:(x^6 - y^6),), 4, "bivariate, homogeneous"),
    ("Factorization", :factor, (:((x + y)^5 * (x - y)^3),), 8, "bivariate, multiplicities"),
    ("Factorization", :factor, (:((x^2 + y*z + 1) * (x*y - z^2 + 3) * (x + y + z)),), 3, "trivariate"),
    ("Factorization", :factor, (:((x^3*y + 2z*w - 1) * (x*y^2 - z^3 + w) * (x + y + z + w + 3)),), 3, "four variables"),

    # ---- gcd
    ("GCD", :gcd, (:(x^2 - 1), :(x^2 + 2x + 1)), :(x + 1), "univariate"),
    ("GCD", :gcd, (:((x + 1)^30 * (x - 2)^20), :((x + 1)^25 * (x + 3)^20)), :((x + 1)^25), "univariate, high degree"),
    ("GCD", :gcd, (:(x^3 + y), :(x^2 - y^3)), 1, "coprime"),
    ("GCD", :gcd, (:((x + y + z + w + 1)^4 * ((x - y + 2z*w + 1)^2 + y)),
                   :((x + y + z + w + 1)^4 * ((y*z - x + 3)^2 + w))), :((x + y + z + w + 1)^4), "dense, four variables"),
    ("GCD", :gcd, (:((x + y + z + w + 1)^6 * ((x - y + 2z*w + 1)^3 + y)),
                   :((x + y + z + w + 1)^6 * ((y*z - x + 3)^3 + w))), :((x + y + z + w + 1)^6), "dense, higher degree"),
    ("GCD", :gcd, (:((3486784401x + 30517578125y*z - 282475249w + 1)^3 * ((x - y + z)^3 + 2)),
                   :((3486784401x + 30517578125y*z - 282475249w + 1)^3 * ((x*w + y - 1)^2 - z))),
                  :((3486784401x + 30517578125y*z - 282475249w + 1)^3), "large coefficients"),
    ("GCD", :gcd, (:((x^10*y^3*z + 2y^7*z^5*w - x^2*w^8 + 3y^9 + x*y*z*w) * (x^9*z^4 + y^8*w^2 + 5)),
                   :((x^10*y^3*z + 2y^7*z^5*w - x^2*w^8 + 3y^9 + x*y*z*w) * (y^10*w^3 - x^6*w^5 + 7z))),
                  :(x^10*y^3*z + 2y^7*z^5*w - x^2*w^8 + 3y^9 + x*y*z*w), "sparse"),
    ("GCD", :gcd, (:(((v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8 + v9 + v10 + v11 + v12)^2 + v1*v2) * (v1 - v3)),
                   :(((v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8 + v9 + v10 + v11 + v12)^2 + v1*v2) * (v2 + v4 - 1))),
                  :((v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8 + v9 + v10 + v11 + v12)^2 + v1*v2), "twelve variables"),

    # ---- rational functions
    ("Rational functions", :cancel, (:((x^2 - 1)/(x + 1)),), 0, "polynomial result"),
    ("Rational functions", :cancel, (:((x^4 - 1)/(x^2 - 1)),), 0, "polynomial result"),
    ("Rational functions", :cancel, (:((x^3 - y^3)/(x^2 - y^2)),), 1, "bivariate"),
    ("Rational functions", :cancel, (:(1/x + 1/x^2),), 2, "sum of fractions"),
    ("Rational functions", :cancel, (:(1/(x - 1) - 1/(x + 1) - 2/(x^2 - 1)),), 0, "sum to zero"),
    ("Rational functions", :cancel, (:((x + y)^5 * (x - z)^3 / ((x + y)^3 * (x + z))),), 1, "common power"),
    ("Rational functions", :cancel, (:((x + y + z + w + 1)^4 * ((x - y + 2z*w + 1)^2 + y) /
                                       ((x + y + z + w + 1)^4 * ((y*z - x + 3)^2 + w))),), 2, "dense, four variables"),
    ("Rational functions", :cancel, (:((x^2 + 2*sqrt(2)*x + 2)/(x + sqrt(2))),), 0, "radical coefficients"),

    ("Partial fractions", :apart, (:(1/(x^2 - 1)), :x), 2, "two linear factors"),
    ("Partial fractions", :apart, (:((x^3 + 1)/(x^2 - 1)), :x), 2, "polynomial part"),
    ("Partial fractions", :apart, (:(1/(x*(x + 1)^2)), :x), 3, "repeated factor"),
    ("Partial fractions", :apart, (:((2x + 3)/((x - 1)*(x^2 + 1))), :x), 2, "quadratic factor"),
    ("Partial fractions", :apart, (:(1/(x^4 - 1)), :x), 3, "three factors"),
    ("Partial fractions", :apart, (:(1/((x - a)*(x - b))), :x), 2, "parameters"),

    # ---- differentiation
    ("Differentiation", :diff, (:(x^x), :x, 1), nothing, "x^x"),
    ("Differentiation", :diff, (:(sin(x)^cos(x)), :x, 1), nothing, "variable exponent"),
    ("Differentiation", :diff, (:(sqrt(x + sqrt(x + sqrt(x)))), :x, 1), nothing, "nested radicals"),
    ("Differentiation", :diff, (:(atan(x/sqrt(1 + x^2))), :x, 1), nothing, "inverse trigonometric"),
    ("Differentiation", :diff, (:(exp(x^2)*log(x)), :x, 3), nothing, "third derivative"),
    ("Differentiation", :diff, (:(sin(cos(tan(x)))), :x, 2), nothing, "composition"),
    ("Differentiation", :diff, (:((x^2 + 1)^20), :x, 5), nothing, "polynomial, fifth derivative"),
    ("Differentiation", :diff, (:(x^10*exp(-x)), :x, 10), nothing, "tenth derivative"),
    ("Differentiation", :diff, (:(a*sin(b*x)^2/(1 + x^2)), :x, 2), nothing, "parameters"),

    # ---- equations
    ("Equations", :solve, (:(x^2 - 5x + 6), 0, :x), 2, "quadratic"),
    ("Equations", :solve, (:(x^3 - 6x^2 + 11x - 6), 0, :x), 3, "cubic, rational roots"),
    ("Equations", :solve, (:(x^3 - 3x + 1), 0, :x), 3, "cubic, three real roots (casus irreducibilis)"),
    ("Equations", :solve, (:(x^4 - 10x^2 + 1), 0, :x), 4, "quartic"),
    ("Equations", :solve, (:(x^5 - x + 1), 0, :x), 1, "quintic, not solvable by radicals"),
    ("Equations", :solve, (:(a*x^2 + b*x + c), 0, :x), 2, "quadratic with parameters"),
    ("Equations", :solve, (:(exp(x)), 2, :x), 1, "exponential"),
    ("Equations", :solve, (:(2^x), 8, :x), 1, "exponential, integer base"),
    ("Equations", :solve, (:(x*exp(x)), 1, :x), 1, "Lambert W"),
    ("Equations", :solve, (:(sqrt(x + 1)), :(x - 1), :x), 1, "radical, an extraneous root"),
    ("Equations", :solve, (:(log(x) + log(x - 1)), :(log(6)), :x), 1, "logarithms, an extraneous root"),
    ("Equations", :linsolve, ([(:(x + y), 3), (:(x - y), 1)], (:x, :y)), nothing, "2 × 2"),
    ("Equations", :linsolve, ([(:(2x + y - z), 8), (:(-3x - y + 2z), -11), (:(-2x + y + 2z), -3)], (:x, :y, :z)), nothing, "3 × 3"),
    ("Equations", :linsolve, ([(:(a*x + b*y), 1), (:(x - a*y), :b)], (:x, :y)), nothing, "2 × 2 with parameters"),
]
