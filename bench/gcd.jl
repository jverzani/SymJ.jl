# Polynomial gcd: SymJ modular vs SymJ PRS vs Nemo (FLINT).
# Run with: julia --project=bench bench/gcd.jl
using BenchmarkTools
import SymJ
using SymJ: Poly, Sym
import Nemo

const x, y, z, w = SymJ.Sym(:x), SymJ.Sym(:y), SymJ.Sym(:z), SymJ.Sym(:w)
P(e) = Poly(e, x, y, z, w)

const R, (nx, ny, nz, nw) = Nemo.polynomial_ring(Nemo.ZZ, ["x", "y", "z", "w"])
# rebuild a SymJ Poly in Nemo from its terms
function tonemo(p)
    s = zero(R)
    for (ev, c) in zip(SymJ.exponent_vectors(p), SymJ.coefficients(p))
        s += Nemo.ZZ(BigInt(c)) * nx^ev[1] * ny^ev[2] * nz^ev[3] * nw^ev[4]
    end
    s
end

ms(t) = lpad(string(round(t * 1e3; digits = 3)), 10) * " ms"

cases = [
    ("g=(x+y+z+w+1)^10, linear cofactors", (x + y + z + w + 1)^10, x - y + 2z + 1, y * z - x + 3),
    ("g=(x+y+z+w+1)^4, deg-4 cofactors", (x + y + z + w + 1)^4, (x - y + 2z * w + 1)^2 + y, (y * z - x + 3)^2 + w),
    ("g=(x+y+z+w+1)^4, deg-6 cofactors", (x + y + z + w + 1)^4, (x - y + 2z * w + 1)^3 + y, (y * z - x + 3)^3 + w),
    ("g=(x+y+z+w+1)^6, deg-6 cofactors", (x + y + z + w + 1)^6, (x - y + 2z * w + 1)^3 + y, (y * z - x + 3)^3 + w),
    ("g=(x+y+z+w+1)^8, deg-6 cofactors", (x + y + z + w + 1)^8, (x - y + 2z * w + 1)^3 + y, (y * z - x + 3)^3 + w),
    ("big coeffs: g=(3^20 x + 5^15 y z - 7^10 w)^3", (big(3)^20 * x + big(5)^15 * y * z - big(7)^10 * w + 1)^3,
        (x - y + z)^3 + 2, (x * w + y - 1)^2 - z),
]

println(rpad("case", 48), lpad("modular", 14), lpad("PRS", 14), lpad("Nemo", 14))
for (name, g, a, b) in cases
    G = P(g)
    A, B = G * P(a), G * P(b)
    @assert gcd(A, B) == G
    nA, nB = tonemo(A), tonemo(B)
    tm = @belapsed gcd($A, $B)
    tp = @elapsed SymJ._gcd(A, B)
    tp = min(tp, @elapsed SymJ._gcd(A, B))
    tn = @belapsed Nemo.gcd($nA, $nB)
    println(rpad(name, 48), ms(tm), ms(tp), ms(tn))
end
