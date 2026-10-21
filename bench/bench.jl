# Compare SymJ against SymEngine on a few standard workloads.
# Run with: julia --project=bench bench/bench.jl
using BenchmarkTools
using SymJ: Sym
import SymJ
import SymEngine

const SJ = (x = SymJ.Sym(:x), y = SymJ.Sym(:y), z = SymJ.Sym(:z), w = SymJ.Sym(:w))
const SE = (x = SymEngine.symbols(:x), y = SymEngine.symbols(:y), z = SymEngine.symbols(:z), w = SymEngine.symbols(:w))

# build sum_{i=1}^n i*x^i*y^(n-i), one term at a time
function build_sum(s, n)
    acc = 0 * s.x
    for i in 1:n
        acc += i * s.x^i * s.y^(n - i)
    end
    acc
end

# repeated multiplication of distinct symbols-ish factors
function build_prod(s, n)
    acc = 1 + 0 * s.x
    for i in 1:n
        acc *= (s.x + i)^i
    end
    acc
end

expand1(s, E) = E.expand((s.x + s.y + s.z + 1)^15)
function expand2(s, E)
    e = E.expand((s.x + s.y + s.z + s.w)^15)
    E.expand(e * (e + s.w))
end
diffsin(s, E) = E.diff(E.expand((s.x + s.y + 1)^8) * sin(s.x * s.y), s.x)
subs1(s, e, E) = E.subs(e, s.x => s.y + 1)

const SJmod = (expand = SymJ.expand, diff = SymJ.diff, subs = SymJ.subs)
const SEmod = (expand = SymEngine.expand, diff = SymEngine.diff, subs = SymEngine.subs)

function report(name, tj, te)
    println(rpad(name, 22), lpad(BenchmarkTools.prettytime(tj), 12), lpad(BenchmarkTools.prettytime(te), 12),
            lpad(string(round(tj / te; digits = 2)), 8))
end

println(rpad("benchmark", 22), lpad("SymJ", 12), lpad("SymEngine", 12), lpad("ratio", 8))
report("x+y (cons)", median(@benchmark $(SJ.x) + $(SJ.y)).time, median(@benchmark $(SE.x) + $(SE.y)).time)
report("build_sum(200)", median(@benchmark build_sum(SJ, 200)).time, median(@benchmark build_sum(SE, 200)).time)
report("build_prod(30)", median(@benchmark build_prod(SJ, 30)).time, median(@benchmark build_prod(SE, 30)).time)
report("expand (x+y+z+1)^15", median(@benchmark expand1(SJ, SJmod)).time, median(@benchmark expand1(SE, SEmod)).time)
report("diff", median(@benchmark diffsin(SJ, SJmod)).time, median(@benchmark diffsin(SE, SEmod)).time)
let ej = SymJ.expand((SJ.x + SJ.y + SJ.z + 1)^10), ee = SymEngine.expand((SE.x + SE.y + SE.z + 1)^10)
    report("subs", median(@benchmark subs1(SJ, $ej, SJmod)).time, median(@benchmark subs1(SE, $ee, SEmod)).time)
end
report("expand2 (big)", (@elapsed expand2(SJ, SJmod)) * 1e9, (@elapsed expand2(SE, SEmod)) * 1e9)
report("expand2 (big)", (@elapsed expand2(SJ, SJmod)) * 1e9, (@elapsed expand2(SE, SEmod)) * 1e9)

# raw polynomial arithmetic, no conversion to/from expressions
let p = SymJ.Poly((SJ.x + SJ.y + SJ.z + SJ.w)^15, SJ.x, SJ.y, SJ.z, SJ.w)
    q = p + SymJ.Poly(SJ.w, SJ.x, SJ.y, SJ.z, SJ.w)
    println(rpad("Poly p*(p+w) only", 22), lpad(BenchmarkTools.prettytime(median(@benchmark $p * $q).time), 12))
    println(rpad("Poly → Sym (6k terms)", 22), lpad(BenchmarkTools.prettytime(median(@benchmark Sym($p * $q)).time), 12))
    println(rpad("Sym → Poly (6k terms)", 22), lpad(BenchmarkTools.prettytime(median(@benchmark SymJ.Poly($(SymJ.Sym(p * q)), SJ.x, SJ.y, SJ.z, SJ.w)).time), 12))
end
