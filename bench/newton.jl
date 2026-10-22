# Newton's method on x^5 - x - 1 from x0 = 1: calling expressions directly vs lambdify.
# Run with: julia --project=bench bench/newton.jl
using BenchmarkTools
using SymJ

function newton(f, fp, x0; tol = 1e-14, maxiter = 100)
    x = x0
    for i in 1:maxiter
        dx = f(x) / fp(x)
        x -= dx
        abs(dx) <= tol * max(1, abs(x)) && return x, i
    end
    error("no convergence")
end

@syms x
u = x^5 - x - 1
du = diff(u, x)

fu, fdu = lambdify(u), lambdify(du)
g(x) = x^5 - x - 1
gp(x) = 5x^4 - 1

root, iters = newton(u, du, 1.0)
@assert newton(fu, fdu, 1.0) == (root, iters) == newton(g, gp, 1.0)
println("root = $root after $iters iterations ($(2iters) function evaluations)\n")

r = Ref(1.0)
t_call = median(@benchmark newton($u, $du, $r[])).time
t_lamb = median(@benchmark newton($fu, $fdu, $r[])).time
t_jl = median(@benchmark newton($g, $gp, $r[])).time

# one-time cost: lambdify and the first call (which compiles). Each trial uses a new
# constant, because compiled code is cached by expression.
t_make = map(1:5) do k
    v = x^5 - x - (1 + k // 1000)
    dv = diff(v, x)
    @elapsed begin
        a, b = lambdify(v), lambdify(dv)
        a(1.0), b(1.0)
    end
end |> minimum
t_make *= 1e9

pt(t) = lpad(BenchmarkTools.prettytime(t), 12)
println(rpad("calling the expressions", 36), pt(t_call))
println(rpad("lambdify'd functions", 36), pt(t_lamb), "   ", round(Int, t_call / t_lamb), "x faster")
println(rpad("hand-written Julia", 36), pt(t_jl))
println(rpad("lambdify + first call (compile)", 36), pt(t_make))
println("\nlambdify pays for itself after about ", round(Int, t_make / (t_call - t_lamb)), " solves (",
        round(Int, t_make / ((t_call - t_lamb) / (2iters))), " evaluations)")
