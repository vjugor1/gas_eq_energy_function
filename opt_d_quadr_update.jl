using JuMP
using Ipopt
using LinearAlgebra
using Random
using Plots
using FiniteDiff
include("jac.jl")
include("init_n_bnd_values.jl")

######## GLOBAL PARAMETERS ########
epsilon_t = 1e-3
epsilon_x = 1e-1
step = 1/10.
#c_s = 300
alpha = 8.57
I_ = 10
M = 1000
######## GLOBAL PARAMETERS ########
d_min_p = []
p_min_p = []
ps      = []
ds      = []
Qs      = []
criterions = []
timing    = []
tmp_out = get_init_p_Q(I_, M, left_bound_const, right_bound_const,
                        initial_value_lin, left_bound_Q_lin_inc, 42)
init_value_p_2 = tmp_out[1]
init_value_Q   = tmp_out[2]
global criterion = 10
global iter = 1
global A  = get_A_block(I_)
while  (criterion > -6)# & (iter < 2)
    t1 = time_ns()
    if iter == 1
        d = init_value_p_2
        global p_prev = sqrt.(init_value_p_2)
        global Q_prev = init_value_Q
        push!(ps, p_prev)
        push!(Qs, Q_prev)
        push!(ds, d)
    end
    println("iter = ", iter)
    #println("making gradient step..")
    if iter > 1
        tmp = (d - reshape(p_prev, I_, M).^2)
    else
        tmp = d
    end
    tmp_r = reshape(tmp, I_ * M, 1)

    #x = rand(10)
    #output = zeros(10,10)
    #FiniteDiff.finite_difference_jacobian!(output,f,x)

    #=jac_fin_diff = zeros(I_ * M, I_ * M)
    p_2(dx, x) = f_d!(x, dx, init_value_p_2, init_value_Q,
                                epsilon_t, epsilon_x, false)
    #p_2(dx, x) = solve_opt!(I_, M, init_value_p_2, init_value_Q,
                                epsilon_t, epsilon_x, false, x, dx)
    FiniteDiff.finite_difference_jacobian!(jac_fin_diff, p_2, d)
    jacie = jac_fin_diff - Diagonal(ones(I_ * M, I_ * M))
    grad = jacie * tmp_r
    global d  = d + step * reshape(grad, I_, M)=#

    if iter > 1
        global d  = d - step* (d - reshape(p_prev, I_, M).^2)
    end
    #println("gradient step done")
    p_prev = zeros(I_, M)
    Q_prev = zeros(I_, M)
    solve_out = [p_prev, Q_prev]
    f_d!(d, solve_out, init_value_p_2, init_value_Q,
                            epsilon_t, epsilon_x, true)
    #solve_opt!(I_, M, init_value_p_2, init_value_Q,
    #                epsilon_t, epsilon_x, true, d, solve_out)
    p_prev = solve_out[1]
    Q_prev = solve_out[2]
    global init_value_p_2 = p_prev.^2
    global init_value_Q   = Q_prev
    #println("logging...")
    push!(p_min_p, norm(p_prev - ps[length(ps)]))
    global criterion = log(10, norm(d - p_prev.^2))
    push!(criterions, criterion)
    push!(d_min_p, criterion)
    push!(ps, p_prev)
    push!(Qs, Q_prev)
    push!(ds, d)
    #println("logging...done")
    t2 = time_ns()
    #println("obj = ", objective_value(model))
    println("verbose...")
    #println("|| d - p || = ", d_min_p[length(d_min_p)])
    #println("|| p_prev - p || = ", p_min_p[length(p_min_p)])
    println("criterion = ", criterion)
    println("verbose...done")
    global iter = iter + 1
    push!(timing, (t2 - t1) / 1e9)
end

println("Total elapsed time = ", sum(timing))

plot_stuff(timing, "Iteration number", "Elapsed time",
            true, "timing.png", 1.0, true)
plot_stuff(d_min_p, "Iteration number", "log||d - p||",
            true, "plot_d_min_p.png", 1.0, true)
plot_stuff(p_prev[:, M], "x", "p(x,T)",
            true, "plot_p_x_T.png", epsilon_x, true)
plot_stuff(p_prev[:, 1], "x", "p(x,0)",
            true, "plot_p_x_0.png", epsilon_x, true)

plot_stuff(Q_prev[:, 1], "x", "Q(x,0)", true,
            "plot_Q_x_0.png", epsilon_x, true)

plot_stuff(Q_prev[:, M], "x", "Q(x,T)", true,
            "plot_Q_x_T.png", epsilon_x, true)


plot_stuff(p_prev[1, :], "t", "p(0,t)", true,
            "plot_p_0_t.png", epsilon_t, true)

plot_stuff(p_prev[I_, :], "t", "p(L,t)", true,
            "plot_p_L_t.png", epsilon_t, true)

for i=1:(trunc(Int, M / 100))
    plot_stuff(p_prev[:, (i * 100)], "x", string("p(x,", string(round(i*100*epsilon_t, digits=3)), ")"),
                true, string("plot_Q_x_", string(round(i*100*epsilon_t, digits=3)), ".png"), epsilon_t, true)
    plot_stuff(Q_prev[:, (i * 100)], "x", string("Q(x,", string(round(i*100*epsilon_t, digits=3)), ")"),
                true, string("plot_Q_x_", string(round(i*100*epsilon_t, digits=3)), ".png"), epsilon_t, true)
end
#


#scatter(Array(1:length(d_min_p)), d_min_p)
using JLD
save("Q_prev.jld", "Q_prev", Q_prev)
save("p_prev.jld", "p_prev", p_prev)
save("ds.jld", "ds", ds)
save("ps.jld", "ps", ps)
save("Qs.jld", "Qs", Qs)
save("criterions.jld", "criterions", criterions)
