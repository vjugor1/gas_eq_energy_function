using JuMP
using Ipopt
using LinearAlgebra
using Random
using Plots
using FiniteDiff
using PlotlyJS
include("jac.jl")
include("init_n_bnd_values.jl")
include("Plots_2d.jl")
include("plotly_3d_graphs.jl")

######## GLOBAL PARAMETERS ########
L = 1.0
T = 0.7
I_ = 10#10
M  = 1200 #200
epsilon_t = T / M
#epsilon_x = L / I_

#=epsilon_t = 1e-3
epsilon_x = 1e-1
L = 1.0
T = 1.0
I_ = 12
M = 1002=#

step = 1e-1
#c_s = 300
alpha_ = 8.57# * 10000

######## GLOBAL PARAMETERS ########
d_min_p = []
p_min_p = []
ps      = []
ds      = []
Qs      = []
criterions = []
timing    = []
#belgian
j_list, e_dict = get_clean_belgian_j_e()

e_dict_p_Q = set_stat_p_Q(j_list, e_dict, alpha_)
j_to_erase = []
j_good = ["1", "4", "5", "81", "11", "14", "17", "16"]
for item in j_list
    if !(item in j_good)
        push!(j_to_erase, item)
    end
end
j_list = setdiff(j_list, j_to_erase)
A_inc = get_inc_matrix(j_list, e_dict_p_Q)
solve_scheme_belgian!(I_, M, j_list,
                                e_dict_p_Q, A_inc, d_min_p, p_min_p, ps, ds,
                                Qs, criterions, timing, f_d_belg!, -7)

#3 pipes
#solve_scheme_3_pipes!(I_, M, d_min_p, p_min_p, ps, ds, Qs, criterions, timing, f_d_3_pipes!, -5.5)

#1 pipe
#solve_scheme!(I_, M, d_min_p, p_min_p, ps, ds, Qs, criterions, timing, f_d!, -6)
#=p_low = ps[length(ps)]
Q_low = Qs[length(Qs)]
eps_x_new = 1e-1 / 2
eps_t_new = 1e-3
step = 1/100.
p_high_sq, Q_high, I_, M = get_higher_resolution_init_vals(L, T, p_low.^2, Q_low, epsilon_x, epsilon_t, eps_x_new, eps_t_new)
solve_scheme!(I_, M, d_min_p, p_min_p, ps, ds, Qs, criterions, timing, f_d!, -6, [p_high_sq, Q_high])=#
println("Total elapsed time = ", sum(timing))
p_prev = ps[length(ps)]
Q_prev = Qs[length(Qs)]
plot_stuff(timing, "Iteration number", "Elapsed time",
            true, "timing.png", 1.0, true)
plot_stuff(d_min_p, "Iteration number", "log||d - p||",
            true, "plot_d_min_p.png", 1.0, true)
plot_stuff(p_prev[:, M-1], "x", "p(x,T)",
            true, "plot_p_x_T.png", epsilon_x, true)
plot_stuff(p_prev[:, 1], "x", "p(x,0)",
            true, "plot_p_x_0.png", epsilon_x, true)

plot_stuff(Q_prev[:, 1], "x", "Q(x,0)", true,
            "plot_Q_x_0.png", epsilon_x, true)

plot_stuff(Q_prev[:, M-1], "x", "Q(x,T)", true,
            "plot_Q_x_T.png", epsilon_x, true)


plot_stuff(p_prev[1, :], "t", "p(0,t)", true,
            "plot_p_0_t.png", epsilon_t, true)

plot_stuff(p_prev[I_, :], "t", "p(L,t)", true,
            "plot_p_L_t.png", epsilon_t, true)

for i=1:(trunc(Int, M / 100))
    plot_stuff(p_prev[:, (i * 100)], "x", string("p(x,", string(round(i*100*epsilon_t, digits=3)), ")"),
                true, string("plot_p_x_", string(round(i*100*epsilon_t, digits=3)), ".png"), epsilon_t, true)
    plot_stuff(Q_prev[:, (i * 100)], "x", string("Q(x,", string(round(i*100*epsilon_t, digits=3)), ")"),
                true, string("plot_Q_x_", string(round(i*100*epsilon_t, digits=3)), ".png"), epsilon_t, true)
end
#


#scatter(Array(1:length(d_min_p)), d_min_p)
#=using JLD
save("Q_prev.jld", "Q_prev", Q_prev)
save("p_prev.jld", "p_prev", p_prev)
save("ds.jld", "ds", ds)
save("ps.jld", "ps", ps)
save("Qs.jld", "Qs", Qs)
save("criterions.jld", "criterions", criterions)=#
