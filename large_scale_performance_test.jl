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
include("gen_simulated.jl")
######## GLOBAL PARAMETERS ########
L = 1.0
T = 0.3
I_ = 5#10
M  = 500 #200
epsilon_t = T / M
epsilon_x = L / I_

#=epsilon_t = 1e-3
epsilon_x = 1e-1
L = 1.0
T = 1.0
I_ = 12
M = 1002=#

step = 1e-1
#c_s = 300
alpha_ = 8.57# * 10000
depth_list = [3, 4, 5, 6, 7, 8, 9, 10]
timing_newty = []
timing_tes   = []
for dee in depth_list
    ######## GLOBAL PARAMETERS ########
    d_min_p = []
    p_min_p = []
    ps      = []
    ds      = []
    Qs      = []
    criterions = []
    timing    = []

    ################################################################simulated
    j_list, e_dict_p_Q = get_e_dict(dee)
    A_inc = get_inc_matrix(j_list, e_dict_p_Q)

    #simulated
    t1 = time_ns()
    solve_scheme_belgian!(I_, M, j_list,
                                    e_dict_p_Q, A_inc, d_min_p, p_min_p, ps, ds,
                                    Qs, criterions, timing, f_d_simul!, -7.230)
    t2 = time_ns()
    println((t2 - t1) / 1e9)
    push!(timing_tes, (t2 - t1) / 1e9)



    ######## GLOBAL PARAMETERS ########
    d_min_p = []
    p_min_p = []
    ps      = []
    ds      = []
    Qs      = []
    criterions = []
    timing    = []

    ################################################################simulated
    j_list, e_dict_p_Q = get_e_dict(dee)
    A_inc = get_inc_matrix(j_list, e_dict_p_Q)

    ################################################################Newton
    t1 = time_ns()
    time_newty = solve_opt_newty(I_, M, j_list,
                                    e_dict_p_Q, A_inc, d_min_p, p_min_p, ps, ds,
                                    Qs, criterions, timing, f_d_simul!, -3.230, true)
    t2 = time_ns()
    println((t2 - t1) / 1e9)
    push!(timing_newty, (t2 - t1) / 1e9)
end

println("newty: ", round.(timing_newty, digits=2))
println("t-ef: ", round.(timing_tes, digits=2))
