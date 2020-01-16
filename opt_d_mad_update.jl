using JuMP
using Ipopt
using LinearAlgebra
using Plots
include("jac.jl")

######## GLOBAL PARAMETERS ########
epsilon = 1e-2
a = zeros(3, 3)
for i = 1:3
    for j = 1:3
        a[i,j] = sqrt(3) / 3 * i * j
    end
end
######## GLOBAL PARAMETERS ########

############ FUNCTIONS ############

foo = p_vec

############ FUNCTIONS ############

model = Model(with_optimizer(Ipopt.Optimizer))
############ VARIABLES ############
@variable(model, d[i=1:3, m=1:3], start=sqrt(3)/3 * i * m)
@variable(model, p_2[i=1:3, m=1:3], start=sqrt(3)/3 * i * m)
############ VARIABLES ############

#global d  =  reshape(p_prev, 3, 3).^2
########### CONSTRAINTS ###########
#aux = foo(reshape(d, 1, 9))
@NLexpression(model, foo_expr[m=1:2], foo(reshape(d, 1, 9))[m])
@NLconstraint(model, im_f[i=2,m=1:2], p_2[i, m] == foo_expr[m])
@NLconstraint(model, ball, sum( (d[i,m] - a[i,m])^2 for i=1:3, m=1:3) <= epsilon)
########### CONSTRAINTS ###########



#@constraint(model, press_bnds[i=1:3, m=1:3], 0.69 <= p[i,m] <= 0.7)
@NLobjective(model, Min, sum( (d[i,m] - p_2[i,m])^2 for i=1:3, m=1:3 ))

optimize!(model)
#length(d_min_p)
#plot(Array(1:length(d_min_p)), d_min_p, seriestype=:scatter,
#         xlabel = "Iteration Number", ylabel="log||d - p^2||")
#savefig("plot_success_log_5000.png")
#scatter(Array(1:length(d_min_p)), d_min_p)
