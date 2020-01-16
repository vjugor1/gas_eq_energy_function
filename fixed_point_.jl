using JuMP
using Ipopt
using LinearAlgebra
using Plots
include("jac.jl")
######## GLOBAL PARAMETERS ########
epsilon = 1
I = 3
M = 3
p_1_i = [sqrt(3)/1; sqrt(3)/2; sqrt(3)/3]
p_m_1 = [sqrt(3)/1; sqrt(3)/5; sqrt(3)/6]
p_m_I = [sqrt(3)/3; sqrt(3)/8; sqrt(3)/9]
global init_val = zeros(3, 3)
#=for i=1:3
        for j=1:3
                init_val[i,j] =  sqrt(3) / 3 * (i + j) * 1e-1
        end
end=#
global init_val = [[0.333333   1.33333   3.0];
                   [5.333      5.33354   5.33346];
                   [3.0       12.0      27.0]]
######## GLOBAL PARAMETERS ########

############ PLOT LISTS ###########
obj_list = []
d_min_p = []
p_min_p = []
ps      = []
ds      = []
rads    = []
cntrs   = []
############ PLOT LISTS ###########


############ FUNCTIONS ############
#=function p_of_d(d1d2, d2d3, d1d2_abs, d2d3_abs, i, m, pi1)
        @assert i > 1
        @assert i < I
        @assert m > 1
        @assert m < M
        @assert size(d1d2) == size(d2d3)
        @assert size(d1d2_abs) == size(d2d3_abs)
        @assert size(d1d2) == size(d2d3_abs)
        @assert size(d1d2) == (m-1,)
        out = pi1 + epsilon * sum( √(d2d3[m_]) * (d2d3[m_] / d2d3_abs[m_]) -
                        √(d1d2[m_]) * (d1d2[m_] / d1d2_abs[m_])  for m_=1:m-1)

        return out
end=#
############ FUNCTIONS ############
global criterion = 1488
#for iter = 1:5000
global iter = 0
#while criterion > -5
while iter < 40
        global iter = iter + 1
        if iter ==1
                global d_par = init_val
                global p_prev = init_val
        end
        if iter > 1
                global d_par = d_opt
        end

        model1 = Model(with_optimizer(Ipopt.Optimizer))
        @variable(model1, Q[i=1:2, m=1:3], start=sqrt(sqrt((d_par[i+1,m] - d_par[i,m])^2)) * sign(d_par[i+1,m] - d_par[i,m]))
        @variable(model1, p[i=1:3, m=1:3], start=init_val[i,m])
        A  = transpose(reshape([[-1.0; 0.0]; [1.0; -1.0]; [0.0, 1.0]], 2, 3))

        #global d  = d -  (1/100)*(d - reshape(p_prev, 3, 3).^2)


        #global d  =  reshape(p_prev, 3, 3).^2

        @constraint(model1, pdes[i=2,m=1:2], (Q[i,m] - Q[i-1,m]) + (p[i,m+1] - p[i,m]) == 0.0)


        obj2 = sum((A' * d_par[:,m])' * Q[:,m] for m=1:3)
        @variable(model1, aux)
        push!(ds, d_par)
        @constraint(model1, aux == obj2)
        #@constraint(model, press_bnds[i=1:3, m=1:3], 0.69 <= p[i,m] <= 0.7)
        @NLobjective(model1, Min, (1. / 3.)*sum(sqrt(Q[i,m]^6) for i=1:2, m=1:3) -  aux)

        optimize!(model1)

        push!(p_min_p, norm(value.(p) - p_prev))
        global p_prev = value.(p)
        #global dual_prev = zeros(1,2)
        #dual_prev[1,1] = dual(pdes[2,1])
        #dual_prev[1,2] = dual(pdes[2,2])
        #global dual_prev_pb = zeros(3,3)
        #=for i=1:3
                for j=1:3
                        dual_prev_pb[i,j] = dual(press_bnds[i,j])
                end
        end=#


        #println("dual(pdes[2,1]) = ", dual(pdes[2,1]))
        #println("dual(pdes[2,2]) = ", dual(pdes[2,2]))
        global criterion = log(10, norm(d_par - p_prev.^2))
        push!(d_min_p, log(10, norm(d_par - p_prev.^2)))
        push!(ps, p_prev)

        #println("obj = ", objective_value(model))
        println("|| d - p || = ", d_min_p[length(d_min_p)])
        println("|| p_prev - p || = ", p_min_p[length(p_min_p)])

        global init_val = p_prev





        model = Model(with_optimizer(Ipopt.Optimizer))
        ############ VARIABLES ############
        @variable(model, d[i=1:3, m=1:3], start=init_val[i,m]^2 )
        @variable(model, psi[i=1:2, m=1:3], start=sqrt((init_val[i+1,m]^2 - init_val[i,m]^2)^2) ) # This is |d_i+1^m - d_i^m|

        @variable(model, aux)
        @variable(model, aux2)
        @variable(model, aux3)
        #@variable(model, aux2[1:2, 1:3])
        ############ VARIABLES ############

        ########### CONSTRAINTS ###########
        #@constraint(model, aux2[i,m] == psi[i,m]^2 for i=1:2, m=1:3)
        @constraint(model, abs_dd[i=1:2, m=1:3], psi[i,m]^2 == (d[i+1,m] - d[i,m])^2)
        @constraint(model, psi_non_neg[i=1:2, m=1:3], psi[i, m] >= 1e-7)
        jacie = get_jac(reshape(init_val, 1, 9))
        svds = svd(jacie).S
        l = svds[argmin(svds)]
        L = svds[argmax(svds)]
        #@constraint(model, init_point1[i=1:3, m=1], d[i,m] == init_val[i,m]^2)
        #@constraint(model, init_point2[i=1, m=1:3], d[i,m] == init_val[i,m]^2)
        @NLconstraint(model, ball, (sum( (d[i,m] - init_val[i,m]^2)^2 for i=1:3, m=1:3)) <= l / (2*L) / iter )
        #@constraint(model, d_non_neg[i=1:3, m=1:3], d[i, m] >= 1e-7)
        #@constraint(model, d_init[i=1:3, m=1], d[i,m] == p_1_i[i]^2)
        #@constraint(model, d_left[i=1, m=1:3], d[i,m] == p_m_1[m]^2)
        #@constraint(model, d_rght[i=1, m=1:3], d[i,m] == p_m_I[m]^2)

        #obj1 = sum((d[i,1] - p_1_i[i]^2)^2 for i=2)
        #obj2 = sum((d[i,m] -
        #        (p_1_i[i] + epsilon * sum( sqrt(psi[i,m_]) * (d[i+1, m_] - d[i, m_])/psi[i,m_]  -
        #        sqrt(psi[i-1,m_]) * (d[i, m_] - d[i-1, m_])/psi[i-1,m_] for m_=1:(m-1)))^2)^2 for i=2, m=2:3)


        #@constraint(model, aux == obj1)
        ########### CONSTRAINTS ###########



        @NLconstraint(model, aux == sum( (d[i,1] - init_val[i,1]^2)^2   for i=1:3))
        @NLconstraint(model, aux2 == sum( (d[3,m] - init_val[3,m]^2)^2   for m=1:3))
        @NLconstraint(model, aux3 == sum( (d[1,m] - init_val[1,m]^2)^2   for m=1:3))
        @NLobjective(model, Min,  aux3 + aux2 + aux + sum((d[i,m] -
                (init_val[i,1] -  sum( sqrt(psi[i,m_]) * (d[i+1, m_] - d[i, m_])/psi[i,m_]  -
                sqrt(psi[i-1,m_]) * (d[i, m_] - d[i-1, m_])/psi[i-1,m_] for m_=1:(m-1)))^2)^2 for i=2, m=2:3)
        )


        optimize!(model)

        #push!(p_min_p, norm(value.(p) - p_prev))
        #global p_prev = value.(p)

        #push!(d_min_p, log(10, norm(d - p_prev.^2)))


        println("obj = ", objective_value(model))
        push!(obj_list, objective_value(model))
        println("d = ", value.(d))
        global d_opt = value.(d)
        push!(rads, l / (2*L))
        push!(cntrs, (sum( (d_opt[i,m] - init_val[i,m]^2)^2 for i=1:3, m=1:3)))
        println("psi = ", value.(psi))
        global dual_to_ball = dual(ball)
        println("dual_to_ball = ", dual_to_ball)
        #global init_val = d_opt




end

        plot(Array(1:length(d_min_p)), (d_min_p), seriestype=:scatter,
                 xlabel = "Iteration Number", ylabel="log||d - p^2||")

        #println("|| d - p || = ", d_min_p[length(d_min_p)])
        #println("|| p_prev - p || = ", p_min_p[length(p_min_p)])

        ##############################################################
        ###################Fixed point search ended###################
        ##############################################################

         #=
        model1 = Model(with_optimizer(Ipopt.Optimizer))
        ############ VARIABLES ############
        @variable(model1, Q[1:2, 1:3], start=sqrt(3)/3)
        @variable(model1, p[i=1:3, m=1:3], start=sqrt(3)/3)
        #for i=1:3
        #        for m=1:3
        #                @variable(model1, p[i,m], start=sqrt(d_opt[i,m]))
        #        end
        #end
        ############ VARIABLES ############
        A  = transpose(reshape([[-1.0; 0.0]; [1.0; -1.0]; [0.0, 1.0]], 2, 3))
        #d_opt = reshape([[1.0, 1.0, 1.0];[1.0, 1.0, 1.0];[1.0, 1.0, 1.0]], 3, 3)
        global d  = d_opt
        #global d  =  reshape(p_prev, 3, 3).^2
        ########### CONSTRAINTS ###########
        @constraint(model1, pdes[i=2,m=1:2], (Q[i,m] - Q[i-1,m]) + (p[i,m+1] - p[i,m]) == 0.0)
        #@constraint(model1, p_1_is[i=1:3,m=1], p[i, m] == p_1_i[i])

        @constraint(model1, p_init[i=1:3, m=1], p[i,m] == sqrt(d_opt[i]))
        #@constraint(model1, p_left[i=1, m=1:3], p[i,m] == sqrt(p_m_1[m]))
        #@constraint(model1, p_rght[i=1, m=1:3], p[i,m] == sqrt(p_m_I[m]))


        obj2 = sum((A' * d[:,m])' * Q[:,m] for m=1:3)
        @variable(model1, aux)
        @constraint(model1, aux == obj2)
        #@constraint(model1, press_bnds[i=1:3, m=1:3], sqrt(d[i,m]) - 1<= p[i,m] <= sqrt(d[i,m])+1)

        ########### CONSTRAINTS ###########
        @NLobjective(model1, Min, (1. / 3.)*sum(sqrt(Q[i,m]^6) for i=1:2, m=1:3) -  aux)

        optimize!(model1)

        #push!(p_min_p, norm(value.(p) - p_prev))
        #global p_prev = value.(p)

        #push!(d_min_p, log(10, norm(d - p_prev.^2)))


        println("obj = ", objective_value(model1))
        p_opt   = value.(p)
        Q_opt   = value.(Q)
        p_opt_2 = p_opt.^2
        println("||p^2 - d|| = ", norm(p_opt_2 - d))
        println("p^2 = ", p_opt_2 )
        println("d = ", d )
        #println("|| d - p || = ", d_min_p[length(d_min_p)])
        #println("|| p_prev - p || = ", p_min_p[length(p_min_p)])
        #for i = 1:3
        #        for m = 1:3
        #                println("i=", i)
        #                println("m=", m)
        #                println("dual=", dual(press_bnds[i,m]))
        #                println("_____")
        #        end
        #end
        for i=1:2
                for m=1:3
                        println("check1=", Q_opt[i,m] * abs(Q_opt[i,m]) - (d_opt[i+1,m] - d_opt[i,m]))


                end
        end
        =#
