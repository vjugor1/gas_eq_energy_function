using JuMP
using Ipopt
using LinearAlgebra
using Random
using Plots
include("jac.jl")

######## GLOBAL PARAMETERS ########
epsilon_t = 1e-5
epsilon_x = 1e-2
step = 1/10.
#c_s = 300
alpha = 8.57
I = 70
M = 70
######## GLOBAL PARAMETERS ########

############ FUNCTIONS ############

# see jac.jl





######## SINUSOIDAL BNDS #########

#=
a1 = 0.2
w1 = 100
b1 = 0.3

a2 = 0.2
w2 = 100
b2 = 0.4

function left_bound(x)
    return a1 * sin(w1 *x) + b1
end

function right_bound(x)
    return a2 * sin.(w2 *x) +  b2
end

a = (left_bound(epsilon_t) - right_bound(epsilon_t)) / (epsilon_x - I * epsilon_x)
b = left_bound(epsilon_t) - a * epsilon_x

function initial_value(x)
    return a * x + b
end
=#
######## SINUSOIDAL BNDS #########

a1 = 0.2
w1 = 100
b1 = 0.3

a2 = -0.2
w2 = 100
b2 = 0.4
function left_bound(x)
    return  b1
end

function right_bound(x)
    return b2 + a2 * x
end

a = (left_bound(epsilon_t) - right_bound(epsilon_t)) / (epsilon_x - I * epsilon_x)
b = left_bound(epsilon_t) - a * epsilon_x


function initial_value(x)
    return a * x + b
end





############ FUNCTIONS ############

#d = reshape([[0.0, 0.0, 1.0];[0.0, 0.0, 0.0];[0.0, 0.0, 0.0]], 3, 3)
d_min_p = []
p_min_p = []
ps      = []
ds      = []
Qs      = []
criterions = []
Random.seed!(42)
init_value_p_2 = (rand(I, M)).^2
#=for i=1:I
        for j=1:M
                if j == 1
                    init_value_p_2[i,j] =  i * 1e-1 + 1.0  #sqrt(3) / 3 * i * j + rand()
                elseif i == 1
                    tmp = ones(M)
                    for m=1:M
                        tmp[m] = tmp[m] + (m-1) * 1e-1
                    end
                    init_value_p_2[i,:] = tmp
                elseif i == I
                    tmp = ones(M)
                    for m=1:M
                        tmp[m] = tmp[m] + (m-1) * 1e-1
                    end
                    init_value_p_2[i,:] = tmp * 2
                else
                    init_value_p_2[i,j] = (init_value_p_2[i-1, j] + 0.1 + (j - 2)*1e-3 )
                end
        end
end
init_value_p_2 = init_value_p_2 + ones(I, M) * 1.5
init_value_p_2[100,:] = init_value_p_2[100,:] + 9 * ones(50)
init_value_p_2[1,1] = init_value_p_2[1,1] + 0.1
#init_value_p_2[2:(I-1),1] = [2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3]
#init_value_p_2[3,1] = 1
#init_value_p_2 = [ [1 1 1]; [1.01 1.01 1.01]; [1.02 1.02 1.02] ]
=#

init_value_p_2[1,:] = left_bound.(collect(1:M) * epsilon_t)
init_value_p_2[I,:] = right_bound.(collect(1:M) * epsilon_t)
init_value_p_2[:,1] = initial_value.(collect(1:I) * epsilon_x)

Random.seed!(43)
init_value_Q = rand(I-1,M)# / alpha * epsilon_x
#init_value_Q[:, 1] = ones(I-1) * a / alpha * epsilon_x

#init_value_Q[2, 1] = -10
#init_value_Q = [ [0.1 0.2]; [0.1 0.3]; [0.1 0.4] ]'
#init_value_Q = init_value_p_2[1:2, :]
#p_prev = reshape([[0.0, 0.0, 1.0];[0.0, 0.0, 0.0];[0.0, 0.0, 0.0]], 3, 3)

global criterion = 10
global iter = 1
while  (criterion > -8)# & (iter < 3000)
    if iter == 1
        d = init_value_p_2 #reshape([[1.0, 1.0, 1.0];[1.0, 1.0, 1.0];[1.0, 1.0, 1.0]], 3, 3)
        p_prev = sqrt.(init_value_p_2)
    end
    println("iter = ", iter)
    println("forming optimization problem...")
    model = Model(with_optimizer(Ipopt.Optimizer))
    @variable(model, Q[i=1:(I-1), m=1:M], start=init_value_Q[i, m] )
    @variable(model, p[i=1:I, m=1:M], start=sqrt(init_value_p_2[i,m]))
    A  = get_A_block(I)'
    println("making gradient step..")
    #=tmp = (d - reshape(p_prev, I, M).^2)
    tmp_r = reshape(tmp, I * M, 1)
    jacie = get_jac(reshape(d, 1, I * M)) - Diagonal(ones(I * M, I * M))
    grad = jacie * tmp_r
    global d  = d + step * reshape(grad, I, M)=#
    if iter > 1
        global d  = d - step* (d - reshape(p_prev, I, M).^2)
    end
    println("gradient step done")

    #global d  =  reshape(p_prev, 3, 3).^2
    println("forming constraints...")
    @constraint(model, pdes[i=2:(I-1),m=1:(M-1)], (Q[i,m + 1] - Q[i-1,m + 1])/ ( 2. * epsilon_x) + (Q[i,m] - Q[i-1,m])/ ( 2. * epsilon_x) + (p[i,m+1] - p[i,m]) / epsilon_t == 0.0)



    obj2 = sum((A' * d[:,m])' * Q[:,m] for m=1:M)
    @variable(model, aux)

    @constraint(model, aux == obj2)
    @constraint(model, feeez[idx=2:(I-1)], p[idx,1] == sqrt(init_value_p_2[idx, 1]))
    #@NLconstraint(model, naive[i=1:3, m=1:3], sqrt(d[i,m]) == p[i,m])
    #@constraint(model, press_bnds[i=1:I, m=1:M], -10.0 <= p[i,m])
    println("forming objective...")
    @NLobjective(model, Min, alpha * (epsilon_x / 3.)*sum(sqrt(Q[i,m]^6) for i=1:(I-1), m=1:M) -  aux)
    println("forming objective...done")
    println("optimize...")
    optimize!(model)
    println("optimize...done")
    #global save_Q =  value.(Q[1,2])
    println("logging...")
    push!(p_min_p, norm(value.(p) - p_prev))
    global p_prev = value.(p)
    global Q_prev = value.(Q)
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

    global criterion = log(10, norm(d - p_prev.^2))#log10(abs(sum( abs(d[i+1,m] - p_prev[i+1,m]^2 - d[i,m] + p_prev[i,m]^2) for i=1:(I-1), m=1:M)))
    push!(criterions, criterion)
    push!(d_min_p, log(10, norm(d - p_prev.^2)))
    push!(ps, p_prev)
    push!(Qs, Q_prev)
    push!(ds, d)
    println("logging...done")
    #println("obj = ", objective_value(model))
    println("|| d - p || = ", d_min_p[length(d_min_p)])
    println("|| p_prev - p || = ", p_min_p[length(p_min_p)])
    println("criterion = ", criterion)

    global iter = iter + 1
    #println("dual(feeez) = ", dual(feeez))
end
#length(d_min_p)

#plot(collect(1:(M))*epsilon_t, p_prev[25, :])
plot(Array(1:length(d_min_p)), d_min_p, seriestype=:scatter,
         xlabel = "Iteration Number", ylabel="log||d - p^2||")
savefig("plot_d_min_p.png")
plot(collect(1:I) * epsilon_x, p_prev[:, M],
         xlabel = "x", ylabel="p(x, T)")
savefig("plot_p_x_T.png")
plot(collect(1:M) * epsilon_t, p_prev[10, :],
         xlabel = "t", ylabel="p(half, t)")
savefig("plot_p_half_t.png")
plot(collect(1:(I-1)) * epsilon_x, Q_prev[:, M],
                xlabel = "x", ylabel="Q(x, T)")
savefig("plot_Q_x_T.png")

for i=1:10
    plot(collect(1:I) * epsilon_x, p_prev[:, (i * 10)],
                    xlabel = "x", ylabel=string("p(x,", string(i*10), ")"))
    savefig(string("plot_p_x_", string(i*10*epsilon_t), ".png"))
end


#scatter(Array(1:length(d_min_p)), d_min_p)
