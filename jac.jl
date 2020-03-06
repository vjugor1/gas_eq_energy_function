using LinearAlgebra
#using ForwardDiff
using JuMP
using Ipopt
using Random
using Plots
include("init_n_bnd_values.jl")

######## GLOBAL PARAMETERS ########
#=epsilon_t = 1e-3
epsilon_x = 1e-1
step = 1/10.
#c_s = 300
alpha_ = 8.57
I_ = 10
M = 1000
p_1 = zeros(I_)
for i=1:I_
        p_1[i] = i
end=#
######## GLOBAL PARAMETERS ########

############ FUNCTIONS ############


#opt_solve_method!(d, solve_out, init_value_p_2, init_value_Q,
#                        epsilon_t, epsilon_x, true)
function solve_opt!(d, out, init_value_p_2, init_value_Q, epsilon_t, epsilon_x, p_Q_out)
    #In-place type
    model = Model(with_optimizer(Ipopt.Optimizer, print_level=0))
    @variable(model, Q[i=1:I_, m=1:M], start=init_value_Q[i, m] )
    @variable(model, p[i=1:I_, m=1:M], start=sqrt(init_value_p_2[i,m]))
    if size(d) != (I_, M)
        d = reshape(d, I_, M)
    end
    #if iter > 1
    #    global d  = d - step* (d - reshape(p_prev, I_, M).^2)
    #end
    #println("gradient step done")
    #global d  =  reshape(p_prev, 3, 3).^2
    #println("forming constraints...")

    obj2 = sum((A' * d[:,m])' * Q[2:I_,m] for m=1:M)
    @variable(model, aux)
    @constraint(model, aux == obj2)
    @constraint(model, feeez[idx=1:(I_)], p[idx,1] == sqrt(init_value_p_2[idx, 1]))
    @constraint(model, feeez_Q[idx=1:(M)], Q[1,idx] == init_value_Q[1, idx])
    @constraint(model, pdes[i=1:I_-1,m=1:(M-1)], (Q[i+1, m] - Q[i, m]) / epsilon_x + (p[i,m+1] - p[i,m]) / epsilon_t == 0.0)
    #println("forming objective...")
    @NLobjective(model, Min, alpha_ * (epsilon_x / 3.)*sum(sqrt(Q[i,m]^6) for i=1:I_, m=1:M) -  aux)
    #println("forming objective...done")
    #println("optimize...")
    optimize!(model)
    #println("optimize...done")
    p_out = value.(p)
    Q_out = value.(Q)
    if p_Q_out == true
        out[1] = p_out
        out[2] = Q_out
        nothing
    else
        out[:] = p_out.^2
        nothing
    end
end

function solve_opt(I_, M, init_value_p_2, init_value_Q, epsilon_t, epsilon_x, p_Q_out, d)
    #Out-of-place type
    model = Model(with_optimizer(Ipopt.Optimizer, print_level=0))
    @variable(model, Q[i=1:I_, m=1:M], start=init_value_Q[i, m] )
    @variable(model, p[i=1:I_, m=1:M], start=sqrt(init_value_p_2[i,m]))
    if size(d) != (I_, M)
        d = reshape(d, I_, M)
    end
    #if iter > 1
    #    global d  = d - step* (d - reshape(p_prev, I_, M).^2)
    #end
    #println("gradient step done")
    #global d  =  reshape(p_prev, 3, 3).^2
    #println("forming constraints...")

    obj2 = sum((A' * d[:,m])' * Q[2:I_,m] for m=1:M)
    @variable(model, aux)
    @constraint(model, aux == obj2)
    @constraint(model, feeez[idx=1:(I_)], p[idx,1] == sqrt(init_value_p_2[idx, 1]))
    @constraint(model, feeez_Q[idx=1:(M)], Q[1,idx] == init_value_Q[1, idx])
    @constraint(model, pdes[i=2:I_,m=1:(M-1)], (Q[i, m] - Q[i-1, m]) / epsilon_x + (p[i,m+1] - p[i,m]) / epsilon_t == 0.0)
    #println("forming objective...")
    @NLobjective(model, Min, alpha_ * (epsilon_x / 3.)*sum(sqrt(Q[i,m]^6) for i=1:I_, m=1:M) -  aux)
    #println("forming objective...done")
    #println("optimize...")
    optimize!(model)
    #println("optimize...done")
    p_out = value.(p)
    Q_out = value.(Q)
    if p_Q_out == true
        return (p_out, Q_out)
    else
        return p_out.^2
    end
end

function f_d_bwd_eul!(d, out, init_value_p_2, init_value_Q, epsilon_t, epsilon_x, p_Q_out)
    #out of place
    p = sqrt.(copy(init_value_p_2))
    Q = copy(init_value_Q)
    I_, M = size(init_value_p_2)
    for i=2:(I_)
        for m=1:M
            Q[i,m] = - sqrt( abs((d[i,m] - d[i-1,m])) / (epsilon_x * alpha_) ) * ( abs((d[i,m] - d[i-1,m])) / ((d[i,m] - d[i-1,m])))
        end
    end

    for i=2:(I_)
        for m=2:M
            p[i,m] = p[i,1] - (epsilon_t) / (epsilon_x) * sum( (Q[i,m_] - Q[i-1,m_]) for m_=2:(m))
        end
    end
    if p_Q_out == false
        out[:] = p.^2
        nothing
    else
        out[1] = p
        out[2] = Q
        nothing
    end
end


function f_d!(d, out, init_value_p_2, init_value_Q, epsilon_t, epsilon_x, p_Q_out)
    #out of place
    p = sqrt.(copy(init_value_p_2))
    Q = copy(init_value_Q)
    for i=2:(I_)
        for m=1:M
            Q[i,m] = - sqrt( abs(d[i,m] - d[i-1,m]) / (epsilon_x * alpha_) ) * sign(d[i,m] - d[i-1,m])
        end
    end

    for i=1:(I_-1)
        for m=2:M
            p[i,m] = p[i,1] - (epsilon_t) / (epsilon_x) * sum( Q[i+1,m_] - Q[i,m_] for m_=1:(m-1))
        end
    end
    if p_Q_out == false
        out[:] = p.^2
        nothing
    else
        out[1] = p
        out[2] = Q
        nothing
    end
end

function f_d_3_pipes!(d, out, init_value_p_2, init_value_Q, epsilon_t, epsilon_x, p_Q_out)
    #out of place
    # (3 I_ - 2) x M
    p = sqrt.(copy(init_value_p_2))
    Q = copy(init_value_Q)

    for m=1:M
        d[2 * I_, m] = d[2 * I_ + 1, m]
        d[I_, m] = d[2 * I_ + 1, m]
        p[2 * I_, m] = p[2 * I_ + 1, m]
        p[I_, m] = p[2 * I_ + 1, m]
        #=d[2 * I_, m] = d[I_, m]
        d[2 * I_ + 1, m] = d[I_, m]
        p[2 * I_, m] = p[I_, m]
        p[2 * I_ + 1, m] = p[I_, m]=#
    end

    p[1 + 2 * I_, 1] = p[I_, 1]


    for m=1:M
        for i=2:(I_)
            Q[i,m] = - sqrt( abs(d[i,m] - d[i-1,m]) / (epsilon_x * alpha_) ) * sign(d[i,m] - d[i-1,m])
        end
    end

    for m=1:M
        for i=(2 + I_):(2 * I_)
            #if i < 2 * I_
                Q[i,m] = - sqrt( abs(d[i,m] - d[i-1,m]) / (epsilon_x * alpha_) ) * sign(d[i,m] - d[i-1,m])
            #=elseif i == 2 * I_
                Q[i,m] = - sqrt( abs(d[I_,m] - d[i-1,m]) / (epsilon_x * alpha_) ) * sign(d[I_,m] - d[i-1,m])
            end=#
        end
    end

    for m=1:M
        Q[1 + 2 * I_, m] = Q[I_, m] + Q[2 * I_, m]
    end
    for m=1:M
        for i=(2 + 2 * I_):(3 * I_)
            #if i == 2 + 2 * I_
                Q[i,m] = - sqrt( abs(d[i,m] - d[i-1,m]) / (epsilon_x * alpha_) ) * sign(d[i,m] - d[i-1,m])
            #=else
                Q[i,m] = - sqrt( abs(d[i,m] - d[i-1,m]) / (epsilon_x * alpha_) ) * sign(d[i,m] - d[i-1,m])
            end=#
        end
    end

    for m=2:M
        for i=1:(I_-1)
            p[i,m] = p[i,1] - (epsilon_t) / (epsilon_x) * sum( Q[i+1,m_] - Q[i,m_] for m_=1:(m-1))
        end
    end
    for m=2:M
        for i=(1 + I_):(2 * I_ - 1)
            p[i,m] = p[i,1] - (epsilon_t) / (epsilon_x) * sum( Q[i+1,m_] - Q[i,m_] for m_=1:(m-1))
        end
    end
    for m=2:M
        for i=(1 + 2 * I_):(3 * I_ - 1)
            #if i == 1 + 2 * I_
                p[i,m] = p[i,1] - (epsilon_t) / (epsilon_x) * sum( Q[i+1,m_] - Q[i,m_] for m_=1:(m-1))
            #=else
                p[i,m] = p[i,1] - (epsilon_t) / (epsilon_x) * sum( Q[i+1,m_] - Q[i,m_] for m_=1:(m-1))
            end=#
        end
    end

    if p_Q_out == false
        out[:] = p.^2
        nothing
    else
        out[1] = p
        out[2] = Q
        nothing
    end
end

function f_d(d, init_value_p_2, init_value_Q, epsilon_t, epsilon_x, p_Q_out)
    #in place
    p = sqrt.(copy(init_value_p_2))
    Q = copy(init_value_Q)
    for i=2:I_
        for m=1:M
            Q[i,m] = sqrt( abs(d[i,m] - d[i-1,m]) / (epsilon_x * alpha_) ) * ((d[i,m] - d[i-1,m]) / abs(d[i,m] - d[i-1,m]))
        end
    end

    for i=2:I_
        for m=2:M
            p[i,m] = p[i,1] - (epsilon_t) / (epsilon_x) * sum( Q[i,m_] - Q[i-1,m_] for m_=1:(m-1))
        end
    end
    if p_Q_out == false
        return p.^2
    else
        return (p, Q)
    end
end

function f_d_e(d, init_value_p_2, init_value_Q, epsilon_t, epsilon_x, p_Q_out)
    #in place
    #p = similar(d)
    d  = reshape(d, I_, M)
    p  = Array{Any, 2}(undef, I_, M)
    #println("kkk")
    #p[1,1] = d[1,1] - init_value_p_2[1,1]
    #println("kkk")
    #initial condition p
    for i=1:I_
        p[i,1] = d[i,1] - d[i,1] + sqrt(init_value_p_2[i,1])
    end

    #left bound p

    for m=1:M
        p[1,m] = d[1,m] - d[1,m] + sqrt(init_value_p_2[1,m])
    end
    #with initial Q
    for m=2:M
        p[2,m] = sqrt(init_value_p_2[2,1]) - (epsilon_t) / (epsilon_x) * sum( sqrt( abs(d[2,m_] - d[1,m_]) / (epsilon_x * alpha_) ) * ((d[2,m_] - d[1,m_]) / abs(d[2,m_] - d[1,m_])) - init_value_Q[1,m_] for m_=1:(m-1))
    end

    for i=3:I_
        for m=2:M
            #Q  = sqrt( abs(d[i,m_] - d[i-1,m_]) / (epsilon_x * alpha_) ) * ((d[i,m_] - d[i-1,m_]) / abs(d[i,m_] - d[i-1,m_]))
            #Q_ = sqrt( abs(d[i-1,m_] - d[i-2,m_]) / (epsilon_x * alpha_) ) * ((d[i-1,m_] - d[i-2,m_]) / abs(d[i-1,m_] - d[i-2,m_]))
            p[i,m] = sqrt(init_value_p_2[i,1]) - (epsilon_t) / (epsilon_x) * sum( sqrt( abs(d[i,m_] - d[i-1,m_]) / (epsilon_x * alpha_) ) * ((d[i,m_] - d[i-1,m_]) / abs(d[i,m_] - d[i-1,m_])) - sqrt( abs(d[i-1,m_] - d[i-2,m_]) / (epsilon_x * alpha_) ) * ((d[i-1,m_] - d[i-2,m_]) / abs(d[i-1,m_] - d[i-2,m_])) for m_=1:(m-1))
        end
    end
    if p_Q_out == false
        return sum((d - p).^2)
    else
        return (p, Q)
    end
end


#wrapper function for solving scheme

function solve_scheme!(I_, M, d_min_p, p_min_p, ps, ds, Qs, criterions, timing, opt_solve_method!, crit, init_values=-112)
    if init_values == -112
        tmp_out = get_init_p_Q(I_, M, right_bound_const, left_bound_const,
                            initial_value_lin, left_bound_Q_const, right_bound_Q_const, 42)
        global init_value_p_2 = tmp_out[1]
        global init_value_Q   = tmp_out[2]
    else
        global init_value_p_2 = init_values[1]
        global init_value_Q   = init_values[2]
    end
    global criterion = 10
    global iter = 1
    global A  = get_A_block(I_)
    while  (criterion > crit)# & (iter < 2)
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
        #=
        if iter > 1
            tmp = (d - reshape(p_prev, I_, M).^2)
        else
            tmp = d
        end
        tmp_r = reshape(tmp, I_ * M, 1)
        jac_fin_diff = zeros(I_ * M, I_ * M)
        p_2(dx, x) = f_d_backward_euler!(x, dx, init_value_p_2, init_value_Q,
                                    epsilon_t, epsilon_x, false)
        #p_2(dx, x) = solve_opt!(I_, M, init_value_p_2, init_value_Q,
        #                            epsilon_t, epsilon_x, false, x, dx)
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
        opt_solve_method!(d, solve_out, init_value_p_2, init_value_Q,
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
    nothing
end

function solve_scheme_3_pipes!(I_, M, d_min_p, p_min_p, ps, ds, Qs, criterions, timing, opt_solve_method!, crit, init_values=-112)
    if init_values == -112
        tmp_out = get_init_p_Q(I_, M, right_bound_const, left_bound_const,
                            initial_value_lin, left_bound_Q_const, right_bound_Q_const, 42)
        tmp_out_third_pipe = get_init_p_Q(I_, M, right_bound_const_3rd, left_bound_const_3rd,
                            initial_value_lin, left_bound_Q_const_3rd, right_bound_Q_const_3rd, 42)
        global init_value_p_2 = [tmp_out[1]; tmp_out[1]; tmp_out_third_pipe[1]]
        global init_value_Q   = [tmp_out[2]; tmp_out[2]; tmp_out_third_pipe[2]]
    else
        global init_value_p_2 = init_values[1]
        global init_value_Q   = init_values[2]
    end
    global criterion = 10
    global iter = 1
    while  (criterion > crit)# & (iter < 2)
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
        #=
        if iter > 1
            tmp = (d - reshape(p_prev, I_, M).^2)
        else
            tmp = d
        end
        tmp_r = reshape(tmp, I_ * M, 1)
        jac_fin_diff = zeros(I_ * M, I_ * M)
        p_2(dx, x) = f_d_backward_euler!(x, dx, init_value_p_2, init_value_Q,
                                    epsilon_t, epsilon_x, false)
        #p_2(dx, x) = solve_opt!(I_, M, init_value_p_2, init_value_Q,
        #                            epsilon_t, epsilon_x, false, x, dx)
        FiniteDiff.finite_difference_jacobian!(jac_fin_diff, p_2, d)
        jacie = jac_fin_diff - Diagonal(ones(I_ * M, I_ * M))
        grad = jacie * tmp_r
        global d  = d + step * reshape(grad, I_, M)=#


        if iter > 1
            global d  = d - step* (d - reshape(p_prev, 3 * I_, M).^2)
        end


        #println("gradient step done")
        p_prev = zeros(3 * I_, M)
        Q_prev = zeros(3 * I_, M)
        solve_out = [p_prev, Q_prev]
        opt_solve_method!(d, solve_out, init_value_p_2, init_value_Q,
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
    nothing
end

# Direct hands solve of pdes


function solve_pde_hands(init_value_p_2, init_value_Q,
                            I_, M, epsilon_x, epsilon_t)
    p = sqrt.(abs.(init_value_p_2))
    Q = init_value_Q
    for m=1:(M-1)
        for i=2:I_
            Q[i,m] = sqrt(abs(p[i,m]^2 - p[i-1,m]^2) / (epsilon_x * alpha_) ) *
                        sign((p[i,m]^2 - p[i-1,m]^2))
            p[i,m+1] = p[i,m] - epsilon_t * (Q[i,m] - Q[i-1,m]) / epsilon_x
        end
    end
    for i=2:I_
        Q[i,M] = sqrt(abs(p[i,M]^2 - p[i-1,M]^2) / (epsilon_x * alpha_) ) *
                    sign((p[i,M]^2 - p[i-1,M]^2))
    end
    return (p, Q)
end


function p_of_d(d1d2, d2d3, d1d2_abs, d2d3_abs, i, m, p_1)
        @assert i > 1
        @assert i < I_
        @assert m > 1
        @assert m <= M
        @assert size(d1d2) == size(d2d3)
        @assert size(d1d2_abs) == size(d2d3_abs)
        @assert size(d1d2) == size(d2d3_abs)
        @assert size(d1d2) == (m-1,)
        out = p_1 + epsilon / alpha_ * sum( √(d2d3_abs[m_]) * sign(d2d3[m_]) -
                        √(d1d2_abs[m_]) * sign(d1d2[m_])   for m_=1:m-1)

        return out
end
function p_vec(d)
        #out = Array{Any, 1}(undef, I_ * M) #zeros(1, (I_ - 2) * (M - 1))
        out = Vector{typeof(zero(d[1])/oneunit(d[1]))}(undef, I_ * M)
        for i=2:(I_-1)
                for m=2:M
                        #if (i in 2:(I_-1)) & (m in 2:M)
                                d1 = d[(M * (i - 2) + 1):(M * (i - 2) + m-1)]
                                d2 = d[(M * (i - 1) + 1):(M * (i - 1) + m-1)]
                                #println("i = ", i)
                                #println("m = ", m)
                                #println("(I_ * (i - 0) + 1) = ", (M * (i - 0) + 1))
                                #println("(I_ * (i - 0) + m-1) = ", (M * (i - 0) + m-1))
                                d3 = d[(M * (i - 0) + 1):(M * (i - 0) + m-1)]
                                d1d2 = d2 - d1
                                d2d3 = d3 - d2
                                d1d2_abs = sqrt.(d1d2.^2)
                                d2d3_abs = sqrt.(d2d3.^2)
                                p_1i = p_1[i]
                                #println("m = ", m)
                                #println("size(d1d2) = ", size(d1d2))
                                out[M * (i - 1) + m] = p_of_d(d1d2, d2d3, d1d2_abs, d2d3_abs, i, m, p_1i)
                        #else
                                #tmp = d[I_ * (i - 1) + m]
                                #out[M * (i - 1) + m] = 0.
                        #end
                end
        end
        #println("out=", size(out))
        return out.^2
end

function get_jac(d)
        foo = p_vec
        return ForwardDiff.jacobian(foo, 10000. * d) / 10000.
end


function gen_fin_diff(n, number)
    @assert number <= n-1
    tmp = zeros(1, n)
    tmp[number]     = -1.0
    tmp[number + 1] = +1.0
    return tmp
end

function get_A_block(I_)
    ### In a A_block ###
    # A ∈ I_-1 x I
    ### In a A_block ###
    A_block = reshape(gen_fin_diff(I_, 1), I_, 1)
    for i=2:(I_-1)
        A_block = [A_block  reshape(gen_fin_diff(I_, i), I_, 1)]
    end
    return A_block
end

function get_A_proj(I_, M)
    ### In total A_proj ###
    # (I_-1) * M -- num of rows
    # I_ * M -- num of cols
    ### In total A_proj ###

    ### In a A_block ###
    # A ∈ I_-1 x I
    ### In a A_block ###

    A_proj = zeros((I_-1) * M, I_ * M)
    for i=0:(M-1)
        A_block = get_A_block(I)
        #print(i)
        #print(i*(I_-1)+1)
        #print((i+1)*(I_-1))
        #print("___")
        A_proj[i*(I_-1)+1:(i+1)*(I_-1), 1+I_*i:I_*(i+1)] = A_block
    end

    return A_proj

end




function Z_c(d1, d2, d3, i, m, p_1)

end

function Z_r(d1, d2, d3, i, m, p_1)


end

function Z_l(d1, d2, d3, i, m, p_1)

end

function jac_string(d, i, m)
        @assert i > 1
        @assert i < I_
        out = zeros(I_ * M)
        p_vec_value = p_vec(d)
        for i_=3:(i)
                for m_=2:(m-1)
                        diff_d_l = (d[M * (i - 1 - 2) + m_] - d[M * (i - 1 - 1) + m_])
                        diff_d_r = (d[M * (i - 1) + m_] - d[M * (i - 1 - 1) + m_])
                        out[M * (i_ - 1 - 2) + m_] = 2 * p_vec_value[M * (i-1) + m] * (- epsilon / 2 * (sign(diff_d_l) / sqrt(abs(diff_d_l))))
                        #=if (sqrt(abs(diff_d_l)) == 0) | (sqrt(abs(diff_d_r)) == 0)
                                println("i_ = ", i_)
                                println("m_ = ", m_)
                        end=#
                        #if (isnan(p_vec_value[M * (i-1) + m]))
                        #        println("i_ = ", i_)
                        #        println("m_ = ", m_)
                        #end
                        #println()

                        out[M * (i_ - 1  -1) + m_] = 2*  p_vec_value[M * (i-1) + m] * ( epsilon / 2 * (sign(diff_d_l) / sqrt(abs(diff_d_l))) - epsilon / 2 * (sign(diff_d_r) / sqrt(abs(diff_d_r))))
                        out[M * (i_ - 1) + m_] = 2 * p_vec_value[M * (i-1) + m] * ( epsilon / 2 * (sign(diff_d_r) / sqrt(abs(diff_d_r))))


                end
        end
        return out
end

function jac_hands(d)
        out = zeros(I_ * M, I_ * M)
        for i=2:(I_-1)
                for m=1:M
                        #println("i = ", i)
                        #println("m = ", m)
                        out[M * (i-1) + m, :] = jac_string(d, i, m)
                end
        end
        return out
end

#function jac_p_sq

#end
############ FUNCTIONS ############
#=
x = collect(1:50)
x = x * epsilon_t
y = 0.2 * sin.(100 *x) + ones(50) * 0.3


plot(x, y, seriestype=:scatter,
         xlabel = "x", ylabel="y")
=#



#sum(abs((p_prev[i, m+1] - p_prev[i,m]) / epsilon_t + (Q_prev[i, m] - Q_prev[i-1, m]) / epsilon_x) for i=2:(I_-1), m=1:(M-1))
#sum(abs(abs(Q_prev[i,m]) * Q_prev[i,m] - (p_prev[i+1, m]^2 - p_prev[i,m]^2 )/(epsilon_x * alpha_)) for i=1:(I_-1), m=1:M )
#sum(abs(abs(Q_prev[i,m]) * Q_prev[i,m] - (d[i+1, m] - d[i,m] )/(epsilon_x * alpha_)) for i=1:(I_-1), m=1:M )
