using LinearAlgebra
#using ForwardDiff
using JuMP
using Ipopt
using Random
using Plots


######## GLOBAL PARAMETERS ########
eps_t = 1e-3
eps_x = 1e-1
epsilon = eps_t / eps_x
alpha = 8.57
I_ = 5
M = 80
p_1 = zeros(I_)
for i=1:I_
        p_1[i] = i
end
######## GLOBAL PARAMETERS ########

############ FUNCTIONS ############



function solve_opt!(I_, M, init_value_p_2, init_value_Q, epsilon_t, epsilon_x, p_Q_out, d, out)
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
    @constraint(model, pdes[i=2:I_,m=1:(M-1)], (Q[i, m] - Q[i-1, m]) / epsilon_x + (p[i,m+1] - p[i,m]) / epsilon_t == 0.0)


    obj2 = sum((A' * d[:,m])' * Q[2:I_,m] for m=1:M)

    @variable(model, aux)

    @constraint(model, aux == obj2)
    @constraint(model, feeez[idx=1:(I_)], p[idx,1] == sqrt(init_value_p_2[idx, 1]))

    @constraint(model, feeez_Q[idx=1:(M)], Q[1,idx] == init_value_Q[1, idx])


    #println("forming objective...")
    @NLobjective(model, Min, alpha * (epsilon_x / 3.)*sum(sqrt(Q[i,m]^6) for i=1:I_, m=1:M) -  aux)
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




function p_of_d(d1d2, d2d3, d1d2_abs, d2d3_abs, i, m, p_1)
        @assert i > 1
        @assert i < I_
        @assert m > 1
        @assert m <= M
        @assert size(d1d2) == size(d2d3)
        @assert size(d1d2_abs) == size(d2d3_abs)
        @assert size(d1d2) == size(d2d3_abs)
        @assert size(d1d2) == (m-1,)
        out = p_1 + epsilon / alpha * sum( √(d2d3_abs[m_]) * sign(d2d3[m_]) -
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
#sum(abs(abs(Q_prev[i,m]) * Q_prev[i,m] - (p_prev[i+1, m]^2 - p_prev[i,m]^2 )/(epsilon_x * alpha)) for i=1:(I_-1), m=1:M )
#sum(abs(abs(Q_prev[i,m]) * Q_prev[i,m] - (d[i+1, m] - d[i,m] )/(epsilon_x * alpha)) for i=1:(I_-1), m=1:M )
