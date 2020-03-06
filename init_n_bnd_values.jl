using Interpolations



######## Bound functions #########



a1 = 0.4
w1 = 100
b1 = 0.35

a2 = -0.2
w2 = 100
b2 = 0.4
function left_bound_const(x)
    return  b1
end
function left_bound_const_3rd(x)
    return  0.15
end

function left_bound_inc(x)
    return  b1 + a1 * x
end

function right_bound_const(x)
    return b2# + a2 * x
end

function right_bound_const_3rd(x)
    return b1# + a2 * x
end


function right_bound_lin(x)
    return b2 + a2 * x
end

function left_bound_Q_lin_dec(x)
    return sqrt(abs(a) / alpha_) * (1 - x * 10)
end
function left_bound_Q_lin_inc(x)
    return sqrt(abs(a) / alpha_) * (1 + (x - epsilon_t) * 0)
end
function left_bound_Q_const(x)
    return  sqrt(abs(a) / alpha_) * 2.0
end
function right_bound_Q_const(x)
    return  sqrt(abs(a) / alpha_) * 1.0
end

function left_bound_Q_const_3rd(x)
    return  sqrt(abs(a) / alpha_) * 2.0
end
function right_bound_Q_const_3rd(x)
    return  sqrt(abs(a) / alpha_) * 2.0
end



function left_bound_Q_right_dec(x)
    return sqrt(abs(a) / alpha_) * (1 - x * 10)
end

function initial_value_lin(x, lb_foo, rb_foo)
    ### consistency ###
    global a = (lb_foo(epsilon_t) - rb_foo(epsilon_t)) / (epsilon_x - I_ * epsilon_x)
    global b = lb_foo(epsilon_t) - a * epsilon_x
    ### consistency ###
    return a * x + b
end



############ Bound functions ############

############ Get initial and boundary conditions ############
function set_lb_pressure!(init_, foo)
    M_ = size(init_)[2]
    init_[1,:] = foo.(collect(1:M_) * epsilon_t)
end
function set_rb_pressure!(init_, foo)
    M_ = size(init_)[2]
    init_[I_,:] = foo.(collect(1:M_) * epsilon_t)
end

function set_init_pressure!(init_, foo, lb_foo, rb_foo)
    I__ = size(init_)[1]
    to_apply(x) = foo(x, lb_foo, rb_foo)
    init_[:,1] = to_apply.(collect(1:I__) * epsilon_x)
end

function set_lb_Q!(init_, foo)
    M_ = size(init_)[2]
    init_[1,:] = foo.(collect(1:M_) * epsilon_t)
end
function set_rb_Q!(init_, foo)
    M_ = size(init_)[2]
    I__ = size(init_)[1]
    init_[I__,:] = foo.(collect(1:M_) * epsilon_t)
end
############ Get initial and boundary conditions ############
function get_init_p_Q(I_, M, lb_foo_p, rb_foo_p, init_foo_p, lb_foo_Q, rb_foo_Q, rnd_seed)
    Random.seed!(rnd_seed)
    init_value_p_2 = (rand(I_, M)).^2
    set_lb_pressure!(init_value_p_2, lb_foo_p)
    set_rb_pressure!(init_value_p_2, rb_foo_p)
    set_init_pressure!(init_value_p_2, init_foo_p, lb_foo_p, rb_foo_p)
    init_value_Q = ones(I_,M) * sqrt(abs(a) / alpha_) # / alpha_ * epsilon_x
    set_lb_Q!(init_value_Q, lb_foo_Q)
    set_rb_Q!(init_value_Q, rb_foo_Q)
    return (init_value_p_2, init_value_Q)
end

function get_higher_resolution_init_vals(L, T, p_low, Q_low, eps_x_low, eps_t_low, eps_x_new, eps_t_new)
    I_new = Int(round(L / eps_x_new))
    M_new = Int(round(T / eps_t_new))
    p_new = zeros(I_new, M_new)
    Q_new = zeros(I_new, M_new)
    I_low, M_low = size(p_low)
    x_low = collect(0:(I_low-1)) * eps_x_low .- eps_x_low
    t_low = collect(0:(M_low-1)) * eps_t_low .- eps_t_low
    itp_p = interpolate((x_low, t_low), p_low, Gridded(Linear()))
    itp_Q = interpolate((x_low, t_low), Q_low, Gridded(Linear()))
    for i=1:I_new
        for m=1:M_new
            #println((i-1) * eps_x_new)
            #println((m-1) * eps_t_new)
            p_new[i,m] = itp_p((i-1) * eps_x_new, (m-1) * eps_t_new)
            Q_new[i,m] = itp_Q((i-1) * eps_x_new, (m-1) * eps_t_new)
        end
    end
    return p_new, Q_new, I_new, M_new
end
