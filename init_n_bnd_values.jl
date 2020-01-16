



######## Bound functions #########



a1 = 0.2
w1 = 100
b1 = 0.3

a2 = -0.2
w2 = 100
b2 = 0.4
function left_bound_const(x)
    return  b1
end

function right_bound_const(x)
    return b2# + a2 * x
end

function right_bound_lin(x)
    return b2 + a2 * x
end

function left_bound_Q_lin_dec(x)
    return sqrt(abs(a) / alpha) * (1 - x * 10)
end
function left_bound_Q_lin_inc(x)
    return sqrt(abs(a) / alpha) * (1 + x * 10)
end
function left_bound_Q_const(x)
    return sqrt(abs(a) / alpha)
end


function left_bound_Q_right_dec(x)
    return sqrt(abs(a) / alpha) * (1 - x * 10)
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
    I__ = size(init_)[1]
    init_[:,1] = foo.(collect(1:I__) * epsilon_x)
end

function get_init_p_Q(I_, M, lb_foo_p, rb_foo_p, init_foo_p, lb_foo_Q, rnd_seed)
    Random.seed!(rnd_seed)
    init_value_p_2 = (rand(I_, M)).^2
    set_lb_pressure!(init_value_p_2, lb_foo_p)
    set_rb_pressure!(init_value_p_2, rb_foo_p)
    set_init_pressure!(init_value_p_2, init_foo_p, lb_foo_p, rb_foo_p)
    init_value_Q = ones(I_,M) * sqrt(abs(a) / alpha) # / alpha * epsilon_x
    set_lb_Q!(init_value_Q[:, 1], lb_foo_Q)
    return (init_value_p_2, init_value_Q)
end




############ Get initial and boundary conditions ############
