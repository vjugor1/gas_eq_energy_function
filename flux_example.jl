using Flux

epsilon_t = 1e-3
epsilon_x = 1e-1
epsilon = epsilon_t / epsilon_x
alpha = 8.57
I_ = 5
M = 20


# Low-level API computation of derivatives.

# Construct some trackable data.
pnt_ = randn(I_ * M) * 100


global A  = get_A_block(I_)
tmp_out = get_init_p_Q(I_, M, left_bound_const, right_bound_const, initial_value_lin, left_bound_Q_const, 42)
init_value_p_2 = tmp_out[1]
init_value_Q   = tmp_out[2]

f(x) = f_d_e(x, init_value_p_2, init_value_Q, epsilon_t, epsilon_x, false)
df(x) = gradient(f, x)

println(df(pnt_))
