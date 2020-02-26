using Nabla

epsilon_t = 1e-3
epsilon_x = 1e-1
epsilon = epsilon_t / epsilon_x
alpha = 8.57
I_ = 5
M = 20


# Low-level API computation of derivatives.

# Construct some trackable data.
pnt_ = randn(I_ * M) * 100

pnt = Leaf(Tape(), pnt_)
println("!!!!!!!!!!")
println(typeof(pnt))
println("!!!!!!!!!!")
# Compute the forward pass.
global A  = get_A_block(I_)
tmp_out = get_init_p_Q(I_, M, left_bound_const, right_bound_const, initial_value_lin, left_bound_Q_const, 42)
init_value_p_2 = tmp_out[1]
init_value_Q   = tmp_out[2]
foo(x) = f_d_e(x, init_value_p_2, init_value_Q, epsilon_t, epsilon_x, false)
z = foo(pnt)

println("Output of the forward pass is:")
println(z)
println()
println("y is $(Nabla.unbox(z)).")
println()

# Get the reverse tape.
z̄ = ∇(z)
println("Output of reverse-pass is")
#println(z̄)
println()

# Index into the reverse tape using x to get the gradient of `y` w.r.t. `x`.
x̄ = z̄[pnt]
println("Gradient of z w.r.t. x at $x_ is $x̄.")
println()

ȳ = z̄[y]
println("Gradient of z w.r.t. y at $y_ is $ȳ")


# (Current) High-Level API computation of derivatives. I will probably maintain this
# interface and extend is significantly as it is currently rather limited. It just returns
# a function which returns the gradient, which isn't really what you want.

# Define the function to be differentiated. Parameters w.r.t. which we want gradients must
# be arguments. Parameters that we don't want gradients w.r.t. should be passed in via a
# closure.
@unionise f(x::AbstractVector, y::AbstractVector) = transpose(x)A * y

# Compute a function `∇f` which computes the derivative of `f` w.r.t. the inputs.
∇f = ∇(f)

# Compute the derivative of `f` w.r.t. `x` at `x_`. Result is currently a 1-Tuple. Might
# introduce a special case for unary functions where it just returns the result.
(x̄, ȳ) = ∇f(x_, y_)

@assert x̄ == z̄[x]
@assert ȳ == z̄[y]
