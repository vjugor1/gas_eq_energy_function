using LinearAlgebra
using StaticArrays
using FiniteDiff
using SparsityDetection, SparseArrays
fcalls = 0
function f(dx,x) # in-place
  global fcalls += 1
  for i in 2:length(x)-1
    dx[i] = x[i-1] - 2x[i] + x[i+1]
  end
  dx[1] = -2x[1] + x[2]
  dx[end] = x[end-1] - 2x[end]
  nothing
end



x = rand(10)
output = zeros(10,10)
FiniteDiff.finite_difference_jacobian!(output,f,x)
println(output)


in = rand(10)
out = similar(in)
sparsity_pattern = sparsity!(f,out,in)
sparsejac = Float64.(sparse(sparsity_pattern))

cache = FiniteDiff.JacobianCache(x)

@time FiniteDiff.finite_difference_jacobian!(output,f,x,cache) # 0.000008 seconds (7 allocations: 224 bytes)

using SparseDiffTools
colors = matrix_colors(sparsejac)

sparsecache = FiniteDiff.JacobianCache(x,colorvec=colors,sparsity=sparsejac)
FiniteDiff.finite_difference_jacobian!(sparsejac,f,x,sparsecache)


fcalls = 0
FiniteDiff.finite_difference_jacobian!(output,f,x,cache)
println(fcalls) #11

fcalls = 0
FiniteDiff.finite_difference_jacobian!(sparsejac,f,x,sparsecache)
println(fcalls) #4
