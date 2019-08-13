push!(LOAD_PATH, pwd())
using fastsum
using LinearAlgebra

println("fastsum test")

# set the parameters:
d = 2
N = 20000
M = 20000
kernel = "multiquadric"
c = 1/sqrt(N)
p = 8
flags = 0
m = p
n = 256
eps_I = p/n
eps_B = max( 1/16, p/n )
nn = 2*n

# create a Plan-Object in Julia
plan = fastsum.Plan( d, N, M, n, p, kernel, c, eps_I, eps_B, nn, m ) 

# generate source nodes in circle of radius 0.25-eps_B/2
r = sqrt.( rand(N) ).*(0.25-eps_B/2)
phi = rand(N).*(2*pi)
X = [ (r.*cos.(phi)) (r.*sin.(phi)) ]
plan.x = X

# generate coefficients alpha_k
alpha = rand(N)+im*rand(N)
plan.alpha = alpha

# generate target nodes in circle of radius 0.25-eps_B/2
r = sqrt.( rand(M) ).*(0.25-eps_B/2)
phi = rand(M).*(2*pi)
Y = [ (r.*cos.(phi)) (r.*sin.(phi)) ]
plan.y = Y

# Start the Transformation
println( "time:" )
@time fastsum.trafo( plan )

f1 = copy( plan.f )
println( "time direct:" )
@time fastsum.trafo_exact( plan )
f2 = copy( plan.f )
error_vector = f1 - f2

E_2 = norm(error_vector)/norm(f1)
E_infty = norm(error_vector, Inf)/norm(plan.alpha,1)
println( "E_2 error:" )
println( E_2 )
println( "E_infty error:" )
println( E_infty );
