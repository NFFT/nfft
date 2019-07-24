push!(LOAD_PATH, pwd())
using fastsum
using LinearAlgebra

println("fastsum test")

# set the parameters:
d = 2
N = 2000
M = 2000
kernel = "multiquadric"
c = 1/sqrt(N)
p = 3
flags = 0
m = p
n = 156
eps_I = p/n
eps_B =  1/16
nn = 2*n

# create a Plan-Object in Julia
p = fastsum.Plan( d, N, M, n, p, kernel, c, eps_I, eps_B, nn, m ) 

# generate source nodes in circle of radius 0.25-eps_B/2
r = sqrt.( rand(N) ).*(0.25-eps_B/2)
phi = rand(N).*(2*pi)
X = vcat( (r.*cos.(phi))', (r.*sin.(phi))' )
p.x = X

# generate coefficients alpha_k
alpha = rand(N)+im*rand(N)
p.alpha = alpha

# generate target nodes in circle of radius 0.25-eps_B/2
r = sqrt.( rand(M) ).*(0.25-eps_B/2)
phi = rand(M).*(2*pi)
Y = vcat( (r.*cos.(phi))', (r.*sin.(phi))' )
p.y = Y

# Start the Transformation
fastsum.trafo( p )

f1 = copy( p.f )
fastsum.trafo_exact( p )
f2 = copy( p.f )
error_vector = f1 - f2

E_2 = norm(error_vector)/norm(f1)
E_infty = norm(error_vector, Inf)/norm(p.alpha,1)
println( "E_2 error:" )
println( E_2 )
println( "E_infty error:" )
println( E_infty );
