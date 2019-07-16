push!(LOAD_PATH, pwd())
using NFST
using LinearAlgebra

println("1d NFST Test")

# bandwidth
N = 50

#number of nodes
M = 100

#create plan
p = NFSTplan((N,),M)

println("Number of Threads: ", p.num_threads)

#generate random nodes
A = 0.5 .* rand(M)

#set nodes
p.x = A

#generate random Fourier coefficients
fhat = rand(N-1)

#set Fourier coefficients
p.fhat = fhat

#transform
println( "trafo time:" )
@time NFST.trafo(p)

#get function values
f2 = p.f

#define Fourier matrix
F = [ sin(2*pi*k_l*x_j) for x_j in A, k_l in 1:N-1 ]

#multiply Fourier matrix with vector of Fourier coefficients
f1 = F*fhat

error_vector = f1-f2
E_2 = norm(error_vector)/norm(f1)
E_infty = norm(error_vector, Inf)/norm(fhat,1)
println( "E_2 error:" )
println( E_2 )
println( "E_infty error:" )
println( E_infty )

if ( E_2 >= 1e-8 ) || ( E_infty >= 1e-8 )
	error( "Errors are too large." )
end

#adjoint
println( "adjoint time:" )
@time NFST.adjoint(p)

#get function values
f2 = p.fhat

#multiply Fourier matrix with vector of Fourier coefficients
f1 = F'*p.f

error_vector = f1-f2
E_2 = norm(error_vector)/norm(f1)
E_infty = norm(error_vector, Inf)/norm(fhat,1)
println( "E_2 error:" )
println( E_2 )
println( "E_infty error:" )
println( E_infty )

if ( E_2 >= 1e-8 ) || ( E_infty >= 1e-8 )
	error( "Errors are too large." )
end
