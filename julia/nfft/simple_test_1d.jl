push!(LOAD_PATH, pwd())
using NFFT
using LinearAlgebra

println("1d NFFT Test")

# bandwidth
N = 10

#number of nodes
M = 10000

#create plan
p = Plan((N,),M)

println("Number of Threads: ", p.num_threads)

#generate random nodes
A = rand(M).-0.5

#set nodes
p.x = A

#generate random Fourier coefficients
fhat = rand(N)+im*rand(N)

#set Fourier coefficients
p.fhat = fhat

#transform
println( "trafo time:" )
@time NFFT.trafo(p)

#get function values
f2 = p.f

#define Fourier matrix
F = [ exp(-2*pi*im*x_j*k_l) for x_j in A, k_l in -N/2:N/2-1 ]

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
@time NFFT.adjoint(p)

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
