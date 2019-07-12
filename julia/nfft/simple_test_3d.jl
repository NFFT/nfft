push!(LOAD_PATH, pwd())
using NFFT
using LinearAlgebra

println( "3d NFFT Test" )

# bandwidth
N = ( 16, 8, 4 )

#number of nodes
M = 10000

#create plan
p = Plan(N,M)

println("Number of Threads: ", p.num_threads)

#generate random nodes
A = rand(3,M).-0.5

#set nodes
p.x = A

#generate random Fourier coefficients
fhat = rand(prod(N))+im*rand(prod(N))

#set Fourier coefficients
p.fhat = fhat

#transform
println( "trafo time:" )
@time NFFT.trafo(p)

#get function values
f2 = p.f

#indices
I = [ [j; i; k] for k in -N[3]/2:N[3]/2-1, i in -N[2]/2:N[2]/2-1, j in -N[1]/2:N[1]/2-1 ]
I = vec(I)

#define Fourier matrix
F = [ exp(-2*pi*im*sum(A[:,j]'*I[l])) for j in 1:M, l in 1:prod(N) ]

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
