push!(LOAD_PATH, pwd())
using NFST
using LinearAlgebra

println( "2d NFST Test" )

#bandwidth
N = ( 16, 8 )

#number of nodes
M = 10000

#create plan
p = NFSTplan( N, M ) 

println("Number of Threads: ", p.num_threads)

#generate random nodes
A = 0.5 .* rand( 2, M )

#set nodes
p.x = A

#generate random Fourier coefficients
fhat = rand( prod(collect(N).-1) )

#set Fourier coefficients
p.fhat = fhat

#transform
println( "trafo time:" )
@time NFST.trafo(p)

#get function values
f2 = p.f

#indices
I = [ [j; i] for i in 1:N[2]-1, j in 1:N[1]-1 ]
I = vec(I)

#define Fourier matrix
F = [ sin(2*pi*A[:,j][1]*I[l][1])*sin(2*pi*A[:,j][2]*I[l][2]) for j in 1:M, l in 1:prod(collect(N).-1) ]

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
