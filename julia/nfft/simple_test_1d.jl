push!(LOAD_PATH, pwd())
using NFFT
using LinearAlgebra

println("1d NFFT Test")

# bandwidth
N = Cint(10)

#number of nodes
M = Cint(10000)

#create plan
p = Plan((N,),M)

#generate random nodes
A = rand(M).-0.5

#set nodes
p.x = A

#node-dependent precomputations
println("precompute time:")
@time NFFT.precompute_psi(p)

#generate random Fourier coefficients
fhat = rand(N)+im*rand(N)

#set Fourier coefficients
p.fhat = fhat

#transform
println("trafo time:")
@time NFFT.trafo(p)

#get function values
f2 = p.f

#define Fourier matrix
F = [ exp(-2*pi*im*x_j*k_l) for x_j in A, k_l in -N/2:N/2-1 ]

#multiply Fourier matrix with vector of Fourier coefficients
f1 = F*fhat

error_vector = f1-f2
println("relative l2 error:")
print( norm(error_vector)/norm(f1) )
println("\nl infinity error:")
print( norm(error_vector, Inf)/norm(fhat,1) )

#adjoint
println("\nadjoint time:")
@time NFFT.adjoint(p)

#get function values
f2 = p.fhat

#multiply Fourier matrix with vector of Fourier coefficients
f1 = F'*p.f

error_vector = f1-f2
println("relative l2 error:")
print( norm(error_vector)/norm(f1) )
println("\nl infinity error:")
print( norm(error_vector, Inf)/norm(p.f,1) )
 
