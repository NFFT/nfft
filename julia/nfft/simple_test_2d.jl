push!(LOAD_PATH, pwd())
using NFFT
using LinearAlgebra

p = Nothing
N = Nothing
M = Nothing
A = Nothing
fhat = Nothing

println("2d NFFT Test")

# bandwidth
N = ( Cint(16), Cint(8) )

#number of nodes
M = Cint(10000)

#create plan
p = Plan(N,M)

#generate random nodes
A = rand(2,M).-0.5

#set nodes
p.x = A

#node-dependent precomputations
println("precompute time:")
@time NFFT.precompute_psi(p)

#generate random Fourier coefficients
fhat = rand(prod(N))+im*rand(prod(N))

#set Fourier coefficients
p.fhat = fhat

#transform
println("trafo time:")
@time NFFT.trafo(p)

#get function values
f2 = p.f

#indices
I = [ [j; i] for i in -N[2]/2:N[2]/2-1, j in -N[1]/2:N[1]/2-1 ]
I = vec(I)

#define Fourier matrix
F = [ exp(-2*pi*im*sum(A[:,j]'*I[l])) for j in 1:M, l in 1:prod(N) ]

#multiply Fourier matrix with vector of Fourier coefficients
f1 = F*fhat

error_vector = f1-f2
println("relative l2 error:")
print( norm(error_vector)/norm(f1) )
println("\nl infinity error:")
print( norm(error_vector, Inf)/norm(fhat,1) ) 
