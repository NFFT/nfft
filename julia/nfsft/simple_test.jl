push!(LOAD_PATH, pwd())
using NFSFT
using LinearAlgebra

println("NFSFT Test")
println(NFSFT.nfsft_nfft_default)
# bandwidth
N = 100

#number of nodes
M = 100

# pseudo-random nodes
X = rand(2,M)
X[1,:] .-= 0.5
X[2,:] .*= 0.5

# test init and setting x
p = NFSFTplan(N, M)
p.x = X

# generate pseudo-random Fourier coefficients

fhat = zeros(p.N_total)+im*zeros(p.N_total)
for k = 0:N
    for n = -k:k
        index = NFSFT.nfsft_index(p, k, n)
        fhat[index+1] = (rand() - 0.5) + im * (rand() - 0.5)
    end
end
println(size(fhat))
p.fhat = fhat

# test trafo direct
NFSFT.trafo_direct(p)
f1 = p.f
#print("Vector f (NDSFT):")
#println(f1)

# test trafo
NFSFT.trafo(p)
f2 = p.f
#print("Vector f (NFSFT):")
#println(f2)

# test adjoint direct
NFSFT.adjoint_direct(p)
f3 = p.fhat
#print("Vector fhat (NDSFT):")
#println(f3)

# test fast approximate adjoint
NFSFT.adjoint(p)
f4 = p.fhat
#print("Vector fhat (NFSFT):")
#println(f4)

# calculate the error vectors
error_vector_traf = f1 - f2
error_vector_adj = f3 - f4

println(norm(error_vector_traf)/norm(f1))
println(norm(error_vector_adj, Inf)/norm(f3,1))
println(p.finalized)
NFSFT.finalize_plan(p)
println(p.finalized)