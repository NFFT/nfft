# Copyright (c) 2018-2019, 2025 The NFFT Authors
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, see <https://www.gnu.org/licenses/>.

push!(LOAD_PATH, pwd())
using NFCT
using LinearAlgebra
using Random

# Set random seed for reproducible results
Random.seed!(1234)

println("1d NFCT Test")

# bandwidth
N = 50

#number of nodes
M = 100

#create plan
p = NFCTplan((N,),M)

println("Number of Threads: ", p.num_threads)

#generate random nodes
A = 0.5 .* rand(M)

#set nodes
p.x = A

#generate random Fourier coefficients
fhat = rand(N)

#set Fourier coefficients
p.fhat = fhat

#transform
println( "trafo time:" )
@time NFCT.trafo(p)

#get function values
f2 = p.f

#define Fourier matrix
F = [ cos(2*pi*k_l*x_j) for x_j in A, k_l in 0:N-1 ]

#multiply Fourier matrix with vector of Fourier coefficients
f1 = F*fhat

error_vector = f1-f2
E_2 = norm(error_vector)/norm(f1)
E_infty = norm(error_vector, Inf)/norm(fhat,1)
println( "E_2 error:" )
println( E_2 )
println( "E_infty error:" )
println( E_infty )

if ( E_2 >= 2.5e-7 ) || ( E_infty >= 2.5e-7 )
	error( "Errors are too large." )
end

#adjoint
println( "adjoint time:" )
@time NFCT.adjoint(p)

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

if ( E_2 >= 2.5e-7 ) || ( E_infty >= 2.5e-7 )
	error( "Errors are too large." )
end
