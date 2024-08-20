module NFSFT

export NFSFTplan,nfsft_plan

# file ending for OS
ending = ".so"

if Sys.iswindows()
	ending = ".dll"
elseif Sys.isapple()
	ending = ".dylib"
end

# path to .so file
const lib_path = string( @__DIR__, "/libnfsftjulia", ending )

# NFFT flags
PRE_PHI_HUT = UInt32(1)<<0
FG_PSI = UInt32(1)<<1
PRE_LIN_PSI = UInt32(1)<<2
PRE_FG_PSI = UInt32(1)<<3
PRE_PSI = UInt32(1)<<4
PRE_FULL_PSI = UInt32(1)<<5
MALLOC_X = UInt32(1)<<6
MALLOC_F_HAT = UInt32(1)<<7
MALLOC_F = UInt32(1)<<8
FFT_OUT_OF_PLACE = UInt32(1)<<9
FFTW_INIT = UInt32(1)<<10
NFFT_SORT_NODES = UInt32(1)<<11
NFFT_OMP_BLOCKWISE_ADJOINT = UInt32(1)<<12
PRE_ONE_PSI = (PRE_LIN_PSI| PRE_FG_PSI| PRE_PSI| PRE_FULL_PSI)

# NFSFT flags
NFSFT_NORMALIZED = UInt32(1)<<0
NFSFT_USE_NDFT = UInt32(1)<<1
NFSFT_USE_DPT = UInt32(1)<<2
NFSFT_MALLOC_X = UInt32(1)<<3
NFSFT_MALLOC_F_HAT = UInt32(1)<<5
NFSFT_MALLOC_F = UInt32(1)<<6
NFSFT_PRESERVE_F_HAT = UInt32(1)<<7
NFSFT_PRESERVE_X = UInt32(1)<<8
NFSFT_PRESERVE_F = UInt32(1)<<9
NFSFT_DESTROY_F_HAT = UInt32(1)<<10
NFSFT_DESTROY_X = UInt32(1)<<11
NFSFT_DESTROY_F = UInt32(1)<<12
NFSFT_NFSFT_NO_DIRECT_ALGORITHM = UInt32(1)<<13
NFSFT_NO_FAST_ALGORITHM = UInt32(1)<<14
NFSFT_ZERO_F_HAT = UInt32(1)<<16
NFSFT_EQUISPACED  = UInt32(1)<<17

# default flag values
nfsft_default = UInt32(NFSFT_MALLOC_X | NFSFT_MALLOC_F | NFSFT_MALLOC_F_HAT)
#nfsft_nfft_default = UInt32(PRE_PHI_HUT | PRE_PSI | FFTW_INIT | NFFT_OMP_BLOCKWISE_ADJOINT)
nfsft_nfft_default = UInt32(PRE_PHI_HUT | PRE_PSI | FFTW_INIT | FFT_OUT_OF_PLACE)

# default window cut off
nfsft_default_nfft_cut_off = 6

# dummy struct for C
mutable struct nfsft_plan
end
# NFFT plan struct

mutable struct NFSFTplan
	N::Int32                # bandwidth tuple
    N_total::Int32          # Fourier coefficients length
	M::Int32                # number of nodes
	flags::UInt32           # NFFT flags
	nfft_flags::UInt32      # FFTW flags
    nfft_cutoff::Int32      # window cut off
	init_done::Bool         # bool for plan init
	finalized::Bool	    	# bool for finalizer
	x::Ref{Float64}         # nodes
	f::Ref{ComplexF64}      # function values
	fhat::Ref{ComplexF64}   # Fourier coefficients
	plan::Ref{nfsft_plan}   # plan (C pointer)
	function NFSFTplan(N::Int32,M::Int32,flags::UInt32,nfft_flags::UInt32,nfft_cutoff::Int32)
	# create plan object
	new(N,(2*N+2)^2,M,flags,nfft_flags,nfft_cutoff,false,false)
	end
end

function NFSFTplan(N::Integer,M::Integer,flags::UInt32=nfsft_default,nfft_flags::UInt32=nfsft_nfft_default,nfft_cutoff::Integer=Int32(nfsft_default_nfft_cut_off))
    # safety checks
	if N <= 0
		error("Invalid N: " + N + ". Argument must be a positive integer")
	end

	if M <= 0
		error("Invalid M: " + M + ". Argument must be a positive integer")
	end

	NFSFTplan(Int32(N),Int32(M),flags,nfft_flags, Int32(nfft_cutoff))
end

# finalizer
function finalize_plan(P::NFSFTplan)
	if !P.init_done
		error("Plan not initialized.")
	end
	
	if !P.finalized
		Core.setfield!(P,:finalized,true)
		ccall(("jnfsft_finalize", lib_path),Nothing,(Ref{nfsft_plan},),P.plan)
	end
end

# allocate plan memory and init
function nfsft_init(p::NFSFTplan)
	# call init for memory allocation
	ptr = ccall(("jnfsft_alloc", lib_path),Ptr{nfsft_plan},())

	# set pointer
	Core.setfield!(p,:plan,ptr)

	# initialize values
	ccall(("jnfsft_init", lib_path),Nothing,(Ref{nfsft_plan},Int32,Int32,UInt32,UInt32,Int32),ptr,p.N,p.M,p.flags,p.nfft_flags,p.nfft_cutoff)
	Core.setfield!(p,:init_done,true)
	finalizer(finalize_plan,p)
end

# overwrite dot notation for plan struct in order to use C memory
function Base.setproperty!(p::NFSFTplan,v::Symbol,val)
	# init plan if not done [usually with setting nodes]
	if !p.init_done
		nfsft_init(p)
	end

	# prevent bad stuff from happening
	if p.finalized
		error("NFSFTplan already finalized")
	end

	# setting nodes, verification of correct size dxM
	if v == :x
        if typeof(val) != Array{Float64,2}
            error("x has to be a Float64 matrix.")
        end
        if size(val)[1] != 2 || size(val)[2] != p.M
            error("x has to be a Float64 matrix of size 2xM.")
        end

		ptr = ccall(("jnfsft_set_x",lib_path),Ptr{Float64},(Ref{nfsft_plan},Ref{Cdouble}),p.plan,val)
		Core.setfield!(p,v,ptr)

    # setting values
    elseif v == :f
		if typeof(val) != Array{ComplexF64,1}
			error("f has to be a ComplexFloat64 vector.")
		end
		if size(val)[1] != p.M
			error("f has to be a ComplexFloat64 vector of size M.")
		end
		ptr = ccall(("jnfsft_set_f",lib_path),Ptr{ComplexF64},(Ref{nfsft_plan},Ref{ComplexF64}),p.plan,val)
		Core.setfield!(p,v,ptr)

	# setting Fourier coefficients
	elseif v == :fhat
		if typeof(val) != Array{ComplexF64,1}
			error("fhat has to be a ComplexFloat64 vector.")
		end
		if size(val)[1] != p.N_total
			error("fhat has to be a ComplexFloat64 vector of size (2*N+2)^2.")
		end
		ptr = ccall(("jnfsft_set_fhat",lib_path),Ptr{ComplexF64},(Ref{nfsft_plan},Ref{ComplexF64}),p.plan,val)
		Core.setfield!(p,v,ptr)

    # prevent modification of NFSFT plan pointer
	elseif v == :plan
		@warn "You can't modify the C pointer to the NFSFT plan."
	elseif v == :num_threads
		@warn "You can't currently modify the number of threads of the NFSFT plan. Use NFSFT.set_num_threads(nthreads) instead."
	elseif v == :init_done
		@warn "You can't modify this flag."
	elseif v == :N
		@warn "You can't modify the bandwidth, please create an additional plan."
	elseif v == :M
		@warn "You can't modify the number of nodes, please create an additional plan."
	elseif v == :flags
		@warn "You can't modify the NFSFT flags, please create an additional plan."
	elseif v == :nfft_flags
		@warn "You can't modify the NFFT flags, please create an additional plan."
	elseif v == :nfft_cutoff
		@warn "You can't modify the nfft_cutoff, please create an additional plan."
	# handle other set operations the default way
	else
		Core.setfield!(p,v,val)
	end
end


# overwrite dot notation for plan struct in order to use C memory
function Base.getproperty(p::NFSFTplan,v::Symbol)
	if v == :x
		if !isdefined(p,:x)
			error("x is not set.")
		end
		ptr = Core.getfield(p,:x)
		return unsafe_wrap(Matrix{Float64},ptr,(2,Int64(p.M)))  # get nodes from C memory and convert to Julia type
	elseif v == :num_threads
		return ccall(("nfft_get_num_threads", lib_path),Int64,())
	elseif v == :f
		if !isdefined(p,:f)
			error("f is not set.")
		end
		ptr = Core.getfield(p,:f)
		return unsafe_wrap(Vector{ComplexF64},ptr,p.M)  # get function values from C memory and convert to Julia type
	elseif v == :fhat
		if !isdefined(p,:fhat)
			error("fhat is not set.")
		end
		ptr = Core.getfield(p,:fhat)
		return unsafe_wrap(Vector{ComplexF64},ptr,p.N_total) # get Fourier coefficients from C memory and convert to Julia type
	else
		return Core.getfield(p,v)
	end
end

function nfsft_index(p::NFSFTplan, k::Integer, n::Integer)::Integer
    return (2*p.N+2)*(p.N-n+1)+(p.N+k+1)
end


# nfsft trafo direct [call with NFSFT.trafo_direct outside module]
function trafo_direct(P::NFSFTplan)
	# prevent bad stuff from happening
	if P.finalized
		error("NFSFTplan already finalized")
	end

	if !isdefined(P, :fhat)
		error("fhat has not been set.")
	end

	if !isdefined(P,:x)
		error("x has not been set.")
	end

	ptr = ccall(("jnfsft_trafo_direct",lib_path),Ptr{ComplexF64},(Ref{nfsft_plan},),P.plan)
	Core.setfield!(P,:f,ptr)
end


# adjoint trafo direct [call with NFSFT.adjoint_direct outside module]
function adjoint_direct(P::NFSFTplan)
	# prevent bad stuff from happening
	if P.finalized
		error("NFSFTplan already finalized")
	end
	if !isdefined(P, :f)
		error("f has not been set.")
	end
	if !isdefined(P,:x)
		error("x has not been set.")
	end
	ptr = ccall(("jnfsft_adjoint_direct",lib_path),Ptr{ComplexF64},(Ref{nfsft_plan},),P.plan)
	Core.setfield!(P,:fhat,ptr)
end

# nfsft trafo [call with NFSFT.trafo outside module]
function trafo(P::NFSFTplan)
	# prevent bad stuff from happening
	if P.finalized
		error("NFSFTplan already finalized")
	end
	if !isdefined(P, :fhat)
		error("fhat has not been set.")
	end
	if !isdefined(P,:x)
		error("x has not been set.")
	end
	ptr = ccall(("jnfsft_trafo",lib_path),Ptr{ComplexF64},(Ref{nfsft_plan},),P.plan)
	Core.setfield!(P,:f,ptr)
end

# adjoint trafo [call with NFSFT.adjoint outside module]
function adjoint(P::NFSFTplan)
	# prevent bad stuff from happening
	if P.finalized
		error("NFSFTplan already finalized")
	end
	if !isdefined(P, :f)
		error("f has not been set.")
 	end
	if !isdefined(P,:x)
		error("x has not been set.")
	end
	ptr = ccall(("jnfsft_adjoint",lib_path),Ptr{ComplexF64},(Ref{nfsft_plan},),P.plan)
	Core.setfield!(P,:fhat,ptr)
end


# module end
end
