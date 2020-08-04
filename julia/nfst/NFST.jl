module NFST

export NFSTplan, nfst_plan

ending = ".so"

if Sys.iswindows()
	ending = ".dll"
end

# path to .so file
const lib_path = string( @__DIR__, "/libnfstjulia", ending )

# NFST flags
PRE_PHI_HUT = UInt32(1)<<0
FG_PSI = UInt32(1)<<1
PRE_LIN_PSI = UInt32(1)<<2 #
PRE_FG_PSI = UInt32(1)<<3 #
PRE_PSI = UInt32(1)<<4 #
PRE_FULL_PSI = UInt32(1)<<5 #
MALLOC_X = UInt32(1)<<6
MALLOC_F_HAT = UInt32(1)<<7
MALLOC_F = UInt32(1)<<8
FFT_OUT_OF_PLACE = UInt32(1)<<9
FFTW_INIT = UInt32(1)<<10
NFST_SORT_NODES = UInt32(1)<<11
NFST_OMP_BLOCKWISE_ADJOINT = UInt32(1)<<12
PRE_ONE_PSI = (PRE_LIN_PSI| PRE_FG_PSI| PRE_PSI| PRE_FULL_PSI)

# FFTW flags
FFTW_MEASURE = UInt32(0)
FFTW_DESTROY_INPUT = UInt32(1)<<0
FFTW_EXHAUSTIVE = UInt32(1)<<3
FFTW_PATIENT = UInt32(1)<<5
FFTW_ESTIMATE = UInt32(1)<<6

#default flag values
f1_default_1d = UInt32(PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE)
f1_default = UInt32(PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE | NFST_SORT_NODES | NFST_OMP_BLOCKWISE_ADJOINT)
f2_default = UInt32(FFTW_ESTIMATE | FFTW_DESTROY_INPUT)

# dummy struct for C
mutable struct nfst_plan
end

# NFST plan struct
mutable struct NFSTplan{D}
	N::NTuple{D,Int32}      # bandwidth tuple
	M::Int32                # number of nodes
	n::NTuple{D,Int32}      # oversampling per dimension
	m::Int32                # windows size
	f1::UInt32              # NFST flags
	f2::UInt32              # FFTW flags
	init_done::Bool         # bool for plan init
	finalized::Bool	    	# bool for finalizer
	x::Ref{Float64}         # nodes
	f::Ref{Float64}      # function values
	fhat::Ref{Float64}   # Fourier coefficients
	plan::Ref{nfst_plan}    # plan (C pointer)
	function NFSTplan{D}(N::NTuple{D,Int32},M::Int32,n::NTuple{D,Int32},m::Int32,f1::UInt32,f2::UInt32) where D
	# create plan object
	new(N,M,n,m,f1,f2,false,false)
	end
end

# additional constructor for easy use [NFSTplan((N,N),M) instead of NFSTplan{2}((N,N),M)]
function NFSTplan(N::NTuple{D,Integer},M::Integer) where {D}
	if any(x->x<=0,N)
		error("Every entry of N has to be an even, positive integer." )
	end

	if sum(N .% 2) != 0
	error("Every entry of N has to be an even, positive integer." )
	end

	if M <= 0
	error("M has to be a positive integer." )
	end

	# convert N to vector for passing it over to C
	Nv = collect(N)

	# default oversampling
	n = Array{Int32}(2 .^(ceil.(log.(Nv)/log(2)).+1))
	n = NTuple{D,Int32}(n)

	# default NFST flags
	f1 = UInt32(0)

	if D > 1
		f1 = f1_default
	else
		f1 = f1_default_1d
	end

	NFSTplan{D}(NTuple{D,Int32}(N),Int32(M),n,Int32(8),f1,f2_default)
end

function NFSTplan(N::NTuple{D,Integer},M::Integer,n::NTuple{D,Integer},m::Integer=Int32(8),f1::UInt32=(D > 1 ? f1_default : f1_default_1d),f2::UInt32=f2_default) where {D}
	@info "You are using the guru interface. Please consult the README if you are having trouble."

    # safety checks
	if any(x->x<=0,N)
		error("Every entry of N has to be an even, positive integer." )
	end

	if M <= 0
		error("M has to be a positive integer." )
	end

	if any(x->x<=0,n)
		error("Every entry of n has to be an even integer." )
	end

	if n <= N
		error("Every entry of n has to be larger than the corresponding entry in N." )
	end

	if sum(n .% 2) != 0
		error("Every entry of n has to be an even integer." )
	end

	if m <= 0
		error("m has to be a positive integer." )
	end

	NFSTplan{D}(NTuple{D,Int32}(N),Int32(M),NTuple{D,Int32}(n),Int32(m),(f1 | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT),f2)
end

# finalizer
function finalize_plan(P::NFSTplan{D}) where {D}
	if !P.init_done
		error("NFSTplan not initialized.")
	end

	if !P.finalized
		Core.setfield!(P,:finalized,true)
		ccall(("jnfst_finalize", lib_path),Nothing,(Ref{nfst_plan},),P.plan)
	end
end

# allocate plan memory and init with D,N,M,n,m,f1,f2
function nfst_init(p::NFSTplan{D}) where {D}
	# convert N and n to vectors for passing them over to C
	Nv = collect(p.N)
	n = collect(p.n)

	# call init for memory allocation
	ptr = ccall(("jnfst_alloc", lib_path),Ptr{nfst_plan},())

	# set pointer
	Core.setfield!(p,:plan,ptr)

	# initialize values
	ccall(("jnfst_init", lib_path),Nothing,(Ref{nfst_plan},Int32,Ref{Int32},Int32,Ref{Int32},Int32,UInt32,UInt32),ptr,D,Nv,p.M,n,p.m,p.f1,p.f2)
	Core.setfield!(p,:init_done,true)
	finalizer(finalize_plan,p)
end

# overwrite dot notation for plan struct in order to use C memory
function Base.setproperty!(p::NFSTplan{D},v::Symbol,val) where {D}
	# init plan if not done [usually with setting nodes]
	if !p.init_done
		nfst_init(p)
	end

	# prevent bad stuff from happening
	if p.finalized
		error("NFSTplan already finalized")
	end

	# setting nodes, verification of correct size dxM
	if v == :x
		if D == 1
			if typeof(val) != Vector{Float64}
				error("x has to be a Float64 vector.")
			end
			if size(val)[1] != p.M
				error("x has to be a Float64 vector of length M.")
			end
		else
			if typeof(val) != Array{Float64,2}
				error("x has to be a Float64 matrix.")
			end
			if size(val)[1] != D || size(val)[2] != p.M
				error("x has to be a Float64 matrix of size dxM.")
			end
		end
		ptr = ccall(("jnfst_set_x",lib_path),Ptr{Float64},(Ref{nfst_plan},Ref{Cdouble}),p.plan,val)
		Core.setfield!(p,v,ptr)

    # setting values
    elseif v == :f
		if typeof(val) != Array{Float64,1}
			error("f has to be a Float64 vector.")
		end
		if size(val)[1] != p.M
			error("f has to be a Float64 vector of size M.")
		end
		ptr = ccall(("jnfst_set_f",lib_path),Ptr{Float64},(Ref{nfst_plan},Ref{Float64}),p.plan,val)
		Core.setfield!(p,v,ptr)

	# setting Fourier coefficients
	elseif v == :fhat
		if typeof(val) != Array{Float64,1}
			error("fhat has to be a Float64 vector.")
		end
		l = prod(collect(p.N).-1)
		if size(val)[1] != l
			error("fhat has to be a Float64 vector of size prod(N-1).")
		end
		ptr = ccall(("jnfst_set_fhat",lib_path),Ptr{Float64},(Ref{nfst_plan},Ref{Float64}),p.plan,val)
		Core.setfield!(p,v,ptr)

    # prevent modification of NFST plan pointer
	elseif v == :plan
		@warn "You can't modify the C pointer to the NFST plan."
	elseif v == :num_threads
		@warn "You can't currently modify the number of threads of the NFST plan. Use NFST.set_num_threads(nthreads) instead."
	elseif v == :init_done
		@warn "You can't modify this flag."
	elseif v == :N
		@warn "You can't modify the bandwidth, please create an additional plan."
	elseif v == :M
		@warn "You can't modify the number of nodes, please create an additional plan."
	elseif v == :n
		@warn "You can't modify the oversampling parameter, please create an additional plan."
	elseif v == :m
		@warn "You can't modify the window size, please create an additional plan."
	elseif v == :f1
		@warn "You can't modify the NFST flags, please create an additional plan."
	elseif v == :f2
		@warn "You can't modify the FFTW flags, please create an additional plan."
	# handle other set operations the default way
	else
		Core.setfield!(p,v,val)
	end
end

# overwrite dot notation for plan struct in order to use C memory
function Base.getproperty(p::NFSTplan{D},v::Symbol) where {D}
	if v == :x
		if !isdefined(p,:x)
			error("x is not set.")
		end
		ptr = Core.getfield(p,:x)
		if D==1
			return unsafe_wrap(Vector{Float64},ptr,p.M)             # get nodes from C memory and convert to Julia type
		else
			return unsafe_wrap(Matrix{Float64},ptr,(D,Int64(p.M)))  # get nodes from C memory and convert to Julia type
		end
	elseif v == :num_threads
                return ccall(("nfft_get_num_threads", lib_path),Int64,())
	elseif v == :f
		if !isdefined(p,:f)
			error("f is not set.")
		end
		ptr = Core.getfield(p,:f)
		return unsafe_wrap(Vector{Float64},ptr,p.M)  # get function values from C memory and convert to Julia type
	elseif v == :fhat
		if !isdefined(p,:fhat)
			error("fhat is not set.")
		end
		ptr = Core.getfield(p,:fhat)
		return unsafe_wrap(Vector{Float64},ptr,prod(collect(p.N).-1)) # get Fourier coefficients from C memory and convert to Julia type
	else
		return Core.getfield(p,v)
	end
end

# nfst trafo direct [call with NFST.trafo_direct outside module]
function trafo_direct(P::NFSTplan{D}) where {D}
	# prevent bad stuff from happening
	if P.finalized
		error("NFSTplan already finalized")
	end

	if !isdefined(P, :fhat)
		error("fhat has not been set.")
	end

	if !isdefined(P,:x)
		error("x has not been set.")
	end

	ptr = ccall(("jnfst_trafo_direct",lib_path),Ptr{Float64},(Ref{nfst_plan},),P.plan)
	Core.setfield!(P,:f,ptr)
end

# adjoint trafo direct [call with NFST.adjoint_direct outside module]
function adjoint_direct(P::NFSTplan{D}) where {D}
	# prevent bad stuff from happening
	if P.finalized
		error("NFSTplan already finalized")
	end
	if !isdefined(P, :f)
		error("f has not been set.")
	end
	if !isdefined(P,:x)
		error("x has not been set.")
	end
	ptr = ccall(("jnfst_adjoint_direct",lib_path),Ptr{Float64},(Ref{nfst_plan},),P.plan)
	Core.setfield!(P,:fhat,ptr)
end

# nfst trafo [call with NFST.trafo outside module]
function trafo(P::NFSTplan{D}) where {D}
	# prevent bad stuff from happening
	if P.finalized
		error("NFSTplan already finalized")
	end
	if !isdefined(P, :fhat)
		error("fhat has not been set.")
	end
	if !isdefined(P,:x)
		error("x has not been set.")
	end
	ptr = ccall(("jnfst_trafo",lib_path),Ptr{Float64},(Ref{nfst_plan},),P.plan)
	Core.setfield!(P,:f,ptr)
end

# adjoint trafo [call with NFST.adjoint outside module]
function adjoint(P::NFSTplan{D}) where {D}
	# prevent bad stuff from happening
	if P.finalized
		error("NFSTplan already finalized")
	end
	if !isdefined(P, :f)
		error("f has not been set.")
 	end
	if !isdefined(P,:x)
		error("x has not been set.")
	end
	ptr = ccall(("jnfst_adjoint",lib_path),Ptr{Float64},(Ref{nfst_plan},),P.plan)
	Core.setfield!(P,:fhat,ptr)
end

function get_num_threads()
	return ccall(("nfft_get_num_threads",lib_path),Int64,())
end

function set_num_threads(nthreads::Number)
	ccall(("nfft_set_num_threads",lib_path),Nothing,(Int64,),convert(Int64,nthreads))
end

# module end
end
