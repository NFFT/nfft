module NFFT

export Plan,nfft_plan

# file ending for OS
ending = ".so"

if Sys.iswindows()
    ending = ".dll"
elseif Sys.isapple()
    ending = ".dylib"
end

# path to .so file
const lib_path = string( "$(pwd())/libnfftjulia", ending )

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

# FFTW flags
FFTW_MEASURE = UInt32(0)
FFTW_DESTROY_INPUT = UInt32(1)<<0
FFTW_UNALIGNED = UInt32(1)<<1
FFTW_CONSERVE_MEMORY = UInt32(1)<<2
FFTW_EXHAUSTIVE = UInt32(1)<<3
FFTW_PRESERVE_INPUT = UInt32(1)<<4
FFTW_PATIENT = UInt32(1)<<5
FFTW_ESTIMATE = UInt32(1)<<6
FFTW_WISDOM_ONLY = UInt32(1)<<21

# dummy struct for C
mutable struct nfft_plan
end

# NFFT plan struct
mutable struct Plan{D}
    N::NTuple{D,Int32}      # bandwidth tuple
    M::Int32                # number of nodes
    n::NTuple{D,Int32}      # oversampling per dimension
    m::Int32                # windows size
    f1::UInt32              # NFFT flags
    f2::UInt32              # FFTW flags
    X::Ref{Float64}         # nodes
    f::Ref{ComplexF64}      # function values
    fhat::Ref{ComplexF64}   # Fourier coefficients
    plan::Ref{nfft_plan}    # plan (C pointer)
    precompute_done::Bool   # bool for precompute
    init_done::Bool         # bool for plan init
    function Plan{D}(N::NTuple{D,Int32},M::Int32) where D
        # convert N to vector for passing it over to C
        Nv = collect(N)
        # default oversampling
        n = Array{Int32}(2 .^(ceil.(log.(Nv)/log(2)).+1))
        n = NTuple{D,Int32}(n)
        # default NFFT flags
        f1 = UInt32(0)
        if D > 1
            f1 = UInt32(PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE | NFFT_SORT_NODES | NFFT_OMP_BLOCKWISE_ADJOINT)
        else
            f1 = UInt32(PRE_PHI_HUT | PRE_PSI | MALLOC_X | MALLOC_F_HAT | MALLOC_F | FFTW_INIT | FFT_OUT_OF_PLACE)
        end
        # default FFTW flags
        f2 = UInt32(FFTW_ESTIMATE | FFTW_DESTROY_INPUT)
        # create plan object
        new(N,M,n,Int32(6),f1,f2)
    end
end

# additional constructor for easy use [Plan((N,N),M) instead of Plan{2}((N,N),M)]
function Plan(N::NTuple{D,Int32},M::Int32) where {D}
    Plan{D}(N,M)
end

# allocate plan memory and init with D,N,M,n,m,f1,f2
function nfft_init(p::Plan{D}) where {D}
    # convert N and n to vectors for passing them over to C
    Nv = collect(p.N)
    n = collect(p.n)
    # call init for memory allocation
    ptr = ccall(("jnfft_alloc", lib_path),Ptr{nfft_plan},())
    # set pointer
    Core.setfield!(p,:plan,ptr)
    # initialize values
    ccall(("jnfft_init", lib_path),Nothing,(Ref{nfft_plan},Int32,Ref{Int32},Int32,Ref{Int32},Int32,UInt32,UInt32),ptr,D,Nv,p.M,n,p.m,p.f1,p.f2)
    Core.setfield!(p,:init_done,true)
end

# overwrite dot notation for plan struct in order to use C memory
function Base.setproperty!(p::Plan{D},v::Symbol,val) where {D}
    # init plan if not done [usually with setting nodes]
    if !p.init_done
        nfft_init(p)
    end
    # setting nodes, verification of correct size Mxd
    if v == :X
        if D == 1
            if typeof(val) != Vector{Float64}
                error("X has to be a vector.")
            end
            if size(val)[1] != p.M
                error("X has to be a vector of length M.")
            end
        else
            if typeof(val) != Array{Float64,2}
                error("X has to be a matrix.")
            end
            if size(val)[1] != D || size(val)[2] != p.M
                error("X has to be a matrix of size Mxd.")
            end
        end
        ptr = ccall(("jnfft_set_x",lib_path),Ptr{Float64},(Ref{nfft_plan},Ref{Cdouble}),p.plan,val)
        Core.setfield!(p,v,ptr)
        # new precomputations necessary if nodes changed
        Core.setfield!(p,:precompute_done,false)
    # setting values
    elseif v == :f
        if typeof(val) != Array{ComplexF64,1}
            error("f has to be a vector.")
        end
        if size(val)[1] != p.M
            error("f has to be a vector of size M.")
        end
        ptr = ccall(("jnfft_set_f",lib_path),Ptr{ComplexF64},(Ref{nfft_plan},Ref{ComplexF64}),p.plan,val)
        Core.setfield!(p,v,ptr)
    # setting Fourier coefficients
    elseif v == :fhat
        if typeof(val) != Array{ComplexF64,1}
            error("fhat has to be a vector.")
        end
        l = prod(p.N)
        if size(val)[1] != l
            error("fhat has to be a vector of size N[1]*N[2]*N[3].")
        end
        ptr = ccall(("jnfft_set_fhat",lib_path),Ptr{ComplexF64},(Ref{nfft_plan},Ref{ComplexF64}),p.plan,val)
        Core.setfield!(p,v,ptr)
    # prevent modification of NFFT plan pointer
    elseif v == :plan
        error("You can't modify the C pointer to the NFFT plan.")
    # handle other set operations the default way
    else
        Core.setfield!(p,v,val)
    end
end

# overwrite dot notation for plan struct in order to use C memory
function Base.getproperty(p::Plan{D},v::Symbol) where {D}
    if v == :X
        if !isdefined(p,:X)
            error("X is not set.")
        end
        ptr = Core.getfield(p,:X)
        if D==1
            return unsafe_wrap(Vector{Float64},ptr,p.M)             # get nodes from C memory and convert to Julia type
        else
            return unsafe_wrap(Matrix{Float64},ptr,(D,Int64(p.M)))  # get nodes from C memory and convert to Julia type
        end
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
        return unsafe_wrap(Vector{ComplexF64},ptr,prod(p.N)) # get Fourier coefficients from C memory and convert to Julia type
    else
        return Core.getfield(p,v)
    end
end

# precomputations [call with NFFT.precompute_psi outside module]
function precompute_psi(P::Plan{D}) where {D}
	if !isdefined(P,:X)
		error("X has not been set.")
	end
    ccall(("jnfft_precompute_psi",lib_path),Nothing,(Ref{nfft_plan},),P.plan)
    P.precompute_done = true
end

# nfft trafo [call with NFFT.trafo outside module]
function trafo(P::Plan{D}) where {D}
	if !isdefined(P, :fhat)
		error("fhat has not been set.")
	end
    if !P.precompute_done
        error("Please do the precomputations with nfft_precompute_psi.")
    end
    ptr = ccall(("jnfft_trafo",lib_path),Ptr{ComplexF64},(Ref{nfft_plan},),P.plan)
    Core.setfield!(P,:f,ptr)
end

# adjoint trafo [call with NFFT.adjoint outside module]
function adjoint(P::Plan{D}) where {D}
    if !isdefined(P, :f)
        error("f has not been set.")
    end
    if !P.precompute_done
        error("Please do the precomputations with nfft_precompute_psi.")
    end
    ptr = ccall(("jnfft_adjoint",lib_path),Ptr{ComplexF64},(Ref{nfft_plan},),P.plan)
    Core.setfield!(P,:fhat,ptr)
end

# module end
end