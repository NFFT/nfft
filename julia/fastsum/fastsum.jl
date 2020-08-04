module fastsum

ending = ".so"

if Sys.iswindows( )
	ending = ".dll"
elseif Sys.isapple()
	ending = ".dylib"
end

const lib_path = string( @__DIR__, "/libfastsumjulia", ending )

mutable struct fastsum_plan
end

mutable struct Plan
	d::Integer					# dimension
	N::Integer					# number of source nodes
	M::Integer					# number of target nodes
	n::Integer					# expansion degree
	p::Integer					# degree of smoothness
	kernel::String				# name of kernel
	c::Vector{Float64}	# kernel parameters
	eps_I::Real					# inner boundary
	eps_B::Real				# outer boundary
	nn_x::Integer				# oversampled nn in x
	nn_y::Integer				# oversampled nn in y
	m_x::Integer				# NFFT-cutoff in x
	m_y::Integer				# NFFT-cutoff in y
	init_done::Bool			# bool for plan init
	finalized::Bool			# bool for finalizer
	flags::UInt32				# flags
	x::Ref{Float64}			# source nodes
	y::Ref{Float64}			# target nodes
	alpha::Ref{ComplexF64}		# source coefficients
	f::Ref{ComplexF64}				# target evaluations
	plan::Ref{fastsum_plan}		# plan (C pointer)
	
	function Plan( d::Integer, N::Integer, M::Integer, n::Integer, p::Integer, kernel::String, c::Vector{<:Real}, eps_I::Real, eps_B::Real, nn_x::Integer, nn_y::Integer, m_x::Integer, m_y::Integer, flags::UInt32 )
		new( d, N, M, n, p, kernel, Vector{Float64}(c), eps_I, eps_B, nn_x, nn_y, m_x, m_y, false, false, flags )
	end
end #struct fastsumplan

function Plan( d::Integer, N::Integer, M::Integer, n::Integer, p::Integer, kernel::String, c::Real, eps_I::Real, eps_B::Real, nn::Integer, m::Integer )

	if N <= 0
		error( "N has to be a positive Integer." )
	end
	
	if M <= 0
		error( "M has to be a positive Integer." )
	end
	
	if n <= 0
		error( "n has to be a positive Integer. ")
	end
	
	if m <= 0
		error( "m has to be a positive Integer." )
	end
	
	cv = Vector{Float64}( undef, 1 )
	cv[1] = Float64(c)
	
	Plan( d, N, M, n, p, kernel, cv, eps_I, eps_B, nn, nn, m, m, UInt32(0) )
	
end #constructor

function fastsum_init( p::Plan ) 

	ptr = ccall( ("jfastsum_alloc", lib_path), Ptr{fastsum_plan}, () )
	Core.setfield!( p, :plan, ptr )
	
	code = ccall( ("jfastsum_init",lib_path), Int64, (Ref{fastsum_plan}, Int32, Cstring, Ref{Float64}, UInt32, Int32, Int32, Float64, Float64, Int32, Int32, Int32, Int32, Int32, Int32), ptr, Int32(p.d), p.kernel, p.c, p.flags, Int32(p.n), Int32(p.p), Float64(p.eps_I), Float64(p.eps_B), Int32(p.N), Int32(p.M), Int32(p.nn_x), Int32(p.nn_y), Int32(p.m_x), Int32(p.m_y) )
	
	if code == 1
		error( "Unkown kernel." )
	end

	Core.setfield!( p, :init_done, true )
	finalizer( finalize_plan, p )
	
end #fastsum_init

function finalize_plan( p::Plan )

	if !p.init_done
		error( "Plan not initialized." )
	end

	if !p.finalized
		ccall( ("jfastsum_finalize",lib_path), Nothing, (Ref{fastsum_plan},), p.plan )
		Core.setfield!( p, :finalized, true )
	end

end #finalize_plan

function Base.setproperty!( p::Plan, v::Symbol, val )

	if !p.init_done
		fastsum_init( p )
	end

	if p.finalized
		error( "Plan already finalized" )
	end

	# edit source nodes
	if v == :x

		if p.d == 1
			if typeof(val) != Vector{Float64}
				error( "x has to be a Float64 vector." )
			end
			if size(val)[1] != p.N
				error( "x has to be a Float64 vector of length N." )
			end
	    else # => D >1
			if typeof(val) != Array{Float64, 2}
				error( "x has to be a Float64 matrix." )
			end
			if size(val)[1] != p.N || size(val)[2] != p.d
				error( "x has to be a Float64 matrix of size N." )
			end
		end

		ptr = ccall( ("jfastsum_set_x", lib_path), Ptr{Float64}, (Ref{fastsum_plan},Ref{Cdouble}), p.plan, val )
		Core.setfield!( p, v, ptr )

	 # edit target nodes
	 elseif v == :y

		if p.d == 1
			if typeof(val) != Vector{Float64}
				error( "y has to be a Float64 vector." )
			end
			if size(val)[1] != p.M
				error( "y has to be a Float64 vector of length M." )
			end
		else # => D > 1
			if typeof(val) != Array{Float64, 2}
				error( "y has to be a Float64 matrix." )
			end
			if size(val)[1] != p.M || size(val)[2] != p.d
				error( "y has to be a Float64 matrix of size M." )
			end
		end

		ptr = ccall( ("jfastsum_set_y", lib_path), Ptr{Float64}, (Ref{fastsum_plan},Ref{Cdouble}), p.plan, val )
		Core.setfield!( p, v, ptr )

   # edit source coefficients
	elseif v == :alpha

		if typeof(val) != Vector{ComplexF64}
			error( "alpha has to be a ComplexF64 vector." )
		end
		if size(val)[1] != p.N
			error( "alpha has to be a ComplexF64 vector of length N." )
		end

		ptr = ccall( ("jfastsum_set_alpha", lib_path), Ptr{ComplexF64}, (Ref{fastsum_plan},Ref{ComplexF64}), p.plan, val )
		Core.setfield!( p, v, ptr )

	elseif v == :M
	  	@warn("You can't modify the number of target nodes.")
  	elseif v == :N
	  	@warn("You can't modify the number of source nodes.")
  	elseif v == :n
	  	@warn("You can't modify the expansion degree.")
  	elseif v == :m
	  	@warn("You can't modify the cut-off parameter.")
  	elseif v == :p
	  	@warn("You can't modify the degree of smoothness.")
  	elseif v == :kernel
    	@warn("You can't modify the kernel.")
  	elseif v == :c
    	@warn("You can't modify the kernel parameters.")
  	elseif v == :eps_I
    	@warn("You can't modify the inner boundary.")
  	elseif v == :eps_B
    	@warn("You can't modify the outer boundary.")
  	elseif v == :plan
    	@warn("You can't modify the pointer to the fastsum plan.")
  	else
    	Core.setfield!( p, v, val )
  	end
  	
end # Base.setproperty!

# overwrite dot notation for plan struct in order to use C memory
function Base.getproperty( p::Plan, v::Symbol ) 

	if v == :x

		if !isdefined( p, :x )
			error("x is not set.")
		end

		ptr = Core.getfield( p, :x )

		if p.d == 1
			return unsafe_wrap( Vector{Float64}, ptr, p.N ) 				# get source nodes from C memory and convert to Julia type
		else
			return unsafe_wrap( Matrix{Float64}, ptr, (p.d, p.N) )  	# get source nodes from C memory and convert to Julia type
		end

	elseif v == :y

		if !isdefined( p, :y )
			error("y is not set.")
		end

		ptr = Core.getfield( p, :y )

		if p.d == 1
			return unsafe_wrap( Vector{Float64}, ptr, p.M )             
		else
			return unsafe_wrap( Matrix{Float64}, ptr, (p.d, p.M) )
		end 

	elseif v == :alpha

		if !isdefined( p, :alpha )
			error("alpha is not set.")
		end

		ptr = Core.getfield( p, :alpha )
		return unsafe_wrap( Vector{ComplexF64}, ptr, p.N ) 			# get coefficients from C memory and convert to Julia type

	elseif v == :f

		if !isdefined( p, :f )
	    error("f is not set.")
		end

		ptr = Core.getfield( p, :f )
		return unsafe_wrap( Vector{ComplexF64}, ptr, p.M )  			# get function values from C memory and convert to Julia type

	else
		return Core.getfield( p, v )
	end
end # Base.getproperty

function trafo( p::Plan )

	if p.finalized
		error( "Plan already finalized." )
	end

	if !isdefined( p, :x )
		error( "x has not been set." )
	end

	if !isdefined( p, :y )
		error( "y has not been set." )
	end

	if !isdefined( p, :alpha )
		error( "alpha has not been set." )
	end

	ptr = ccall( ( "jfastsum_trafo", lib_path), Ptr{ComplexF64}, (Ref{fastsum_plan},), p.plan )
	Core.setfield!( p, :f, ptr )

end #trafo

function trafo_exact( p::Plan )

	if p.finalized
		error( "Plan already finalized." )
	end

	if !isdefined( p, :x )
		error( "x has not been set." )
	end

	if !isdefined( p, :y )
		error( "y has not been set." )
	end

	if !isdefined( p, :alpha )
		error( "alpha has not been set." )
	end

	ptr = ccall( ("jfastsum_exact", lib_path), Ptr{ComplexF64}, (Ref{fastsum_plan},), p.plan )
	Core.setfield!( p, :f, ptr )

end #trafo

function get_num_threads()
	return ccall(("nfft_get_num_threads",lib_path),Int64,())
end

function set_num_threads(nthreads::Number)
	ccall(("nfft_set_num_threads",lib_path),Nothing,(Int64,),convert(Int64,nthreads))
end

end #module
