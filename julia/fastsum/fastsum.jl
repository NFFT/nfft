module fastsum

export fastsumplan, fastsum_plan

ending = ".so"

if Sys.iswindows()
	ending = ".dll"
elseif Sys.isapple()
	ending = ".dylib"
end

const lib_path = string( @__DIR__, "/libfastsumjulia", ending )

#dummy-struct for C
mutable struct fastsum_plan
  end


mutable struct fastsumplan{D}
	N::Int32                # number of source nodes
	M::Int32                # number of target nodes
	n::Int32                # expansion degree
	m::Int32                # cut-off parameter
	p::Int32                # degree of smoothness
	kernel::String          # name of kernel
	c::NTuple{D,Float64}    # kernel parameters
	eps_I::Float64          # inner boundary
	eps_B::Float64          # outer boundary
	init_done::Bool         # bool for plan init
	finalized::Bool	      	# bool for finalizer
	x::Ref{Float64}         # source nodes
	y::Ref{Float64}         # target nodes
  alpha::Ref{ComplexF64}  # source coefficients
	f::Ref{ComplexF64}      # target evaluations
	plan::Ref{fastsum_plan} # plan (C pointer)
	function fastsumplan{D}(N::Int32,M::Int32,n::Int32,m::Int32,p::Int32,kernel::String,c::NTuple{D,Float64},eps_I::Float64,eps_B::Float64) where {D}
	# create plan object
	new(N,M,n,m,p,kernel,c,eps_I,eps_B,false,false)
	end

end #struct fastsumplan

function fastsumplan(N::Integer,M::Integer,n::Integer,m::Integer,p::Integer,kernel::String,c::NTuple{D,Float64},eps_I::Float64,eps_B::Float64) where{D}


  if N <= 0
    error("N has to be a positive Integer.")
  end
  if M <= 0
    error("M has to be a positive Integer.")
  end
  if n <= 0
    error("n has to be a positive Integer.")
  end
  if m <= 0
    error("m has to be a positive Integer.")
  end


  fastsumplan{D}(Int32(N),Int32(M),Int32(n),Int32(m),Int32(p),kernel,c,Float64(eps_I),Float64(eps_B))
end #constructor


function fastsum_init(fp::fastsumplan{D}) where {D}

  ptr = ccall(("jfastsum_alloc", lib_path), Ptr{fastsum_plan}, ())

  Core.setfield!(fp, :plan, ptr)
	# c noch in pointer umwandeln
	c = collect(fp.c)
  ccall(("jfastsum_init",lib_path),Nothing,(Ref{fastsum_plan},Int32,Int32,Int32,Cstring,Ref{Float64},Int32,Int32,Int32,Float64,Float64),ptr,D,fp.N,fp.M,fp.kernel,c,fp.n,fp.m,fp.p,fp.eps_I,fp.eps_B)

	Core.setfield!(fp,:init_done, true)
  finalizer(finalize_plan, fp)
end #fastsum_init

function finalize_plan(fp::fastsumplan{D}) where{D}
  if !fp.init_done
    error("Plan not initialized.")
  end

  if !fp.finalized
    ccall(("jfastsum_finalize",lib_path),Nothing,(Ref{fastsum_plan},),fp.plan)
    Core.setfield!(fp,:finalized,true)
  end
end #finalize_plan

function Base.setproperty!(fp::fastsumplan{D},v::Symbol,val) where {D}
  if !fp.init_done
	  fastsum_init(fp)
  end

  if fp.finalized
	  error("Plan already finalized")
  end

  # edit source nodes
  if v == :x
      if D==1
				 if typeof(val) != Vector{Float64}
           error("x has to be a Float64 vector.")
         end
         if (size(val)[1]) != fp.N
           error("x has to be a Float64 vector of length N.")
         end
	    else # => D >1
         if typeof(val) != Array{Float64, 2}
           error("x has to be a Float64 matrix.")
         end
         if size(val)[1] != D || size(val)[2] != fp.N
           error("x has to be a Float64 matrix of size N.")
         end
			end
      ptr = ccall(("jfastsum_set_x", lib_path), Ptr{Float64}, (Ref{fastsum_plan},Ref{Cdouble}), fp.plan,val)
      Core.setfield!(fp,v,ptr)

	 # edit target nodes
	 elseif v == :y
		  if D==1
			  if typeof(val) != Vector{Float64}
				  error("y has to be a Float64 vector.")
			  end
			  if (size(val)[1]) != fp.M
				  error("y has to be a Float64 vector of length M.")
			  end
	    else # => D > 1
			  if typeof(val) != Array{Float64, 2}
				  error("y has to be a Float64 matrix.")
			  end
			  if size(val)[1] != D || size(val)[2] != fp.M
				  error("y has to be a Float64 matrix of size M.")
			  end
	    end
              
			ptr = ccall(("jfastsum_set_y", lib_path), Ptr{Float64}, (Ref{fastsum_plan},Ref{Cdouble}), fp.plan, val)
      Core.setfield!(fp,v,ptr)

   # edit source coefficients
	 elseif v == :alpha
	    	if typeof(val) != Vector{ComplexF64}
		  	  error("alpha has to be a ComplexF64 vector.")
		    end
		    if (size(val)[1]) != fp.N
			    error("alpha has to be a ComplexF64 vector of length N.")
		    end
			ptr = ccall(("jfastsum_set_alpha", lib_path), Ptr{ComplexF64}, (Ref{fastsum_plan},Ref{ComplexF64}), fp.plan, val)
      Core.setfield!(fp,v,ptr)

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
    	Core.setfield!(fp,v,val)
  	end
end # Base.setproperty!

# overwrite dot notation for plan struct in order to use C memory
function Base.getproperty(fp::fastsumplan{D},v::Symbol) where {D}
	if v == :x
		if !isdefined(fp,:x)
			error("x is not set.")
		end
		ptr = Core.getfield(fp,:x)
		if D==1
			return unsafe_wrap(Vector{Float64},ptr,fp.N)             # get source nodes from C memory and convert to Julia type
		else
			return unsafe_wrap(Matrix{Float64},ptr,(D,Int64(fp.N)))  # get source nodes from C memory and convert to Julia type
		end

	elseif v == :y
		if !isdefined(fp,:y)
			error("y is not set.")
		end
		ptr = Core.getfield(fp,:y)
		if D==1
			return unsafe_wrap(Vector{Float64},ptr,fp.M)             # get target nodes from C memory and convert to Julia type
		else
			return unsafe_wrap(Matrix{Float64},ptr,(D,Int64(fp.M)))
		end 

	elseif v == :alpha
		if !isdefined(fp,:alpha)
			error("alpha is not set.")
		end
		ptr = Core.getfield(fp,:alpha)
		return unsafe_wrap(Vector{ComplexF64},ptr,fp.N)             # get coefficients from C memory and convert to Julia type

	elseif v == :f
		if !isdefined(fp,:f)
			error("f is not set.")
		end
		ptr = Core.getfield(fp,:f)
		return unsafe_wrap(Vector{ComplexF64},ptr,fp.M)  # get function values from C memory and convert to Julia type
	elseif v == :c
		if !isdefined(fp,:c)
			error("c is not set.")
		end
		c_temp = Core.getfield(fp,:c)
		return c_temp #unsafe_wrap(Vector{Float64},ptr,D)
	else
		return Core.getfield(fp,v)
	end
end # Base.getproperty

function trafo(fp::fastsumplan{D}) where {D}
  if fp.finalized
    error("Plan already finalized.")
  end
  if !isdefined(fp, :x)
    error("x has not been set.")
  end
  if !isdefined(fp, :y)
    error("y has not been set.")
  end
  if !isdefined(fp, :alpha)
    error("alpha has not been set.")
	end
	ptr = ccall(("jfastsum_trafo", lib_path), Ptr{ComplexF64}, (Ref{fastsum_plan},), fp.plan)
  Core.setfield!(fp,:f,ptr)
end #trafo

function trafo_exact(fp::fastsumplan{D}) where {D}
  if fp.finalized
    error("Plan already finalized.")
  end
  if !isdefined(fp, :x)
    error("x has not been set.")
  end
  if !isdefined(fp, :y)
    error("y has not been set.")
  end
  if !isdefined(fp, :alpha)
    error("alpha has not been set.")
	end
	ptr = ccall(("jfastsum_exact", lib_path), Ptr{ComplexF64}, (Ref{fastsum_plan},), fp.plan)
  Core.setfield!(fp,:f,ptr)
end #trafo
end #module
