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
	c::Ref{Float64}         # kernel parameters
	eps_I::Float64          # inner boundary
	eps_B::Float64          # outer boundary
	init_done::Bool         # bool for plan init
	finalized::Bool	    	# bool for finalizer
	x::Ref{ComplexF64}      # source nodes
	y::Ref{ComplexF64}      # target nodes
    alpha::Ref{ComplexF64}   # source coefficients
	f::Ref{ComplexF64}      # target evaluations
	plan::Ref{fastsum_plan} # plan (C pointer)
	function fastsumplan{D}(N::Int32,M::Int32,n::Int32,m::Int32,p::Int32,kernel::String,c::Ref{Float64},eps_I::Float64,eps_B::Float64) where D
	# create plan object
	new(N,M,n,m,p,kernel,c,eps_I,eps_B,false,false)
	end

end #struct fastsumplan

function fastsum_init(fp::fastsumplan{D}) where {D}

  ptr = ccall(("jfastsum_alloc", lib_path), Ptr{fastsum_plan}, ())

  Core.setfield!(fp, :plan, ptr)
  # c noch in pointer umwandeln
  ccall(("jfastsum_init",lib_path),Nothing,(Ref{fastsum_plan},Int32,Int32,Int32,Cstring,Ref{Float64},Int32,Int32,Int32,Float64,Float64),
  ptr,D,fp.N,fp.M,kernel,c,fp.n,fp.m,fp.p,fp.eps_I,fp.eps_B)

  Core.setfield!(fp,:init_done, true)
end #fastsum_init

function Base.setproperty!(fp::fastsumplan{D},v::Symbol,val) where {D}
  if !fp.init_done
	  fastsum_init(fp)
  end

  if !fp.finalized
	  error("Plan already finalized")
  end

  if v == :M
	  @warn("You can't modify the number of target nodes.")
  elseif v == :N
	  @warn("You can't modify the number of source nodes.")
  elseif v == :n
	  @warn("You can't modify the expansion degree.")
  elseif v == :m
	  @warn("You can't modify the cut-off parameter.")
  elseif v == :p
	  @warn("You can't modify the degree of smoothness.")

  else
    Core.setfield!(fp,v,val)
  end


end # Base.setproperty!


end #module
