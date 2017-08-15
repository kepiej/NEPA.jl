abstract type AbstractDataEnvelopment end

# Abstract DEA model
abstract type AbstractDEA{S<:AbstractDataEnvelopment,T<:RS} <: AbstractArray{Any,1} end

function getdata(A::AbstractDEA)
	# Concrete types should overload this method to get array-like indexing features
end

# Return parametric types as a tuple
Base.eltype{S,T}(::Type{AbstractDEA{S,T}}) = (S, T)

Base.size(A::AbstractDEA) = size(getdata(A))
Base.IndexStyle(::Type{AbstractDEA}) = IndexLinear()

function Base.getindex(A::AbstractDEA, I)
	all(1 .<= I .<= size(A,1)) || throw(BoundsError(A, I))
	return getdata(A)[I]
end

# Broken since Julia 0.5! Will it ever come back?
# function Base.call{T<:AbstractDEA}(A::T)
# 	Data = getdata(A)
# 	res = Array(DEAResult,getnrdmu(Data))
#
# 	@sync @parallel for k in eachindex(Data)
# 		res[k] = Base.call(A,Data[k]...)
# 	end
#
# 	return res
# end
