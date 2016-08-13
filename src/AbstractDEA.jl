abstract AbstractDataEnvelopment

# Abstract DEA model
abstract AbstractDEA{S<:AbstractDataEnvelopment,T<:RS} <: AbstractArray{Any,1}

function getdata(A::AbstractDEA)
	# Concrete types should overload this method to get array-like indexing features
end

# Return parametric types as a tuple
Base.eltype{S,T}(::Type{AbstractDEA{S,T}}) = (S, T)

Base.size(A::AbstractDEA) = size(getdata(A))
Base.linearindexing(::Type{AbstractDEA}) = Base.LinearFast()

function Base.getindex(A::AbstractDEA, i::Int)
	1 <= i <= size(A,1) || throw(BoundsError(A, i))
	return getdata(A)[i]
end

function Base.call(A::AbstractDEA)
	Data = getdata(A)
	res = Array(DEAResult,getnrdmu(Data))

	@sync @parallel for k in eachindex(Data)
		res[k] = Base.call(A,Data[k]...)
	end

	return res
end
