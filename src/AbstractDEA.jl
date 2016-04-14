abstract AbstractDataEnvelopment

# Abstract DEA model
abstract AbstractDEA <: AbstractArray{Any,1}

function getData(A::AbstractDEA)
	# Concrete types should overload this method to get array-like indexing features
end

Base.size(A::AbstractDEA) = size(getData(A))
Base.linearindexing(::Type{AbstractDEA}) = Base.LinearFast()

function Base.getindex(A::AbstractDEA, i::Int)
	1 <= i <= size(A,1) || throw(BoundsError(A, i))
	return getData(A)[i]
end

function Base.call(A::AbstractDEA)
	Data = getData(A)
	beta = Array(Float64,getNrDMU(Data))

	@sync @parallel for k in eachindex(Data)
		beta[k] = Base.call(A,Data[k]...)
	end

	return beta
end
