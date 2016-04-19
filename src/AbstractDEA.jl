abstract AbstractDataEnvelopment

# Abstract DEA model
abstract AbstractDEA <: AbstractArray{Any,1}

function getdata(A::AbstractDEA)
	# Concrete types should overload this method to get array-like indexing features
end

Base.size(A::AbstractDEA) = size(getdata(A))
Base.linearindexing(::Type{AbstractDEA}) = Base.LinearFast()

function Base.getindex(A::AbstractDEA, i::Int)
	1 <= i <= size(A,1) || throw(BoundsError(A, i))
	return getdata(A)[i]
end

function Base.call(A::AbstractDEA)
	Data = getdata(A)
	beta = Array(Float64,getnrdmu(Data))

	@sync @parallel for k in eachindex(Data)
		beta[k] = Base.call(A,Data[k]...)
	end

	return beta
end
