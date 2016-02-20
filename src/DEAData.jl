immutable DEAData <: AbstractArray{Any,1}
  X::Array
  Y::Array
  gx::Array
  gy::Array
  indexes::Vector{Int}

  function DEAData(X::Array,Y::Array,gx::Array,gy::Array,indexes::Vector{Int})
    if(size(X,1) != size(Y,1))
			error("Number of DMUs in X and Y are not the same: X has $(size(X,1)) DMUs and Y has $(size(Y,1)) DMUs!")
		end
		if(size(gx) != size(X) || size(gy) != size(Y))
			error("Dimensions of gx and X or gy and Y are not the same: size(gx) = $(size(gx)), size(X) = $(size(X)), size(gy) = $(size(gy)) and size(Y) = $(size(Y)) DMUs!")
		end
		if(gx == zeros(size(gx)) && gy == zeros(size(gy)))
			error("Direction vectors cannot both be zero!")
		end
    new(X,Y,gx,gy,indexes)
  end
end

DEAData(X::Array,Y::Array,gx::Array,gy::Array) = DEAData(X,Y,gx,gy,Vector{Int}(collect(1:size(X,1))))
DEAData(X,Y) = DEAData(X,Y,X,Y)

function getNrDMU(Data::DEAData)
  return size(Data.X,1)
end

function getIODim(Data::DEAData)
  return size(Data.X,2), size(Data.Y,2)
end

function getIndexes(Data::DEAData)
  return Data.indexes
end

Base.size(Data::DEAData) = (length(Data.indexes),)
Base.linearindexing(::Type{DEAData}) = Base.LinearFast()

function Base.getindex(Data::DEAData,i::Int)
  1 <= i <= length(Data.indexes) || throw(BoundsError(Data, i))
  return Data.X[Data.indexes[i],:],Data.Y[Data.indexes[i],:],Data.gx[Data.indexes[i],:],Data.gy[Data.indexes[i],:]
end

function Base.getindex(Data::DEAData, I)
  N,M = getIODim(Data)
  K = length(I)
  X = Array(eltype(Data.X),K,N)
  Y = Array(eltype(Data.Y),K,M)
  gx = Array(eltype(Data.gx),K,N)
  gy = Array(eltype(Data.gy),K,M)
  for i in eachindex(I)
    X[i,:],Y[i,:],gx[i,:],gy[i,:] = Data[I[i]]
  end
  return X,Y,gx,gy
end
