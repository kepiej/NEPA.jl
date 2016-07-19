type DEAData <: AbstractArray{Any,1}
  X::Array
  Y::Array
  gx::Array
  gy::Array
  indexes::Vector{Int}

  function DEAData(X::Array, Y::Array, gx::Array = [], gy::Array = []; indexes::Vector{Int} = Vector{Int}(collect(1:size(X,1))))
    if(size(X,1) != size(Y,1))
			error("Number of DMUs in X and Y are not the same: X has $(size(X,1)) DMUs and Y has $(size(Y,1)) DMUs!")
		end
    if(!isempty(gx) && size(gx) != size(X))
      error("Dimensions of gx and X are not the same: size(gx) = $(size(gx)), size(X) = $(size(X))!")
    end
    if(!isempty(gy) && size(gy) != size(Y))
      error("Dimensions of gy and Y are not the same: size(gy) = $(size(gy)) and size(Y) = $(size(Y))!")
    end
    if(isempty(gx) != isempty(gy))
      error("If gx is empty then gy must also be empty (and vice versa)!")
    end
    if(!isempty(gx) && !isempty(gy) && gx == zeros(size(gx)) && gy == zeros(size(gy)))
      error("Direction vectors cannot both be zero!")
    end
    new(X,Y,gx,gy,indexes)
  end
end

function getnrdmu(Data::DEAData)
  return size(Data,1)
end

function getiodim(Data::DEAData)
  return size(Data.X,2), size(Data.Y,2)
end

function getindexes(Data::DEAData)
  return copy(Data.indexes)
end

function setindexes!(Data::DEAData, newindexes::Vector{Int})
  if(length(newindexes) > size(Data.X,1))
    error("$(length(newindexes)) > $(size(Data.X,1))!")
  end
  Data.indexes = copy(newindexes)
end

Base.size(Data::DEAData) = (length(Data.indexes),)
Base.linearindexing(::Type{DEAData}) = Base.LinearFast()

function Base.getindex(Data::DEAData,i::Int)
  1 <= i <= size(Data,1) || throw(BoundsError(Data, i))
  if(!isempty(Data.gx))
    return Data.X[Data.indexes[i],:],Data.Y[Data.indexes[i],:],Data.gx[Data.indexes[i],:],Data.gy[Data.indexes[i],:]
  else
    return Data.X[Data.indexes[i],:],Data.Y[Data.indexes[i],:]
  end
end

function Base.getindex(Data::DEAData, I)
  N,M = getiodim(Data)
  K = length(I)
  X = Array(eltype(Data.X),K,N)
  Y = Array(eltype(Data.Y),K,M)
  if(!isempty(Data.gx))
    gx = Array(eltype(Data.gx),K,N)
    gy = Array(eltype(Data.gy),K,M)
    for i in eachindex(I)
      X[i,:],Y[i,:],gx[i,:],gy[i,:] = Data[I[i]]
    end
    return X,Y,gx,gy
  else
    for i in eachindex(I)
      X[i,:],Y[i,:] = Data[I[i]]
    end
    return X,Y
  end
end
