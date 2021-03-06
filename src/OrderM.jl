
# Compute order-m efficiency score
function orderm{T<:AbstractDEA}(A::T,M::Int,B::Int)
  Data = getdata(A)
  K = getnrdmu(Data)
  # getdata() returns a shallow copy of the data in A, but we need a deep copy!
  Datacopy = deepcopy(Data)
  if(M > K)
    error("M cannot be larger than the number of DMUs!")
  end
  theta = Array(Float64,K,B)
  @sync @parallel for i=1:B
    setindexes!(Data,rand(1:K,M))
    for j in eachindex(Datacopy)
      theta[j,i] = geteff(A(Datacopy[j]...))
    end
  end
  setindexes!(Data,collect(1:K))
  return mean(theta,2)
end

# Compute order-m efficiency score for FDH
function orderm{T<:RS}(A::FDH{T},M::Int,B::Int)
  Data = getdata(A)
  K = getnrdmu(Data)
  # getdata() returns a shallow copy of the data in A, but we need a deep copy!
  Datacopy = deepcopy(Data)
  if(M > K)
    error("M cannot be larger than the number of DMUs!")
  end
  X,Y = Data[1:end]
  theta = Array(Float64,K,B)
  for j in eachindex(Datacopy)
    # Condition on Y >= y for input-oriented and X <= x for output-oriented
    Xk,Yk = Datacopy[j]
    if(A.input)
      dom = find(Y .>= Yk)
    else
      dom = find(X .<= Xk)
    end
    @sync @parallel for i=1:B
      setindexes!(Data,rand(1:K,M))
      theta[j,i] = geteff(A(Xk,Yk))
    end
  end
  setindexes!(Data,collect(1:K))
  return mean(theta,2)
end
