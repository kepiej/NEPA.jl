
# Compute order-m efficiency score
function OrderM{T <: AbstractDEA}(A::T,M::Int,B::Int)
  Data = getData(A)
  if(M > getNrDMU(Data))
    error("M cannot be larger than the number of DMUs!")
  end
  theta = Array(Float64,getNrDMU(Data),B)
  @sync @parallel for k=1:B
    Am = T(Data[rand(1:getNrDMU(Data),M)]...)
    for j in eachindex(Data)
      theta[j,k] = Am(Data[j]...)
    end
  end
  return mean(theta,2)
end

# Compute order-m efficiency score for FDH
function OrderM{T<:RS}(A::FDH{T},M::Int,B::Int)
  Data = getData(A)
  if(M > getNrDMU(Data))
    error("M cannot be larger than the number of DMUs!")
  end
  X,Y = Data[1:end]
  theta = Array(Float64,getNrDMU(Data),B)
  for j in eachindex(Data)
    # Condition on Y >= y for input-oriented and X <= x for output-oriented
    Xk,Yk = Data[j]
    if(A.input)
      dom = find(Y .>= Yk)
    else
      dom = find(X .<= Xk)
    end
    @sync @parallel for k=1:B
      Am = FDH{T}(Data[rand(dom,M)]...,A.input)
      theta[j,k] = Am(Xk,Yk)
    end
  end
  return mean(theta,2)
end
