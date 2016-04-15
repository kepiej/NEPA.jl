# Define shorthand functions depending on the returns to scale
FDH_CRS(X,Y,input) = FDH{CRS}(X,Y,input)
FDH_VRS(X,Y,input) = FDH{VRS}(X,Y,input)
FDH_NIRS(X,Y,input) = FDH{NIRS}(X,Y,input)
FDH_NDRS(X,Y,input) = FDH{NDRS}(X,Y,input)

immutable FDH{T<:RS} <: AbstractDEA
  Data::DEAData
  input::Bool

  function FDH(X::Array,Y::Array,input::Bool)
    new(DEAData(X,Y),input)
  end

  function FDH(Data::DEAData)
    new(Data)
  end
end

function getData(F::FDH)
  return F.Data
end

# Solve FDH program under VRS with (Xk,Yk) as evaluation point
function Base.call(F::FDH{VRS},Xk::Array,Yk::Array)
  X,Y = getData(F)[1:end]
  if(F.input) #Input-oriented
    theta = Inf
    domy = find(all(Y .>= Yk,2))
    for k in domy
      Ixk = X[k,:] .> 0.0
      curmin = maximum(X[k,Ixk]./Xk[:,Ixk])
      if(curmin < theta)
        theta = curmin
      end
    end
  else #Output-oriented
    theta = -Inf
    domx = find(all(X .<= Xk,2))
    for k in domx
      Iyk = Y[k,:] .> 0.0
      curmax = minimum(Y[k,Iyk]./Yk[:,Iyk])
      if(curmax > theta)
        theta = curmax
      end
    end
  end
  return theta
end

# Solve FDH program under NDRS with (Xk,Yk) as evaluation point
function Base.call(F::FDH{NDRS},Xk::Array,Yk::Array)
  X,Y = getData(F)[1:end]
  K = size(X,1)
  if(F.input) #Input-oriented
    theta = Inf
    for k=1:K
      Ixk = X[k,:] .> 0.0
      Jyk = Y[k,:] .> 0.0
      curmin = max(maximum(Yk[:,Jyk]./Y[k,Jyk]),1.0)*maximum(X[k,Ixk]./Xk[:,Ixk])
      if(curmin < theta)
        theta = curmin
      end
    end
  else #Output-oriented
    theta = -Inf
    domx = find(all(X .<= Xk,2))
    for k in domx
      Ixk = X[k,:] .> 0.0
      Jyk = Y[k,:] .> 0.0
      curmax = minimum(Y[k,Jyk]./Yk[:,Jyk])*minimum(Xk[:,Ixk]./X[k,Ixk])
      if(curmax > theta)
        theta = curmax
      end
    end
  end
  return theta
end

# Solve FDH program under CRS with (Xk,Yk) as evaluation point
function Base.call(F::FDH{CRS},Xk::Array,Yk::Array)
  X,Y = getData(F)[1:end]
  K = size(X,1)
  if(F.input) #Input-oriented
    theta = Inf
    for k=1:K
      Ixk = X[k,:] .> 0.0
      Jyk = Y[k,:] .> 0.0
      curmin = maximum(Yk[:,Jyk]./Y[k,Jyk])*maximum(X[k,Ixk]./Xk[:,Ixk])
      if(curmin < theta)
        theta = curmin
      end
    end
  else #Output-oriented
    theta = -Inf
    for k=1:K
      Ixk = X[k,:] .> 0.0
      Jyk = Y[k,:] .> 0.0
      curmax = minimum(Y[k,Jyk]./Yk[:,Jyk])*minimum(Xk[:,Ixk]./X[k,Ixk])
      if(curmax > theta)
        theta = curmax
      end
    end
  end
  return theta
end

# Solve FDH program under NIRS with (Xk,Yk) as evaluation point
function Base.call(F::FDH{NIRS},Xk::Array,Yk::Array)
  X,Y = getData(F)[1:end]
  K = size(X,1)
  if(F.input) #Input-oriented
    theta = Inf
    domy = find(all(Y .>= Yk,2))
    for k in domy
      Ixk = X[k,:] .> 0.0
      Jyk = Y[k,:] .> 0.0
      curmin = maximum(Yk[:,Jyk]./Y[k,Jyk])*maximum(X[k,Ixk]./Xk[:,Ixk])
      if(curmin < theta)
        theta = curmin
      end
    end
  else #Output-oriented
    theta = -Inf
    for k=1:K
      Ixk = X[k,:] .> 0.0
      Jyk = Y[k,:] .> 0.0
      curmax = minimum(Y[k,Jyk]./Yk[:,Jyk])*min(minimum(Xk[:,Ixk]./X[k,Ixk]),1.0)
      if(curmax > theta)
        theta = curmax
      end
    end
  end
  return theta
end
