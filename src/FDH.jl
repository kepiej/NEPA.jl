# Define shorthand functions depending on the returns to scale
FDH_CRS(X,Y,input) = FDH{CRS}(X,Y,input)
FDH_VRS(X,Y,input) = FDH{VRS}(X,Y,input)
FDH_NIRS(X,Y,input) = FDH{NIRS}(X,Y,input)
FDH_NDRS(X,Y,input) = FDH{NDRS}(X,Y,input)

immutable FDH{T<:RS} <: AbstractDEA{FreeDisposal,T}
  Data::DEAData
  input::Bool

  function FDH{T}(X::Array,Y::Array,input::Bool) where T<:RS
    new(DEAData(X,Y),input)
  end
end

function getdata(F::FDH)
  return F.Data
end

# Solve FDH program under VRS with (Xk,Yk) as evaluation point
function (F::FDH{VRS})(Xk::Array,Yk::Array)
  X,Y = F[1:end]
  dompeer = NaN
  dom = find(all(Y .>= Yk,2) & all(X .<= Xk,2))
  if(F.input) #Input-oriented
    theta = Inf
    for k in dom
      Ixk = find(X[k,:] .> 0.0)
      curmin = maximum(X[k,Ixk]./Xk[Ixk])
      if(curmin < theta)
        theta = curmin
        dompeer = k
      end
    end
  else #Output-oriented
    theta = -Inf
    for k in dom
      Iyk = find(Y[k,:] .> 0.0)
      curmax = minimum(Y[k,Iyk]./Yk[Iyk])
      if(curmax > theta)
        theta = curmax
        dompeer = k
      end
    end
  end
  return DEAResult(theta,[],[],[k == dompeer ? 1.0 : 0.0 for k=1:size(X,1)])
end

# Solve FDH program under NDRS with (Xk,Yk) as evaluation point
function (F::FDH{NDRS})(Xk::Array,Yk::Array)
  X,Y = F[1:end]
  K = size(X,1)
  dompeer = NaN
  domx = find(all(X .<= Xk,2))
  if(F.input) #Input-oriented
    theta = Inf
    for k in domx
      Ixk = find(X[k,:] .> 0.0)
      Jyk = find(Y[k,:] .> 0.0)
      curmin = max(maximum(Yk[Jyk]./Y[k,Jyk]),1.0)*maximum(X[k,Ixk]./Xk[Ixk])
      if(curmin < theta)
        theta = curmin
        dompeer = k
      end
    end
  else #Output-oriented
    theta = -Inf
    for k in domx
      Ixk = find(X[k,:] .> 0.0)
      Jyk = find(Y[k,:] .> 0.0)
      curmax = minimum(Y[k,Jyk]./Yk[Jyk])*minimum(Xk[Ixk]./X[k,Ixk])
      if(curmax > theta)
        theta = curmax
        dompeer = k
      end
    end
  end
  return DEAResult(theta,[],[],[k == dompeer ? 1.0 : 0.0 for k=1:K])
end

# Solve FDH program under CRS with (Xk,Yk) as evaluation point
function (F::FDH{CRS})(Xk::Array,Yk::Array)
  X,Y = F[1:end]
  K = size(X,1)
  dompeer = NaN
  dom = find(maximum(reshape(Yk,1,:)./Y,2) .<= minimum(reshape(Xk,1,:)./X,2))
  if(F.input) #Input-oriented
    theta = Inf
    for k in dom
      Ixk = find(X[k,:] .> 0.0)
      Jyk = find(Y[k,:] .> 0.0)
      curmin = maximum(Yk[Jyk]./Y[k,Jyk])*maximum(X[k,Ixk]./Xk[Ixk])
      if(curmin < theta)
        theta = curmin
        dompeer = k
      end
    end
  else #Output-oriented
    theta = -Inf
    for k in dom
      Ixk = find(X[k,:] .> 0.0)
      Jyk = find(Y[k,:] .> 0.0)
      curmax = minimum(Y[k,Jyk]./Yk[Jyk])*minimum(Xk[Ixk]./X[k,Ixk])
      if(curmax > theta)
        theta = curmax
        dompeer = k
      end
    end
  end
  return DEAResult(theta,[],[],[k == dompeer ? 1.0 : 0.0 for k=1:K])
end

# Solve FDH program under NIRS with (Xk,Yk) as evaluation point
function (F::FDH{NIRS})(Xk::Array,Yk::Array)
  X,Y = F[1:end]
  K = size(X,1)
  dompeer = NaN
  domy = find(all(Y .>= Yk,2))
  if(F.input) #Input-oriented
    theta = Inf
    for k in domy
      Ixk = X[k,:] .> 0.0
      Jyk = Y[k,:] .> 0.0
      curmin = maximum(Yk[:,Jyk]./Y[k,Jyk])*maximum(X[k,Ixk]./Xk[:,Ixk])
      if(curmin < theta)
        theta = curmin
        dompeer = k
      end
    end
  else #Output-oriented
    theta = -Inf
    for k in domy
      Ixk = X[k,:] .> 0.0
      Jyk = Y[k,:] .> 0.0
      curmax = minimum(Y[k,Jyk]./Yk[:,Jyk])*min(minimum(Xk[:,Ixk]./X[k,Ixk]),1.0)
      if(curmax > theta)
        theta = curmax
        dompeer = k
      end
    end
  end
  return DEAResult(theta,[],[],[k == dompeer ? 1.0 : 0.0 for k=1:K])
end

function (F::FDH{T})() where T<:RS
  	Data = getdata(F)
  	res = Array{DEAResult}(getnrdmu(Data))

  	@sync @parallel for k in eachindex(Data)
  		res[k] = F(Data[k]...)
  	end

  	return res
end
