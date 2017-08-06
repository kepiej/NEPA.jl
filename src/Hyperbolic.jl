immutable Hyperbolic{S<:AbstractDataEnvelopment,T<:RS} <: AbstractDEA{S,T}
  Data::DEAData

  function Hyperbolic{S,T}(X,Y) where {S<:AbstractDataEnvelopment,T<:RS}
    new(DEAData(X,Y))
  end
end

function getdata(DMU::Hyperbolic)
  return DMU.Data
end

#Solve hyperbolic distance function under CRS with (Xk,Yk) as evaluation point
function (H::Hyperbolic{FreeDisposal,CRS})(Xk::Array,Yk::Array)
  X,Y = getdata(H)[1:end]
  K = size(X,1)
  dompeer = NaN
  dom = find(maximum(Yk./Y,2) .<= minimum(Xk./X,2))
  gamma = Inf
  for k in dom
    Ixk = X[k,:] .> 0.0
    Iyk = Y[k,:] .> 0.0
    domgamma = sqrt(maximum(Yk[:,Iyk]./Y[k,Iyk])./minimum(Xk[:,Ixk]./X[k,Ixk]))
    if(domgamma < gamma)
      gamma = domgamma
      dompeer = k
    end
  end
  return DEAResult(gamma,[],[],[k == dompeer ? 1.0 : 0.0 for k=1:K])
end

#T::Union{NIRS,NDRS}??
#Solve hyperbolic distance function under NIRS/NDRS with (Xk,Yk) as evaluation point
function (H::Hyperbolic{FreeDisposal,T}){T<:RS}(Xk::Array,Yk::Array)
  X,Y = getdata(H)[1:end]
  K = size(X,1)
  gamma = Inf
  dompeer = NaN
  dom = find(all(Y .>= Yk,2))
  for k in dom
    Ixk = X[k,:] .> 0.0
    Iyk = Y[k,:] .> 0.0
    alphak = sqrt(minimum(Xk[:,Ixk]./X[k,Ixk])*maximum(Yk[:,Iyk]./Y[k,Iyk]))
    if(T == NIRS)
      if(0.0 <= alphak <= 1)
        domgamma = sqrt(maximum(Yk[:,Iyk]./Y[k,Iyk])./minimum(Xk[:,Ixk]./X[k,Ixk]))
      else
        domgamma = max((1/minimum(Xk[:,Ixk]./X[k,Ixk])),maximum(Yk[:,Iyk]./Y[k,Iyk]))
      end
    elseif(T == NDRS)
      if(alphak >= 1)
        domgamma = sqrt(maximum(Yk[:,Iyk]./Y[k,Iyk])./minimum(Xk[:,Ixk]./X[k,Ixk]))
      else
        domgamma = max((1/minimum(Xk[:,Ixk]./X[k,Ixk])),maximum(Yk[:,Iyk]./Y[k,Iyk]))
      end
    end
    if(domgamma < gamma)
      gamma = domgamma
      dompeer = k
    end
  end
  return DEAResult(gamma,[],[],[k == dompeer ? 1.0 : 0.0 for k=1:K])
end

#Solve hyperbolic distance function under VRS with (Xk,Yk) as evaluation point
function (H::Hyperbolic{FreeDisposal,VRS})(Xk::Array,Yk::Array)
  X,Y = getdata(H)[1:end]

  dom = find(all(Y .>= Yk,2) .& all(X .<= Xk,2))
  gamma = Inf
  dompeer = NaN
  for k in dom
    Ixk = X[k,:] .> 0.0
    Iyk = Y[k,:] .> 0.0
    curmin = max(maximum(X[k,Ixk]./Xk[:,Ixk]),1/minimum(Y[k,Iyk]./Yk[:,Iyk]))
    if(curmin < gamma)
      gamma = curmin
      dompeer = k
    end
  end
  return DEAResult(gamma,[],[],[k == dompeer ? 1.0 : 0.0 for k=1:size(X,1)])
end
