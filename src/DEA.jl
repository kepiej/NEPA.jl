# Define shorthand functions depending on the returns to scale
DEA_CRS(X,Y,input) = DEA{CRS}(X,Y,input)
DEA_VRS(X,Y,input) = DEA{VRS}(X,Y,input)
DEA_NIRS(X,Y,input) = DEA{NIRS}(X,Y,input)
DEA_NDRS(X,Y,input) = DEA{NDRS}(X,Y,input)

# DEA is a special case of the DDF where gx or gy is zero depending on the orientation
immutable DEA{T<:RS} <: AbstractDEA{Convex,T}
  D::DDF
  input::Bool

  function DEA(X,Y,input)
    input ? new(DDF{Convex,T}(X,Y,X,zeros(size(Y))),input) : new(DDF{Convex,T}(X,Y,zeros(size(X)),Y),input)
  end
end

function getdata(DMU::DEA)
	return getdata(DMU.D)
end

function Base.call(DMU::DEA,Xk::Array,Yk::Array)
  if DMU.input
    res = DMU.D(Xk,Yk,Xk,zeros(size(Yk)))
    res.eff = 1 - res.eff
  else
    res = DMU.D(Xk,Yk,zeros(size(Xk)),Yk)
    res.eff = 1 + res.eff
  end
  return res
end

function Base.call(DMU::DEA)
  res = DMU.D()
  for k in eachindex(res)
    if DMU.input
      res[k].eff = 1 - res[k].eff
    else
      res[k].eff = 1 + res[k].eff
    end
  end
  return res
end
