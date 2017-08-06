# Define shorthand functions depending on the returns to scale
DEA_CRS(X,Y,input) = DEA{CRS}(X,Y,input)
DEA_VRS(X,Y,input) = DEA{VRS}(X,Y,input)
DEA_NIRS(X,Y,input) = DEA{NIRS}(X,Y,input)
DEA_NDRS(X,Y,input) = DEA{NDRS}(X,Y,input)

# DEA is a special case of the DDF where gx or gy is zero depending on the orientation
immutable DEA{T<:RS} <: AbstractDEA{Convex,T}
  D::DDF
  input::Bool

  function DEA{T}(X::Array,Y::Array,input::Bool) where T<:RS
    input ? new(DDF{Convex,T}(X,Y,X,zeros(size(Y))),input) : new(DDF{Convex,T}(X,Y,zeros(size(X)),Y),input)
  end
end

function getdata(DMU::DEA)
	return getdata(DMU.D)
end

function (DMU::DEA)(Xk::Array,Yk::Array)
  if DMU.input
    res = DMU.D(Xk,Yk,Xk,zeros(size(Yk)))
    res.eff = 1 - res.eff
  else
    res = DMU.D(Xk,Yk,zeros(size(Xk)),Yk)
    res.eff = 1 + res.eff
  end
  res.sx = []
  res.sy = []
  return res
end

function (DMU::DEA)()
  res = DMU.D()
  for k in eachindex(res)
    if DMU.input
      res[k].eff = 1 - res[k].eff
    else
      res[k].eff = 1 + res[k].eff
    end
    res[k].sx = []
    res[k].sy = []
  end
  return res
end
