# Define shorthand functions depending on the returns to scale
DEA_CRS(X,Y,input) = DEA{CRS}(X,Y,input)
DEA_VRS(X,Y,input) = DEA{VRS}(X,Y,input)
DEA_NIRS(X,Y,input) = DEA{NIRS}(X,Y,input)
DEA_NDRS(X,Y,input) = DEA{NDRS}(X,Y,input)

# DEA is a special case of the DDF where gx or gy is zero depending on the orientation
immutable DEA{T<:RS} <: AbstractDEA
  D::DDF
  input::Bool

  function DEA(X,Y,input)
    input ? new(DDF{Tuple{Convex,T}}(X,Y,X,zeros(size(Y))),input) : new(DDF{Tuple{Convex,T}}(X,Y,zeros(size(X)),Y),input)
  end
end

function getdata(DMU::DEA)
	return getdata(DMU.D)
end

function Base.call(DMU::DEA,Xk::Array,Yk::Array)
  return DMU.input ? (1 - DMU.D(Xk,Yk,Xk,zeros(size(Yk)))) : (1 + DMU.D(Xk,Yk,zeros(size(Xk)),Yk))
end

function Base.call(DMU::DEA)
  return DMU.input ? (1 - DMU.D()) : (1 + DMU.D())
end
