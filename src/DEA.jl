# Define shorthand functions depending on the returns to scale
DEA_CRS(X,Y,input) = DEA(X,Y,input,CRS())
DEA_VRS(X,Y,input) = DEA(X,Y,input,VRS())
DEA_NIRS(X,Y,input) = DEA(X,Y,input,NIRS())
DEA_NDRS(X,Y,input) = DEA(X,Y,input,NDRS())

# DEA is a special case of the DDF where gx or gy is zero depending on the orientation
immutable DEA <: AbstractDEA
  D::DDF
  input::Bool

  function DEA(X,Y,input,RStype::RS)
    input ? new(DDF(X,Y,X,zeros(size(Y)),RStype),input) : new(DDF(X,Y,zeros(size(X)),Y,RStype),input)
  end
end

function Base.call(DMU::DEA,Xk::Array,Yk::Array)
  return DMU.input ? DMU.D(Xk,Yk,Xk,zeros(size(Yk))) : DMU.D(Xk,Yk,zeros(size(Xk)),Yk)
end

function Base.call(DMU::DEA)
  return DMU.D()
end
