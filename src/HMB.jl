# Hicks-Moorsteen-Bjurek (HMB) total factor productivity index as defined in Bjurek (SJE, 1996)
immutable HMB
  In0::DEA
  In1::DEA
  Out0::DEA
  Out1::DEA
  K::Int64

  function HMB(X0::Array,Y0::Array,X1::Array,Y1::Array)
    if(size(X0) != size(X1) || size(Y0) != size(Y1))
      error("The number of DMUs must be the same in both periods!")
    end
    new(DEA(X0,Y0,true,VRS()),DEA(X1,Y1,true,VRS()),DEA(X0,Y0,false,VRS()),DEA(X1,Y1,false,VRS()),size(X0,1))
  end
end

# Compute the HMB index
function Base.call(H::HMB)
  TFP = Array(Float64,H.K)
  MO0 = Array(Float64,H.K)
  MO1 = Array(Float64,H.K)
  MI0 = Array(Float64,H.K)
  MI1 = Array(Float64,H.K)
  for k=1:H.K
    MO0[k] = H.Out0(H.Out0.D.X[k,:],H.Out1.D.Y[k,:])/H.Out0(H.Out0.D.X[k,:],H.Out0.D.Y[k,:])
    MO1[k] = H.Out1(H.Out1.D.X[k,:],H.Out1.D.Y[k,:])/H.Out1(H.Out1.D.X[k,:],H.Out0.D.Y[k,:])

    MI0[k] = H.In0(H.In0.D.X[k,:],H.In0.D.Y[k,:])/H.In0(H.In1.D.X[k,:],H.In0.D.Y[k,:])
    MI1[k] = H.In1(H.In0.D.X[k,:],H.In1.D.Y[k,:])/H.In1(H.In1.D.X[k,:],H.In1.D.Y[k,:])
  end
  return sqrt((MO0./MI0).*(MO1./MI1))
end

# Technical efficiency change
function TEI(H::HMB)
  return H.Out1()./H.Out0()
end

# Technical change
function TC(H::HMB)
  TC0 = Array(Float64,H.K)
  TC1 = Array(Float64,H.K)
  for k=1:H.K
    TC0[k] = H.Out1(H.Out0.D.X[k,:],H.Out0.D.Y[k,:])
    TC1[k] = H.Out0(H.Out1.D.X[k,:],H.Out1.D.Y[k,:])
  end
  TC0 = H.Out0()./TC0
  TC1 = TC1./H.Out1()
  return sqrt(TC0.*TC1)
end

# Returns to scale measure
function SEC(H::HMB)
  return H()./(TC(H).*TEI(H))
end
