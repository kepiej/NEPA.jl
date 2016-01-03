# Hicks-Moorsteen-Bjurek (HMB) total factor productivity index as defined in Bjurek (SJE, 1996)
immutable HMB
  In0::DDF
  In1::DDF
  Out0::DDF
  Out1::DDF
  K::Int64

  function HMB(X0::Array,Y0::Array,X1::Array,Y1::Array)
    if(size(X0) != size(X1) || size(Y0) != size(Y1))
      error("The number of DMUs must be the same in both periods!")
    end
    new(DDF(X0,Y0,X0,zeros(size(Y0)),VRS()),DDF(X1,Y1,X1,zeros(size(Y1)),VRS()),DDF(X0,Y0,zeros(size(X0)),Y0,VRS()),DDF(X1,Y1,zeros(size(X1)),Y1,VRS()),size(X0,1))
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
    MO0[k] = H.Out0(H.Out0.X[k,:],H.Out1.Y[k,:],zeros(size(H.Out0.X[k,:])),H.Out1.Y[k,:])/H.Out0(H.Out0.X[k,:],H.Out0.Y[k,:],zeros(size(H.Out0.X[k,:])),H.Out0.Y[k,:])
    MO1[k] = H.Out1(H.Out1.X[k,:],H.Out1.Y[k,:],zeros(size(H.Out1.X[k,:])),H.Out1.Y[k,:])/H.Out1(H.Out1.X[k,:],H.Out0.Y[k,:],zeros(size(H.Out1.X[k,:])),H.Out0.Y[k,:])

    MI0[k] = H.In0(H.In0.X[k,:],H.In0.Y[k,:],H.In0.X[k,:],zeros(size(H.In0.Y[k,:])))/H.In0(H.In1.X[k,:],H.In0.Y[k,:],H.In1.X[k,:],zeros(size(H.In0.Y[k,:])))
    MI1[k] = H.In1(H.In0.X[k,:],H.In1.Y[k,:],H.In0.X[k,:],zeros(size(H.In1.Y[k,:])))/H.In1(H.In1.X[k,:],H.In1.Y[k,:],H.In1.X[k,:],zeros(size(H.In1.Y[k,:])))
  end
  return sqrt((MO0./MI0).*(MO1./MI1))
end
