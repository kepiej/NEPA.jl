# Implements the Luenberger productivity growth indicator
immutable Luenberger
  D0::DDF
  D1::DDF
  K::Int64

  function Luenberger(X0::Array,Y0::Array,gx0::Array,gy0::Array,X1::Array,Y1::Array,gx1::Array,gy1::Array)
    if(size(X0) != size(X1) || size(Y0) != size(Y1))
      error("The number of DMUs must be the same in both periods!")
    end
    new(DDF(X0,Y0,gx0,gy0,VRS()),DDF(X1,Y1,gx1,gy1,VRS()),size(X0,1))
  end
end

function Base.call(L::Luenberger)
  val = Array(Float64,L.K)
  for k=1:L.K
    val[k] = -L.D0(L.D1.X[k,:],L.D1.Y[k,:],L.D1.gx[k,:],L.D1.gy[k,:]) + L.D1(L.D0.X[k,:],L.D0.Y[k,:],L.D0.gx[k,:],L.D0.gy[k,:])
  end
  return (val + L.D0() - L.D1())./2
end

function TEI(L::Luenberger)
  return L.D0() - L.D1()
end

function TC(L::Luenberger)
  # Compute technical change (TC0) from period 0 to period 1 using observations at time 0
  TC0 = Array(Float64,L.K)
  for k=1:L.K
      TC0[k] = L.D1(L.D0.X[k,:],L.D0.Y[k,:],L.D0.gx[k,:],L.D0.gy[k,:]) - L.D0(L.D0.X[k,:],L.D0.Y[k,:],L.D0.gx[k,:],L.D0.gy[k,:])
  end

  # Compute technical change (TC1) from period 0 to period 1 using observations at time 1
  TC1 = Array(Float64,L.K)
  for k=1:L.K
      TC1[k] = L.D1(L.D1.X[k,:],L.D1.Y[k,:],L.D1.gx[k,:],L.D1.gy[k,:]) - L.D0(L.D1.X[k,:],L.D1.Y[k,:],L.D1.gx[k,:],L.D1.gy[k,:])
  end

  # Technical change (TC) is an arithmetic average of TC0 and TC1
  return (TC0+TC1)./2
end
