# Implements the Luenberger productivity growth indicator
immutable Luenberger{T<:Tuple{AbstractDataEnvelopment,RS}}
  D0::DDF
  D1::DDF
  K::Int64

  function Luenberger(X0::Array,Y0::Array,gx0::Array,gy0::Array,X1::Array,Y1::Array,gx1::Array,gy1::Array)
    if(size(X0) != size(X1) || size(Y0) != size(Y1))
      error("The number of DMUs must be the same in both periods!")
    end
    new(DDF{T}(X0,Y0,gx0,gy0),DDF{T}(X1,Y1,gx1,gy1),size(X0,1))
  end

  function Luenberger(DDF0::DDF,DDF1::DDF)
    new(DDF0,DDF1)
  end
end

function Base.call(L::Luenberger)
  val = Array(Float64,L.K)
  for k=1:L.K
    val[k] = -L.D0(L.D1[k]...) + L.D1(L.D0[k]...)
  end
  return (val + L.D0() - L.D1())./2
end

function TEI(L::Luenberger)
  return L.D0() - L.D1()
end

function TC(L::Luenberger)
  TC0 = Array(Float64,L.K)
  TC1 = Array(Float64,L.K)
  for k=1:L.K
      # Compute technical change (TC0) from period 0 to period 1 using observations at time 0
      TC0[k] = L.D1(L.D0[k]...) - L.D0(L.D0[k]...)
      # Compute technical change (TC1) from period 0 to period 1 using observations at time 1
      TC1[k] = L.D1(L.D1[k]...) - L.D0(L.D1[k]...)
  end

  # Technical change (TC) is an arithmetic average of TC0 and TC1
  return (TC0+TC1)./2
end