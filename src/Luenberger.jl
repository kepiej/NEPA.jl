# Implements the Luenberger productivity growth indicator
immutable Luenberger{S<:AbstractDataEnvelopment,T<:RS}
  D0::DDF
  D1::DDF
  K::Int64

  function Luenberger{S,T}(X0::Array,Y0::Array,gx0::Array,gy0::Array,X1::Array,Y1::Array,gx1::Array,gy1::Array) where {S<:AbstractDataEnvelopment,T<:RS}
    if(size(X0) != size(X1) || size(Y0) != size(Y1))
      error("The number of DMUs must be the same in both periods!")
    end
    new(DDF{S,T}(X0,Y0,gx0,gy0),DDF{S,T}(X1,Y1,gx1,gy1),size(X0,1))
  end

  function Luenberger{S,T}(DDF0::DDF,DDF1::DDF) where {S<:AbstractDataEnvelopment,T<:RS}
    new(DDF0,DDF1)
  end
end

function (L::Luenberger)()
  val = Array{Float64}(L.K)
  for k=1:L.K
    val[k] = -geteff(L.D0(L.D1[k]...)) + geteff(L.D1(L.D0[k]...))
  end
  return (val + geteff(L.D0()) - geteff(L.D1()))./2
end

function TEI(L::Luenberger)
  return geteff(L.D0()) - geteff(L.D1())
end

function TC(L::Luenberger)
  TC0 = Array{Float64}(L.K)
  TC1 = Array{Float64}(L.K)
  for k=1:L.K
      # Compute technical change (TC0) from period 0 to period 1 using observations at time 0
      TC0[k] = geteff(L.D1(L.D0[k]...)) - geteff(L.D0(L.D0[k]...))
      # Compute technical change (TC1) from period 0 to period 1 using observations at time 1
      TC1[k] = geteff(L.D1(L.D1[k]...)) - geteff(L.D0(L.D1[k]...))
  end

  # Technical change (TC) is an arithmetic average of TC0 and TC1
  return (TC0+TC1)./2
end
