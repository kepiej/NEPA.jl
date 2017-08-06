# Hicks-Moorsteen-Bjurek (HMB) total factor productivity index as defined in Bjurek (SJE, 1996)
immutable HMB{S<:AbstractDataEnvelopment,T<:RS}
  In0::Union{FDH,DEA}
  In1::Union{FDH,DEA}
  Out0::Union{FDH,DEA}
  Out1::Union{FDH,DEA}
  K::Int64

  function HMB{S,T}(X0::Array,Y0::Array,X1::Array,Y1::Array) where {S<:AbstractDataEnvelopment,T<:RS}
    if(size(X0) != size(X1) || size(Y0) != size(Y1))
      error("The number of DMUs must be the same in both periods!")
    end
    if(eltype(HMB{S,T})[1] == Convex)
      new(DEA{T}(X0,Y0,true),DEA{T}(X1,Y1,true),DEA{T}(X0,Y0,false),DEA{T}(X1,Y1,false),size(X0,1))
    elseif(eltype(HMB{S,T})[1] == FreeDisposal)
      new(FDH{T}(X0,Y0,true),FDH{T}(X1,Y1,true),FDH{T}(X0,Y0,false),FDH{T}(X1,Y1,false),size(X0,1))
    else
      # Not defined yet
    end
  end

end

# Return parametric types as a tuple
eltype{S,T}(::Type{HMB{S,T}}) = (S, T)

# Compute the HMB index
function (H::HMB)()
  MO0 = Array{Float64}(H.K)
  MO1 = Array{Float64}(H.K)
  MI0 = Array{Float64}(H.K)
  MI1 = Array{Float64}(H.K)
  for k=1:H.K
    X0,Y0 = getdata(H.Out0)[k]
    X1,Y1 = getdata(H.Out1)[k]

    MO0[k] = geteff(H.Out0(X0,Y1))/geteff(H.Out0(X0,Y0))
    MO1[k] = geteff(H.Out1(X1,Y1))/geteff(H.Out1(X1,Y0))

    X0,Y0 = getdata(H.In0)[k]
    X1,Y1 = getdata(H.In1)[k]

    MI0[k] = geteff(H.In0(X0,Y0))/geteff(H.In0(X1,Y0))
    MI1[k] = geteff(H.In1(X0,Y1))/geteff(H.In1(X1,Y1))
  end
  return sqrt((MO0./MI0).*(MO1./MI1))
end

# Technical efficiency change
function TEI(H::HMB)
  return geteff(H.Out1())./geteff(H.Out0())
end

# Technical change
function TC(H::HMB)
  TC0 = Array{Float64}(H.K)
  TC1 = Array{Float64}(H.K)
  for k=1:H.K
    X0,Y0 = getdata(H.Out0)[k]
    TC0[k] = geteff(H.Out1(X0,Y0))
    X1,Y1 = getdata(H.Out1)[k]
    TC1[k] = geteff(H.Out0(X1,Y1))
  end
  TC0 = geteff(H.Out0())./TC0
  TC1 = TC1./geteff(H.Out1())
  return real(sqrt(complex(TC0.*TC1)))
end

# Scale efficiency change
function SEC(H::HMB)
  return H()./(TC(H).*TEI(H))
end
