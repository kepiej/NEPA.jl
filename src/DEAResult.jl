import Base.show

type DEAResult
  eff::Float64
  sx::Array
  sy::Array
  peers::Array

  function DEAResult(eff::Float64 = NaN,sx::Array = [],sy::Array = [],peers::Array = [])
    new(eff,sx,sy,peers)
  end
end

function geteff(res::DEAResult)
  return res.eff
end

function geteff(res::Array{DEAResult})
  eff = Array{Float64}(length(res))
  for k in eachindex(res)
    eff[k] = geteff(res[k])
  end
  return eff
end

function getpeers(res::DEAResult)
  return res.peers
end

function getpeers(res::Array{DEAResult})
  peers = Array{Float64}(length(res),length(res))
  for k in eachindex(res)
    peers[k,:] = getpeers(res[k])
  end
  return peers
end

function getsx(res::DEAResult)
  return res.sx
end

function getsx(res::Array{DEAResult})
  sx = Array{Float64}(length(res),length(getsx(res[1])))
  for k in eachindex(res)
    sx[k,:] = getsx(res[k])
  end
  return sx
end

function getsy(res::DEAResult)
  return res.sy
end

function getsy(res::Array{DEAResult})
  sy = Array{Float64}(length(res),length(getsy(res[1])))
  for k in eachindex(res)
    sy[k,:] = getsy(res[k])
  end
  return sy
end

function Base.show(io::IO,res::DEAResult)
  write(io,"eff = $(geteff(res)), sx = $(getsx(res)), sy = $(getsy(res))\n, peers = $(getpeers(res))")
end
