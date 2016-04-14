# Define shorthand functions depending on the returns to scale
FDH_CRS(X,Y,input) = FDH(X,Y,input,CRS())
FDH_VRS(X,Y,input) = FDH(X,Y,input,VRS())
FDH_NIRS(X,Y,input) = FDH(X,Y,input,NIRS())
FDH_NDRS(X,Y,input) = FDH(X,Y,input,NDRS())

# Solve FDH program under VRS with (Xk,Yk) as evaluation point
function FDH(X::Array,Y::Array,input::Bool,RStype::VRS,Xk::Array,Yk::Array)
  if(input) #Input-oriented
    theta = Inf
    domy = find(all(Y .>= Yk,2))
    for k in domy
      Ixk = X[k,:] .> 0.0
      curmin = maximum(X[k,Ixk]./Xk[:,Ixk])
      if(curmin < theta)
        theta = curmin
      end
    end
  else #Output-oriented
    theta = -Inf
    domx = find(all(X .<= Xk,2))
    for k in domx
      Iyk = Y[k,:] .> 0.0
      curmax = minimum(Y[k,Iyk]./Yk[:,Iyk])
      if(curmax > theta)
        theta = curmax
      end
    end
  end
  return theta
end

# Solve FDH program under NDRS with (Xk,Yk) as evaluation point
function FDH(X::Array,Y::Array,input::Bool,RStype::NDRS,Xk::Array,Yk::Array)
  K = size(X,1)
  if(input) #Input-oriented
    theta = Inf
    for k=1:K
      Ixk = X[k,:] .> 0.0
      Jyk = Y[k,:] .> 0.0
      curmin = max(maximum(Yk[:,Jyk]./Y[k,Jyk]),1.0)*maximum(X[k,Ixk]./Xk[:,Ixk])
      if(curmin < theta)
        theta = curmin
      end
    end
  else #Output-oriented
    theta = -Inf
    domx = find(all(X .<= Xk,2))
    for k in domx
      Ixk = X[k,:] .> 0.0
      Jyk = Y[k,:] .> 0.0
      curmax = minimum(Y[k,Jyk]./Yk[:,Jyk])*minimum(Xk[:,Ixk]./X[k,Ixk])
      if(curmax > theta)
        theta = curmax
      end
    end
  end
  return theta
end

# Solve FDH program under CRS with (Xk,Yk) as evaluation point
function FDH(X::Array,Y::Array,input::Bool,RStype::CRS,Xk::Array,Yk::Array)
  K = size(X,1)
  if(input) #Input-oriented
    theta = Inf
    for k=1:K
      Ixk = X[k,:] .> 0.0
      Jyk = Y[k,:] .> 0.0
      curmin = maximum(Yk[:,Jyk]./Y[k,Jyk])*maximum(X[k,Ixk]./Xk[:,Ixk])
      if(curmin < theta)
        theta = curmin
      end
    end
  else #Output-oriented
    theta = -Inf
    for k=1:K
      Ixk = X[k,:] .> 0.0
      Jyk = Y[k,:] .> 0.0
      curmax = minimum(Y[k,Jyk]./Yk[:,Jyk])*minimum(Xk[:,Ixk]./X[k,Ixk])
      if(curmax > theta)
        theta = curmax
      end
    end
  end
  return theta
end

# Solve FDH program under NIRS with (Xk,Yk) as evaluation point
function FDH(X::Array,Y::Array,input::Bool,RStype::NIRS,Xk::Array,Yk::Array)
  K = size(X,1)
  if(input) #Input-oriented
    theta = Inf
    domy = find(all(Y .>= Yk,2))
    for k in domy
      Ixk = X[k,:] .> 0.0
      Jyk = Y[k,:] .> 0.0
      curmin = maximum(Yk[:,Jyk]./Y[k,Jyk])*maximum(X[k,Ixk]./Xk[:,Ixk])
      if(curmin < theta)
        theta = curmin
      end
    end
  else #Output-oriented
    theta = -Inf
    for k=1:K
      Ixk = X[k,:] .> 0.0
      Jyk = Y[k,:] .> 0.0
      curmax = minimum(Y[k,Jyk]./Yk[:,Jyk])*min(minimum(Xk[:,Ixk]./X[k,Ixk]),1.0)
      if(curmax > theta)
        theta = curmax
      end
    end
  end
  return theta
end

# Compute efficiency score for all DMUs
function FDH(X::Array,Y::Array,input::Bool,RStype::RS)
	K = size(X,1)
	theta = Array(Float64,K)
	for k=1:K
		theta[k] = FDH(X,Y,input,RStype,X[k,:],Y[k,:])
	end
	return theta
end
