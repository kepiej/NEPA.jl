type Convex <: AbstractDataEnvelopment
end

type FreeDisposal <: AbstractDataEnvelopment
end

# Directional distance function
immutable DDF{T<:Tuple{AbstractDataEnvelopment,RS}} <: AbstractDEA
	Data::DEAData

	function DDF(X::Array,Y::Array,gx::Array,gy::Array)
		new(DEAData(X,Y,gx,gy))
	end

	function DDF(Data::DEAData)
		new(Data)
	end
end

function getdata(D::DDF)
	return D.Data
end

# Solve convex DDF program with (Xk,Yk) as evaluation point in the direction of (gxk,gyk)
function Base.call{T<:RS}(D::DDF{Tuple{Convex,T}},Xk::Array,Yk::Array,gxk::Array,gyk::Array)
	K = getnrdmu(D.Data)
	N,M = getiodim(D.Data)

	# Set the appropriate RTS constraint depending on T
	RSA,RSb,RSsense = RSconstraint(K,T())

  f = zeros(K+1)
  f[1] = -1
  l = [-Inf;zeros(K)]
  u = [Inf for i=1:K+1]

  sense = Array(Char,N+M)
  sense[1:N+M] = '>'

	X,Y,gx,gy = D.Data[1:end]

	A = [-gxk' -X'; -gyk' Y'; RSA]
	b = [-Xk[1,:]'; Yk[1,:]'; RSb]
	b = b[:]

  # Solve linear program
  sol = linprog(f,A,[sense;RSsense],b,l,u)

  if sol.status == :Optimal
    beta = sol.sol[1]
  else
    beta = -Inf
    println("Error: solution status $(sol.status)")
  end

	return beta
end

# Solve free disposal DDF program with (Xk,Yk) as evaluation point in the direction of (gxk,gyk)
#FIXME Check that this is correct!
function Base.call(D::DDF{Tuple{FreeDisposal,VRS}},Xk::Array,Yk::Array,gxk::Array,gyk::Array)
  #TODO Implement for (gxk,gyk) < 0
  yind = gyk .> 0.0
  xind = gxk .> 0.0
  beta = 0.0
  for i in eachindex(D.Data)
    Xi,Yi,gxi,gyi = D.Data[i]
    curmin = minimum([((Yi[yind]-Yk[yind])./gyk[yind]); ((Xk[xind]-Xi[xind])./gxk[xind])])
    if(curmin > beta)
      beta = curmin
    end
  end
	
	return beta
end
