abstract AbstractDataEnvelopment

type Convex <: AbstractDataEnvelopment
end

type FreeDisposal <: AbstractDataEnvelopment
end

# Abstract DEA model
abstract AbstractDEA <: AbstractArray{Any,1}

function getData(A::AbstractDEA)
	# Concrete types should overload this method
end

Base.size(A::AbstractDEA) = size(getData(A))
Base.linearindexing(::Type{AbstractDEA}) = Base.LinearFast()

function Base.getindex(A::AbstractDEA, i::Int)
	1 <= i <= size(A,1) || throw(BoundsError(A, i))
	return getData(A)[i]
end

function Base.getindex(A::AbstractDEA, I)
	return getData(A)[I]
end

function Base.call(A::AbstractDEA)
	Data = getData(A)
	beta = Array(Float64,getNrDMU(Data))

	@sync @parallel for k in eachindex(Data)
		beta[k] = Base.call(A,Data[k]...)
	end

	return beta
end

# Directional distance function
immutable DDF{T<:Tuple{AbstractDataEnvelopment,RS}} <: AbstractDEA
	Data::DEAData

	function DDF(X,Y,gx,gy)
		new(DEAData(X,Y,gx,gy))
	end

	function DDF(Data::DEAData)
		new(Data)
	end
end

function getData(D::DDF)
	return D.Data
end

# Solve convex DDF program with (Xk,Yk) as evaluation point in the direction of (gxk,gyk)
function Base.call{T<:RS}(D::DDF{Tuple{Convex,T}},Xk::Array,Yk::Array,gxk::Array,gyk::Array)
	K = getNrDMU(D.Data)
	N,M = getIODim(D.Data)

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

	#Check for special cases where one of the direction vectors is zero
	if gxk == zeros(size(gxk))
		return 1.0 + beta
	elseif gyk == zeros(size(gyk))
		return 1.0 - beta
	else
		return beta
	end
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

	#Check for special cases where one of the direction vectors is zero
	if gxk == zeros(size(gxk))
		return 1.0 + beta
	elseif gyk == zeros(size(gyk))
		return 1.0 - beta
	else
		return beta
	end
end
