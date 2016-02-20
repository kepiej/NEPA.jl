# The default RTS constraint returns an all zero vector.
# Constant Returns to Scale (CRS) only imposes \lambda_i >= 0
function RSconstraint{T <: RS}(K,RStype::T)
	return zeros(Float64,1,K+1),0,'='
end

# Variable Returns to Scale (VRS) imposes sum_i{\lambda_i} = 1
function RSconstraint(K,RStype::VRS)
	return [0 ones(1,K)],1,'='
end

# Non-Increasing Returns to Scale (NIRS) imposes sum_i{\lambda_i} <= 1
function RSconstraint(K,RStype::NIRS)
	return [0 ones(1,K)],1,'<'
end

# Non-Decreasing Returns to Scale (NDRS) imposes sum_i{\lambda_i} >= 1
function RSconstraint(K,RStype::NDRS)
	return [0 ones(1,K)],1,'>'
end

# Abstract DEA model
abstract AbstractDEA

# Convex directional distance function
immutable DDF <: AbstractDEA
	Data::DEAData
	RSA::Array
	RSb::Int64
	RSsense::Char

	function DDF(X,Y,gx,gy,RStype::RS)
		# Set the appropriate RTS constraint depending on RStype
		RSA,RSb,RSsense = RSconstraint(size(X,1),RStype)
		new(DEAData(X,Y,gx,gy),RSA,RSb,RSsense)
	end
end

# Solve DDF program with (Xk,Yk) as evaluation point in the direction of (gxk,gyk)
function Base.call(D::DDF,Xk::Array,Yk::Array,gxk::Array,gyk::Array)
	K = getNrDMU(D.Data)
	N,M = getIODim(D.Data)

  f = zeros(K+1)
  f[1] = -1
  l = [-Inf;zeros(K)]
  u = [Inf for i=1:K+1]

  sense = Array(Char,N+M)
  sense[1:N+M] = '>'

	X,Y,gx,gy = D.Data[1:end]

	A = [-gxk' -X'; -gyk' Y'; D.RSA]
	b = [-Xk[1,:]'; Yk[1,:]'; D.RSb]
	b = b[:]

  # Solve linear program
  sol = linprog(f,A,[sense;D.RSsense],b,l,u)

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

function Base.call(D::DDF)
	beta = Array(Float64,getNrDMU(D.Data))

	@sync @parallel for k in eachindex(D.Data)
		beta[k] = Base.call(D,D.Data[k]...)
	end

	return beta
end
