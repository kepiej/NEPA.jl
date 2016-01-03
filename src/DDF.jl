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

immutable DDF
	X::Array
  Y::Array
  gx::Array
	gy::Array
	RSA::Array
	RSb::Int64
	RSsense::Char

	function DDF(X,Y,gx,gy,RStype::RS)
		if(size(X,1) != size(Y,1))
			error("Number of DMUs in X and Y are not the same: X has $(size(X,1)) DMUs and Y has $(size(Y,1)) DMUs!")
		end
		if(size(gx) != size(X) || size(gy) != size(Y))
			error("Dimensions of gx and X or gy and Y are not the same: size(gx) = $(size(gx)), size(X) = $(size(X)), size(gy) = $(size(gy)) and size(Y) = $(size(Y)) DMUs!")
		end
		if(gx == zeros(size(gx)) && gy == zeros(size(gy)))
			error("Direction vectors cannot both be zero!")
		end
		# Set the appropriate RTS constraint depending on RStype
		RSA,RSb,RSsense = RSconstraint(size(X,1),RStype)
		new(X,Y,gx,gy,RSA,RSb,RSsense)
	end
end

# Solve DDF program with (Xk,Yk) as evaluation point
function Base.call(D::DDF,Xk::Array,Yk::Array,gxk::Array,gyk::Array)
  K,N = size(D.X)
  M = size(D.Y,2)

  f = zeros(K+1)
  f[1] = -1
  l = [0;zeros(K)]
  u = [Inf for i=1:K+1]

  sense = Array(Char,N+M)
  sense[1:N+M] = '>'

	A = [-gxk' -D.X'; -gyk' D.Y'; D.RSA]
	b = [-Xk[1,:]'; Yk[1,:]'; D.RSb]
	b = b[:]

  # Solve linear program
  sol = linprog(f,A,[sense;D.RSsense],b,l,u)

  if sol.status == :Optimal
    beta = sol.sol[1]
  else
    beta = Inf
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
	K = size(D.X,1)
	beta = Array(Float64,K)

	@sync @parallel for k=1:K
		beta[k] = Base.call(D,D.X[k,:],D.Y[k,:],D.gx[k,:],D.gy[k,:])
	end

	return beta
end
