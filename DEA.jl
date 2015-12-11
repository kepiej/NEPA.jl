# Define shorthand functions depending on the returns to scale
DEA_CRS = (X,Y,input) -> DEA(X,Y,input,CRS())
DEA_VRS = (X,Y,input) -> DEA(X,Y,input,VRS())
DEA_NIRS = (X,Y,input) -> DEA(X,Y,input,NIRS())
DEA_NDRS = (X,Y,input) -> DEA(X,Y,input,NDRS())

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

immutable DEA
  X::Array
  Y::Array
  input::Bool
	RSA::Array
	RSb::Int64
	RSsense::Char

	function DEA(X,Y,input,RStype::RS)
		if(size(X,1) != size(Y,1))
			error("Number of DMUs in X and Y are not the same: X has $(size(X,1)) DMUs and Y has $(size(Y,1)) DMUs!")
		end
		# Set the appropriate RTS constraint depending on RStype
		RSA,RSb,RSsense = RSconstraint(size(X,1),RStype)
		new(X,Y,input,RSA,RSb,RSsense)
	end
end

# Solve DEA program with (Xk,Yk) as evaluation point
function Base.call(D::DEA,Xk::Array,Yk::Array)
  K,N = size(D.X)
  M = size(D.Y,2)

  f = zeros(K+1)

  if(D.input) # Input oriented
    f[1] = 1
    l = zeros(K+1)
    u = [1;[Inf for i=1:K]]
  else # Output oriented
    f[1] = -1
    l = [1;zeros(K)]
    u = [Inf for i=1:K+1]
  end

  sense = Array(Char,N+M)
  sense[1:N+M] = '>'

  if(D.input) # Input oriented
    A = [Xk[1,:]' -D.X'; zeros(M,1) D.Y'; D.RSA]
    b = [zeros(N+M); D.RSb]
    b[N+1:N+M] = Yk[1,:]'
  else # Output oriented
    A = [zeros(N,1) -D.X'; -Yk[1,:]' D.Y'; D.RSA]
    b = [zeros(N+M); D.RSb]
    b[1:N] = -Xk[1,:]'
  end

  # Solve linear program
  sol = linprog(f,A,[sense;D.RSsense],b,l,u)

  if sol.status == :Optimal
    theta = sol.sol[1]
  else
    theta = NaN
    println("Error: solution status $(sol.status)")
  end

  return theta
end

# Compute DEA efficiency score for all DMUs
function Base.call(D::DEA)
	K = size(D.X,1)
	theta = Array(Float64,K)
	@parallel for k=1:K
		theta[k] = Base.call(D,D.X[k,:],D.Y[k,:])
	end
	return theta
end
