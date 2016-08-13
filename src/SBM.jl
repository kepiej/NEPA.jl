immutable SBM{T<:RS} <: AbstractDEA{Convex,T}
  Data::DEAData

  function SBM(X,Y)
    new(DEAData(X,Y))
  end
end

function getdata(DMU::SBM)
  return DMU.Data
end

function Base.call{T<:RS}(DMU::SBM{T},Xk::Array,Yk::Array; sx::Array = [],sy::Array = [],lambda::Array = [])
  K = getnrdmu(DMU.Data)
	N,M = getiodim(DMU.Data)

	# Set the appropriate RTS constraint depending on T
	RSA,RSb,RSsense = RSconstraint(K,T())

  f = zeros(K+1+N+M)
  f[1] = 1
  f[K+2:K+N+1] = -1./(N.*Xk)
  l = zeros(K+1+N+M)
  l[1] = eps()
  u = [Inf for i=1:K+1+N+M]

  sense = Array(Char,N+M+1)
  sense[1:N+M+1] = '='

	X,Y = DMU.Data[1:end]

	A = [Xk[1,:]' -X' -eye(N) zeros(N,M); Yk[1,:]' -Y' zeros(M,N) eye(M); 1 zeros(1,K+N) (1./(M.*Yk)).*ones(1,M); RSA zeros(1,N+M)]
	b = [zeros(N); zeros(M); 1; RSb]
	b = b[:]
  # Solve linear program
  sol = linprog(f,A,[sense;RSsense],b,l,u)

  if sol.status == :Optimal
    rho = sol.objval
    t = sol.sol[1]
    lambda = sol.sol[2:K+1]./t
    sx = sol.sol[K+2:K+1+N]./t
    sy = sol.sol[K+1+N+1:K+1+N+M]./t
    res = DEAResult(rho,sx,sy,lambda)
  else
    res = DEAResult(NaN)
    println("Error: solution status $(sol.status)")
  end

	return res
end
