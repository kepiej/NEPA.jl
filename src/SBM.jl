immutable SBM{T<:RS} <: AbstractDEA{Convex,T}
  Data::DEAData

  function SBM{T}(X,Y) where T<:RS
    new(DEAData(X,Y))
  end
end

function getdata(DMU::SBM)
  return DMU.Data
end

function (DMU::SBM{T})(Xk::Array,Yk::Array; sx::Array = [],sy::Array = [],lambda::Array = []) where T<:RS
  K = getnrdmu(DMU.Data)
  N,M = getiodim(DMU.Data)

  # Set the appropriate RTS constraint depending on T
  RSA,RSb,RSsense = RSconstraint(K,T())

  f = zeros(K+1+N+M)
  f[1] = 1
  f[K+2:K+N+1] = -1./(N.*reshape(Xk,:,1))
  l = zeros(K+1+N+M)
  l[1] = eps()
  u = [Inf for i=1:K+1+N+M]

  sense = Array{Char}(N+M+1)
  sense[1:N+M+1] = '='

  X,Y = DMU[1:end]

  A = [reshape(Xk',:,1) -X' -eye(N) zeros(N,M); reshape(Yk',:,1) -Y' zeros(M,N) eye(M); 1 zeros(1,K+N) reshape(1./(M.*Yk),1,:); RSA zeros(1,N+M)]
  b = [zeros(N); zeros(M); 1; RSb]
  b = b[:]

  # Solve linear program
  sol = linprog(f,A,[sense;RSsense],b,l,u,ClpSolver())

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

function (DMU::SBM{T})() where T<:RS
  	Data = getdata(DMU)
  	res = Array{DEAResult}(getnrdmu(Data))

  	@sync @parallel for k in eachindex(Data)
  		res[k] = DMU(Data[k]...)
  	end

  	return res
end
