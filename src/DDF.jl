type Convex <: AbstractDataEnvelopment
end

type FreeDisposal <: AbstractDataEnvelopment
end

# Directional distance function
immutable DDF{S<:AbstractDataEnvelopment,T<:RS} <: AbstractDEA{S,T}
	Data::DEAData

	function DDF(X::Array,Y::Array,gx::Array,gy::Array)
		new(DEAData(X,Y,gx,gy))
	end
end

function getdata(D::DDF)
	return D.Data
end

# Solve convex DDF program with (Xk,Yk) as evaluation point in the direction of (gxk,gyk)
function Base.call{T<:RS}(D::DDF{Convex,T},Xk::Array,Yk::Array,gxk::Array,gyk::Array)
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
		res = DEAResult(sol.sol[1],[],[],sol.sol[2:end])
  else
		res = DEAResult(-Inf)
    println("Error: solution status $(sol.status)")
  end

	return res
end

# Solve free disposal DDF program with (Xk,Yk) as evaluation point in the direction of (gxk,gyk)
#FIXME Check that this is correct!
function Base.call(D::DDF{FreeDisposal,VRS},Xk::Array,Yk::Array,gxk::Array,gyk::Array)
	gxk0 = map(isapprox,gxk,0.0.*ones(size(gxk)))
	gyk0 = map(isapprox,gyk,0.0.*ones(size(gyk)))
	if(any(gxk0) || any(gyk0))
		Q = []
		for i in eachindex(D.Data)
			Xi,Yi,gxi,gyi = D.Data[i]
			if(all(Yi[gyk0] .>= Yk[gyk0]) && all(Xi[gxk0] .<= Xk[gxk0]))
				Q = [Q;i]
			end
		end
	else
		Q = D.Data
	end

  ypos = gyk .> 0.0
  xpos = gxk .> 0.0
	yneg = gyk .< 0.0
	xneg = gxk .< 0.0
  alpha = -Inf
	for i in eachindex(Q)
    Xi,Yi,gxi,gyi = D.Data[i]

		if(any(xneg) || any(yneg))
			gamma = maximum([((Yi[yneg]-Yk[yneg])./gyk[yneg]); ((Xk[xneg]-Xi[xneg])./gxk[xneg])])
		else
			gamma = 0.0
		end

		beta = minimum([((Yi[ypos]-Yk[ypos])./gyk[ypos]); ((Xk[xpos]-Xi[xpos])./gxk[xpos])])

		if(beta < gamma)
			beta = 0.0
		end

    if(beta > alpha)
      alpha = beta
    end
  end

	return DEAResult(alpha)
end
