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
		res = DEAResult(sol.sol[1],sol.sol[1].*gxk,sol.sol[1].*gyk,sol.sol[2:end])
  else
		res = DEAResult(-Inf)
    println("Error: solution status $(sol.status)")
  end

	return res
end

# Solve free disposal DDF program under VRS with (Xk,Yk) as evaluation point in the direction of (gxk,gyk)
#FIXME Check that this is correct!
# function Base.call(D::DDF{FreeDisposal,VRS},Xk::Array,Yk::Array,gxk::Array,gyk::Array)
# 	gxk0 = map(isapprox,gxk,0.0.*ones(size(gxk)))
# 	gyk0 = map(isapprox,gyk,0.0.*ones(size(gyk)))
# 	if(any(gxk0) || any(gyk0))
# 		Q = []
# 		for i in eachindex(D.Data)
# 			Xi,Yi,gxi,gyi = D.Data[i]
# 			if(all(Yi[gyk0] .>= Yk[gyk0]) && all(Xi[gxk0] .<= Xk[gxk0]))
# 				Q = [Q;i]
# 			end
# 		end
# 	else
# 		Q = D.Data
# 	end
#
#   ypos = gyk .> 0.0
#   xpos = gxk .> 0.0
# 	yneg = gyk .< 0.0
# 	xneg = gxk .< 0.0
#   alpha = -Inf
# 	alphapeer = 0
# 	gamma = zeros(getnrdmu(D.Data))
# 	beta = zeros(getnrdmu(D.Data))
# 	alpha = zeros(getnrdmu(D.Data))
# 	for i in eachindex(Q)
#     Xi,Yi,gxi,gyi = D.Data[i]
#
# 		if(any(xneg) || any(yneg))
# 			gamma[i] = maximum([((Yi[yneg]-Yk[yneg])./gyk[yneg]); ((Xk[xneg]-Xi[xneg])./gxk[xneg])])
# 		else
# 			gamma[i] = 0.0
# 		end
#
# 		beta[i] = minimum([((Yi[ypos]-Yk[ypos])./gyk[ypos]); ((Xk[xpos]-Xi[xpos])./gxk[xpos])])
#
# 		if(beta[i] < gamma[i])
# 			alpha[i] = 0.0
# 		else
# 			alpha[i] = beta[i]
# 		end
#
#   end
# 	alphapeer = indmax(alpha)
# 	alpha = maximum(alpha)
#
# 	peers = zeros(getnrdmu(D.Data))
# 	peers[alphapeer] = 1.0
# 	return DEAResult(alpha,alpha.*gxk,alpha.*gyk,peers)
# end

function Base.call{T<:RS}(D::DDF{FreeDisposal,T},Xk::Array,Yk::Array,gxk::Array,gyk::Array)
	X,Y,gx,gy = getdata(D)[1:end]

	#Find zero entries in gxk,gyk
	gxk0 = map(isapprox,gxk,0.0.*ones(size(gxk)))
	gyk0 = map(isapprox,gyk,0.0.*ones(size(gyk)))

	#Simar and Vanhems(2012) trick using transformed dataset and computing hyperbolic distance function on the transformed data
	#Note: This assumes gxk,gyk >= 0!
	DGR = Hyperbolic{FreeDisposal,T}(exp(hcat(X[:,!gxk0]./gxk[:,!gxk0],X[:,gxk0])),exp(hcat(Y[:,!gyk0]./gyk[:,!gyk0],Y[:,gyk0])))
	DGRres = DGR(exp(hcat(Xk[:,!gxk0]./gxk[:,!gxk0],Xk[:,gxk0])),exp(hcat(Yk[:,!gyk0]./gyk[:,!gyk0],Yk[:,gyk0])))

	# xpos = gxk .> 0.0
	# ypos = gyk .> 0.0
	# xneg = gxk .< 0.0
	# yneg = gyk .< 0.0

	# Xtr = hcat(exp(X[:,xpos]./gxk[:,xpos]),exp(Y[:,yneg]./gyk[:,yneg]))
	# if(isempty(Xtr))
	# 	Xtr = zeros(size(X))
	# end
	# Xktr = hcat(exp(Xk[:,xpos]./gxk[:,xpos]),exp(Yk[:,yneg]./gyk[:,yneg]))
	# if(isempty(Xktr))
	# 	Xktr = zeros(size(Xk))
	# end
	# Ytr = hcat(exp(Y[:,ypos]./gyk[:,ypos]),exp(X[:,xneg]./gxk[:,xneg]))
	# if(isempty(Ytr))
	# 	Ytr = zeros(size(Y))
	# end
	# Yktr = hcat(exp(Yk[:,ypos]./gyk[:,ypos]),exp(Xk[:,xneg]./gxk[:,xneg]))
	# if(isempty(Yktr))
	# 	Yktr = zeros(size(Yk))
	# end
	# println(Xtr)
	# println(Xktr)
	# DGR = Hyperbolic{FreeDisposal,T}(Xtr,Ytr)
	# DGRres = DGR(Xktr,Yktr)

	beta = log(1./geteff(DGRres))
	return DEAResult(beta,beta.*gxk,beta.*gyk,getpeers(DGRres))
end
