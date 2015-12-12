function WACM(W,X,Y)
	K = size(W,1)
	eff = zeros(K,1)
	for k=1:K
		eff[k,:] = [minimum(X[all(Y .>= Y[k,:],2),:]*W[k,:]')]/(W[k,:]*X[k,:]')
	end
	return eff
end
