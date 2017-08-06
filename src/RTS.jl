# We define several Returns to Scale (RS) types. The RS constraint used in DDF depends on these RTS type
abstract type RS end

type CRS <: RS
end

type VRS <: RS
end

type NIRS <: RS
end

type NDRS <: RS
end

# The default RTS constraint returns an all zero vector.
# Constant Returns to Scale (CRS) only imposes \lambda_i >= 0
function RSconstraint{T<:RS}(K,RStype::T)
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
