# We define several Returns to Scale (RS) types. The RS constraint in DEA depends on these RTS type
abstract RS

type CRS <: RS
end

type VRS <: RS
end

type NIRS <: RS
end

type NDRS <: RS
end