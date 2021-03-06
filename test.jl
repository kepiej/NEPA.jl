cd("$(homedir())/Documents/GitHub")

# Push Documents/GitHub on the search path
push!(LOAD_PATH, "$(homedir())/Documents/GitHub")

using FileIO, ExcelFiles, DataFrames
#using ExcelReaders
#using MAT
using NEPA

#mat = matread("C:/Users/u0093191/Documents/MATLAB/imputeLCM.mat")
#println(mat)

X = DataFrame(load("C:/Users/u0093191/Documents/MATLAB/WACM/Atkinson-RSanalyis-11May2011.xls","Data!F4:H195"))
W = DataFrame(load("C:/Users/u0093191/Documents/MATLAB/WACM/Atkinson-RSanalyis-11May2011.xls","Data!J4:L195"))
Y = DataFrame(load("C:/Users/u0093191/Documents/MATLAB/WACM/Atkinson-RSanalyis-11May2011.xls","Data!E4:E195"))

#X = readxl("C:/Users/u0093191/Documents/MATLAB/WACM/Atkinson-RSanalyis-11May2011.xls","Data!F4:H195")
#W = readxl("C:/Users/u0093191/Documents/MATLAB/WACM/Atkinson-RSanalyis-11May2011.xls","Data!J4:L195")
#Y = readxl("C:/Users/u0093191/Documents/MATLAB/WACM/Atkinson-RSanalyis-11May2011.xls","Data!E4:E195")

# Convert DataFrame or DataArray to Array so that the  function WACM doesn't crash
X = convert(Array{Float64},X)
W = convert(Array{Float64},W)
Y = convert(Array{Float64},Y)

# Input-oriented efficiencies computed using Maple. We use these to check our computed results
#MapleEff = readxl("C:/Users/u0093191/Documents/MATLAB/WACM/Atkinson-RSanalyis-11May2011.xls","MapleEfficiencies!B2:O193")
MapleEff = DataFrame(load("C:/Users/u0093191/Documents/MATLAB/WACM/Atkinson-RSanalyis-11May2011.xls","MapleEfficiencies!B2:O193"))
MapleEff = convert(Array{Float64},MapleEff)

#eff = WACM(W,X,Y)
#Time execution using @time
#@time eff = WACM(w,x,y)
#@code_warntype WACM(w,x,y)
#eff = WACM(w,x,y)
#print(eff)

#Input- or output-oriented?
input = true

D = FDH_VRS(X,Y,input)
println(maximum(abs(geteff(D()) - MapleEff[:,8])))
D = FDH_CRS(X,Y,input)
println(maximum(abs(geteff(D()) - MapleEff[:,9])))
D = FDH_NIRS(X,Y,input)
println(maximum(abs(geteff(D()) - MapleEff[:,10])))
D = FDH_NDRS(X,Y,input)
println(maximum(abs(geteff(D()) - MapleEff[:,11])))

#M = convert(Int,round(3*size(D,1)/4))
#thetaM = orderm(D,M,100)
#println(thetaM - D())

# This should be equal to the FDH_VRS efficiency!
D = DDF{FreeDisposal,VRS}(X,Y,X,zeros(size(Y)))
println(maximum(abs((1 - geteff(D())) - MapleEff[:,8])))
#println(hcat((1 - geteff(D())), MapleEff[:,8]))

# Hyperbolic equals square root of input oriented FDH under CRS
GDR = Hyperbolic{FreeDisposal,CRS}(X,Y)
D = FDH_CRS(X,Y,input)
println(maximum(abs(geteff(GDR()) - sqrt(geteff(D())))))

D = DEA_VRS(X,Y,input)
println(maximum(abs(geteff(D()) - MapleEff[:,1])))
D = DEA_CRS(X,Y,input)
println(maximum(abs(geteff(D()) - MapleEff[:,2])))
D = DEA_NIRS(X,Y,input)
println(maximum(abs(geteff(D()) - MapleEff[:,3])))
D = DEA_NDRS(X,Y,input)
println(maximum(abs(geteff(D()) - MapleEff[:,4])))
#@time theta = D() #4.389059 seconds
#@time theta = D() #0.123012 seconds

D = DDF{Convex,VRS}(X,Y,X,zeros(size(Y)))
println(maximum(abs((1 - geteff(D())) - MapleEff[:,1])))

# Test indexing of DDF object
println(D[1])
println(X[1,:],Y[1,:],X[1,:],zeros(1,size(Y,2)))

# Do order-m efficiency
#M = convert(Int,round(0.9*size(D,1)))
#thetaM = orderm(D,M,100)
#println(thetaM - geteff(D()))

#data = readdlm("./NEPA/data/GriffellTatjéLovell.txt")
#X0 = data[:,1]
#Y0 = data[:,2]
#X1 = data[:,3]
#Y1 = data[:,4]

#LTFP = Luenberger{Convex,VRS}(X,Y,X,Y,X,Y,X,Y)
#println(LTFP())
#println(TEI(LTFP))
#println(TC(LTFP))

# HMBTFP = HMB{Convex,VRS}(X,Y,X,Y)
# println(HMBTFP())
# println(TEI(HMBTFP))
# println(TC(HMBTFP))
# println(SEC(HMBTFP))
