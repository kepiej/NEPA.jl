cd("$(homedir())/Documents/Julia")

# Push Documents/Julia on the search path
push!(LOAD_PATH, "$(homedir())/Documents/Julia")

using ExcelReaders
#using MAT
using NEPA

#mat = matread("C:/Users/u0093191/Documents/MATLAB/imputeLCM.mat")
#println(mat)

X = readxl("C:/Users/u0093191/Documents/MATLAB/WACM/Atkinson-RSanalyis-11May2011.xls","Data!F4:H195")
W = readxl("C:/Users/u0093191/Documents/MATLAB/WACM/Atkinson-RSanalyis-11May2011.xls","Data!J4:L195")
Y = readxl("C:/Users/u0093191/Documents/MATLAB/WACM/Atkinson-RSanalyis-11May2011.xls","Data!E4:E195")

# Convert DataFrame or DataArray to Array so that the  function WACM doesn't crash
X = convert(Array{Float64},X)
W = convert(Array{Float64},W)
Y = convert(Array{Float64},Y)

#eff = WACM(W,X,Y)
#Time execution using @time
#@time eff = WACM(w,x,y)
#@code_warntype WACM(w,x,y)
#eff = WACM(w,x,y)
#print(eff)

#Input- or output-oriented?
input = true

#theta = FDH_VRS(X,Y,input)
#println(theta)
#theta = FDH_CRS(X,Y,input)
#println(theta)
#theta = FDH_NIRS(X,Y,input)
#println(theta)
#theta = FDH_NDRS(X,Y,input)
#println(theta)

#D = DEA(X,Y,input,VRS())
D = DEA_VRS(X,Y,input)
#D = DEA_CRS(X,Y,input)
#D = DEA_NIRS(X,Y,input)
#D = DEA_NDRS(X,Y,input)
@time theta = D() #0.187408 seconds
@time theta = D() #0.000068 seconds
println(theta)
