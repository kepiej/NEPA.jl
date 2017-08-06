#cd("$(homedir())/Documents/GitHub")

# Push Documents/GitHub on the search path
push!(LOAD_PATH, "$(homedir())/Documents/GitHub")

using NEPA
using Base.Test

# Test DEAData.jl
X = [1 2 3;4 5 6; 7 8 9]
Y = [2;5;8]
Data = DEAData(X,Y,X,Y)

@test getnrdmu(Data) == 3
@test getiodim(Data) == (3,1)
@test getindexes(Data) == [1,2,3]
@test size(Data) == (3,)

for i in eachindex(Data)
  Xk,Yk,gxk,gyk = Data[i]
  @test (Xk,Yk,gxk,gyk) == (X[i,:],Y[i,:],X[i,:],Y[i,:])
  @test Xk == gxk
  @test Yk == gyk
end

ind = [1,3]
@test Data[ind] == (X[ind,:],Y[ind,:],X[ind,:],Y[ind,:])
@test Data[1:end] == (X[1:end,:],Y[1:end,:],X[1:end,:],Y[1:end,:])

Data = DEAData(X,Y,X,Y,indexes = ind)
@test Data[1:end] == (X[ind,:],Y[ind,:],X[ind,:],Y[ind,:])

# Test AbstractDEA
@test eltype(AbstractDEA{Convex,VRS}) == (Convex,VRS)
@test eltype(AbstractDEA{FreeDisposal,VRS}) != (Convex,VRS)
@test eltype(AbstractDEA{FreeDisposal,CRS}) == (FreeDisposal,CRS)
@test eltype(AbstractDEA{FreeDisposal,VRS}) != (FreeDisposal,CRS)

#Test SBM model using Tone(2001) example
X = [4.0 3.0; 6.0 3.0; 8.0 1.0; 8.0 1.0; 2.0 4.0]
Y = [2.0 3.0; 2.0 3.0; 6.0 2.0; 6.0 1.0; 1.0 4.0]

D = SBM{CRS}(X,Y)
sbmres = D()
@test geteff(sbmres) ≈ [0.798 0.568 1.0 0.667 1.0]' atol=1e-3
D = DEA_CRS(X,Y,true)
deares = D()
@test geteff(deares) ≈ [0.9 0.833 1.0 1.0 1.0]' atol=1e-3

#Base.show(deares[1])
#Base.show(sbmres)
#Base.show(deares)

#Test DDF model using Cherchye et al. (2001) numerical example
X = [3.0;6.0;4.0;6.0;5.0;8.0;12.0;14.0;18.0]
Y = [4.0;5.0;6.0;7.0;8.0;9.0;11.0;13.0;14.0]

DDFMcFadden = DDF{FreeDisposal,VRS}(X,Y,-X,Y)
@test geteff(DDFMcFadden(X[2,:],Y[2,:],-X[2,:],Y[2,:])) ≈ 1.6 atol=1e-3
@test find(getpeers(DDFMcFadden(X[2,:],Y[2,:],-X[2,:],Y[2,:])) .> 0.0) == [8] # peer H

Ysens = zeros(9)
Ysens[9,:] = 1.0
DDFMcFadden = DDF{FreeDisposal,VRS}(X,Y+Ysens,-X,Y+Ysens)
@test_approx_eq_eps geteff(DDFMcFadden(X[2,:],Y[2,:],-X[2,:],Y[2,:])) 2.0 1e-3
@test find(getpeers(DDFMcFadden(X[2,:],Y[2,:],-X[2,:],Y[2,:])) .> 0.0) == [9] # peer I

#FIXME For some reason incorrect results: it should be 9.5 but we find 1.8...
Xsens = zeros(9)
Xsens[9,:] = 1.2
DDFMcFadden = DDF{FreeDisposal,VRS}(X-Xsens,Y,-X+Xsens,Y)
@test geteff(DDFMcFadden(X[2,:],Y[2,:],-X[2,:],Y[2,:])) ≈ 9.5 atol=1e-3
@test find(getpeers(DDFMcFadden(X[2,:],Y[2,:],-X[2,:],Y[2,:])) .> 0.0) == [9] # peer I

Graph = Hyperbolic{FreeDisposal,VRS}(X,Y)
println(geteff(Graph()))
