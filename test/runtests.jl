#cd("$(homedir())/Documents/GitHub")

# Push Documents/GitHub on the search path
push!(LOAD_PATH, "$(homedir())/Documents/GitHub")

using NEPA
using Base.Test
using Clp

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
@test_approx_eq_eps geteff(sbmres) [0.798 0.568 1.0 0.667 1.0] 1e-3
D = DEA_CRS(X,Y,true)
deares = D()
@test_approx_eq_eps geteff(deares) [0.9 0.833 1.0 1.0 1.0] 1e-3

Base.show(deares[1])
Base.show(sbmres)
Base.show(deares)
