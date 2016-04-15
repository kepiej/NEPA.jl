#cd("$(homedir())/Documents/GitHub")

# Push Documents/GitHub on the search path
push!(LOAD_PATH, "$(homedir())/Documents/GitHub")

using NEPA
using Base.Test

# Test DEAData.jl
X = [1 2 3;4 5 6; 7 8 9]
Y = [2;5;8]
Data = DEAData(X,Y,X,Y)

@test getNrDMU(Data) == 3
@test getIODim(Data) == (3,1)
@test getIndexes(Data) == [1,2,3]
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
