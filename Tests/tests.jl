#author: Tymoteusz Surynt
include("../Module/blocksys.jl")
using blocksys
using matrixgen
#Function responsible for calculating b vector for some random matrixes
function cVector(ml,vl)
    r=importMatrix(ml)
    n::Int64=r[1]
    l::Int64=r[2]
    A=r[3]
    v::Int64=floor(r[1]/r[2])
    it::Int64=2
    sum::Float64=0.0
    b=Array{Float64}(n+1)
    b[1]=n
    for i in 1:v
        for j in 1:l
            sum=0.0
            for k in 1:2
                sum+=A[i,1][j,k]
            end
            for k in 1:l
                sum+=A[i,2][j,k]
                sum+=A[i,3][j,k]
            end
            b[it]=sum
            it+=1
        end
    end
    exportVectorNoError(vl,b)
end
cVector("./TestMatrix/ck10/A.txt","./TestMatrix/ck10/b.txt")
cVector("./TestMatrix/ck10_2/A.txt","./TestMatrix/ck10_2/b.txt")
cVector("./Tests/TestMatrix/ck10_5/A.txt","./TestMatrix/ck10_5/b.txt")
cVector("./TestMatrix/ck10_10/A.txt","./TestMatrix/ck10_10/b.txt")
cVector("./TestMatrix/ck10_13/A.txt","./TestMatrix/ck10_13/b.txt")
