#author: Tymoteusz Surynt
include("../Module/blocksys.jl")
using blocksys

result=importMatrix("./TestMatrix/ck10_10/A.txt")
result2=importVector("./TestMatrix/ck10_10/b.txt")
if result[1]==result2[1]
                A=result[3]
                b=result2[2]
                w=gaussElimination2(result[1],result[2],A,b)
                if w[4]==0
                                exportVectorError("./TestMatrix/ck10_10/x_gauss2.txt", w[3], getError(w[3],1.0,result[1]))
                else
                                println("Error");
                end
else
                println("Error")
end
