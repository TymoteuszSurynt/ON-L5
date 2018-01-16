#author: Tymoteusz Surynt
include("../Module/blocksys.jl")
using blocksys

result=importMatrix("/home/timmi/Desktop/ON/Dane16_1_1/A.txt")
result2=importVector("/home/timmi/Desktop/ON/Dane16_1_1/b.txt")
if result[1]==result2[1]
                A=result[3]
                b=result2[2]
                w=gaussElimination(result[1],result[2],A,b)
                if w[4]==0
                                exportVectorError("/home/timmi/Desktop/ON/Tymoteusz_Surynt_Lista5/Tests/test16_gauss1.txt", w[3], getError(w[3],1.0,result[1]))
                else
                                printf("Error");
                end
                w=gaussElimination2(result[1],result[2],A,b)
                if w[4]==0
                                exportVectorError("/home/timmi/Desktop/ON/Tymoteusz_Surynt_Lista5/Tests/test16_gauss2.txt", w[3], getError(w[3],1.0,result[1]))
                else
                                printf("Error");
                end
else
                printf("Error")
end
