#author: Tymoteusz Surynt
include("../Module/blocksys.jl")
using blocksys

result=importMatrix("/home/timmi/Desktop/ON/Tymoteusz_Surynt_Lista5/Tests/TestMatrix/ck10_10/A.txt")
result2=importVector("/home/timmi/Desktop/ON/Tymoteusz_Surynt_Lista5/Tests/TestMatrix/ck10_10/b.txt")
A=result[3]
#b=ones(Float64,10000)
b=result2[2]
w=gaussElimination2(result[1],result[2],A,b)
if w[4]==0
                exportVectorNoError("/home/timmi/Desktop/ON/Tymoteusz_Surynt_Lista5/Tests/TestMatrix/ck10_10/c.txt", w[3])
else
                printf("Error");
end
