#author: Tymoteusz Surynt
include("../Module/blocksys.jl")
using blocksys
b=importVector("/home/timmi/Desktop/ON/Dane50000_1_1/b.txt")
A=importMatrix("/home/timmi/Desktop/ON/Dane50000_1_1/A.txt")
w=gaussElimination2(A[1],A[2],A[3],b[2])

#printfMatrix(A[3],b[2],A[2],A[1])
# result2=importVector("/home/timmi/Desktop/ON/Tymoteusz_Surynt_Lista5/Tests/TestMatrix/ck10_10/b.txt")
# A=result[3]
# #b=ones(Float64,10000)
# b=result2[2]
# w=gaussElimination2(result[1],result[2],A,b)
# if w[4]==0
#                 exportVectorNoError("/home/timmi/Desktop/ON/Tymoteusz_Surynt_Lista5/Tests/TestMatrix/ck10_10/c.txt", w[3])
# else
#                 printf("Error");
# end
