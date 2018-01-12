#author: Tymoteusz Surynt

using ObliczeniaNaukowe5
result=importMatrix("/home/timmi/Desktop/ON/Dane16_1_1/A.txt")
A=result[3]

abc=importVector("/home/timmi/Desktop/ON/Dane16_1_1/b.txt")
b=abc[2]
printfMatrix(A,b,result[2],result[1])
w=gaussElimination(result[1],result[2],A,b)
printfMatrix(w[1],w[2],result[2],result[1])
eps=1/10^8
function gaussElimination(n,l,AB,v)
    A=deepcopy(AB)
    b=deepcopy(v)
    it::Int64=2
    #First iteration
    for i in 1:(l-1)
        for j in (i+1):l
            #if abs(A[1,2][i,i])<eps
            #    return (A,1)
            #end
            m=A[1,2][j,i]/A[1,2][i,i]
            for k in i:l
                A[1,2][j,k]=A[1,2][j,k]-m*A[1,2][i,k]
            end
            for k in 1:l
                A[1,3][j,k]=A[1,3][j,k]-m*A[1,3][i,k]
            end
            b[j]=b[j]-m*b[j-1]
        end
    end
    #Whole matrix
    for x in 2:(trunc(Int64,n/l)-1)
        #First col of B
        for j in 1: l
            #if abs(A[x-1,2][l-1,l-1])<eps
            #    return (A,1)
            #end
            m=A[x,1][j,1]/A[x-1,2][l-1,l-1]

            A[x,1][j,1]=A[x,1][j,1]-m*A[x-1,2][l-1,l-1]
            for k in 1:l
                A[x,2][j,k]=A[x,2][j,k]-m*A[x-1,3][l-1,k]
            end
        end
        #Second col of B
        for j in 1: l
            #if abs(A[x-1,2][l,l])<eps
            #    return (A,1)
            #end
            m=A[x,1][j,2]/A[x-1,2][l,l]
            A[x,1][j,2]=A[x,1][j,2]-m*A[x-1,2][l,l]
            for k in 1:l
                A[x,2][j,k]=A[x,2][j,k]-m*A[x-1,3][l,k]
            end
        end
        #Rest
        for i in 1:(l-1)
            for j in (i+1):l
                #if abs(A[x,2][i,i])<eps
                #    return (A,1)
                #end
                m=A[x,2][j,i]/A[x,2][i,i]
                for k in i:l
                    A[x,2][j,k]=A[x,2][j,k]-m*A[x,2][i,k]
                end
                for k in 1:l
                    A[x,3][j,k]=A[x,3][j,k]-m*A[x,3][i,k]
                end
                b[trunc(Int64,x*l+j)]=b[trunc(Int64,x*l+j)]-m*b[trunc(Int64,x*l+j)-1]
            end
            it=it+1
        end
    end
    #Last iteration
    #First col of B
    for j in 1: l
        m=A[l,1][j,1]/A[l-1,2][l-1,l-1]
        #if abs(A[l-1,2][l-1,l-1])<eps
        #    return (A,1)
        #end
        A[l,1][j,1]=A[l,1][j,1]-m*A[l-1,2][l-1,l-1]
        for k in 1:l
            A[l,2][j,k]=A[l,2][j,k]-m*A[l-1,3][l-1,k]
        end
    end
    #Second col of B
    for j in 1: l
        #if abs(A[l-1,2][l,l])<eps
        #    return (A,1)
        #end
        m=A[l,1][j,2]/A[l-1,2][l,l]
        A[l,1][j,2]=A[l,1][j,2]-m*A[l-1,2][l,l]
        for k in 1:l
            A[l,2][j,k]=A[l,2][j,k]-m*A[l-1,3][l,k]
        end
    end
    #Rest
    for i in 1:(l-1)
        for j in (i+1):l
            #if abs(A[l,2][i,i])<eps
            #    return (A,1)
            #end
            m=A[l,2][j,i]/A[l,2][i,i]
            for k in i:l
                A[l,2][j,k]=A[l,2][j,k]-m*A[l,2][i,k]
            end
            b[trunc(Int64,(n/l-1)+j)]=b[trunc(Int64,(n/l-1)+j)]-m*b[trunc(Int64,(n/l-1)+j-1)]
        end
    end
    it=n
    println(it)
    for i in trunc(Int64,n/l):-1:1
        for j in l:-1:1
            s=b[it]
            for k in l:-1:1
                s=s-A[i,2][k,j]*b[k]
            end
            for k in l:-1:1
                s=s-A[i,3][k,j]*b[k]
            end
            b[it]=s/A[i,2][j,j]
            it=it-1
        end
    end
    return (A,b,0)
end
