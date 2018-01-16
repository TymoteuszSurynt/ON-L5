#author: Tymoteusz Surynt

module blocksys
    export importMatrix, printfMatrix, printfMatrix2, importVector, gaussElimination, gaussElimination2, exportVectorNoError, exportVectorError, getError
    epsilon=1/10^8
    #Function responsible for importing the A matrix
    #fileLocation - string containing the file from which we want to import the matrix
    function importMatrix(fileLocation)
        open(fileLocation) do f
            input = split(readline(f))
            n::Int64=parse(Int64, input[1])
            l::Int64=parse(Int64, input[2])
            #@printf "n: %d l: %d\n" n l
            v::Int64=floor(n/l)
            data=Array{Any}(v,3);
            for i in 1:v
                data[i,1]=zeros(Float64,l,2)
                data[i,2]=Array{Float64}(l,l)
                data[i,3]=zeros(Float64,l,l)
            end
            offsetA::Int64=0
            offsetB::Int64=0
            offsetC::Int64=l
            line::Int64=0
            x::Int64=0
            y::Int64=0

            while !eof(f)
                input = split(readline(f))
                y=parse(Int64, input[1])
                x=parse(Int64, input[2])
                var=parse(Float64, input[3])
                line=floor((y-1)/l)
                if x>(line+1)*l
                    #@printf "x-offsetC: %d y-line*l: %d x: %d y: %d line: %d \n" x-offsetC y-line*l x y line
                    data[line+1,3][y-line*l,y-line*l]=var
                elseif x>line*l
                    #@printf "x-offsetA: %d y-line*l: %d x: %d y: %d line: %d \n" x-offsetA y-line*l x y line
                    data[line+1,2][y-line*l,x-(line*l)]=var
                elseif x>(line-1)*l
                    #@printf "x-offsetB: %d y-line*l: %d x: %d y: %d line: %d \n" x-((line-1)*l)-(l-2) y-line*l x y line
                    data[line+1,1][y-line*l,x-((line-1)*l)-(l-2)]=var
                end

            end

            return n, l, data
        end
    end
    #Function responsible for calculating the error
    #x- vector x from Ax=b
    #value::Float64 - desired value
    #n::Int64 - size of vector x
    function getError(x, value::Float64, n::Int64)
        sum::Float64=0;
        for i in x
            sum=sum+abs(i-value)
        end
        sum=sum/n
        return sum
    end
    #Function responsible for exporting the X vector
    #fileLocation - string containing the file from where you want to save the vector
    function exportVectorNoError(fileLocation, x)
        open(fileLocation, "w") do f

            for i in x
                write(f,"$i\n")
            end
        end
    end
    #Function responsible for exporting the X vector
    #fileLocation - string containing the file from where you want to save the vector
    function exportVectorError(fileLocation, x, err)
        open(fileLocation, "w") do f
            write(f,"$err\n")
            for i in x
                write(f,"$i\n")
            end
        end
    end
    #Function responsible for importing the b vector
    #fileLocation - string containing the file from which we want to import the vector
    function importVector(fileLocation)
        open(fileLocation) do f
            n::Int64=parse(Int64,readline(f))
            data=Array{Float64}(n)
            i::Int64=1
            while !eof(f)
                var=parse(Float64, readline(f))
                data[i]=var
                i=i+1
            end

            return n, data
        end
    end

    #Function responsible for printing out the A matrix
    #A- matrix
    #l::Int64 - the size of the block matrix
    #n::Int64 - the size of A matrix
    function printfMatrix(A,b, l::Int64, n::Int64)
        count=1
        @printf "\n%s\n" "-"^(l*27)
        for j::Int64 in 1:n
            @printf("| ")
            @printf " %8.2f" A[trunc(Int64, floor((j-1)/l)+1),1][trunc(Int64, j-floor((j-1)/l)*l),1]
            @printf " %8.2f" A[trunc(Int64, floor((j-1)/l)+1),1][trunc(Int64, j-floor((j-1)/l)*l),2]
            @printf(" |")
            for i::Int64 in 1:l
                @printf " %8.2f" A[trunc(Int64, floor((j-1)/l)+1),2][trunc(Int64, j-floor((j-1)/l)*l),i]
            end
            @printf(" |")
            for i::Int64 in 1:l
                @printf " %8.2f" A[trunc(Int64, floor((j-1)/l)+1),3][trunc(Int64, j-floor((j-1)/l)*l),i]
            end
            @printf(" |")
            @printf " %8.2f" b[j]
            @printf(" |")
            if count>=l
                @printf "\n%s\n" "-"^(l*27)
                count=1
            else
                @printf "\n"
                count=count+1
            end
        end
    end
    #Function responsible for printing out the A matrix
    #A- matrix
    #l::Int64 - the size of the block matrix
    #n::Int64 - the size of A matrix
    function printfMatrix2(A,b, l::Int64, n::Int64)
        count=1
        @printf "\n%s\n" "-"^(l*27)
        for j::Int64 in 1:n
            @printf("| ")
            @printf " %8.2f" A[trunc(Int64, floor((j-1)/l)+1),1][trunc(Int64, j-floor((j-1)/l)*l),1]
            @printf " %8.2f" A[trunc(Int64, floor((j-1)/l)+1),1][trunc(Int64, j-floor((j-1)/l)*l),2]
            @printf(" |")
            for i::Int64 in 1:l
                @printf " %8.2f" A[trunc(Int64, floor((j-1)/l)+1),2][trunc(Int64, j-floor((j-1)/l)*l),i]
            end
            @printf(" |")
            for i::Int64 in 1:l
                @printf " %8.2f" A[trunc(Int64, floor((j-1)/l)+1),3][trunc(Int64, j-floor((j-1)/l)*l),i]
            end
            @printf(" |")
            for i::Int64 in 1:l
                @printf " %8.2f" A[trunc(Int64, floor((j-1)/l)+1),4][trunc(Int64, j-floor((j-1)/l)*l),i]
            end
            @printf(" |")
            @printf " %8.2f" b[j]
            @printf(" |")
            if count>=l
                @printf "\n%s\n" "-"^(l*27)
                count=1
            else
                @printf "\n"
                count=count+1
            end
        end
    end
    function swapRow(ra::Int64, a::Int64,rb::Int64, b::Int64,A,v,l::Int64)

        if ra==rb
            for i in 1:l
                temp=A[ra,2][a,i]
                A[ra,2][a,i]=A[rb,2][b,i]
                A[rb,2][b,i]=temp

                temp=A[ra,3][a,i]
                A[ra,3][a,i]=A[rb,3][b,i]
                A[rb,3][b,i]=temp

                temp=A[ra,4][a,i]
                A[ra,4][a,i]=A[rb,4][b,i]
                A[rb,4][b,i]=temp

            end
            temp = v[(ra-1)*l+a]
            v[(ra-1)*l+a]=v[(rb-1)*l+b]
            v[(rb-1)*l+b]=temp
        else
            for i in 1:l-2
                A[ra,2][a,i]=0
            end
            for i in (l-1):l
                temp=A[ra,2][a,i]
                A[ra,2][a,i]=A[rb,1][b,i-(l-2)]
                A[rb,1][b,i-(l-2)]=temp
            end

            for i in 1:l
                temp=A[ra,3][a,i]
                A[ra,3][a,i]=A[rb,2][b,i]
                A[rb,2][b,i]=temp

                A[ra,4][a,i]=A[rb,3][b,i]
                A[rb,3][b,i]=0


            end
            temp = v[(ra-1)*l+a]
            v[(ra-1)*l+a]=v[(rb-1)*l+b]
            v[(rb-1)*l+b]=temp
        end
    end
    #Funkcja odpowiedzialna za liczenie wektora X z równania Ax=b
    #n::Int64 - wielkość macierzy
    #l::Int64 - wielkość pojedynczego bloku
    #AB - macierz A z równania Ax=b
    #v - wektor b z równania Ax=b
    function gaussElimination(n::Int64,l::Int64,AB,v)
        A=deepcopy(AB)
        b=deepcopy(v)
        X=Array{Float64}(n)
        it::Int64=2
        #First iteration
        for i in 1:(l-1)
            for j in (i+1):l
                if abs(A[1,2][i,i]) <epsilon
                    return (A,b,X,1)
                end
                m=A[1,2][j,i]/A[1,2][i,i]
                for k in i:l
                    A[1,2][j,k]=A[1,2][j,k]-m*A[1,2][i,k]
                end
                for k in 1:l
                    A[1,3][j,k]=A[1,3][j,k]-m*A[1,3][i,k]
                end
                b[j]=b[j]-m*b[i]
            end
        end
        #Whole matrix
        for x in 2:(trunc(Int64,n/l)-1)
            #First col of B
            for j in 1: l
                if abs(A[x-1,2][l-1,l-1]) <epsilon
                    return (A,b,X,1)
                end
                m=A[x,1][j,1]/A[x-1,2][l-1,l-1]
                A[x,1][j,2]=A[x,1][j,2]-m*A[x-1,2][l-1,l]
                A[x,1][j,1]=A[x,1][j,1]-m*A[x-1,2][l-1,l-1]
                for k in 1:l
                    A[x,2][j,k]=A[x,2][j,k]-m*A[x-1,3][l-1,k]
                end
                b[(x-1)*l+j]=b[(x-1)*l+j]-m*b[(x-2)*l+l-1]
            end
            #Second col of B
            for j in 1: l
                if abs(A[x-1,2][l,l]) <epsilon
                    return (A,b,X,1)
                end
                m=A[x,1][j,2]/A[x-1,2][l,l]
                A[x,1][j,2]=A[x,1][j,2]-m*A[x-1,2][l,l]
                for k in 1:l
                    A[x,2][j,k]=A[x,2][j,k]-m*A[x-1,3][l,k]
                end
                b[(x-1)*l+j]=b[(x-1)*l+j]-m*b[(x-2)*l+l]
            end
            #Rest
            for i in 1:(l-1)
                for j in (i+1):l
                    if abs(A[x,2][i,i]) <epsilon
                        return (A,b,X,1)
                    end
                    m=A[x,2][j,i]/A[x,2][i,i]
                    for k in i:l
                        A[x,2][j,k]=A[x,2][j,k]-m*A[x,2][i,k]
                    end
                    for k in 1:l
                        A[x,3][j,k]=A[x,3][j,k]-m*A[x,3][i,k]
                    end
                    b[(x-1)*l+j]=b[(x-1)*l+j]-m*b[(x-1)*l+i]
                end
                it=it+1
            end
        end
        #Last iteration
        #First col of B
        for j in 1: l
            if abs(A[trunc(Int64,n/l)-1,2][l-1,l-1]) <epsilon
                return (A,b,X,1)
            end
            m=A[trunc(Int64,n/l),1][j,1]/A[trunc(Int64,n/l)-1,2][l-1,l-1]
            A[trunc(Int64,n/l),1][j,1]=A[trunc(Int64,n/l),1][j,1]-m*A[trunc(Int64,n/l)-1,2][l-1,l-1]
            A[trunc(Int64,n/l),1][j,2]=A[trunc(Int64,n/l),1][j,2]-m*A[trunc(Int64,n/l)-1,2][l-1,l]
            for k in 1:l
                A[trunc(Int64,n/l),2][j,k]=A[trunc(Int64,n/l),2][j,k]-m*A[trunc(Int64,n/l)-1,3][l-1,k]
            end
            b[n-l+j]=b[n-l+j]-m*b[n-l-1]
        end
        #Second col of B
        for j in 1: l
            if abs(A[trunc(Int64,n/l)-1,2][l,l]) <epsilon
                return (A,b,X,1)
            end
            m=A[trunc(Int64,n/l),1][j,2]/A[trunc(Int64,n/l)-1,2][l,l]
            A[trunc(Int64,n/l),1][j,2]=A[trunc(Int64,n/l),1][j,2]-m*A[trunc(Int64,n/l)-1,2][l,l]
            for k in 1:l
                A[trunc(Int64,n/l),2][j,k]=A[trunc(Int64,n/l),2][j,k]-m*A[trunc(Int64,n/l)-1,3][l,k]
            end
            b[n-l+j]=b[n-l+j]-m*b[n-l]
        end
        #Rest
        for i in 1:(l-1)
            for j in (i+1):l
                if abs(A[trunc(Int64,n/l),2][i,i]) <epsilon
                    return (A,b,X,1)
                end
                m=A[trunc(Int64,n/l),2][j,i]/A[trunc(Int64,n/l),2][i,i]
                for k in i:l
                    A[trunc(Int64,n/l),2][j,k]=A[trunc(Int64,n/l),2][j,k]-m*A[trunc(Int64,n/l),2][i,k]
                end
                b[n-l+j]=b[n-l+j]-m*b[n-l+i]
            end
        end
        it=n
        #First iteration
        for i in l:-1:1
            s=b[it]
            if abs(A[trunc(Int64,n/l),2][i,i]) <epsilon
                return (A,b,X,1)
            end
            for j in l:-1:(i+1)
                s=s-A[trunc(Int64,n/l),2][i,j]*X[n-l+j]
            end
            X[it]=s/A[trunc(Int64,n/l),2][i,i]
            it=it-1
        end

        #Rest
        offsetA::Int64=n-2*l
        offsetC::Int64=n-l
        for i in trunc(Int64,n/l-1):-1:1
            for j in l:-1:1
                s=b[it]
                if abs(A[i,2][j,j]) <epsilon
                    return (A,b,X,1)
                end
                for k in l:-1:(j+1)
                    s=s-A[i,2][j,k]*X[offsetA+k]
                end
                for k in j:-1:1
                    s=s-A[i,3][j,k]*X[offsetC+k]
                end
                X[it]=s/A[i,2][j,j]
                it=it-1
            end
            offsetA=offsetA-l
            offsetC=offsetC-l
        end
        return (A,b,X,0)
    end

    #Funkcja odpowiedzialna za liczenie wektora X z równania Ax=b
    #n::Int64 - wielkość macierzy
    #l::Int64 - wielkość pojedynczego bloku
    #AB - macierz A z równania Ax=b
    #v - wektor b z równania Ax=b
    function gaussElimination2(n::Int64,l::Int64,AB,v)
        A=Array{Any}(trunc(Int64,floor(n/l)),4);
        for i in 1:trunc(Int64,floor(n/l))
            A[i,1]=deepcopy(AB[i,1])
            A[i,2]=deepcopy(AB[i,2])
            A[i,3]=deepcopy(AB[i,3])
            A[i,4]=zeros(Float64,l,l)
        end
        b=deepcopy(v)
        X=Array{Float64}(n)
        #Whole matrix
        for x in 1:(trunc(Int64,n/l)-1)
            for i in 1:l
                max=(i,i,2,x,0)
                if i== l-1
                    if abs(A[x,2][l,i]) > abs(A[x,2][i,i])
                        max=(l,i,2,x,1)
                    end
                    for j in 1:l
                        if abs(A[x+1,1][j,1]) > abs(A[max[4],max[3]][max[1],max[2]])
                            max=(j,1,1,x+1,1)
                        end
                    end
                    if max[5]!=0
                        swapRow(x,i,max[4],max[1],A,b,l)
                    end
                    for j in (i+1):l
                        if abs(A[x,2][i,i]) <epsilon
                            return (A,b,X,1)
                        end
                        m=A[x,2][j,i]/A[x,2][i,i]
                        for k in i:l
                            A[x,2][j,k]=A[x,2][j,k]-m*A[x,2][i,k]
                        end
                        for k in 1:l
                            A[x,3][j,k]=A[x,3][j,k]-m*A[x,3][i,k]
                        end
                        for k in 1:l
                            A[x,4][j,k]=A[x,4][j,k]-m*A[x,4][i,k]
                        end
                        b[(x-1)*l+j]=b[(x-1)*l+j]-m*b[(x-1)*l+i]
                    end
                    #First col of B
                    for j in 1: l
                        if abs(A[x,2][l-1,l-1]) <epsilon
                            return (A,b,X,1)
                        end
                        m=A[x+1,1][j,1]/A[x,2][l-1,l-1]
                        A[x+1,1][j,1]=A[x+1,1][j,1]-m*A[x,2][l-1,l-1]
                        A[x+1,1][j,2]=A[x+1,1][j,2]-m*A[x,2][l-1,l]
                        for k in 1:l
                            A[x+1,2][j,k]=A[x+1,2][j,k]-m*A[x,3][l-1,k]
                        end
                        for k in 1:l
                            A[x+1,3][j,k]=A[x+1,3][j,k]-m*A[x,4][l-1,k]
                        end
                        b[x*l+j]=b[x*l+j]-m*b[x*l-1]
                    end
                elseif i==l
                    for j in 1:l
                        if abs(A[x+1,1][j,2]) > abs(A[max[4],max[3]][max[1],max[2]])
                            max=(j,2,1,x+1,1)
                        end
                    end
                    if max[5]!=0
                        swapRow(x,i,max[4],max[1],A,b,l)
                    end

                    #Second col of B
                    for j in 1: l
                        if abs(A[x,2][l,l]) <epsilon
                            return (A,b,X,1)
                        end
                        m=A[x+1,1][j,2]/A[x,2][l,l]
                        A[x+1,1][j,2]=A[x+1,1][j,2]-m*A[x,2][l,l]
                        for k in 1:l
                            A[x+1,2][j,k]=A[x+1,2][j,k]-m*A[x,3][l,k]
                        end
                        for k in 1:l
                            A[x+1,3][j,k]=A[x+1,3][j,k]-m*A[x,4][l,k]
                        end
                        b[x*l+j]=b[x*l+j]-m*b[x*l]
                    end
                else
                    for j in (i+1):l
                        if abs(A[x,2][j,i]) > abs(A[x,2][max[1],max[2]])
                            max=(j,i,2,x,1)
                        end
                    end
                    if max[5]!=0
                        swapRow(x,i,max[4],max[1],A,b,l)
                    end
                    for j in (i+1):l
                        if abs(A[x,2][i,i]) <epsilon
                            return (A,b,X,1)
                        end
                        m=A[x,2][j,i]/A[x,2][i,i]
                        for k in i:l
                            A[x,2][j,k]=A[x,2][j,k]-m*A[x,2][i,k]
                        end
                        for k in 1:l
                            A[x,3][j,k]=A[x,3][j,k]-m*A[x,3][i,k]
                        end
                        for k in 1:l
                            A[x,4][j,k]=A[x,4][j,k]-m*A[x,4][i,k]
                        end
                        b[(x-1)*l+j]=b[(x-1)*l+j]-m*b[(x-1)*l+i]
                    end
                end
            end
        end

        x=(trunc(Int64,n/l))
        #Last iteration
        for i in 1:(l-1)
            max=(i,i,2,x,0)
            if i== l-1
                if abs(A[x,2][i,i]) > abs(A[x,2][max[1],max[2]])
                    max=(i,l,2,x,1)
                end
                if max[5]!=0
                    swapRow(x,i,max[4],max[1],A,b,l)
                end
                for j in (i+1):l
                    if abs(A[x,2][i,i]) <epsilon
                        return (A,b,X,1)
                    end
                    m=A[x,2][j,i]/A[x,2][i,i]
                    for k in i:l
                        A[x,2][j,k]=A[x,2][j,k]-m*A[x,2][i,k]
                    end
                    for k in 1:l
                        A[x,3][j,k]=A[x,3][j,k]-m*A[x,3][i,k]
                    end
                    b[(x-1)*l+j]=b[(x-1)*l+j]-m*b[(x-1)*l+i]
                end
            else
                for j in (i+1):l
                    if abs(A[x,2][j,i]) > abs(A[x,2][max[1],max[2]])
                        max=(j,i,2,x,1)
                    end
                end
                if max[5]!=0
                    swapRow(x,i,max[4],max[1],A,b,l)
                end
                for j in (i+1):l
                    if abs(A[x,2][i,i]) <epsilon
                        return (A,b,X,1)
                    end
                    m=A[x,2][j,i]/A[x,2][i,i]
                    for k in i:l
                        A[x,2][j,k]=A[x,2][j,k]-m*A[x,2][i,k]
                    end
                    for k in 1:l
                        A[x,3][j,k]=A[x,3][j,k]-m*A[x,3][i,k]
                    end
                    b[(x-1)*l+j]=b[(x-1)*l+j]-m*b[(x-1)*l+i]
                end
            end
        end

        it=n
        #First iteration
        for i in l:-1:1
            s=b[it]
            if abs(A[trunc(Int64,n/l),2][i,i]) <epsilon
                return (A,b,X,1)
            end
            for j in l:-1:(i+1)
                s=s-A[trunc(Int64,n/l),2][i,j]*X[n-l+j]
            end
            X[it]=s/A[trunc(Int64,n/l),2][i,i]
            it=it-1
        end

        #Rest
        offsetA::Int64=n-2*l
        offsetC::Int64=n-l
        offsetD::Int64=n
        for i in trunc(Int64,n/l-1):-1:1
            for j in l:-1:1
                s=b[it]
                if abs(A[i,2][j,j]) <epsilon
                    return (A,b,X,1)
                end
                for k in l:-1:(j+1)
                    s=s-A[i,2][j,k]*X[offsetA+k]
                end
                for k in l:-1:1
                    s=s-A[i,3][j,k]*X[offsetC+k]
                end
                if offsetD<n
                    for k in l:-1:1
                        s=s-A[i,4][j,k]*X[offsetD+k]
                    end
                end
                X[it]=s/A[i,2][j,j]
                it=it-1
            end
            offsetA=offsetA-l
            offsetC=offsetC-l
            offsetC=offsetD-l
        end
        return (A,b,X,0)
    end
end
