#author: Tymoteusz Surynt

module ObliczeniaNaukowe5
    export importMatrix, printfMatrix, importVector, gaussianElimination, simpleGauss
    eps = 1/(10^5)
    #Function responsible for importing the A matrix
    #fileLocation - string containing the file from which we want to import the matrix
    function importMatrix(fileLocation)
        open(fileLocation) do f
            input = split(readline(f))
            n::Int64=parse(Int64, input[1])
            l::Int64=parse(Int64, input[2])
            #@printf "n: %d l: %d\n" n l
            v::Int64=floor(n/l)
            data=Array(Any,v,3);
            for i in 1:v
                data[i,1]=zeros(Float64,l,2)
                data[i,2]=Array(Float64,l,l)
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
    #Function responsible for importing the b vector
    #fileLocation - string containing the file from which we want to import the vector
    function importVector(fileLocation)
        open(fileLocation) do f
            n::Int64=parse(Int64,readline(f))
            data=Array(Float64, n)
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
        for j::Int64 in 1:n
            @printf " %6.2f" A[trunc(Int64, floor((j-1)/l)+1),1][trunc(Int64, j-floor((j-1)/l)*l),1]
            @printf " %6.2f" A[trunc(Int64, floor((j-1)/l)+1),1][trunc(Int64, j-floor((j-1)/l)*l),2]
            @printf(" |")
            for i::Int64 in 1:l
                @printf " %6.2f" A[trunc(Int64, floor((j-1)/l)+1),2][trunc(Int64, j-floor((j-1)/l)*l),i]
            end
            @printf(" |")
            for i::Int64 in 1:l
                @printf " %6.2f" A[trunc(Int64, floor((j-1)/l)+1),3][trunc(Int64, j-floor((j-1)/l)*l),i]
            end
            @printf(" |")
            @printf " %6.2f" b[j]
            if count>=l
                @printf "\n%s\n" "-"^(l*14)
                count=1
            else
                @printf "\n"
                count=count+1
            end
        end
    end
    function gaussianElimination(A, l::Int64, n::Int64)
        v=n/l;
    end
    function simpleGauss(n,A)
        for i in 1:(n-1)
            for j in (i+1):n
                m=A[j,i]/A[i,i]
                for k in i:n
                    A[j,k]=A[j,k]-m*A[i,k]
                end
            end
        end
        return A
    end

end
