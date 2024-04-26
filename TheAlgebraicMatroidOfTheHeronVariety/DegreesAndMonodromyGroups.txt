using Pandora
#Find Pandora at 
#https://github.com/tbrysiewicz/Pandora
#This code was executed on 4/26/2024 after commit ffd6080
#The other dependency info is below
#=
  [aaaa29a8] Clustering v0.15.7
  [861a8166] Combinatorics v1.0.2
  [c863536a] GAP v0.9.8
  [f213a82b] HomotopyContinuation v2.9.2
  [f1435218] Oscar v0.13.0
  [91a5bcdd] Plots v1.39.0
  [49802e3a] ProgressBars v1.5.1
  [295af30f] Revise v3.5.14
  [37e2e46d] LinearAlgebra
  [44cfe95a] Pkg v1.10.0
=#



using LinearAlgebra
using HomotopyContinuation
using Oscar

#= 
This code snippet corresponds to the paper
[1] "The Algebraic Matroid of the Heron Variety"
by Seth K. Asante, Taylor Brysiewicz, and Michelle Hatzel
https://arxiv.org/pdf/2401.06286

Please contact Taylor Brysiewicz with any questions or concerns at tbrysiew@uwo.ca
=#

function GenericEDMatrix(n)
    @var x[1:n+1,1:n+1]
    M = [x[min(i,j),max(i,j)]*(i!=j) for i in 1:n+1, j in 1:n+1]
end


function CayleyMengerDeterminant(D;flag="nothing")
    n=size(D)[1]
    newRow = [1 for i in 1:n]
    newCol = vcat([1 for i in 1:n],0)
    M = hcat(vcat(D,newRow'),newCol)
    C = (-1)^(n)//((factorial(n-1))^2*2^(n-1))
    if flag=="NoConstant"
    	return(det(M))
    else
    	return(C*det(M))
    end
end


function CMDeterminants(D,C; flag = "nothing")
    n = size(D)[1]
    #C = collect(filter(x->length(x)>1,collect(combinations(1:n))))
    DetBucket=[]
    for c in C
        push!(DetBucket,CayleyMengerDeterminant(D[c,c];flag=flag))
    end
    return(DetBucket)
end


function ConstructHeronSystemGeneral(B,n)
      X = GenericEDMatrix(n)

      Vals = CMDeterminants(X,B)
      Vars = Vector{Variable}(unique(vcat(X...))[2:end])

      @var y[1:length(B)]
      F = System([Vals[i]-y[i] for i in 1:length(B)],variables = Vars,parameters = y)

      M = GenericEDMatrix(n)
      SVars = Vector{Variable}(unique(vcat(M...))[2:end])
      Msubbed = subs(M,SVars=>randn(ComplexF64,length(SVars)))
      StartSol = Vector{ComplexF64}(unique(vcat(Msubbed...))[2:end])
      StartParams = CMDeterminants(Msubbed,B)
      return(F,StartSol,StartParams)
end

CandidateBases = [[[1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 3], [1, 4], [2, 3], [2, 4]],
[[1, 2], [1, 2, 3], [1, 3], [1, 4], [2, 4], [3, 4]],
[[1, 2], [1, 2, 3, 4], [1, 3], [1, 4], [2, 3], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [1, 4], [2, 3]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [2, 3], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [1, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 3], [1, 4], [2, 3, 4], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [2, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 3], [2, 3, 4], [2, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [1, 3, 4], [2, 3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3], [1, 4], [2, 3]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3], [1, 4], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 4], [2, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3], [2, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3, 4], [2, 3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [1, 3, 4], [1, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [1, 4], [2, 3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [1, 3, 4], [2, 3]],
[[1, 2], [1, 2, 3], [1, 3, 4], [1, 4], [2, 3, 4], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [1, 3, 4], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 4], [1, 3], [2, 3, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3], [1, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3], [1, 4], [2, 3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3], [2, 3]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3, 4], [1, 4], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3], [2, 3, 4], [2, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3, 4], [2, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3], [1, 3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3], [2, 3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 4], [1, 3, 4], [3, 4]],
[[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 3, 4], [1, 4], [2, 3, 4]]]

BI = [1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]

HeronProblems = []
for i in 1:35
    if BI[i]==1
        E = EnumerativeProblem(ConstructHeronSystemGeneral(CandidateBases[i],3)[1])
        degree(E)
        push!(HeronProblems,E)
    else
        push!(HeronProblems,nothing)
    end
end

#Collects Galois Groups
GaloisGroups = [] 
for i in 1:35
    if BI[i]==1
        E = HeronProblems[i]
        push!(GaloisGroups,galois_group(E;nloops = 20))
    else
        push!(GaloisGroups,nothing)
    end
end

#Used to produce all volumes
function EDfrom3Sol(s)
    M = [0 s[1] s[2] s[3]; s[1] 0 s[4] s[5]; s[2] s[4] 0 s[6]; s[3] s[5] s[6] 0]
end

allvols = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4],[1,2,3],[1,2,4],[1,3,4],[2,3,4],[1,2,3,4]]

function VolVectorFromSols(s)
    D = EDfrom3Sol(s)
    CM = CMDeterminants(D,allvols; flag = "nothing")
    return(CM)
end

#Collects Coordinate Symmetry groups
CSG = []
for i in 1:35
    if BI[i]==1
        E = HeronProblems[i]
        push!(CSG,coordinate_symmetry_group(E;F = VolVectorFromSols))
    else
        push!(CSG,nothing)
    end
end


#Prints the cases where the Galois group numerically computed
#  is not equal to the Coordinate Symmetry Group
for i in 1:35
    if BI[i]==1
        if GaloisGroups[i]!=CSG[i]
            println(i,"    ",describe(GaloisGroups[i]),"     ",describe(CSG[i]))
        end
    end
end


