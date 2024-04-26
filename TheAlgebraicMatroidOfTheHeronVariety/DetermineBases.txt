using HomotopyContinuation
using Combinatorics
using LinearAlgebra
using LinearAlgebraX
using Oscar

#= 
This code snippet corresponds to the paper
[1] "The Algebraic Matroid of the Heron Variety"
by Seth K. Asante, Taylor Brysiewicz, and Michelle Hatzel
https://arxiv.org/pdf/2401.06286

Please contact Taylor Brysiewicz with any questions or concerns at tbrysiew@uwo.ca
=#

#######################################
#### For extracting known results #####
#######################################

function LoadOrbits(n;verbose=false)
	loc1 = pwd()*"/Results/CandidateBasesX"*string(n)*".txt"
	R = readlines(loc1)
	CandidateOrbits=[]
	counter=0
	portion=0
	for r in R
		counter=counter+1
		if counter>length(R)/(10)
			portion=portion+10
			counter=0
			if verbose
				println(string(portion)*"% done")
			end
		end
		push!(CandidateOrbits,eval(Meta.parse(r)))
	end
	return(CandidateOrbits)
end

function LoadBasisIndicator(n;verbose=false)
    loc2 = pwd()*"/Results/BasisIndicatorX"*string(n)*".txt"
    basisindicator=[]
    counter=0
    portion=0
    R = readlines(loc2)
    for r in R
        counter=counter+1
        if counter>length(R)/(10)
            portion=portion+10
            counter=0
            if verbose
                println(string(portion)*"% done")
            end
        end

        push!(basisindicator,eval(Meta.parse(r)))
    end
    return(basisindicator)
end

function LoadBasisOrbits(n)
    O = LoadOrbits(n)
    B = LoadBasisIndicator(n)
    return(O[findall(x->x==1,B)])
end


###########################################
#### For constructing Heron Systems   #####
###########################################

#Constructs Generic Euclidean Distance Matrix
function GenericEDMatrix(n)
	@var x[1:n+1,1:n+1]
	M = [x[min(i,j),max(i,j)]*(i!=j) for i in 1:n+1, j in 1:n+1]
end

#Pads an ED matrix with appropriate rows/columns so the determinant
#  is a Cayley Menger determinant. Can be performed with numbers or 
#  with indeterminants. The constant which makes the determinant a
#  Euclidean volume is optional.
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
       
#For each subset of indices collected in C, produce
#the CM determinants of D indexed by those.     
function CMDeterminants(D,C; flag = "nothing")
	n = size(D)[1]
	DetBucket=[]
	for c in C
		push!(DetBucket,CayleyMengerDeterminant(D[c,c];flag=flag))
	end
	return(DetBucket)
end

#Construct the system which parametrizes the Heron variety as 
# a HomotopyContinuation.jl system (so we can easily access jacobian)
function HeronParametrization(n; flag = "nothing")
      D = GenericEDMatrix(n+1)
      C = collect(filter(x->length(x)>1,collect(combinations(1:n+1))))
      F = System(CMDeterminants(D,C; flag = flag))
      return(F)
end

#Consruct the Jacobian of the Heron parametrization
function HeronJacobian(n; flag = "nothing")
      F = HeronParametrization(n; flag = flag)
      return(jacobian(F))
end


#Produce the parametrized polynomial system (enumerative problem)
#  with starting parameter and solution, corresponding to the 
#  branched cover associated to B from X_n
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


function BasisBKK(B,n)
	F = ConstructHeronSystemGeneral(B,n)[1]
	HomotopyContinuation.mixed_volume(F)
end


# IsBasis

function isBasis(B,n;J=nothing)
	#If the Jacobian is not already passed, construct it
	if J == nothing
		J = HeronJacobian(n; flag = "NoConstant")
	end
	#This is the order of the volumes
	C = collect(filter(x->length(x)>1,collect(combinations(1:n+1))))
	I = findall(x->x∈B,C)
	#So now I gives the indices of the columns of the jacobian associated to B
	SubMat = J[I,1:end] #SubMat is the correpsonding submatrix
	D = LinearAlgebraX.detx(SubMat) #D is the corresponding determinant
	if D==0
		return(false) #This occurs if D is obviously zero
	else #otherwise, D is in a non-expanded form which is not obviously zero
		ExpandedD = expand(D)
		if ExpandedD ==0
			return(false)
		else
			return(true)
		end
	end
end
	


# IsBasisMC - Evaluates the Jacobian at random integers. If it produces independent columns indexed by B,
#  we quickly determine the columns are independent when left as indeterminants, but the converse is not true
#  We amplify this one-sided monte-carlo method once. 
function isBasisMC(B,n; J=nothing, verbose=false)
	#If the Jacobian is not already passed, construct it
	if J == nothing
		J = HeronJacobian(n; flag = "NoConstant")
	end
	
	#This is the order of the volumes
	C = collect(filter(x->length(x)>1,collect(combinations(1:n+1))))
	I = findall(x->x∈B,C)
	#So now I gives the indices of the columns of the jacobian associated to B
	SubMat = J[I,1:end] #SubMat is the correpsonding submatrix
	
	for i in 1:2
		DQ = LinearAlgebraX.detx(subs(SubMat,variables(SubMat)=>rand(Int64,length(variables(SubMat)))))
		if DQ != 0
			return(true)
		end
	end
end


# IsBasisBKK

function isBasisBKK(B,n)
	BKK = BasisBKK(B,n)
	if BKK==0
		return(false)
	end
end


B10 = [[1, 2], [1, 2, 3], [1, 3], [2, 3, 4], [2, 4], [3, 4]]
G = symmetric_group(5)
AllB10InX4 = [sort([sort(a) for a in C10]) for C10 in collect(orbit(G,on_sets_sets,B10))]


function Restrict(B,S)
	filter(x->issubset(x,S),B)
end

# InductiveNonBasis
# Warning: this only works for X4->X3

function InductiveNonBasis(B)
	#First we check if there are any full triangles
	for T in collect(combinations(1:5,3))
		R = Restrict(B,T)
		if length(R)==4
			println("Too much support on triangle")
			return(true)
		end
	end
	#The only other way to be dependent in X3 is 
	#  to have more than rank many elements or to
	#  contain a copy of Orbit 10 from [1]
	
	for F in collect(combinations(1:5,4))
		R = Restrict(B,F)
		if length(R)>6
			println("Too much support on facet")
			return(true)
		end
		R = sort([sort(r) for r in R])
		if in(R,AllB10InX4)
			println("B10 copy")
			return(true)
		end
	end
	print(".")
end



function DetermineBases(n; SaveFile=false)
	O = LoadOrbits(n; verbose=true)
	BI = [2 for i in 1:length(O)]
	println("Determining many bases via one-sided evaluation test")
	Bases = []
	MCcount=0
	J = HeronJacobian(n; flag = "NoConstant")
	PotentialNonBasisIndices=[]
	for i in 1:length(O)
		println(i)
		B = O[i]
		if isBasisMC(B,n; J=J, verbose=false)==true
			push!(Bases,B)
			BI[i]=1
		else
			push!(PotentialNonBasisIndices,i)
		end
	end
	MCcount = length(Bases)
	println("This has established ",length(Bases)," many bases")
	NonBases = []
	FullCheckList = []
	println("Now checking BKK nonbases")
	for i in PotentialNonBasisIndices
		println(i)
		B = O[i]
		if isBasisBKK(B,n) == false
			push!(NonBases,B)
			BI[i]=0
		else
			push!(FullCheckList,i)
		end
	end
	BKKcount = length(NonBases)
	println("This has established ",length(NonBases), " many non-bases")
	println("That leaves ",length(FullCheckList), " potential bases/nonbases which must be checked explicitly")
	
	for i in FullCheckList
		B = O[i]
		if isBasis(B,n;J=J) == true
			push!(Bases,B)
			BI[i]=1
		else
			push!(NonBases,B)
			BI[i]=0
		end
	end
	if SaveFile==true
		if isdir(pwd()*"/Results")==false
			mkdir(pwd()*"/Results")
		end
		loc = pwd()*"/Results/BasisIndicatorX"*string(n)*".txt"
		touch(loc)
		io = open(loc,"w");
		for i in 1:length(O)
			write(io,string(BI[i]))
			write(io,"\n")
		end
		close(io)
	end
	println("Bases recognized by Monte Carlo:",MCcount)
	println("Nonbases recognized by BKK:",BKKcount)
	println("Bases:",length(Bases))
	println("NonBases:",length(NonBases))
	
		
end


#Theorems 4.1,4.2,4.3, 4.6 and 4.8
DetermineBases(2; SaveFile=true)
DetermineBases(3; SaveFile=true)
DetermineBases(4; SaveFile=true)
#Cor 4.9
O4 = LoadOrbits(4; verbose=true)
INBrecognized = [InductiveNonBasis(O4[i]) for i in 1:length(O4)]
length(findall(x->x==true,INBrecognized))
