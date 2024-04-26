using GAP
GAP.Globals.LoadPackage(GapObj("images"))
using Combinatorics
using Oscar

#= 
This code snippet corresponds to the paper
[1] "The Algebraic Matroid of the Heron Variety"
by Seth K. Asante, Taylor Brysiewicz, and Michelle Hatzel
https://arxiv.org/pdf/2401.06286

Please contact Taylor Brysiewicz with any questions or concerns at tbrysiew@uwo.ca
=#

#Returns number of edges on an n-simplex
function e(n)
	binomial(n+1,2)
end

#Returns number of positive dimensional faces of an n-simplex
function N(n)
	2^(n+1)-(n+1)-1
end

function orb(G,A)
    A = sort([sort(a) for a in A])
    O = sort(collect(orbit(G,on_sets_sets,A)))
    return(O)
end


#Given a collection of faces of an n-simplex (each indexed by its vertices)
#the symmetric group S_(n+1) acts on the collection by permuting the vertices
#of the n-simplex. This function is a well-defined map from orbits of this 
#action to orbit-representatives. 
function MinimalImage(CB,n)
	G = symmetric_group(n+1)
	m=GAP.Globals.MinimalImage(GAP.julia_to_gap(G),GapObj(Array([GapObj(Array(q)) for q in CB])),GAP.Globals.OnSetsSets)
	K=[GAP.gap_to_julia(k) for k in GAP.gap_to_julia(m)]
	return(Vector{Vector{Int64}}(K))
end

#Returns all N(n) positive dimensional faces of an n-simplex
function AllPosFaces(n)
	F = []
	for i in 2:n+1
		for f in combinations(1:n+1,i)
			push!(F,f)
		end
	end
	return(F)
end

#Returns an iterator for all k-subsets of the positive dimensional faces of an n-simplex
function AllSizeKOrbits(n,k)
	return(combinations(AllPosFaces(n),k))
end

#This returns true if and only if each element of a is less than or equal to its
# counterpart in b
function StrongLessEq(a,b)
	possiblyLess=true
	counter=1
	while possiblyLess==true && counter<=length(a)
		if a[counter]>b[counter]
			return(false)
		else
			counter=counter+1
		end
	end
	return(true)
end


#Define the f-vector of a collection of positive dimensional faces as the counts
# of the dimensions of faces in that collection. 
function FVector(S,n)
	return([count(x->length(x)==i,S) for i in 2:n+1])
end


# This function returns the possible f-vectors of bases of the n-th Heron variety
function CandidateFVectors(n)
	AllSignatures=[]
	FacesByDimension = [combinations(1:n+1,i) for i in 2:n+1]
	#This vector bounds the nubmer of faces in each dimension
	fvec = [length(a) for a in FacesByDimension] 
	#This collects all candidate partitions (i<->i-1 so we partition e(n)+n)
	P = Oscar.partitions(e(n)+n,n)
	for p in P #for each partition
		for q in unique(collect(Combinatorics.permutations(p))) #and each permutation of it
			#check if it is bounded by fvec
			if StrongLessEq(q.-1,fvec) #if so, include it as possible
				push!(AllSignatures, q.-1)
			end
		end
	end
	return(AllSignatures)
end

#Given a partition, we will compute minimal images for all 
#  faces of a particular dimension, and then build up from there
# e.g. two edges of a simplex can either touch or not: there are two such orbits.
#The following function determines which dimension to perform this reduction in
# by a heuristic.
function DetermineTrimDim(P,n)
	fv = [binomial(n+1,i) for i in 2:n+1]
	Fv = [binomial(fv[i],P[i]) for i in 1:n]
	return(findfirst(x->x==max(Fv...),Fv))
end


#####################################################
####		Algorithm 1 in [1]		 ####
#####################################################
#Input: n and an f-vector P of dim(X_n) faces of n
#Output: a representative of each S_n+1 orbit of such a subset of faces
function OrbitsFromFVector(P,n;verbose = false)
	if verbose
		println("\n Computing all ways of choosing ",sum(P)," faces according to the dimension partition ",P)
	end
	G = symmetric_group(n+1)
	#This the dimension which has the most faces involved
	maxDim = findfirst(x->x==max(P...),P)
	#The i-th entry of this list is the set of all i-faces of the n-simplex
	#  represented by the indices of vertices (1...n+1) included in them
	FacesByDimension = [collect(combinations(1:n+1,i)) for i in 2:n+1]
	MainDimensionChoices = combinations(FacesByDimension[maxDim],P[maxDim])
	StartOrbits = [orb(G,S)[1] for S in MainDimensionChoices]
	unique!(StartOrbits)
	if verbose
		println("There are ",length(MainDimensionChoices), " ways to choose ",P[maxDim]," faces of the main dimension ",maxDim)
		println("  These occur in ",length(StartOrbits), " orbit(s).")
	end
	FullCombinations = StartOrbits
	RemainingDimensions = filter(x->x!=maxDim,1:n)
	#We now exhaustively choose P[i] i-dimensional faces of Delta_n for each remaining
	#  dimension.
	for d in RemainingDimensions
		NewCombinations=[]
		C = collect(combinations(FacesByDimension[d],P[d]))
		for f in FullCombinations
			for c in C
				push!(NewCombinations,vcat(f,c))
			end
		end
		FullCombinations = NewCombinations
	end
	FullOrbits = [orb(G,S)[1] for S in FullCombinations]
	unique!(FullOrbits)
	if verbose
		println("  Including all choices for the remaining faces yields a total of ",length(FullOrbits), " orbit(s)")
	end
	return(FullOrbits)
end


#This function produces all candidate f-vectors for bases of the algebraic
#  matroid of the n-ther Heron variety and calls Algorithm 1 in [1] on each
#  to collect orbits for all candidate bases.
#Optionally, one may save these orbits to a file. 
function ProduceHeronOrbits(n; SaveFile=false, verbose=false)
	P = CandidateFVectors(n)
	HeronianOrbits=[]
	if verbose==true
		println("Must check ",length(P), " candidate f-vectors")
	end
	for i in 1:length(P)
		if verbose==true
			println("f-vector ",i,": ",P[i])
		end
		FVecOrbit = OrbitsFromFVector(P[i],n;)
		for O in FVecOrbit
			push!(HeronianOrbits,O)
		end
	end
	if SaveFile==true
		if isdir(pwd()*"/Results")==false
			mkdir(pwd()*"/Results")
		end
		loc = pwd()*"/Results/CandidateBasesX"*string(n)*".txt"
		touch(loc)
		io = open(loc,"w");
		for i in 1:length(HeronianOrbits)
			write(io,string(HeronianOrbits[i]))
			write(io,"\n")
		end
		close(io)
	end
	return(HeronianOrbits)
end



ProduceHeronOrbits(3; SaveFile=true, verbose = true)
ProduceHeronOrbits(4; SaveFile=true, verbose = true)

