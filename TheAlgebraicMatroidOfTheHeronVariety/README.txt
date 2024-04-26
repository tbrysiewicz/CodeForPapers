The code in this directory corresponds to the paper
[1] "The Algebraic Matroid of the Heron Variety"
by Seth K. Asante, Taylor Brysiewicz, and Michelle Hatzel
https://arxiv.org/pdf/2401.06286

Please contact Taylor Brysiewicz with any questions or concerns at tbrysiew@uwo.ca



ProduceOrbits.jl
	Dependencies
		GAP's images package
		OSCAR
		Combinatorics
	This file contains code to produce the \beta_n candidate basis
	orbits of the n-th Heron variety. This code is an implementation
	of Algorithm 1 in [1], along with an evaluation of Algorithm 1
	on each candidate f-vector, which produces the exhaustive list of
	candidate bases.
	
	The results of this computation are saved in /Results as 
		/Results/CandidateBasesX3.txt
	and
		/Results/CandidateBasesX4.txt
		
	These computations support
		Theorem 4.1

DetermineBases.jl
	Dependencies
		HomotopyContinuation.jl
		Combinatorics.jl
	This file contains implementations of 
		Algorithm 2 - IsBasis
			Determines symbolically, whether a subset of 
			coordinates is a basis for the algebraic
			matroid of a parametrized variety.
			No false positives
			No false negatives
		Algorithm 3 - IsBasisMC
			Attempts to do the same thing as IsBasis by
			evaluating the Jacobian and a random rational
			number. This is a one-sided randomized Monte-
			Carlo algorithm. 
			Returns true: input is indeed a basis
				or
			Returns nothing: no conclusion drawn
		Algorithm 4 - IsBasisBKK
			Also a symbolic computation, but is one-sided
			in the opposite direction. Ascertains whether
			the input is a basis by checking whether the
			(affine) mixed volume of the system whose 
			solutions are the fibres of the projection map
			is zero. If zero, there can be no solutions.
			This is the BKK theorem (adapted to handle
			non-torus solutions).
			Returns false: input is not a basis
				or
			Returns nothing: no conclusion drawn
		Lemma 4.4 - InductiveNonBasis
			Returns true if a given candidate basis orbit
			is a non-basis because it contains a dependent
			subset of positive-dimensional faces which are
			supported on some face of the n-simplex
	The following function puts them all together
		DetermineBases(n,CandidateBases)
			Reads the candidate bases from the file which
			should have been produced by ProduceOrbits.jl
			and first performs Algorithm 3 on each twice.
			Those for which Algorithm 3 draws no conclusion (twice),
			Algorithm 4 is evaluated. Those which still
			have no conclusion are decided to be bases or
			not bases via Algorithm 2.
			It optionally writes an indicator file in /Results
			which is necessary for subsequent computations.
	
	The results of these computations are saved in /Results as 
		/Results/BasisIndicatorX3.txt
		/Results/BasisIndicatorX4.txt
	These computations support
		Theorem 4.2
		Theorem 4.3
		Proposition 4.6
		Theorem 4.8
		Corollary 4.9
		

DegreesAndMonodromyGroups
	Dependencies
		Pandora.jl
			HomotopyContinuation.jl
			OSCAR
	The software system Pandora.jl contains implementations of 
	Algorithms 5 and 6 in [1]. This file creates, for each
	basis orbit of the Heron variety, an "EnumerativeProblem"
	which is the main data type of Pandora.jl. It solves these
	systems using HomotopyContinuation.jl so counting the 
	solutions gives the degree of the system. It then applies
	Pandora's versions of Algorithms 5 and 6 (of [1]) to find
	the Galois/monodromy groups as well as the coordinate
	symmetry groups of each basis of the third Heron variety. 
	
	These computations support
		Theorem 4.12*
		Theorem 4.13*
		Corollary 4.15*
	
Experiments in the paper can be reproduced using Pandora.jl's "Explore" function
		
