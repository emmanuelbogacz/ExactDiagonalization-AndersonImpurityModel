# NOTE: parameters can be edited in the main, at the bottom of this file

using LinearAlgebra # for matrix operations
using Plots # plotting of the spectral function

# functions are located in this file:
include("ed_functions.jl")


# main function to perform the full ED calculation
function do_ed(Nf::Int, eps::Float64, U::Float64, V::Float64, t::Float64, 
		wmax::Float64, dw::Float64, nu::Float64, T::Float64)

	n::Int = 2Nf + 2 # total number of sites (unfolded rpz)
	impsiteL::Int = Nf + 1 # impurity site left
	impsiteR::Int = Nf + 2 # impurity site right

	basis = BitArray(zeros(Bool, n, 2^n)) # basis states
	basis_state = BitArray(zeros(Bool, size(basis,1))) # basis state


	# STEP 1: generate the basis states 
	counter::Ref{Int} = Ref(0) # can pass in global recursive func
	gen_basis(1, basis, basis_state, counter)
	println("Basis states generated")


	# STEP 2: construct the Hamiltonian
	H = setup_ham(basis, impsiteL, impsiteR, eps, U, t, V)

	println("Hamiltonian matrix generated")


	# STEP 3: diagonalize the Hamiltonian
	eigens = eigen(H)
	println("ground state energy: ", eigens.values[1])


	# STEP 4: calculate the impurity occupation at T=0
	# CRITICAL: rescale eigenvalues wrt to the ground state
	eigens.values .-= eigens.values[1]

	# partition function
	Z = sum(exp.(-eigens.values/T))

	# impurity occupation (nd) in occupation basis matrix
	nd_bs = BitArray(zeros(Bool, size(H,1), size(H,1))) 
	for i in 1:size(nd_bs,1)
		if basis[impsiteR,i] # if up spin present on right impurity
			nd_bs[i,i] = true
		end
	end

	# matrix multiplication method
	nd_eig = eigens.vectors' * nd_bs * eigens.vectors
	nd = (1/Z) * tr(nd_eig .* exp.(-eigens.values/T))
	println("impurity occupation: ", 2*nd) # = nd(down) + nd(up), as no magnetic field


	# STEP 5: compute the matrix elements for Lehmann
	# for impurity = d:
	cdag_bs = BitArray(zeros(Bool, size(H,1), size(H,1))) # creation(d) occupation basis

	for j in 1:size(H,1)
		for k in 1:size(H,1)

			if cdag_condition(basis[:,j], basis[:,k], impsiteR)
				cdag_bs[j,k] = true
			end
		end
	end
	println("creation basis elements computed")
	

	# convert to eigens basis
	cdag_eig = eigens.vectors' * cdag_bs * eigens.vectors
	println("cdag eigen elements computed")
	

	# STEP 6: using Lehmann's representation, compute the spectral function
	# frequency grid
	w::Vector{ComplexF64} = collect(-wmax:dw:wmax) .+ nu*im
	# retarded GF
	G = zeros(ComplexF64, length(w))
	G .= (1/Z) * lehmann_finiteT(w, eigens.values, cdag_eig, T)

	# spectral function
	A = zeros(Float64, length(w))
	A .= -imag(G) ./ pi 

	println("normalization of spectral function: ", sum(A)*dw)

	# saving spectral function to file (optional) as w (column 1) and A (column 2)
	open("spectral_function.dat", "w") do file
		[println(file, "$(real(w[i]))\t$(A[i])") for i in 1:length(A)]
	end


	# plot the spectral function (optional)
	plot(real(w), A, dpi=300, legend=false)
	xlabel!("w")
	ylabel!("A(w)")
	title!("(Nf, eps, U, V, T) = ($Nf, $eps, $U, $V, $T)")
	savefig("spectral_function.png")

end


		### MAIN ###
let

	# model params
	Nf::Int = 3 # bath sites (folded representation)
	T::Float64 = 1e-4 # temperature

	eps::Float64 = 0.0 # impurity potential
	U::Float64 = 1.0 # interaction
	V::Float64 = 0.2 # hybridization
	t::Float64 = 0.5 # bath hopping

	# spectral function parameters
	wmax::Float64 = 1.5 # max frequency to compute
	dw::Float64 = 1e-3 # frequency step
	nu::Float64 = 1e-2 # imaginary linear broadening


	# main compute function
	do_ed(Nf, eps, U, V, t, wmax, dw, nu, T)

end # let
