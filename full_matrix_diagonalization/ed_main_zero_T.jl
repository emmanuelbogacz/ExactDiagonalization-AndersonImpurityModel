# NOTE: parameters can be edited in the main, at the bottom of this file

using LinearAlgebra # for matrix operations
using Plots # plotting of the spectral function

include("ed_functions.jl")


# annihilation operator Hamiltonian term condition
function c_condition(arr1::BitArray, arr2::BitArray, k::Int) # k=impsite
        # if electron is present in bra or absent in ket
        if  arr1[k] || !arr2[k]
                return false
        end

        mask = trues(length(arr1))
        mask[k] = false

        # if the rest of the states are the same then true
        return all(arr1[mask] .== arr2[mask])
end


# partial (lesser or greater) spectral function, T=0
function lehmann_T0(w::Vector{ComplexF64}, en::Vector{Float64}, elem::Vector{Float64})
	# returns Gle or gr, without +-i factor
	G = zeros(ComplexF64, length(w))

	for i in 1:length(w)
		for j in 1:length(en)
			G[i] += elem[j] ^2 / (w[i] + en[j])
		end
	end
	return G
end


# main function to perform the full ED calculation
function do_ed(Nf::Int, eps::Float64, U::Float64, V::Float64, t::Float64, 
		wmax::Float64, dw::Float64, nu::Float64)

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
	# CRITICAL for Z: rescale eigenvalues wrt to the ground state, avoid numerical issues
	eigens.values .-= eigens.values[1]

	nd::Float64 = 0.0 # <gs|n_d|gs>
	for i in 1:size(H,1) 
		if basis[impsiteR,i] 
			nd += (eigens.vectors[i,1]) .^2
		end
	end
	println("impurity occupation: ", 2*nd) # = nd(down) + nd(up), as no magnetic field


	# STEP 5: compute the matrix elements for Lehmann
	# for impurity = d:
	c_bs = BitArray(zeros(Bool, size(H,1), size(H,1))) # annihilation(d) occupation basis
	cdag_bs = BitArray(zeros(Bool, size(H,1), size(H,1))) # creation(d) occupation basis

	for j in 1:size(H,1)
		for k in 1:size(H,1)

			if c_condition(basis[:,j], basis[:,k], impsiteR)
				c_bs[j,k] = true
			end

			if cdag_condition(basis[:,j], basis[:,k], impsiteR)
				cdag_bs[j,k] = true
			end
		end
	end
	println("creation/annihilation basis elements computed")
	
	
	# convert to eigens basis, only for <n|cdag|gs> and <n|cd|gs>, at T=0
	c_eig = eigens.vectors' * c_bs * eigens.vectors[:,1]
	cdag_eig = eigens.vectors' * cdag_bs * eigens.vectors[:,1]
	println("creation/annihilation eigen elements computed")
	

	# STEP 6: using Lehmann's representation, compute the spectral function
	# frequency grid
	w::Vector{ComplexF64} = collect(-wmax:dw:wmax) .+ nu*im
	# greater GF
	Ggr = lehmann_T0(w, -(eigens.values.-eigens.values[1]), cdag_eig)
	# lesser GF
	Gle = -lehmann_T0(w, eigens.values.-eigens.values[1], c_eig)

	# retarded GF:
	G = similar(Ggr)
	G .= Ggr .- Gle

	# spectral function
	A = zeros(Float64, length(G))
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
	title!("(Nf, eps, U, V) = ($Nf, $eps, $U, $V)")
	savefig("spectral_function.png")

end


		### MAIN ###
let

	# model params
	Nf::Int = 3 # bath sites (folded representation)

	eps::Float64 = 0.0 # impurity potential
	U::Float64 = 1.0 # interaction
	V::Float64 = 0.2 # hybridization
	t::Float64 = 0.5 # bath hopping

	# spectral function parameters
	wmax::Float64 = 1.5 # max frequency to compute
	dw::Float64 = 1e-3 # frequency step
	nu::Float64 = 1e-2 # imaginary linear broadening


	# main compute function
	do_ed(Nf, eps, U, V, t, wmax, dw, nu)

end
