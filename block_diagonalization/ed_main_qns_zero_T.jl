# NOTE: parameters can be edited in the main, at the bottom of this file

using LinearAlgebra # for matrix operations
using Plots # plotting of the spectral function

include("ed_qns_functions.jl")


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
		wmax::Float64, dw::Float64, nu::Float64, qmax::Int, szmax::Int, qnmax::Int)

	n::Int = 2Nf + 2 # total number of sites (unfolded rpz)
	impsiteL::Int = Nf + 1 # impurity site left
	impsiteR::Int = Nf + 2 # impurity site right


	# Define a dictionary to qn_store the quantum numbers and their positions
	qn_store = Dict{Tuple{Int, Int}, Int}()
	size_bl = zeros(Int, qnmax) # degeneracy in each block

	basis = BitArray(zeros(Bool, n, length(size_bl), length(size_bl)))
	basis_state = BitArray(zeros(Bool, size(basis,1)))


	# STEP 1: generate the basis states 
	gen_basis(1, basis, basis_state, qn_store, size_bl)
	println("Basis states generated")
	

	# STEP 2: construct the Hamiltonian + diagonalize
	# all tunneling elements
	tun::Vector{Float64} = fill(t, size(basis,1)-1)
	tun[impsiteL-1] = V
	tun[impsiteR] = V
	tun[impsiteL] = 0.0 # here interaction

	# construct the Hamiltonian by block, Define a dictionary to store matrices
	H_set = Dict{Symbol, Matrix{Float64}}()
	eigen_set = Dict{Symbol, Tuple{Vector{Float64}, Matrix{Float64}}}()

	for q in -qmax:qmax
	for sz in -szmax:szmax
		qnind::Int = qn(q, sz, qn_store)
		if qnind == 0 # no qn number associated
			continue
		end

		eigen_set, H = setup_ham(basis, impsiteL, impsiteR, eps, U, tun, size_bl, qnind,
				eigen_set)

	end # sz
	end # q

	println("Hamiltonian  set up and diagonalized")


	# STEP 4: find the ground state
	gs_en::Float64 = Inf # large ansatz
	qnind_gs::Int = 0
	q_gs::Int = 0
	sz_gs::Int = 0
	
	for q in -qmax:qmax
	for sz in -szmax:szmax
		qnind::Int = qn(q, sz, qn_store)
		if qnind == 0 # no qn number associated
			continue
		end

		# get the eigensystem for this block
		eigenvalues = get_eigensystem(eigen_set, qnind)[1]
		if eigenvalues[1] < gs_en
			gs_en = eigenvalues[1]
			qnind_gs = qnind
			q_gs = q
			sz_gs = sz
		end
	end # sz
	end # q

	println("ground state energy: $gs_en")
	println("Quantum numbers of the ground state: QNind=$qnind_gs, Q=$q_gs, S_z=$sz_gs")


	# shift all eigenvalues by the ground state energy
	shift_eigenvalues!(eigen_set, gs_en)
	gs_en = get_eigensystem(eigen_set, qnind_gs)[1][1] # =0 by construction
	eigvec_gs = get_eigensystem(eigen_set, qnind_gs)[2]


	# STEP 5: compute the occupation number
	nd::Float64 = 0.0
	for i in 1:size(eigvec_gs,2) # loop over all basis states
		if basis[impsiteL,i,qnind_gs]
			nd += (eigvec_gs[i,1]) .^2
		end
	end
	println("total occupation number (spin up + down): $(2*nd)")

	
	# STEP 6: compute the GF using Lehmann
	w::Vector{ComplexF64} = collect(-wmax:dw:wmax) .+ nu*im

	qnind_cdag::Int = qn(q_gs+1, sz_gs+1, qn_store) # cdag bra
	qnind_c::Int = qn(q_gs-1, sz_gs-1, qn_store) # c bra
	
	# creation/annihilation operators in occupation basis
	cdag_bs = BitArray(zeros(Bool, size_bl[qnind_cdag], size_bl[qnind_gs])) 
	c_bs = BitArray(zeros(Bool, size_bl[qnind_c], size_bl[qnind_gs]))
	

	# find matrix elements for GF:
	for j in 1:size_bl[qnind_cdag] # bra
	for k in 1:size_bl[qnind_gs] # ket

		if cdag_condition(basis[:,j,qnind_cdag], basis[:,k,qnind_gs], impsiteL)
			cdag_bs[j,k] = true
		end
	end
	end
	
	for j in 1:size_bl[qnind_c] # bra
	for k in 1:size_bl[qnind_gs] # ket

		if c_condition(basis[:,j,qnind_c], basis[:,k,qnind_gs], impsiteL)
			c_bs[j,k] = true
		end
	end
	end


	# diagonalize them: <n|cdag|gs> and <n|cd|gs>
	eigval_cdag, eigvec_cdag = get_eigensystem(eigen_set, qnind_cdag)
	cdag_eig = eigvec_cdag' * cdag_bs * eigvec_gs[:,1]

	eigval_c, eigvec_c = get_eigensystem(eigen_set, qnind_c)
	c_eig = eigvec_c' * c_bs * eigvec_gs[:,1]

	# greater GF
	Ggr = lehmann_T0(w, -(eigval_cdag.-gs_en), cdag_eig)
	# lesser GF
	Gle = -lehmann_T0(w, eigval_c.-gs_en, c_eig)

	# retarded GF
	G = similar(Ggr)
	G .= Ggr .- Gle

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
	title!("(Nf, eps, U, V) = ($Nf, $eps, $U, $V)")
	savefig("spectral_function.png")

end # do_ed


		### MAIN ###
let

	# model parameters
	Nf::Int = 3 # bath sites (folded representation)

	eps::Float64 = 0.0 # impurity potential
	U::Float64 = 1.0 # interaction
	V::Float64 = 0.2 # hybridization
	t::Float64 = 0.5 # bath hopping

	# spectral function parameters
	wmax::Float64 = 1.5 # max frequency to compute
	dw::Float64 = 1e-3 # frequency step
	nu::Float64 = 1e-2 # imaginary linear broadening

	# quantum numbers config (need to increase size if Nf is large)
	qmax::Int = 15 # max charge at half filling
	szmax::Int = 15 # max Sz
	qnmax::Int = 6000 # max number of qn


	# main compute function
	do_ed(Nf, eps, U, V, t, wmax, dw, nu, qmax, szmax, qnmax)

end
