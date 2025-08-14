# NOTE: parameters can be edited in the main, at the bottom of this file

using LinearAlgebra # for matrix operations
using Plots # plotting of the spectral function

include("ed_qns_functions.jl")


# partial (greater or lesser) spectral function, finite T
function Gpartial(w::Vector{ComplexF64}, en::Vector{Float64}, 
		enp::Vector{Float64}, elem::Matrix{Float64}, T::Float64, sign::Int)
	# returns Gle or gr, without pre-factor
	G = zeros(ComplexF64, length(w))

	for i in 1:length(w)
		for n in 1:length(en)
		for np in 1:length(enp)
			G[i] += exp(-en[n]/T)* elem[n,np]^2 /
				(w[i] + sign*(enp[np] -en[n]))
		end
		end
	end
	return G
end


# main function to perform the full ED calculation
function do_ed(Nf::Int, eps::Float64, U::Float64, V::Float64, t::Float64, 
		wmax::Float64, dw::Float64, nu::Float64, T::Float64, 
		qmax::Int, szmax::Int, qnmax::Int)

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

	println("Hamiltonian setup and diagonalized")


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


	# CRUCIAL for Z: shift all eigenvalues by the ground state energy
	shift_eigenvalues!(eigen_set, gs_en)

	Z::Float64 = 0.0 # partition function
	for q in -qmax:qmax
	for sz in -szmax:szmax
		qnind::Int = qn(q, sz, qn_store)
		if qnind == 0 # no qn number associated
			continue
		end

		# get the eigensystem for this block
		eigenvalues = get_eigensystem(eigen_set, qnind)[1]
		Z += sum(exp.(-eigenvalues/T))
	end
	end
	println("Z: $Z")


	# STEP 5: compute the impurity occupation number
	nd::Float64 = 0.0

	for q in -qmax:qmax
	for sz in -szmax:szmax
		qnind::Int = qn(q, sz, qn_store)
		if qnind == 0 # no qn number associated
			continue
		end

		# nd in occ basis matrix
		nd_bs = BitArray(zeros(Bool, size_bl[qnind], size_bl[qnind]))
		for i in 1:size(nd_bs,1)
			if basis[impsiteL,i,qnind]
				nd_bs[i,i] = true
			end
		end

		# diagonalize
		eigval, eigvec = get_eigensystem(eigen_set, qnind)
		nd_eig = eigvec' * nd_bs * eigvec
		nd += (1/Z) * tr(nd_eig .* exp.(-eigval/T))

	end # sz
	end # q

	println("total occupation number (spin up + down): $(2*nd)")

	
	# STEP 6: compute the GF using Lehmann
	w::Vector{ComplexF64} = collect(-wmax:dw:wmax) .+ nu*im
	Gle = zeros(ComplexF64, length(w))
	Ggr = zeros(ComplexF64, length(w))


	# CDAG (creation matrix element)
	for q in -qmax:qmax
	for sz in -szmax:szmax
		qnind::Int = qn(q, sz, qn_store)
		qnind_cdag::Int = qn(q+1, sz+1, qn_store) # cdag bra

		if qnind == 0 || qnind_cdag == 0
			continue # no qn number associated
		end
		
		eigval, eigvec = get_eigensystem(eigen_set, qnind)

		# CDAG
		# c_dag in occupation basis
		cdag_bs = BitArray(zeros(Bool, size_bl[qnind_cdag], size_bl[qnind])) 
		
		for j in 1:size_bl[qnind_cdag] # bra
		for k in 1:size_bl[qnind] # ket

			if cdag_condition(basis[:,j,qnind_cdag], basis[:,k,qnind], 
								  impsiteL)
				cdag_bs[j,k] = true
			end
		end
		end
		
		# diagonalize 
		eigval_cdag, eigvec_cdag = get_eigensystem(eigen_set, qnind_cdag)
		cdag_eig = eigvec_cdag' * cdag_bs * eigvec

		Gle += (-1/Z)*Gpartial(w, eigval_cdag, eigval, cdag_eig, T, 1)
	end # sz
	end # q


	# C (annihilation matrix element)
	for q in -qmax:qmax
	for sz in -szmax:szmax
		qnind::Int = qn(q, sz, qn_store)
		qnind_c::Int = qn(q-1, sz-1, qn_store) # c bra

		if qnind == 0 || qnind_c == 0
			continue # no qn number associated
		end

		eigval, eigvec = get_eigensystem(eigen_set, qnind)

		# c in occupation basis
		c_bs = BitArray(zeros(Bool, size_bl[qnind_c], size_bl[qnind]))

		for j in 1:size_bl[qnind_c] # bra
		for k in 1:size_bl[qnind] # ket

			if c_condition(basis[:,j,qnind_c], basis[:,k,qnind], 
								  impsiteL)
				c_bs[j,k] = true
			end
		end
		end

		# diagonalize
		eigval_c, eigvec_c = get_eigensystem(eigen_set, qnind_c)
		c_eig = eigvec_c' * c_bs * eigvec

		Ggr += (1/Z)*Gpartial(w, eigval_c, eigval, c_eig, T, -1)
	end # sz
	end # q


	# retarded GF
	G = similar(Gle)
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
	T::Float64 = 1e-4 # temperature

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
	do_ed(Nf, eps, U, V, t, wmax, dw, nu, T, qmax, szmax, qnmax)

end # let
