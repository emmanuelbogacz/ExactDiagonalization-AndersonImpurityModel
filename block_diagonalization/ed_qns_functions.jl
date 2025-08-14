# NOTE: functions in this file are called in both finite and zero temperature codes

# Function to insert new quantum numbers in a dictionary
# generates a single quantum number from a given pair of quantum numbers
function newqn(q1::Int, q2::Int, qn_store::Dict{Tuple{Int, Int}, Int})
    qnum = (q1, q2)
    pos = get(qn_store, qnum, length(qn_store)+1)
    qn_store[qnum] = pos
    return pos
end


# Function to lookup referenced quantum numbers
function qn(q1::Int, q2::Int, qn_store::Dict{Tuple{Int, Int}, Int})
    qnum = (q1, q2)
    return get(qn_store, qnum, 0)  # Return 0 if qnum is not found
end


# generate the basis states recursively
function gen_basis(level::Int, basis::BitArray, basis_state::BitArray,
		qn_store::Dict{Tuple{Int, Int}, Int}, size_bl::Vector{Int})

		# outer loop is the base (2)
        for b in [true, false] # electron present or absent
			basis_state[level] = b

			# inner loop is the exponent (1,2,...,2(N+1))
			if level == size(basis,1) # construct the basis state
				impsiteL = div(size(basis,1),2)

				# charge at half-filling
				q::Int = sum(basis_state) - impsiteL
				# 2*Sz: spin up left, spin down right
				sz::Int = sum(basis_state[1:impsiteL]) -
						sum(basis_state[impsiteL+1:end])

				qnind::Int = qn(q, sz, qn_store) # quantum number block
				if qnind == 0 # if qn not referenced
					qnind = newqn(q, sz, qn_store)
				end
				size_bl[qnind] += 1 # size of block with label qnind

				for i in 1:size(basis,1)
					basis[i, size_bl[qnind], qnind] = basis_state[i]
				end

			else
				gen_basis(level+1, basis, basis_state, qn_store, size_bl)

			end # if
        end # b
end # function


# Function to add an eigensystem to the dictionary
function add_eigensystem!(eigensystems::Dict{Symbol, Tuple{Vector{Float64},
			Matrix{Float64}}}, index::Int, matrix::Matrix{Float64})
    sym = Symbol(string(index))
    eigenvalues, eigenvectors = eigen(matrix)
    eigensystems[sym] = (eigenvalues, eigenvectors)
end


# Function to get an eigensystem from the dictionary
function get_eigensystem(eigensystems::Dict{Symbol, Tuple{Vector{Float64},
			Matrix{Float64}}}, index::Int)
    sym = Symbol(string(index))
    return get(eigensystems, sym, nothing)  # Return `nothing` if the eigensystem is not found
end


# remove the ground state energy from all eigenvalues (for numerical stability)
function shift_eigenvalues!(eigen_set::Dict{Symbol, Tuple{Vector{Float64},
			Matrix{Float64}}}, gs_en::Float64)
    for (key, value) in eigen_set
		value[1] .-= gs_en # shift all eigenvalues by gs
    end
end


# setup Hamiltonian matrix
function setup_ham(basis::BitArray, impsiteL::Int, impsiteR::Int, eps::Float64,
				   U::Float64, tun::Array, size_bl::Array, qnind::Int,
				   eigen_set::Dict)

		# construct the Hamiltonian for this block
		H = zeros(Float64, size_bl[qnind], size_bl[qnind])

		for i in 1:size(H,1)
			# potential, interaction
			if basis[impsiteL,i,qnind] || basis[impsiteR,i,qnind]
				H[i,i] += eps * (basis[impsiteL,i,qnind] +
						basis[impsiteR,i,qnind]) +
						U * basis[impsiteL,i,qnind] *
						basis[impsiteR,i,qnind]
			end
		end # i

		# tunneling elements
		for i in 1:size(H,1) # bra basis state
		for j in 1:size(H,2) # ket basis state
			if i == j
				continue # diagonal elements already done
			end

			# filling convention: |123> = cdag1 cdag2 cdag3 |vac>
			for k in 1:length(tun)
				if tunneling_condition(basis[:,i,qnind],
						       basis[:,j,qnind], k)
					H[i,j] += tun[k]
				end
			end # k
		end # j
		end # i

        # diagonalize the Hamiltonian
		add_eigensystem!(eigen_set, qnind, H)
		return eigen_set, H
end # func


# tunneling terms in the Hamiltonian, check if the tunneling condition is satisfied
function tunneling_condition(arr1::BitArray, arr2::BitArray, k::Int)
	# Check if the bits at k and k+1 are different, in both arrays
	if arr1[k] == arr1[k+1] || arr2[k] == arr2[k+1]
		return false
	end

	# Check if the bits at k in 1 arrays is the same at k+1 in the other array
	if arr1[k] == arr2[k] || arr1[k+1] == arr2[k+1]
		return false
	end

	# Create a mask to ignore elements at positions k and k+1
	mask = trues(length(arr1))
	mask[k] = false
	mask[k+1] = false

	# Compare the arrays using the mask
	return all(arr1[mask] .== arr2[mask])
end


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


# creation operator Hamiltonian term condition
function cdag_condition(arr1::BitArray, arr2::BitArray, k::Int) # k=impsite
	# if electron is absent in bra or present in ket
	if  !arr1[k] || arr2[k]
		return false
	end

	mask = trues(length(arr1))
	mask[k] = false

	# if the rest of the states are the same then true
	return all(arr1[mask] .== arr2[mask])
end


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
