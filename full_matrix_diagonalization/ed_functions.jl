# NOTE: functions in this file are called in both finite and zero temperature codes

# generate the basis states recursively
function gen_basis(level::Int, basis::BitArray, basis_state::BitArray, counter::Ref{Int})
        # counter is passed by reference, removes the need of a global variable

	# outer loop is the base (2)
        for b in [true, false] # loop over all possible bits at level
                basis_state[level] = b

		# inner loop is the exponent (1,2,...,2(N+1))
                if level == size(basis,1)  # if the last level

                        counter[] += 1
                        for i in 1:size(basis,1)
				# construct the basis state
                                basis[i,counter[]] = basis_state[i]
                        end

                else
                        gen_basis(level+1, basis, basis_state, counter)
                end # if

        end # b
end # function


# setup Hamiltonian matrix
function setup_ham(basis::BitArray, impsiteL::Int, impsiteR::Int, eps::Float64,
				   U::Float64, t::Float64, V::Float64)

	H = zeros(Float64, size(basis,2), size(basis,2))

	# potential epsilon and interaction U
	for i in 1:size(H,1)
		# potential term
		if basis[impsiteL,i] || basis[impsiteR,i]
			H[i,i] += eps * (basis[impsiteL,i] + basis[impsiteR,i])
		end

		# interaction term
		if basis[impsiteL,i] && basis[impsiteR,i]
			H[i,i] += U
		end
	end # i


	# tunneling terms setup: V and t
	tun::Vector{Float64} = fill(t, size(basis,1)-1)
	tun[impsiteL-1] = V
	tun[impsiteR] = V
	tun[impsiteL] = 0.0 # no tunneling between left- and right-impurity (ficticious)


	# tunneling terms in the Hamiltonian
	for i in 1:size(H,1) # bra basis state
		for j in 1:size(H,2) # ket basis state
			if i == j
				continue # diagonal elements already done
			end

			# filling convention: |123> = cdag1 cdag2 cdag3 |vac>
			for k in 1:length(tun)
				if tunneling_condition(basis[:,i], basis[:,j], k)
					H[i,j] += tun[k] # ket
				end
			end # k
		end # j
	end # i

	return H
end


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


# retarded Green's function, finite T
function lehmann_finiteT(w::Vector{ComplexF64}, en::Vector{Float64},
		elem::Matrix{Float64}, T::Float64)
	# returns G_retarded, without 1/Z factor
	G = zeros(ComplexF64, length(w))

	for i in 1:length(w)
		for n in 1:length(en)
		for np in 1:length(en)
			G[i] += (elem[np,n])^2 *
			(exp(-en[n]/T) + exp(-en[np]/T)) / (w[i] +en[n] -en[np])
		end
		end
	end
	return G
end
