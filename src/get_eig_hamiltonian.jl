function print_eigenvector(vec, precision)
	s = size(vec, 1)
	for i=1:s
		print(round(vec[i], precision), "\t")
	end
	println()
end

function eigenvalues_ham(hamiltonian, precision, print, Braak_form)
	size_ham = size(hamiltonian)
	λ = eigfact!(view(hamiltonian, 1:size_ham[1]-2, 1:size_ham[2]-2))
	states = Array{String}(length(λ[:values]))
	#print_matrix(hamiltonian, 2)
	if print
		println("\n Eigenvalues:")
		for i=1:size(λ[:values], 1)
			@printf("#%d |     ", i)
			@printf("%.2f \t current state |", round(λ[:values][i], precision))
			str = ""
			if Braak_form
				if i%2 == 0
					str = str * "+ "* join(convert(Int, floor((i-1)/2))) * ">"
				else
					str = str * "- " * join(convert(Int, floor(i/2))) * ">"
				end
			else
				if i%2 == 0
					str = str * "e, "* convert(Int, floor((i-1)/2))*">"
				else
					str = str * "g, " * convert(Int, floor(i/2)) *">"
				end
			end
			println(str)
			states[i] = str
		end
		#println("Extra λ |  ", round(hamiltonian[size_ham[1], size_ham[2]], precision))

		println("\n Eigenvectors:")
		for i=1:size(λ[:vectors], 2)
			@printf("#%d |     ", i)
			print_eigenvector(λ[:vectors][:, i], precision)
		end
	end
	return λ#(λ, states)
end

function sort_eigs(eigs, unsorted_eigs, ε)
	incr_ratio = eigs[:, 2] - eigs[:, 1]
	new_incr = unsorted_eigs - eigs[:, 2]
	for i=1:size(new_incr, 1)
		if abs(new_incr[i]-incr_ratio[i]) > ε ### out acceptable bounds...
			#info("qui")
			for j=i+1:size(new_incr, 1)
				if abs(unsorted_eigs[j] - eigs[i,2]-incr_ratio[i]) <= ε ##each other very close
					temp = unsorted_eigs[j]
					unsorted_eigs[j] = unsorted_eigs[i]
					unsorted_eigs[i] = temp
					new_incr = unsorted_eigs - eigs[:, 2]
				end
			end
		end
	end
	return unsorted_eigs
end

function get_corresponding_eigenvalue(hamiltonian, λ, str)
	size_ham = size(hamiltonian)
	level = eval(parse(str[search(str,',')+1:end-1]))
  # Grabbig qubit state
  index = -1
  if str[search(str,'|')+1] == 'g'
    index = level*2+1
  else
   index = level*2+2
  end
	test_hamiltonian = deepcopy(hamiltonian)
	test_hamiltonian[index, index] = 0.0
	λ_test = eigfact!(sub(test_hamiltonian, 1:size_ham[1]-1, 1:size_ham[2]-1))
	for x=1: size(λ, 1)
		found = false
		for y =1: size(λ_test[:values], 1)
			if round(λ[x], 5)==round(λ_test[:values][y, 1], 5)
				found = true
			end
		end
		found == false && return λ[x]
	end
	return NaN
end
