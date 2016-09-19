function print_matrix(M, precision)
  # This function prints a 2D matrix
  a, b = size(M)
  println("Hamiltonian:")
  for i=1:b
  	for j=1:a
  		print(round(M[j, i], precision), "\t")
  	end
  	println()
  end
end
