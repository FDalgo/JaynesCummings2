function gen_timeevoarray(hamiltonian,finalTime,samples)
  # Generate the time evolution array for a time independent Hamiltonian
  time_vec=((1:samples)-1)*finalTime/samples
  time_evo_array = Array(typeof(hamiltonian),samples)

  for i=1:samples 
    time_evo_array[i] = expm(-1im .* hamiltonian .* time_vec[i] ./ h_bar)
  end
	
  return time_vec, time_evo_array
end
