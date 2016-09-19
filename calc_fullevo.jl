function calc_fullevo(initialstate,time_evo_array)
  state_prob = zeros(Complex128,length(time_evo_array), size(initialstate, 1),size(initialstate, 2))
  # Time evolve the system density matrix and "measure" the qubit at each time sample taken. This is done by
  # tracing out the resonator and keeping the excited state entry of the density matrix (ie, entry 2,2 in the matrix)
  for i=1:length(time_evo_array)
    state_prob[i, :, :] = partialtrace(time_evo_array[i] * initialstate * time_evo_array[i]',size(initialstate, 1))
  end
  return state_prob
end
