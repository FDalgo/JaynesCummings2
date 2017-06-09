function gen_hamiltonian(omegaQ,omegaR,g,cutoffN, strong_interaction, Braak)
  # This function generates a times invariant matrix for the Jaynes-Cummings
  # Hamiltonian. The parameters omegaR, omegaQ and g should be given in
  # rad/s. The cutoffN parameter is the highest energy level of the harmonic
  # oscillator that is considered. Higher levels are truncated.

  # Some constant parameters for the Jaynes-Cummings
  # Some constant parameters for the Jaynes-Cummings
  identity = [1 0 ; 0 1]; sigmaZ = 0.5*[1 0 ; 0 -1];
  sigmaMinus = [0 1 ; 0 0]; sigmaPlus = [0 0 ; 1 0];
  global h_bar
	#info(a_dagger(cutoffN) * a(cutoffN))
	half_matrix = Braak ? zeros(cutoffN, cutoffN) : eye(cutoffN)*0.5
  # The Jaynes-Cummings Hamiltonian is written with the following
  # conventions: omegaR is the resonator frequency, omegaQ is the qubit
  # frequency and g is the coupling.
  jaynesCummings = zeros(Complex128,2*cutoffN,2*cutoffN)
  jaynesCummings = jaynesCummings +
     kron(h_bar .* omegaR .* ((a_dagger(cutoffN) * a(cutoffN))+half_matrix), identity) -	# harmonic oscillator term
     h_bar .*omegaQ .* kron(eye(cutoffN), sigmaZ) +						# qubit term
     h_bar .* g .* ( kron(a(cutoffN), sigmaPlus) + kron(a_dagger(cutoffN), sigmaMinus) )	# interaction term
	if Braak
		jaynesCummings-=h_bar .*omegaQ .* kron(eye(cutoffN), sigmaZ)
	end
  if strong_interaction
  	jaynesCummings+= h_bar .* g .* (kron(a(cutoffN), sigmaMinus) + kron(a_dagger(cutoffN), sigmaPlus))# strong interaction term
 	end
  return jaynesCummings
end

#This function returns the desired B block
function get_B_block(hamiltonian, index)
	low = (index+1)*2
	return sub(hamiltonian, low:low+1, low:low+1)
end
