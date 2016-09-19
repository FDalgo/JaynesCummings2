module JaynesCummings

global const h_bar = 1.0#1.05457173E-34 # in Joules*seconds		

export a, a_dagger

include("gen_initialstate.jl")
include("rabi_hamiltonian.jl")
include("gen_timeevoarray.jl")
include("calc_fullevo.jl")
include("partialtrace.jl")
include("print_matrix.jl")
include("get_eig_hamiltonian.jl")
include("get_prob.jl")

function a(cutoffN) # Truncated matrix for the lowering operator
  out = zeros(Complex128,cutoffN,cutoffN);
  for i = 1:cutoffN-1
      out[i,i+1] = sqrt(i);
  end
  return out
end

function a_dagger(cutoffN) # Truncated matrix for the raising operator
  out = zeros(Complex128,cutoffN,cutoffN)
  for i = 1:cutoffN-1
    out[i+1,i] = sqrt(i);
  end
  return out
end

end
