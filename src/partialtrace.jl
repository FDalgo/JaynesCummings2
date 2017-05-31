function partialtrace(rho,dimOfSubsystemToKeep)
  # First check if the density matrix is Hermitian
  if maximum(abs(rho'-rho))>8E-15
    error("This density matrix is not Hermitian!")
  end

  dimOfSubsystemToTraceOut = div(size(rho)[1] , dimOfSubsystemToKeep);

  matrix = zeros(Complex128,dimOfSubsystemToKeep,dimOfSubsystemToKeep);

  for i=1:dimOfSubsystemToTraceOut
    vector = zeros(Complex128,dimOfSubsystemToTraceOut,1);
    vector[i] = 1;
    fullVector = kron(eye(dimOfSubsystemToKeep),vector);
    matrix = matrix + fullVector' * rho * fullVector;
  end

  return matrix
end
