function get_prob(matrix,str)
	# Grabbing resonator state
  level = eval(parse(str[search(str,',')+1:end-1]))
  # Grabbig qubit state
  if str[search(str,'|')+1] == 'g'
    return(matrix[:, level*2+1, level*2+1])
  else
   return(matrix[:, level*2+2, level*2+2])
  end
end
