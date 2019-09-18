function error = CheckAdmissiblePath(path)
% a L2 sanity check on the admissible path condition.
% returns the L2 norm of f_x v + f_t - z. 

pathError = ComputeZFromFV(path) - path.z;

error = mean(mean(pathError .* pathError));
end