% A function that computes the action of a path. 

% INPUTS:
% Path - a stuct with fields f, v, z, representing a discretized admissibe 
% path. 

% OUTPUTS:
% Action - a number, the action associated with the path inputted. 

function Action = ComputeAction(Path)

% compute vx, and take L2 norms. 

% Use balanced derivatives. 
% From some testing, circshift seems to have better performance than 
% matrix multiplication. 

vdiff = circshift(Path.v, -1) - circshift(Path.v, 1);

vx = vdiff * (size(Path.v, 1) - 1) / 2;

Action = L2Squared(Path.v) + L2Squared(Path.z) + L2Squared(vx);
%disp(L2(Path.v))
%disp(L2(Path.z))
%disp(L2(vx))

end

% a very simple way to compute the square of the l2 norm of a function. 
function y = L2Squared(f)
y = mean( f .* f, 'all');
end