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
% This is most likely a negligible amout of time cost, and the clarity of
% the finite difference is most likely much more valuable. 

%vdiff = circshift(Path.v, -1) - circshift(Path.v, 1);

% Changing to one directional finite differences because of suspected
% parasitic modes

vdiff = circshift(Path.v, -1) - Path.v;
%vdiff = Path.v - circshift(Path.v, 1);

vx = vdiff * size(Path.v, 1);% / 2;

newZ = ComputeZFromFV(Path);

%Action = L2Squared(Path.v) + L2Squared(Path.z) + L2Squared(vx);
Action = L2Squared(Path.v) + L2Squared(newZ) + L2Squared(vx);
%disp(L2(Path.v))
%disp(L2(Path.z))
%disp(L2(vx))

end

% a very simple way to compute the square of the l2 norm of a function. 
% here assumes 
function y = L2Squared(f)
f2 = f .* f;
y = (mean(mean(f2))) / 2;
end