% A function that computes the action of a path. 

% INPUTS:
% Path - a stuct with fields f, v, z, representing a discretized admissibe 
% path. 

% OUTPUTS:
% Action - a number, the action associated with the path inputted. 

function Action = ComputeAction(Path, config)

% compute vx, and take L2 norms. 

% Use balanced derivatives. 
% From some testing, circshift seems to have better performance than 
% matrix multiplication. 
% This is most likely a negligible amout of time cost, and the clarity of
% the finite difference is most likely much more valuable. 

%vdiff = circshift(Path.v, -1) - circshift(Path.v, 1);

% Changing to one directional finite differences because of suspected
% parasitic modes

n = size(Path.f, 1);
DeltaX = 1 / n;

Dx = circshift( eye(n), -1) - eye(n);
Dx = Dx / DeltaX ;
Dxx = Dx * (-Dx');

vx = Dx * Path.v;
vxx = Dxx * Path.v;

newZ = ComputeZFromFV(Path);

%Action = L2Squared(Path.v) + L2Squared(Path.z) + L2Squared(vx);
Action = L2Squared(Path.v) + L2Squared(newZ);
Action = Action + config.lambda * L2Squared(vx);
Action = Action + config.epsilon * L2Squared(vxx);
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