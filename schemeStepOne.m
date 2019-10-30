% Breaking down the parts of the scheme for ease of use in the Multiscale
% case. 

function [v, ft, fx] = schemeStepOne(oldPath, config)

% size is NOT taken from the config because the number of timesteps can be
% increased throughout the numerical scheme. 
n = size(oldPath.f, 1);
m = size(oldPath.f, 2) - 1;
DeltaX = 1 / n;
DeltaT = 1 / m;

% Dx is taken as a one sided derivative since it does not affect the second
% derivative for v. 
% May shift to center of grid for computation of fx to preserve symmetry in
% the future. 
Dx = circshift( eye(n), -1) - eye(n);
Dx = Dx / DeltaX ;

% Dt needs to be a m + 1 by m matrix. It is on the right when multiplied
% with f, and f is a n by m + 1 matrix. 
% Then manually set the non-zero values of Dt since 
Dt = zeros(m + 1 , m);
for i = 1:m
    Dt(i, i) = -1;
    Dt(i + 1, i) = 1;
end
Dt = Dt / DeltaT;

% Step 1: Finding the optimal v from the old f. 
% - v_xx + (1 + f_x^2) v = - f_x f_t 
% need to do this separately for each time step, since the f_x^2 term 
% is different for every time slice. 

% fx needs to be computed on the grid for v, which is offset in time. 
fxtemp = Dx * oldPath.f; % values on the grid.
fx = 1 / 2 * ( fxtemp(:, 1:m) + fxtemp(:, 2:m+1) ); 
ft = oldPath.f * Dt;
%disp(size(fx))
%disp(size(ft))

v = zeros(n, m);
for j = 1:m
    D2 = config.lambda * Dx * (- Dx');
    D4 = config.epsilon * D2 * D2;
    A = - D2 - D4 + eye(n) + eye(n) .* ( fx(:, j) .* fx(:, j) ) ;
    b = - fx(:, j) .* ft(:, j);
    v(:, j) = A \ b;
    %[result, ~, ~, ~] = pcg(A, b);
    %v(:, j) = result;
end

end