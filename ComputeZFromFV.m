function z = ComputeZFromFV(path)

n = size(path.f, 1);
m = size(path.f, 2) - 1;
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
fxtemp = Dx * path.f; % values on the grid.
fx = 1 / 2 * ( fxtemp(:, 1:m) + fxtemp(:, 2:m+1) ); 
ft = path.f * Dt;

% now need to check for the size of v. This could be bad.
newV = path.v;

if size(path.v, 2) ~= m
    oldM = size(path.v, 2);
    Yq = (0:m-1) / (m-1) * (oldM - 1);
    Yq = Yq + 1;
    Yq = ones(n, 1) * Yq;
    Xq = 1:n;
    Xq = Xq';
    Xq = Xq * ones(1, m);
    %disp(size(path.v))
    %disp(size(Xq))
    %disp(size(Yq))
    newV = interp2( path.v, Xq, Yq);
end

z = ft + fx .* newV;
end