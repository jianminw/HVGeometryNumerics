function z = ComputeZFromFV(path)

% start by computing f_x and f_t. 
n = size(path.f, 1);
m = size(path.f, 2);
DeltaX = 1 / n;
DeltaT = 1 / (m-1);

%Dx = circshift( eye(n), -1 ) - circshift( eye(n), 1) ;
%Dx = Dx / 2;
%Switching to one sided derivatives
Dx = circshift( eye(n), -1) - eye(n);
% Dx = - Dx';
Dx = Dx / DeltaX ;

Dtplus = circshift( eye(m), -1 ) - eye (m) ;
% make sure the wrap around in time doesn't happen. 
Dtplus(m, 1) = 0;
Dtplus(m, m-1) = -1;
Dtplus(m, m) = 1; 
Dtplus = Dtplus / DeltaT;

fx = Dx * path.f;
ft = path.f * Dtplus';

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