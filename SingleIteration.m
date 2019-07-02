% This file contains a function that execute a single iteration of the 
% numerical scheme. 
% This file should be split into a number of smaller files/functions for
% ease of reading/debugging. 

% INPUTS: f 
% OUTPUTS: f, v, z

function [f, v, z] = SingleIteration(oldPath, f0, f1)

n = size(oldPath.f, 1);
m = size(oldPath.f, 2);
DeltaX = 1 / (n - 1);
DeltaT = 1 / (m - 1);

Dx = circshift( eye(n), -1 ) - circshift( eye(n), 1) ;
Dx = Dx / ( 2 * DeltaX );

Dtplus = circshift( eye(m), -1 ) - eye (m) ;
% make sure the wrap around in time doesn't happen. 
Dtplus(m, 1) = 0;
Dtplus(m, m-1) = -1;
Dtplus(m, m) = 1; 
Dtplus = Dtplus / DeltaT;

% Step 1: Finding the optimal v from the old f. 
% - v_xx + (1 + f_x^2) v = - f_x f_t 
% need to do this separately for each time step, since the f_x^2 term 
% is different for every time slice. 

fx = Dx * oldPath.f;
ft = oldPath.f * Dtplus';

v = zeros(n, m);
for j = 1:m
    A = - Dx^2 + eye(n) + eye(n) .* ( fx(:, j) .* fx(:, j) ) ;
    b = - fx(:, j) .* ft(:, j);
    v(:, j) = A \ b;
    %[result, ~, ~, ~] = pcg(A, b);
    %v(:, j) = result;
end

% Step 2: Computing rho, which follows the continuity equation with v, 
% and has initial condition 1. 

rho = ones(n-1, m);

for j = 2:m
    rhoPrev = rho(:, j-1);
    rhoLeft = circshift(rhoPrev, 1);
    rhoRight = circshift(rhoPrev, -1);
    flux = ((rhoLeft + rhoPrev) / 2) .* v(1:n-1, j-1) - ((rhoRight + rhoPrev) / 2 ) .* v(2:n, j-1);
    rho(:, j) = rhoPrev + DeltaT / DeltaX * flux;
end
%disp(rho)

% Step 3: Finding the flow map of v. 
% this is a mess. I think I should do this without the validity check
% first. 

validityCheck = false;
timeSteps = m-1;
phi = ones(n, timeSteps + 1);

while (~validityCheck)
    validityCheck = true;
    phi = ones(n, timeSteps + 1);
    phi(:, 1) = (0:n-1) / (n - 1);

    for j = 2:(timeSteps + 1)
        % need to mess with this to make it work for smaller time steps. 
        phi(:, j) = phi(:, j-1) + interp2( (0:m-1) / (m-1), phi(:, 1), v, (j-1)/timeSteps, phi(:, j-1) ) / timeSteps;
        phi(:, j) = phi(:, j) + (phi(:, j) < 0) - (phi(:, j) > 1);
        % check that none of the intervals shrink too much. 
        prevInterval = phi(:, j-1) - circshift( phi(:, j-1), 1);
        prevInterval = prevInterval + (prevInterval < 0);
        prevInterval = prevInterval(2:n);
        interval = phi(:, j) - circshift( phi(:, j-1), 1);
        interval = interval + (interval < 0);
        interval = interval(2:n);
        validityCheck = sum( (2 * interval) < prevInterval ) == 0;
        if (~validityCheck)
            timeSteps = timeSteps * 2;
            break;
        end
    end
end

% Step 4: Finding the scaling factors for z, and finish. 

% Probably start by interpolating rho. 
% Keep the time step in step 3, since we would expect vx to be about the
% same for each interation, and then the next iteration would need to
% refine the grid anyways. 

rho_Phi = zeros(n, timeSteps + 1);
rho_wrap =  cat(1, rho, rho, rho);
for j = 1:(timeSteps + 1)
    x = (1/2:1:n-3/2) / (n - 1);
    x = cat(2, x - ones(size(x)), x, x + ones(size(x)));
    y = (0:m-1) / (m-1);
    times = phi(:, j);
    temp = interp2( y, x, rho_wrap , (j-1)/timeSteps,  times ) ;
    rho_Phi(:, j) = temp; %interp2( y, x, rho_wrap , (j-1)/timeSteps,  phi(:, j) ) ;
    %rho_Phi(:, j) = interp2( (1/2:1:n-3/2) / (n - 1), (0:m-1) / (m-1), rho, phi(:, j-1), (j-1)/timeSteps ) ;
end
%disp(rho)
%disp(phi)
%disp(rho_Phi)
figure('Name', 'phi')
hold on
plot(phi(:, 1), phi(:, timeSteps + 1) - phi(:, 1) );
%plot(phi(:, 1), f0);
%plot(phi(:, timeSteps + 1), f1);

C = eye(n);
for i = 1:n
    s = sum( rho_Phi(i, :) ) - rho_Phi(i, 1) / 2 - rho_Phi(i, timeSteps + 1) / 2;
    %disp(s)
    s = s / timeSteps;
    %disp(s)
    C(i, i) = (interp1( (0:n-1) / (n-1), f1, phi(i, timeSteps + 1) ) - f0(i) ) / s;
end
z_Phi = C * rho_Phi;

f_Phi = zeros(n, timeSteps + 1);
f_Phi(:, 1) = f0;

for j = 2:(timeSteps + 1)
    f_Phi(:, j) = f_Phi(:, j-1) + (z_Phi(:, j-1) + z_Phi(:, j) ) / 2 / timeSteps ;
end

% interpolate back to the Eulerian coordinates. 
f = zeros(n, timeSteps + 1);
f(:, 1) = f0;
f(:, timeSteps + 1) = f1;
for j = 2:timeSteps
    %size(f_Phi(:,j))
    %size(phi(:,j))
    x = phi(:, j);
    y = f_Phi(:, j);
    %disp(x)
    %disp(y)
    %disp(x)
    f(:, j) = interpOnS1( x, y, (0:n-1)/(n-1)');
end

z = zeros(n, timeSteps + 1);
for j = 2:timeSteps
    z(:, j) = interpOnS1( phi(:, j), z_Phi(:, j), (0:n-1)/(n-1));
end

end


function Vq = interpOnS1(X, V, Xq)
    newX = cat(1, X - ones(size(X)), X, X + ones(size(X)));
    newV = cat(1, V, V, V);
    Vq = interp1( newX, newV, Xq );
end

