

function [f, phi] = schemeStepTwo(v, f0, f1, config)

% Step 2: Computing rho, which follows the continuity equation with v, 
% and has initial condition 1. 

n = size(v, 1);
m = size(v, 2);
DeltaX = 1 / n;
DeltaT = 1 / m;

rho = ones(n, m);

for j = 2:m
    rhoPrev = rho(:, j-1);
    rhoLeft = circshift(rhoPrev, 1);
    rhoRight = circshift(rhoPrev, -1);
    vLeft = circshift( v(:, j-1), 1 );
    vRight = v(:, j-1);
    flux = ((rhoLeft + rhoPrev) / 2) .* vLeft - ((rhoRight + rhoPrev) / 2 ) .* vRight;
    rho(:, j) = rhoPrev + DeltaT / DeltaX * flux;
end
%disp(rho)

% Step 3: Finding the flow map of v. 
% this is a mess. I think I should do this without the validity check
% first. 

phi = ODEScheme( v, config);
%disp(v)

% Step 4: Finding the scaling factors for z, and finish. 

% Probably start by interpolating rho. 
% Keep the time step in step 3, since we would expect vx to be about the
% same for each interation, and then the next iteration would need to
% refine the grid anyways. 

rho_Phi = zeros(n, m + 1);
for j = 1:(m + 1)
    x = (1/2:1:n-1/2) / n;
    t = (0:m-1) / (m-1);
    xquery = phi(:, j);
    temp = interpOnS1andTime( t, x', rho, (j-1)/m,  xquery ) ;
    rho_Phi(:, j) = temp; 
end
%figure('Name', 'phi')
%hold on
%plot(phi(:, 1), phi(:, timeSteps + 1) - phi(:, 1) );
%plot(phi(:, 1), f0);
%plot(phi(:, timeSteps + 1), f1);

C = eye(n);
for i = 1:n
    s = sum( rho_Phi(i, :) ) - rho_Phi(i, 1) / 2 - rho_Phi(i, m + 1) / 2;
    %disp(s)
    s = s / m;
    %disp(s)
    %disp(size(f1))
    %disp(size(phi(:, 1)))
    C(i, i) = (interpOnS1( phi(:, 1) , f1', phi(i, m + 1) ) - f0(i) ) / s;
end
z_Phi = C * rho_Phi;

f_Phi = zeros(n, m + 1);
f_Phi(:, 1) = f0;

for j = 2:(m + 1)
    f_Phi(:, j) = f_Phi(:, j-1) + (z_Phi(:, j-1) + z_Phi(:, j) ) / 2 / m ;
end

% interpolate back to the Eulerian coordinates. 
f = zeros(n, m + 1);
f(:, 1) = f0;
%disp(size(f_Phi(:, 1)))
%disp(size(f0))
%figure('Name', 'f_Phi')
%plot(phi(:, 1), f_Phi(:, 1) - f0')

f(:, m + 1) = f1;
for j = 2:m
    %size(f_Phi(:,j))
    %size(phi(:,j))
    x = phi(:, j);
    y = f_Phi(:, j);
    %disp(x)
    %disp(y)
    f(:, j) = interpOnS1( x, y, phi(:, 1));
end
%figure('Name', 'f')
%hold on
%plot( phi(:, timeSteps + 1), f_Phi(:, timeSteps + 1), 'r')
%plot( phi(:, 1), f0, 'g')
%plot( phi(:, 1), f1, 'b')

end