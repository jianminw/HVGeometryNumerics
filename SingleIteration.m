% This file contains a function that execute a single iteration of the 
% numerical scheme. 
% This file should be split into a number of smaller files/functions for
% ease of reading/debugging. 

% INPUTS: 
% oldPath - a stuct containing an admissible path. Only f is necessary. 
% f0 - the starting function. Should correspond to the 0 time slice in f
% provided in oldPath. 
% f1 - the ending function. Should correspond to the 1 time slice in f
% provided in oldPath. 

% OUTPUTS: 
% f, v, z - the new path obtained by the interative scheme. At one point
% there were some errors that seem to be caused by passing the outputs out
% as a struct, so they are passed out as a vector for now. 

function newPath = SingleIteration(oldPath, f0, f1, iteration)

config = options();

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
    A = - Dx * (- Dx') + eye(n) + eye(n) .* ( fx(:, j) .* fx(:, j) ) ;
    b = - fx(:, j) .* ft(:, j);
    v(:, j) = A \ b;
    %[result, ~, ~, ~] = pcg(A, b);
    %v(:, j) = result;
end

% If set in config, compute z for the above f and v to track the action
% within the iteration. 
% z = f_t + f_x v
% should be a simple computation. 

if config.computeActionMidIteration
    z = ft + fx .* v;
    midIterationPath = struct;
    midIterationPath.f = oldPath.f;
    midIterationPath.v = v;
    midIterationPath.z = z;
    newPath.midIterationAction = ComputeAction(midIterationPath);
end



% Step 2: Computing rho, which follows the continuity equation with v, 
% and has initial condition 1. 

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

validityCheck = false;
timeSteps = m;
phi = ones(n, timeSteps + 1);

while (~validityCheck)
    validityCheck = true;
    phi = ones(n, timeSteps + 1);
    phi(:, 1) = (1:n) / n;

    for j = 2:(timeSteps + 1)
        % need to mess with this to make it work for smaller time steps. 
        Dphi = interpOnS1andTime( (0:m-1) / (m-1), phi(:, 1), v, (j-1)/timeSteps, phi(:, j-1) );
        %Dphi = interpOnS1( phi(:, 1), v(:, j), phi(:, j-1) );
        phi(:, j) = phi(:, j-1) + Dphi / timeSteps;
        phi(:, j) = phi(:, j) + (phi(:, j) <= 0) - (phi(:, j) > 1);
        % check that none of the intervals shrink too much. 
        prevInterval = phi(:, j-1) - circshift( phi(:, j-1), 1);
        prevInterval = prevInterval + (prevInterval < 0);
        interval = phi(:, j) - circshift( phi(:, j-1), 1);
        interval = interval + (interval < 0);
        % check for both shrinking and growing intervals for symmetry. 
        scaleFactor = 2;
        shrinkingIntervals = sum( ( scaleFactor * interval) < prevInterval );
        growingIntervals = sum( interval > (scaleFactor * prevInterval) );
        if shrinkingIntervals > 0 || growingIntervals > 0
            timeSteps = timeSteps * 2;
            %disp(j)
            %disp(timeSteps)
            validityCheck = false;
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
for j = 1:(timeSteps + 1)
    x = (1/2:1:n-1/2) / n;
    t = (0:m-1) / (m-1);
    xquery = phi(:, j);
    temp = interpOnS1andTime( t, x', rho, (j-1)/timeSteps,  xquery ) ;
    rho_Phi(:, j) = temp; 
end
%figure('Name', 'phi')
%hold on
%plot(phi(:, 1), phi(:, timeSteps + 1) - phi(:, 1) );
%plot(phi(:, 1), f0);
%plot(phi(:, timeSteps + 1), f1);

C = eye(n);
for i = 1:n
    s = sum( rho_Phi(i, :) ) - rho_Phi(i, 1) / 2 - rho_Phi(i, timeSteps + 1) / 2;
    %disp(s)
    s = s / timeSteps;
    %disp(s)
    %disp(size(f1))
    %disp(size(phi(:, 1)))
    C(i, i) = (interpOnS1( phi(:, 1) , f1', phi(i, timeSteps + 1) ) - f0(i) ) / s;
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
%disp(size(f_Phi(:, 1)))
%disp(size(f0))
%figure('Name', 'f_Phi')
%plot(phi(:, 1), f_Phi(:, 1) - f0')

f(:, timeSteps + 1) = f1;
for j = 2:timeSteps
    %size(f_Phi(:,j))
    %size(phi(:,j))
    x = phi(:, j);
    y = f_Phi(:, j);
    %disp(x)
    %disp(y)
    %disp(x)
    f(:, j) = interpOnS1( x, y, phi(:, 1));
end
%figure('Name', 'f')
%hold on
%plot( phi(:, timeSteps + 1), f_Phi(:, timeSteps + 1), 'r')
%plot( phi(:, 1), f0, 'g')
%plot( phi(:, 1), f1, 'b')

z = zeros(n, timeSteps + 1);
for j = 1:(timeSteps + 1)
    z(:, j) = interpOnS1( phi(:, j), z_Phi(:, j), (0:n-1)/(n-1));
end

newPath.f = f;
newPath.v = v;
newPath.z = z;
newPath.phi = phi;


plotName = ['Iteration: ' int2str(iteration)];
plotting( newPath, config.plotFOnIteration, config.plotVOnIteration, config.plotZOnIteration, plotName);

if config.plotTrajectoriesOnIteration
    figure('Name', [plotName, ' Trajectories'])
    a = size(phi, 1);
    b = size(phi, 2);
    hold on
    for i = 1:a
        plot(phi(i, :), 1:b)
    end
end

end

% A pair of functions to handle interpolation with the wrap around in time.
% The point at 0 should be interpolated between the last point and the 
% first point. 

function Vq = interpOnS1(X, V, Xq)
    newX = cat(1, X - ones(size(X)), X, X + ones(size(X)));
    newV = cat(1, V, V, V);
    Vq = interp1( newX, newV, Xq );
end

function Vq = interpOnS1andTime(T, X, V, Tq, Xq)
    newX = cat(1, X - ones(size(X)), X, X + ones(size(X)));
    newV = cat(1, V, V, V);
    %disp(X)
    Vq = interp2( T, newX, newV, Tq, Xq );
end

