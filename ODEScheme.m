% Separate out the part where I solve the ODE for better clarity. 

function phi = ODEScheme( v, config )
if config.odeScheme == 1
    phi = FirstOrderScheme(v);
end
if config.odeScheme == 2
    phi = SecondOrderScheme(v);
end
end

function phi = FirstOrderScheme( v )
m = size(v, 2);
n = size(v, 1);

validityCheck = false;
timeSteps = m;
phi = ones(n, timeSteps + 1);

while (~validityCheck)
    validityCheck = true;
    phi = ones(n, timeSteps + 1);
    phi(:, 1) = (1:n) / n;
    T_v = ( (0:m-1) + 1/2) / (m-1);
    X_v =  phi(:, 1);

    for j = 2:(timeSteps + 1)
        % need to mess with this to make it work for smaller time steps. 
        Dphi = interpOnS1andTime( T_v, X_v, v, (j-1)/timeSteps, phi(:, j-1) );
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
            %timeSteps = timeSteps * 2;
            %disp(j)
            %disp(timeSteps)
            %validityCheck = false;
            %break;
        end
    end
end
end

function phi = SecondOrderScheme( v )
m = size(v, 2);
n = size(v, 1);

validityCheck = false;
timeSteps = m;
phi = ones(n, timeSteps + 1);

while (~validityCheck)
    validityCheck = true;
    phi = ones(n, timeSteps + 1);
    phi(:, 1) = (1:n) / n;
    T_v = ( (0:m-1) + 1/2) / (m-1);
    X_v =  phi(:, 1);

    for j = 2:(timeSteps + 1)
        % need to mess with this to make it work for smaller time steps. 
        Dphi = interpOnS1andTime( T_v, X_v, v, (j-1)/timeSteps, phi(:, j-1) );
        midpoints = phi(:, j-1) + Dphi / timeSteps / 2;
        Dphi2 = interpOnS1andTime( T_v, X_v, v, (j-1/2)/timeSteps, midpoints );
        %Dphi = interpOnS1( phi(:, 1), v(:, j), phi(:, j-1) );
        phi(:, j) = midpoints + Dphi2 / timeSteps / 2;
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
            %timeSteps = timeSteps * 2;
            %disp(j)
            %disp(timeSteps)
            %validityCheck = false;
            %break;
        end
    end
end
end

%{

% A pair of functions to handle interpolation with the wrap around in time.
% The point at 0 should be interpolated between the last point and the 
% first point. 

function Vq = interpOnS1(X, V, Xq)
    newX = cat(1, X - ones(size(X)), X, X + ones(size(X)));
    newV = cat(1, V, V, V);
    Vq = interp1( newX, newV, Xq, 'spline');
end

function Vq = interpOnS1andTime(T, X, V, Tq, Xq)
    newX = cat(1, X - ones(size(X)), X, X + ones(size(X)));
    newV = cat(1, V, V, V);
    %disp(X)
    Vq = interp2( T, newX, newV, Tq, Xq, 'spline');
end
%}