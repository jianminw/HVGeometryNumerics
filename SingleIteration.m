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

function newPath = SingleIteration(oldPath, f0, f1, iteration, config)

[v, ft, fx] = schemeStepOne(oldPath);

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

[f, phi] = schemeStepTwo(v, f0, f1, config);

%{
z = zeros(n, m + 1);
for j = 1:(m + 1)
    z(:, j) = interpOnS1( phi(:, j), z_Phi(:, j), (0:n-1)/(n-1));
end
%}

newPath.f = f;
newPath.v = v;
newPath.z = ComputeZFromFV(newPath);
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