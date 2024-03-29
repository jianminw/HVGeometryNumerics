% This is a function to be called when running the numerical scheme. 

% INPUTS:
% f0 - vector that represents the starting function for the admissible
% paths.
% f1 - vector that represents that ending function for the admissible
% paths. Should be the same size as f0. 

% OUTPUTS:
% finalPath - a struct with the fields f, v, z, that represent the last
% path that is computed in the numerical scheme. 

% PARAMETERS:
% maxIterations - the maximum number of iterations that the scheme should
% run through. Currently initialized on line 37. 

function finalPath = NumericalScheme( f0, f1, config )

% Load up config struct from options.m
%config = options();
% Depreciated: config is now passed in to facilitate use of this function
% in the multi scale method.

timeSteps = config.timeSteps;
maxIterations = config.maxIterations;

% Starting path is linear interpolation in time.
% Only f is necessary here for the iterative scheme, but v and z are
% included so that the action of this first path can be computed and
% compared to that of later paths. 

path = struct;
path.f = LinearInterpolation( f0, f1, timeSteps);
path.v = zeros(size(path.f, 1), timeSteps);
path.z = (f1 - f0)' * ones(1, timeSteps); 
disp(ComputeAction(path, config))
bestPath = path;

% Now iterate the scheme. 
% For a start, the stopping condition will just be a set number of 
% iterations. 

action = zeros(maxIterations, 1);
error = zeros(maxIterations, 1);
if config.computeActionMidIteration
    midIterationAction = zeros(maxIterations, 1);
end

for i = 1:maxIterations
    %disp('iteration:')
    %disp(i)
    %[f, v, z] 
    newPath = SingleIteration(path, f0, f1, i, config);
    %newPath.f = f;
    %newPath.v = v;
    %newPath.z = z;
    %disp(f)
    %disp(z)
    
    % check that the new Path is actually better than the old path
    disp(ComputeAction(newPath, config))
    action(i) = ComputeAction(newPath, config);
    %error(i) = CheckAdmissiblePath(newPath);
    
    if config.computeActionMidIteration
        midIterationAction(i) = newPath.midIterationAction;
    end
    
    if ComputeAction(newPath, config) < ComputeAction(bestPath, config)
        bestPath = newPath;
    end
    
    if ComputeAction(newPath, config) >= ComputeAction(path, config)
        % Then the process will have "converged". No more accuracy can
        % be gained due to the error in the discretization/scheme. 
        
        %disp("Sequence has converged.")
        %break;
    end
    path = newPath;
end
if i == maxIterations
    disp("Maximum iterations have been reached.")
end

%{
finalPath.f = bestPath.f;
finalPath.v = bestPath.v;
finalPath.z = bestPath.z;
finalPath.phi = bestPath.phi;
finalPath.action = action;
finalPath.error = error;
%}
finalPath = path;
finalPath.action = action;
finalPath.error = error;


if config.computeActionMidIteration
    finalPath.midIterationAction = midIterationAction;
end
save('OptimalPath.mat', 'finalPath')
end