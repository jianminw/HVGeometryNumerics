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

function finalPath = NumericalScheme( f0, f1 )

% Intervals in space are given by f0 and f1
% Time step needs to be set. 
% Idea 1: Use the same as steps in space. 

timeSteps = length(f0);

% Starting path is linear interpolation in time.
% Only f is necessary here for the iterative scheme, but v and z are
% included so that the action of this first path can be computed and
% compared to that of later paths. 

path = struct;
path.f = LinearInterpolation( f0, f1, timeSteps);
path.v = zeros(size(path.f));
path.z = (f1 - f0)' * ones(1, timeSteps+1); 
disp(ComputeAction(path))

% Now iterate the scheme. 
% For a start, the stopping condition will just be a set number of 
% iterations. 

maxIterations = 2;

for i = 1:maxIterations
    %disp('iteration:')
    %disp(i)
    [f, v, z] = SingleIteration(path, f0, f1);
    newPath.f = f;
    newPath.v = v;
    newPath.z = z;
    %disp(f)
    %disp(z)
    
    % check that the new Path is actually better than the old path
    format long
    %disp(ComputeAction(path))
    disp(ComputeAction(newPath))
    
    if ComputeAction(newPath) >= ComputeAction(path)
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

save('OptimalPath.mat', 'path')
finalPath.f = path.f;
finalPath.v = path.v;
finalPath.z = path.z;

end