function startWithDefinedPath()

% This is mainly a script to analyse what the numerics does for the paths
% between constant functions. 

config = options();

timeSteps = config.timeSteps;
spaceIntervals = config.spaceIntervals;
maxIterations = config.maxIterations;

H = 10;
L = 0.3;

% Only need to define f for the scheme to start going. 
path = struct;
path.f = zeros(spaceIntervals, timeSteps);
path.v = zeros(spaceIntervals, timeSteps);
path.z = zeros(spaceIntervals, timeSteps);
for i = 1:spaceIntervals
    for j = 1:timeSteps
        x = i / spaceIntervals;
        t = j / timeSteps;
        if t < 1/3
            height = 3 * H * t;
            if x < L
                path.f(i, j) = height;
                path.z(i, j) = 3 * H;
            end
        elseif t > 2/3
            height = 3 * H * (t - 2/3);
            if x > 1 - L
                path.f(i, j) = height;
                path.z(i, j) = 3 * H;
            else
                path.f(i, j) = H;
            end
        else
            cutoff = L + (1 - 2 * L) * (t - 1/3) * 3;
            if x < cutoff
                path.f(i, j) = H;
                path.v(i, j) = (1 - 2 * L) * x / cutoff * 3;
            else
                path.v(i, j) = (1 - 2 * L) * (1 - x) / (1 - cutoff) * 3;
            end
        end
    end
end
n = size(path.f, 1);
DeltaX = 1 / n;

Dx = circshift( eye(n), -1) - eye(n);
Dx = Dx / DeltaX ;
Dxx = Dx * (-Dx');

vx = Dx * path.v;
vxx = Dxx * path.v;
Action = L2Squared(path.v) + L2Squared(path.z);
Action = Action + config.lambda * L2Squared(vx);
Action = Action + config.epsilon * L2Squared(vxx);

disp(Action)

bestPath = path;
plotName = 'Start';
figure('Name', ['f ' plotName])
mesh( path.f )

figure('Name', ['v ' plotName])
mesh( path.v )

figure('Name', ['z ' plotName])
mesh( path.z )

f0 = zeros(1, spaceIntervals);
f1 = H * ones(1, spaceIntervals);

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

path = finalPath;

plotting(path, true, true, true, '');

if config.plotTrajectory
    figure('Name', 'Trajectories')
    a = size(path.phi, 1);
    b = size(path.phi, 2);
    hold on
    for i = 1:a
        plot( path.phi(i, :), 1:b)
    end
end

end

function y = L2Squared(f)
f2 = f .* f;
y = (mean(mean(f2))) / 2;
end