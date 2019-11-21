% a dedicated file to keep track of various options and configurations
% while running the numerical scheme

function config = options()
    config = struct;
    % every option should be a field in the struct. 
    
    % Options regarding how the scheme runs. 
    config.plotFOnIteration = false;
    config.plotVOnIteration = false;
    config.plotZOnIteration = false;
    config.plotTrajectoriesOnIteration = false;
    
    config.computeActionMidIteration = false;
    
    config.plotTrajectory = true;
    
    % switch between first and second order schemes. 
    config.odeScheme = 2;
    
    % Constants of the problem/action
    config.lambda = 1;
    config.epsilon = 0.1;
    
    % Constants related to the scheme. 
    config.maxIterations = 1;
    config.spaceIntervals = 50;
    config.timeSteps = 100;
    
end