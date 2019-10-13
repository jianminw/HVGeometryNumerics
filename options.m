% a dedicated file to keep track of various options and configurations
% while running the numerical scheme

function config = options()
    config = struct;
    % every option should be a field in the struct. 
    config.plotFOnIteration = false;
    config.plotVOnIteration = false;
    config.plotZOnIteration = false;
    config.plotTrajectoriesOnIteration = false;
    
    % switch between first and second order schemes. 
    config.odeScheme = 2;
    
    config.maxIterations = 10;
    config.spaceIntervals = 50;
    config.timeSteps = 200;
    
    config.computeActionMidIteration = false;
    
    config.plotTrajectory = true;
end