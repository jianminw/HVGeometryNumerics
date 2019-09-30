% a dedicated file to keep track of various options and configurations
% while running the numerical scheme

function config = options()
    config = struct;
    % every option should be a field in the struct. 
    config.plotFOnIteration = false;
    config.plotVOnIteration = false;
    config.plotZOnIteration = false;
    config.plotTrajectoriesOnIteration = false;
    
    config.maxIterations = 1000;
    config.spaceIntervals = 100;
    config.timeSteps = 200;
    
    config.computeActionMidIteration = false;
    
    config.plotTrajectory = false;
end