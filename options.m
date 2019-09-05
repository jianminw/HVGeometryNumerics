% a dedicated file to keep track of various options and configurations
% while running the numerical scheme

function config = options()
    config = struct;
    % every option should be a field in the struct. 
    config.plotFOnIteration = false;
    config.plotVOnIteration = true;
    config.plotZOnIteration = false;
    
    config.maxIterations = 3;
    config.spaceIntervals = 97;
    config.timeSteps = config.spaceIntervals;
    
    config.computeActionMidIteration = true;
end