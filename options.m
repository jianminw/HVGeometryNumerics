% a dedicated file to keep track of various options and configurations
% while running the numerical scheme

function config = options()
    config = struct;
    % every option should be a field in the struct. 
    config.plotFOnIteration = false;
    config.plotVOnIteration = false;
    config.plotZOnIteration = false;
    
    config.maxIterations = 20;
    config.spaceIntervals = 200;
    config.timeSteps = 100;
    
    config.computeActionMidIteration = true;
end