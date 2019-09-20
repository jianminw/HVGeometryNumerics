% a dedicated file to keep track of various options and configurations
% while running the numerical scheme

function config = options()
    config = struct;
    % every option should be a field in the struct. 
    config.plotFOnIteration = false;
    config.plotVOnIteration = true;
    config.plotZOnIteration = false;
    
    config.maxIterations = 20;
    config.spaceIntervals = 50;
    config.timeSteps = 800;
    
    config.computeActionMidIteration = true;
end