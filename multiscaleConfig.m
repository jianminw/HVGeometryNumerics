function mConfig = multiscaleConfig()
% A function to provide configuration options for the multiscale scheme. 
% Almost all options here should be arrays, representing options for each
% round of the multi scale scheme. 
mConfig = struct;

mConfig.maxIterations = [5, 5, 2];
mConfig.spaceIntervals = [25, 100, 200];
mConfig.timeSteps = [50, 200, 400];

mConfig.odeScheme = [1, 1, 1];
end