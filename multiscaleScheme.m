% A function to be called for running a multiscale scheme. 
% From each level, the final v is preserved/interpolated for the start of
% the new level. 
% Options for this scheme are found in multiscaleConfig.


function finalPath = multiscaleScheme(f0, f1, mConfig)
% Basically just go through the levels of the multiscale scheme. 
% Assuming good formatting in the config file. 
schemeLevels = length(mConfig.maxIterations);
        
path = struct;

for level = 1:schemeLevels
    disp(['Level: ', int2str(level)])
    % Read the config for this level and set up the initial and final data
    timeSteps = mConfig.timeSteps(level);
    spaceIntervals = mConfig.spaceIntervals(level);
    maxIterations = mConfig.maxIterations(level);
    
    nOld = length(f0);
    nNew = spaceIntervals;
    X = (1:nOld)/nOld;
    Xq = (1:nNew)/nNew;
    f0Level = interpOnS1(X', f0', Xq')';
    f1Level = interpOnS1(X', f1', Xq')';
    
    % Initiate config for the level. 
    config = options();
    config.timeSteps = timeSteps;
    config.spaceIntervals = spaceIntervals;
    config.maxIterations = maxIterations;
    config.odeScheme = mConfig.odeScheme(level);
    
    if level == 1
        path.f = LinearInterpolation( f0Level, f1Level, timeSteps);
        path.v = zeros(size(path.f, 1), timeSteps);
        path.z = (f1Level - f0Level)' * ones(1, timeSteps); 
    end
    if level > 1
        v = pathInterpolation(path.v, spaceIntervals, timeSteps);

        [f, ~] = schemeStepTwo(v, f0Level, f1Level, config);

        path.f = f;
        path.v = v;
        path.z = ComputeZFromFV(path);
    end
    
    % Now run through the scheme. 
    for i = 1:maxIterations
        disp(['Iteration: ', int2str(i)])
        newPath = SingleIteration(path, f0Level, f1Level, i, config);
        
        path = newPath;
    end
end

finalPath = path;
end