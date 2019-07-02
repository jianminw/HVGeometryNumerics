% A simple function that creates a linear interpolation between two
% functions with the specified number of time steps. 

% INPUTS: 
% f0 - the function that will be the starting point of the linear
% interpolation. Should be a 1 by n vector. 
% f1 - the function that will be the ending point of the linear
% interpolation. Should have the same size as f0. 
% timeSteps - the number of intervals in time to split the linear
% interpolation by. The total number of time slices will be timeSteps + 1. 

% OUTPUTS:
% f - The function that represents the linear interpolation between f0 and
% f1. It will be a n by (timeSteps + 1) matrix, where the first column is
% f0 and the last column is f1. 

function f = LinearInterpolation(f0, f1, timeSteps)
theta1 = (0:timeSteps) / timeSteps;
theta0 = (timeSteps:-1:0) / timeSteps; 
f = f1' * theta1 + f0' * theta0;
end