% A function that takes in a path, and interpolates onto a different grid. 
% Considering replaceing the start and end values with interpolations from
% f0 and f1 originally, so that details lost in lower resolutions can be
% recovered. 
% I think that the best idea would be to just interpolate v and the start
% with step 2. 

function v = pathInterpolation(oldV, spaceIntervals, timeSteps)
n = size(oldV, 1);
m = size(oldV, 2);
% First set out the grid for the old v. 
Xold = (1:n)' / n;
Told = ( (0:m-1) + 1/2) / (m-1);

n = spaceIntervals;
m = timeSteps;
Xnew = (1:n)' / n;
Tnew = ( (0:m-1) + 1/2) / (m-1);

%disp(size(Told))
%disp(size(Xold))
%disp(size(oldV))
v = interpOnS1andTime( Told, Xold, oldV, Tnew, Xnew);
end