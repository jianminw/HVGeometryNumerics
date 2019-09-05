% A file to plot the various functions (f, v, z, etc.) so it can be used in
% multiple files. 

% First input is the path stuct, containing f, v, z, data, etc)
% Next three inputs are whether to plot f, v, or z
% Final input is a string that will be attached to the name of each figure.
%       This should contain an iteration number, or something similar. 

function plotting(path, plotF, plotV, plotZ, plotName)
[C, C_v] = FlowMapColoring( path );
if plotF
    figure('Name', ['f ' plotName])
    mesh( path.f, C )
end

if plotV
    figure('Name', ['v ' plotName])
    mesh( path.v, C_v )
end

if plotZ
    figure('Name', ['z ' plotName])
    mesh( path.z, C )
end

end


function [C, C_v] = FlowMapColoring( path )
    C = zeros(size(path.f));
    for j = 1:size(path.f, 2)
        C( :, j) = interpOnS1( path.phi( :, j) , path.phi( :, 1), path.phi(:, 1));
    end
    
    C_v = zeros(size(path.v));
    scaleFactor = ( size(C, 2) - 1 ) / ( size(C_v, 2) - 1) ;
    for j = 1:( (size(path.f, 2) - 1) / scaleFactor + 1 )
        C_v( :, j) = C( :, (j - 1) * scaleFactor + 1 );
    end
end

function Vq = interpOnS1(X, V, Xq)
    newX = cat(1, X - ones(size(X)), X, X + ones(size(X)));
    newV = cat(1, V, V, V);
    Vq = interp1( newX, newV, Xq );
end