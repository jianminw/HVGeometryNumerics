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
    
    m_C = size(C, 2) - 1;
    m_v = size(path.v, 2);
    T = (0:m_C)/m_C;
    Tq = (-1/2 + 1:m_v ) / m_v;
    %disp(size(T))
    %disp(size(path.phi(:, 1)))
    %disp(size(C))
    C_v = interpOnS1andTime( T, path.phi( :, 1), C, Tq, path.phi( :, 1));
end

function Vq = interpOnS1(X, V, Xq)
    newX = cat(1, X - ones(size(X)), X, X + ones(size(X)));
    newV = cat(1, V, V, V);
    Vq = interp1( newX, newV, Xq );
end

function Vq = interpOnS1andTime(T, X, V, Tq, Xq)
    newX = cat(1, X - ones(size(X)), X, X + ones(size(X)));
    newV = cat(1, V, V, V);
    %disp(X)
    Vq = interp2( T, newX, newV, Tq, Xq );
end