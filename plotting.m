% A file to plot the various functions (f, v, z, etc.) so it can be used in
% multiple files. 

% First input is the path stuct, containing f, v, z, data, etc)
% Next three inputs are whether to plot f, v, or z
% Final input is a string that will be attached to the name of each figure.
%       This should contain an iteration number, or something similar. 

function plotting(path, plotF, plotV, plotZ, plotName)
[C_f, C_v, C_z] = FlowMapColoring( path );
if plotF
    figure('Name', ['f ' plotName])
    mesh( path.f, C_f )
end

if plotV
    figure('Name', ['v ' plotName])
    mesh( path.v, C_v )
end

if plotZ
    figure('Name', ['z ' plotName])
    mesh( path.z, C_z )
end

end


function [C_f, C_v, C_z] = FlowMapColoring( path )
    C_f = zeros(size(path.f));
    for j = 1:size(path.f, 2)
        C_f( :, j) = interpOnS1( path.phi( :, j) , path.phi( :, 1), path.phi(:, 1));
    end
    
    m_C = size(C_f, 2) - 1;
    m_v = size(path.v, 2);
    T = (0:m_C)/m_C;
    Tq = (-1/2 + 1:m_v ) / m_v;
    %disp(size(T))
    %disp(size(path.phi(:, 1)))
    %disp(size(C))
    C_v = interpOnS1andTime( T, path.phi( :, 1), C_f, Tq, path.phi( :, 1));
    
    m_z = size(path.z, 2);
    Tq = (-1/2 + 1:m_z ) / m_z;
    C_z = interpOnS1andTime( T, path.phi( :, 1), C_f, Tq, path.phi( :, 1));
end

%{
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
%}