% A script wrapper for the rest of the functions, for ease of use. 
% n is the number of intervals, both in time and space. 
% f0 and f1 are the functions passed into the numerical scheme. 
% The output a struct, with the field f, v, z corresponding to the past
% path that was computed in the scheme.

n = 200;
% Sine to cosine
%f0 = sin( (1:n) / n * 2 * pi ) / 2;
%f1 = cos( (1:n) / n * 2 * pi ) / 2;
% narrow bump to wide bump
f0 = ( sin( (1:n) / n * pi ) ).^32;
f1 = ( sin( ( (1:n) / n + 0.25 ) * pi ) ).^4;


path = NumericalScheme( f0, f1 );

% plotting
[n, m] = size(path.f);
% a coloring by time. 
% C = ones(n, 1) * (1:m);
% a coloring by flow maps.
C = FlowMapColoring( path );
figure('Name', 'f')
mesh(path.f, C)
figure('Name', 'v')
mesh(path.v, C)
figure('Name', 'z')
mesh(path.z, C)

function C = FlowMapColoring( path )
    C = zeros(size(path.f));
    for j = 1:size(path.f, 2)
        C( :, j) = interpOnS1( path.phi( :, j) , path.phi( :, 1), path.phi(:, 1));
    end
end

function Vq = interpOnS1(X, V, Xq)
    newX = cat(1, X - ones(size(X)), X, X + ones(size(X)));
    newV = cat(1, V, V, V);
    Vq = interp1( newX, newV, Xq );
end