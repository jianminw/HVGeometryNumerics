% A script wrapper for the rest of the functions, for ease of use. 
% n is the number of intervals, both in time and space. 
% f0 and f1 are the functions passed into the numerical scheme. 
% The output a struct, with the field f, v, z corresponding to the past
% path that was computed in the scheme.

n = 500;
%f0 = sin( (1:n) / n * 2 * pi ) / 2;
%f1 = cos( (1:n) / n * 2 * pi ) / 2;
f0 = ( sin( (1:n) / n * pi ) ).^32;
f1 = ( sin( ( (1:n) / n + 0.05 ) * pi ) ).^4;
path = NumericalScheme( f0, f1 );

% plotting
figure('Name', 'f')
mesh(path.f)
figure('Name', 'v')
mesh(path.v)
figure('Name', 'z')
mesh(path.z)