% A script wrapper for the rest of the functions, for ease of use. 
% n is the number of intervals, both in time and space. 
% f0 and f1 are the functions passed into the numerical scheme. 
% The output a struct, with the field f, v, z corresponding to the past
% path that was computed in the scheme.

config = options();

n = config.spaceIntervals;
% Sine to cosine
%f0 = sin( (1:n) / n * 2 * pi ) / 2;
%f1 = cos( (1:n) / n * 2 * pi ) / 2;

% Sine to negative sine
f0 = sin( (1:n) / n * 4 * pi );
f1 = sin( (1:n) / n * 4 * pi + pi / 2);
scaling = 1;
f1 = scaling * f1;
f0 = scaling * f0;


% narrow bump to wide bump
%f0 = 5 * ( sin( (1:n) / n * pi ) ).^32;
%f1 = 5 * ( sin( ( (1:n) / n + 0.25 ) * pi ) ).^4;

% small and large bumps 
% similar in size
%f0 = 1.2 * exp( - ((1:n) / n - 0.25).^2 * 100 ) + 0.8 * exp( - ((1:n) / n - 0.75).^2 * 100);
%f1 = 0.8 * exp( - ((1:n) / n - 0.25).^2 * 100 ) + 1.2 * exp( - ((1:n) / n - 0.75).^2 * 100);

% small and large bumps 
% dissimilar in size
%f0 = 1 * exp( - ((1:n) / n - 0.25).^2 * 1600 ) + 0.2 * exp( - ((1:n) / n - 0.35).^2 * 1600);
%f1 = 0.2 * exp( - ((1:n) / n - 0.25).^2 * 1600 ) + 1 * exp( - ((1:n) / n - 0.35).^2 * 1600);

path = NumericalScheme( f0, f1 );

% plotting
[n, m] = size(path.f);
% a coloring by time. 
% C = ones(n, 1) * (1:m);
% a coloring by flow maps.
plotting(path, true, true, true, '');

%{ 
Depreciated code.
[C, C_v] = FlowMapColoring( path );
figure('Name', 'f')
mesh(path.f, C)
figure('Name', 'v')
mesh(path.v, C_v)
figure('Name', 'z')
mesh(path.z, C)
%}

figure('Name', 'Action over Iterations')
if ~config.computeActionMidIteration
    plot(path.action, 'DisplayName', 'Action after iteration');
end
if config.computeActionMidIteration
    k = length(path.action);
    y = cat(2, path.midIterationAction, path.action);
    x = cat(2, (1:k)', (1:k)');
    y = reshape(y', [], 1);
    x = reshape(x', [], 1);
    hold on
    plot(x, y);
    scatter(1:k, path.action);
end

if config.plotTrajectory
    figure('Name', 'Trajectories')
    a = size(path.phi, 1);
    b = size(path.phi, 2);
    hold on
    for i = 1:a
        plot( path.phi(i, :), 1:b)
    end
end
%}

%{
figure('Name', 'Error of end of iteration paths')
plot(path.error)
%}