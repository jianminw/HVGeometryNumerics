% A script wrapper for the rest of the functions, for ease of use. 
% n is the number of intervals, both in time and space. 
% f0 and f1 are the functions passed into the numerical scheme. 
% The output a struct, with the field f, v, z corresponding to the past
% path that was computed in the scheme.

config = options();
mConfig = multiscaleConfig;

format long

n = config.spaceIntervals;

% Sine to cosine
%f0 = sin( (1:n) / n * 2 * pi ) / 2;
%f1 = cos( (1:n) / n * 2 * pi ) / 2;

% Sine to negative sine
%f1 = sin( (1:n) / n * 4 * pi );
%f0 = sin( (1:n) / n * 4 * pi + pi / 2);


% narrow bump to wide bump
%f0 = ( sin( (1:n) / n * pi ) ).^32;
%f1 = ( sin( ( (1:n) / n + 0.25 ) * pi ) ).^4;

% bumps with high-frequency perturbation
%f1 = ( sin( (1:n) / n * pi ) ).^6;
%f0 = ( sin( ( (1:n) / n + 0.25 ) * pi ) ).^6+ 0.15*sin((1:n) / n * 24* pi );

% narrow bumps phase shift
%space=0.15;
%f0 = ( sin( (1:n) / n * pi-space ) ).^16;
%f1 = ( sin( ( (1:n) / n + space) * pi ) ).^16;

% narrow bump to wide bump phase shift
% space=0.15;
% f0 = ( sin( (1:n) / n * pi -space -0.3) ).^64;
% f1 = 0.6*( sin( ( (1:n) / n + space ) * pi ) ).^6;

% narrow bump to wide bump phase shift
% space=0.1;
% f1 = (sin(((1:n) / n  -space -0.1)*pi )).^64;
% f0 = 0.6*( sin(  ((1:n) / n + space)  * pi )).^6;

% for Irina
space=0.1;
f1 = (sin(((1:n) / n  -space -0.05)*pi )).^70+(sin(((1:n) / n  -space -0.1)*pi )).^170;
f0 = 0.9*( sin(  ((1:n) / n + space)  * pi )).^6 -0.25*(sin(((1:n) / n  +space +0.25)*pi )).^70 ;

% narrow bumps phase shift
% space=0.15;
% f0 = ( sin( (1:n) / n * pi-space ) ).^16;
% f1 = ( sin( ( (1:n) / n + space) * pi ) ).^16;


% small and large bumps 
% similar in size
%f0 = 1.2 * exp( - ((1:n) / n - 0.25).^2 * 100 ) + 0.8 * exp( - ((1:n) / n - 0.75).^2 * 100);
%f1 = 0.8 * exp( - ((1:n) / n - 0.25).^2 * 100 ) + 1.2 * exp( - ((1:n) / n - 0.75).^2 * 100);

% small and large bumps 
% dissimilar in size
%f0 = 1 * exp( - ((1:n) / n - 0.25).^2 * 1600 ) + 0.2 * exp( - ((1:n) / n - 0.35).^2 * 1600);
%f1 = 0.2 * exp( - ((1:n) / n - 0.25).^2 * 1600 ) + 1 * exp( - ((1:n) / n - 0.35).^2 * 1600);

% f00 = zeros( 1, n );
% delx=1/n;
% x=0:delx:0.999999;
% f0 = 3*exp(-64*(x-0.5).^4);
% f1 = 2*(0.4-2*abs(x-0.5));
% f1=max(f00,f1);

% narrow bumps phase shift
% space=0.2;
% f0 = 0.5*( sin( (1:n) / n * pi-space ) ).^32;
%  f00 = zeros( 1, n );
%  delx=1/n;
%  x=0:delx:0.999999;
% f1 = (1-5*abs((x-0.5+space)));
%  f1=max(f00,f1);

scaling = 1;
f1 = scaling * f1;
f0 = scaling * f0;

%path = NumericalScheme( f0, f1, config );
path = multiscaleScheme( f0, f1, mConfig );

% plotting
[n, m] = size(path.f);
plotting(path, true, true, true, '');

%{ 
Depreciated code. Moved to plotting.m
[C, C_v] = FlowMapColoring( path );
figure('Name', 'f')
mesh(path.f, C)
figure('Name', 'v')
mesh(path.v, C_v)
figure('Name', 'z')
mesh(path.z, C)
%}

% Some plots only really make sense over lots of iterations, so it would
% seemt that the only place where these plots are created should be the end
% of the entire scheme. These are not moved to plotting.m

%{
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
%}

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