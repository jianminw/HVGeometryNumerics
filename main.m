n = 200;
%f0 = sin( (0:n-1) / (n-1) * 2 * pi ) / 2;
%f1 = cos( (0:n-1) / (n-1) * 2 * pi ) / 2;
f0 = 10 * ( sin( (0:n-1) / (n - 1) * pi ) ).^32;
f1 = 10 * ( sin( ((0:n-1) / (n - 1) + 0.05 ) * pi ) ).^4;
path = NumericalScheme( f0, f1 );
figure('Name', 'f')
mesh(path.f)
figure('Name', 'v')
mesh(path.v)
figure('Name', 'z')
mesh(path.z)