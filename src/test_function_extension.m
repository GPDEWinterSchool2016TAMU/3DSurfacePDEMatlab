%  A test to ensure that our function_extension.m is working properly as a
%  vectorized function.
%
%  It maps from a point in space to the closest point on the torus using
%  the distance function and its gradient.
%
%  Spencer Patty
%  Dec 29, 2015


% pick a random point on the surface of the torus
R=1; r=0.6;
n = 4; % number of points
phi=2*pi*rand(n,1)-pi;
theta=2*pi*rand(n,1)-pi;

% The closest point on torus
xtrue = [(R+r*cos(theta)).*cos(phi), (R+r*cos(theta)).*sin(phi), r*sin(theta)];

% We will test in a neighborhood of the torus in particular, we will check
% in the band r \pm 0.3 which is guaranteed to have a unique closest point.
r_ext = r + 2*0.3*rand(n,1)-0.3;

% The point in space
p = [(R+r_ext.*cos(theta)).*cos(phi), (R+r_ext.*cos(theta)).*sin(phi), r_ext.*sin(theta)];

% The computed closest point on torus to p
x = function_extension(p);

if ( norm(x-xtrue) > 1e-10 )
    fprintf('  Error: function_extension.m is not working correctly!\n');
else
    fprintf('  function_extension.m is working correctly.\n');
end

