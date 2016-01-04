% Test the get_parameterization() function which maps from x on torus to y
% in parameters.  In particular y = [theta, phi] with each in the range 
% [-pi, pi].
%
%  We test the vectorization of this function as well. 

r=0.6; R=1.0;
n = 4000; % number of points
theta = 2*pi*rand(n,1)-pi; % in [-pi,pi]
phi   = 2*pi*rand(n,1)-pi; % in [-pi,pi]
x = [(R+r*cos(theta)).*cos(phi), (R+r*cos(theta)).*sin(phi), r*sin(theta)];
y = get_parameterization(x);

errorphi = norm(y(:,2) - phi);
errortheta = norm(y(:,1) - theta);

if (errorphi > 1e-10)
    fprintf('  Error: phi inverse parameterization is incorrect!\n');
end

if (errortheta > 1e-10)
    fprintf('  Error: theta inverse parameterization is incorrect!\n');
end

if (errorphi < 1e-10 && errortheta < 1e-10)
   fprintf('  get_parameterization() is working as desired!\n'); 
end

