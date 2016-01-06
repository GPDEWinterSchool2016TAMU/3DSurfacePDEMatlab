function [ y ] = get_parameterization( p )
% Return the parametric coordinates for 3D point p on the surface of Torus
%
% Torus parameterization:
%  y = [theta, phi]       with each defined on [-pi,pi]
%  p = [(R+r*cos(theta))*cos(phi), (R+r*cos(theta))*sin(phi), r*sin(theta)]
%
%  so the inverse mapping from p to y is
%
%  y(2) = phi = atan2(p(2),p(1))
%
%                 { pi - asin(x(3)/r)  if p(3) > 0 and p(1)^2 + p(2)^2 < R^2
%  y(1) = theta = {-pi - asin(x(3)/r)  if p(3) < 0 and p(1)^2 + p(2)^2 < R^2
%                 { asin(x(3)/r)       if p(1)^2 + p(2)^2 > R^2
%
%
%  input:
%    p = [nx3] vector of n points on torus
%  output:
%    y = [nx2] vector of n parameteric coordinates [theta(:), phi(:)]
%

r=0.6; R=1.0; % should we make these global parameters or pass them in?

[n,~] = size(p);

% initialize size
y=zeros(n,2);  % nx2

% solve for phi
y(:,2) = atan2(p(:,2),p(:,1));

% solve for theta
rad = sqrt(p(:,1).^2+p(:,2).^2);
asinzr = real(asin(p(:,3)/r));
for j=1:n % I wish we could vectorize this completely
    if (rad(j)<R && p(j,3)>0)
        y(j,1) = pi-asinzr(j);
    elseif (rad(j)<R && p(j,3)<0);
        y(j,1) = -pi-asinzr(j);
    else
        y(j,1) = asinzr(j);
    end
end

end

