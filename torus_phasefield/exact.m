function val = exact ( x )
% Exact solution for all of R^3
%
%  input:
%    x = [nx3] vector of points in R^3
%  output:
%    val = [nx1] vector of solution values
%
% Since the exact function is only defined on the surface, we extend it to be
% constant in the normal direction.  Thus we first project from general
% point x in R^3 to p on torus, then we convert it to parametric coordinate
% and then solve for the solution using those parameteric coordinates.

p   = function_extension( x ); % [nx3]
y   = get_parameterization( p );  %[nx2]
val = sin(3*y(:,2)).*cos(3*y(:,1)+y(:,2));  %[nx1]

end