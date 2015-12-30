function res = exact (x)
% Exact solution
% Here x is a point on the surface
% y is the corresponding parametric coordinates

 y  = get_parameterization( x );
res=sin(3*y(2))*cos(3*y(1)+y(2));
end