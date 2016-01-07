function val = exact ( x )
% Exact solution for all of R^3

val = cos(pi*x(:,1)).*cos(pi*x(:,2)).*cos(pi*x(:,3));  %[nx1]

end