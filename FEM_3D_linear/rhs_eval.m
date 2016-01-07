function [ val ] = rhs_eval( x )
% Right hand side evaluation function

val = (3*pi^2+1)*cos(pi*x(:,1)).*cos(pi*x(:,2)).*cos(pi*x(:,3));
end

