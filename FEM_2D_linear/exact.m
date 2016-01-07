function [ val ] = exact( x )
% Exact solution
val = cos(pi*x(:,1)).*cos(2*pi*x(:,2));

end

