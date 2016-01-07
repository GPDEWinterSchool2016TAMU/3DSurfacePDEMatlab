function [ val ] = rhs_eval( x )
% Right hand side evaluation
val=(5*pi*pi+2)*cos(pi*x(:,1)).*cos(2.*pi*x(:,2));

end

