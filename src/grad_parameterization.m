function [ grad ] = grad_parameterization( y )
% Gradient of parameterization.
% returns [3x2] so don't really want to vectorize since it would add third
% dimension complication.
%
R=1.0; r=0.6;
grad=[-r*sin(y(1)).*cos(y(2)), -(R+r*cos(y(1))).*sin(y(2));...
      -r*sin(y(1)).*sin(y(2)), (R+r*cos(y(1))).*cos(y(2));...
       r*cos(y(1)),0];

end

