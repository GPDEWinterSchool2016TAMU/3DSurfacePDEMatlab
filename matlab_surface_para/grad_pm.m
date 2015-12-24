function [ grad_chi ] = grad_pm( y )
% GRAD_PM 

R=1.0; r=0.6;
grad_chi=[-r*sin(y(1))*cos(y(2)), -(R+r*cos(y(1)))*sin(y(2));...
          -r*sin(y(1))*sin(y(2)), (R+r*cos(y(1)))*cos(y(2));...
          r*cos(y(1)),0];

end

