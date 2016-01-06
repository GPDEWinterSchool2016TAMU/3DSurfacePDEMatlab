function [ grad_chi ] = grad_pm( chi )
% Gradient of parametrizaion.

R=1.0; r=0.6;
grad_chi=[-r*sin(chi(1))*cos(chi(2)), -(R+r*cos(chi(1)))*sin(chi(2));...
          -r*sin(chi(1))*sin(chi(2)), (R+r*cos(chi(1)))*cos(chi(2));...
          r*cos(chi(1)),0];

end

