function [ hat_phi,hat_phix,hat_phiy ] = FEEVAL( hatx,haty,nq )
% Return shape value and shape gradient on quadrature points
% in the reference triangle
% The following is the linear case. 
% input:
%   hatx: x-component of the quadrature points.
%   haty: y-component of the quadrature points.
%   nq: number of the quadrature points.
% output:
%   hat_phi: shape value 
%   hat_phix: x-component of the shape gradient.
%   hat_phiy: y-component of the shape gradient.
%
%   Wenyu Lei
%   Dec 30, 2015

hat_phi=zeros(nq,3);
hat_phix=zeros(nq,3);
hat_phiy=zeros(nq,3);

for i=1:nq
        hat_phi(i,1) = 1-hatx(i)-haty(i);
        hat_phi(i,2) = hatx(i);
        hat_phi(i,3) = haty(i);
        
        hat_phix(i,1) = -1; hat_phix(i,2) = 1; hat_phix(i,3) = 0;
        hat_phiy(i,1) = -1; hat_phiy(i,2) = 0; hat_phiy(i,3) = 1;
end

end

