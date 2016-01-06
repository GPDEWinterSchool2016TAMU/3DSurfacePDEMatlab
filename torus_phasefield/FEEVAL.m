function [ hat_phi,hat_phix,hat_phiy,hat_phiz ] = FEEVAL( hatx,haty,hatz,nq )
% Finite element infomation on the reference cell
% input:
%   hatx,haty,hatz: [1xnq] x,y,z components of quadrature points
%   nq: number of quadrature points
% output:
%   hat_phi: [nqx4] shape values on quadrature points
%   hat_phix, hat_phiy, hat_phiz: [nqx4] x,y,z component of shape gradient on 
%                                        quadrature points
hat_phi=zeros(nq,4);
hat_phix=zeros(nq,4);
hat_phiy=zeros(nq,4);
hat_phiz=zeros(nq,4);

for i=1:nq
        hat_phi(i,1) = 1-hatx(i)-haty(i)-hatz(i);
        hat_phi(i,2) = hatx(i);
        hat_phi(i,3) = haty(i);
        hat_phi(i,4) = hatz(i);
        
        hat_phix(i,1) = -1; hat_phix(i,2) = 1; hat_phix(i,3) = 0; hat_phix(i,4) = 0;
        hat_phiy(i,1) = -1; hat_phiy(i,2) = 0; hat_phiy(i,3) = 1; hat_phiy(i,4) = 0;
        hat_phiz(i,1) = -1; hat_phiz(i,2) = 0; hat_phiz(i,3) = 0; hat_phiz(i,4) = 1;
end

end

