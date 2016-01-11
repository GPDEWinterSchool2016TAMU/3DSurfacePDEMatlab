function [ hat_phi,hat_phix,hat_phiy ] = FEEVAL( yhat,nq )
% Return shape value and shape gradient on quadrature points
% in the reference triangle
% The following is the linear case. 
% input:
%   yhat: [nqx2] list of nq quadrature points [x y]
%   nq: number of the quadrature points.
% output:
%   hat_phi: shape value 
%   hat_phix: x-component of the shape gradient.
%   hat_phiy: y-component of the shape gradient.
%
%   Wenyu Lei
%   Dec 30, 2015

hat_phi = [1-yhat(:,1) - yhat(:,2), yhat(:,1), yhat(:,2)];

hat_phix = [-1*ones(nq,1), 1*ones(nq,1), zeros(nq,1)];
hat_phiy = [-1*ones(nq,1), zeros(nq,1), 1*ones(nq,1)];

end

