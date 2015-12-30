function [ dist ] = distfunc( x )
% Signed distance function from torus at point x in R^3. 
% This function is vectorized.
%
% input:
%    x = [nx3] list of n points in R^3
% output:
%    dist = [1xn] list of signed distance values

R=1.0; r=0.6;
dist = sqrt((sqrt(x(:,1).^2 + x(:,2).^2)-R).^2+x(:,3).^2)-r;

% want it to return as a 1xn vector
dist = dist';

end

