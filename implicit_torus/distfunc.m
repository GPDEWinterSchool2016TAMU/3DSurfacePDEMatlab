function [ dist ] = distfunc( x )
% Distance function for torus at point x (This function is vectorized)
%
% input:
%    x = [nx3] list of n points in 3D
% output:
%    dist = [1xn] list of distance values

R=1.0; r=0.6;
dist = sqrt((sqrt(x(:,1).^2 + x(:,2).^2)-R).^2+x(:,3).^2)-r;

% want it to return as a 1xn vector
dist = dist';

end

