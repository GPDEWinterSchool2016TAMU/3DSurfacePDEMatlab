function [ res ] = distfunc( p )
% Distance function for torus
% Here R = 1.0 and r=0.6

R=1.0; r=0.6;

res = sqrt((sqrt(p(1)^2+p(2)^2)-R)^2+p(3)^2)-r;

end

