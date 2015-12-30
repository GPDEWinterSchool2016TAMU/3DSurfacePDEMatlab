function [xij] = find_xij_gamma(xi, xj, di, dj)
% find_xij determines the coordinate xij along line between xi, xj where 
% d(xij) = \pm h. We have assumed that the solution is linear along this
% line and we have values di and dj opposite sides of \pm h.  Thus xij is 
% between xi and xj.  The same formula works for xij outside, but we don't 
% want this to happen in the application since then the proposed 
% tetrahedron is not a subset of the original tetrahedron.
%
% Assert that either h in interval [di, dj]  or  -h in interval [di,dj] but
% not both.
%
% Spencer Patty
% Dec 23, 2015

t = abs(dj)/(abs(di)+abs(dj));
% Assert 0 <= t <= 1  so that xij lies between xi and xj
xij = xj + t*(xi-xj);
return;

end