function [xij] = find_xij(xi, xj, di, dj, h)
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

if (di>-h && dj>-h ) % find d(xij) = h
    t = abs(dj-h)/(abs(di-h)+abs(dj-h));
else % find d(xij) = -h
    t = abs(dj+h)/(abs(di+h)+abs(dj+h));
end

% Assert 0 <= t <= 1  so that xij lies between xi and xj
xij = xj + t*(xi-xj);
return;

end