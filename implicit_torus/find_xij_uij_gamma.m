function [xij, uij] = find_xij_uij_gamma(xi, xj, di, dj, ui, uj)
% [xij, uij] = find_xij_uij_gamma(xi, xj, di, dj, ui, uj) 
%
% This function finds the xij coordinate on the line between xi and xj
% with di = d_h(xi) and dj = d_h(xj) with d_h(xij) = 0,  ie finds the
% intersection of the line segment [xi, xj] with Gamma_h and returns the 
% coordinate of the point on Gamma_h. It then evaluates the linearly
% interpolated solution uij as the interpolation between ui and uj at xij.
%
% It is necessary that di and dj have opposite sign values or have one of 
% them be zero.  If both are zero, then we have an issue since there is no 
% unique intersection point.  It is left to the calling function to avoid
% this issue.
%
%  input: 
%    xi = [1x3] coordinate point in R^3
%    xj = [1x3] coordinate point in R^3
%    di = [1x1] d_h(xi) distance function value at xi
%    dj = [1x1] d_h(xj) distance function value at xj
%    ui = [1x1]  u(xi)  solution function value at xi
%    uj = [1x1]  u(xj)  solution function value at xj
%  output:
%    xij = [1x3] coordinate point on Gamma_h in R^3
%    uij = [1x1] interpolated solution function value at xij
%
% Spencer Patty
% Dec 30, 2015

t = abs(dj)/(abs(di)+abs(dj));
% Assert 0 <= t <= 1  so that xij lies between xi and xj
xij = xj + t*(xi-xj);
uij = uj + t*(ui-uj);

end