function [ val ] = rhs_eval( y )
% Right hand side evaluation function
%
% f = -surflaplacian uexact  (formula computed with Maple)
%
%  Since our exact solution is a function of theta and phi, our rhs
%  function will be also.
%
%  input:
%    y = [nx2] vector of n coordinates in parametric domain in R^2
%  output:
%    val = [nx1] vector of n values of rhs()
%


r=0.6; R=1.0;

ct=cos(y(:,1)); st=sin(y(:,1));
c3p=cos(3*y(:,2)); s3p=sin(3*y(:,2));
c3tp=cos(3*y(:,1)+y(:,2)); s3tp=sin(3*y(:,1)+y(:,2));


val=s3p.*c3tp + ...
      ( 9*r^2*ct.*ct.*s3p.*c3tp ...
       -3*r^2*ct.*st.*s3p.*s3tp ...
      +18*r*R*ct.*s3p.*c3tp ...
       -3*r*R*st.*s3p.*s3tp ...
       +6*r^2*c3p.*s3tp ...
       +9*R^2*s3p.*c3tp ...
      +10*r^2*s3p.*c3tp ) ./(r^2*( r^2*ct.^2 ...
                                   +2*r*R*ct ...
                                   +R^2      ) );

% STILL NEED TO CHECK

end

