function [ res ] = rhs_eval( x )
% Right hand side evaluation function
%
% f = -surflaplacian uexact  (formula computed with Maple)
%
%  Since our exact solution is a function of theta and phi, our rhs
%  function will be also.
%
%  input:
%    x = [nx3] vector of n coordinates in R^3
%  output:
%    val = [nx1] vector of n values of rhs()
%
%  Since the rhs is derived from the exact solution which is only defined
%  on the torus but extended constant in the normal direction, the rhs is
%  also extended constant in the normal direction, using
%  function_extension().
%  
%  p(x) = closest point to torus = x - d(x) * grad d(x)
%

r=0.6; R=1.0;

p = function_extension(x); % [nx3]
y = get_parameterization( p ); % [nx2] = [theta(:), phi(:)]
ct=cos(y(:,1)); st=sin(y(:,1));
c3p=cos(3*y(:,2)); s3p=sin(3*y(:,2));
c3tp=cos(3*y(:,1)+y(:,2)); s3tp=sin(3*y(:,1)+y(:,2));


res=s3p.*c3tp + ...
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

