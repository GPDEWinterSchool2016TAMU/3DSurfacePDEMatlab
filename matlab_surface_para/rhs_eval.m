function [ res ] = rhs_eval( y )
% RHS_EVAL 
% ALSO MODIFY r AND R HERE

r=0.6; R=1.0;


ct=cos(y(1)); st=sin(y(1));
c3p=cos(3*y(2)); s3p=sin(3*y(2));
c3tp=cos(3*y(1)+y(2));
s3tp=sin(3*y(1)+y(2));


res=s3p*c3tp+...
    (9*r^2*ct^2*c3tp*s3p-3*ct*st*s3p*s3tp*r^2+...
    18*r*ct*c3tp*s3p*R-3*st*s3p*s3tp*R*r+...
    6*r^2*c3p*s3tp+9*c3tp*s3p*R^2+10*r^2*c3tp*s3p)...
    /(r^2*(r^2*ct^2+2*ct*r*R+R^2));

% STILL NEED TO CHECK

end

