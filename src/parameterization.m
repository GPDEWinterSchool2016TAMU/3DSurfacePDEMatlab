function [ x ] = parameterization( y )
% Parameterization 
%  chi(y) = x \]in R^3 on torus
%
R=1; r=0.6;

[n,~] = size(y);

x=zeros(n,3);
x(:,1)=(R+r*cos(y(:,1))).*cos(y(:,2));
x(:,2)=(R+r*cos(y(:,1))).*sin(y(:,2));
x(:,3)=r*sin(y(:,1));

end
