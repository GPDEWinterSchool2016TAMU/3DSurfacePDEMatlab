function [ x ] = parameterization( chi )
% Parameterization 
R=1; r=0.6;


x=zeros(1,3);
x(1)=(R+r*cos(chi(1)))*cos(chi(2));
x(2)=(R+r*cos(chi(1)))*sin(chi(2));
x(3)=r*sin(chi(1));

end

