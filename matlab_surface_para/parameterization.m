function [ x ] = parameterization( y )
% PARAMETERIZATION OF TORUS
% MODIFY r AND R AT HERE
R=1; r=0.6;


x=zeros(1,3);
x(1)=(R+r*cos(y(1)))*cos(y(2));
x(2)=(R+r*cos(y(1)))*sin(y(2));
x(3)=r*sin(y(1));

end

