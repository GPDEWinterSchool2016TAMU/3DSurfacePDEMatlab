function [nq,q,w ] = volquad( x,deg )
% Quadrature rule on a tetrahedron.
% input: 
%   x: 4 by 3 matrix. Each row indicates the coordinates of the vertiex
%      of the triangle.
%   deg: degree of the quadrature rule.
% output:
%   nq: number of quadrature points.
%   q: quadrature point list.
%   w: weights for cooresponding quadrature points.
% ref: Theory and practice of finite elements 
%      By A. Ern and JL Guermond
%      Page 361
%
% Wenyu Lei
% Dec 29, 2015

vec1=x(4,:)-x(1,:);
vec2=x(3,:)-x(1,:);
vec3=x(2,:)-x(1,:);
vol = abs(cross(vec1,vec2)*vec3')/6.;

if (deg == 1)
    nq=1;
    q=zeros(nq,3);
    w=zeros(1,nq);
    q(1,:)=(x(4,:)+x(2,:)+x(3,:)+x(1,:))*0.25; w(1)=vol;    
    return;
end

if (deg == 2)
    nq=4;
    q=zeros(nq,3);
    w=zeros(1,nq);
    a = (5-sqrt(5))/20;
    q(1,:)=a*(x(4,:)+x(2,:)+x(3,:))+(1-3*a)*x(1,:); w(1)=vol/4;
    q(2,:)=a*(x(1,:)+x(4,:)+x(3,:))+(1-3*a)*x(2,:); w(2)=vol/4;
    q(3,:)=a*(x(1,:)+x(2,:)+x(4,:))+(1-3*a)*x(3,:); w(3)=vol/4;
    q(4,:)=a*(x(1,:)+x(2,:)+x(3,:))+(1-3*a)*x(4,:); w(4)=vol/4;
    return;
end

if (deg == 3)
    nq=5;
    q=zeros(nq,3);
    w=zeros(1,nq);
    q(1,:)=(x(4,:)+x(2,:)+x(3,:))/6.+0.5*x(1,:); w(1)=vol*0.45;
    q(2,:)=(x(1,:)+x(4,:)+x(3,:))/6.+0.5*x(2,:); w(2)=vol*0.45;
    q(3,:)=(x(1,:)+x(2,:)+x(4,:))/6.+0.5*x(3,:); w(3)=vol*0.45;
    q(4,:)=(x(1,:)+x(2,:)+x(3,:))/6.+0.5*x(4,:); w(4)=vol*0.45;
    q(5,:)=(x(4,:)+x(2,:)+x(3,:)+x(1,:))*0.25; w(5)=-vol*0.8;    
    return;
end

end

