function [nq,q,w ] = facequad( x,deg )
% Quadrature rule on a triangle.
% input: 
%   x: 3 by 3 matrix. Each row indicates the coordinates of the vertiex
%      of the triangle.
%   deg: degree of the quadrature rule.
% output:
%   nq: number of quadrature points.
%   q: quadrature point list.
%   w: weights for cooresponding quadrature points.
% ref: Theory and practice of finite elements 
%      By A. Ern and JL Guermond
%      Page 360
%
% Wenyu Lei
% Dec 29, 2015

vol=norm(cross(x(3,:)-x(1,:),x(3,:)-x(2,:)))/2.;

if (deg==1)
    nq=1;
    q=zeros(nq,3);
    w=zeros(1,nq);
    q(1,:)=x(1,:)/3.+x(2,:)/3.+x(3,:)/3.; w(1)=vol;
    return;
end

if (deg==2)
    nq=3;
    q=zeros(nq,3);
    w=zeros(1,nq);
    q(1,:)=x(1,:)/6.+x(2,:)/6.+2*x(3,:)/3.; w(1)=vol/3.;
    q(2,:)=x(1,:)/6.+x(3,:)/6.+2*x(2,:)/3.; w(2)=vol/3.;
    q(3,:)=x(3,:)/6.+x(2,:)/6.+2*x(1,:)/3.; w(3)=vol/3.;
    return;
end

if (deg==3)
    nq=4;
    q=zeros(nq,3);
    w=zeros(1,nq);
    q(1,:)=x(1,:)*0.2+x(2,:)*0.2+x(3,:)*0.6; w(1)=25*vol/48;
    q(2,:)=x(3,:)*0.2+x(1,:)*0.2+x(2,:)*0.6; w(2)=25*vol/48;
    q(3,:)=x(3,:)*0.2+x(2,:)*0.2+x(1,:)*0.6; w(3)=25*vol/48;
    q(4,:)=x(1,:)/3.+x(2,:)/3.+x(3,:)/3.; w(4)=-9*vol/16;
    return;
end

end

