function [ n_node,n_ele,node,ele ] = uniform_box_mesh( n,h )
%UNIFORM_BOX_MESH Summary of this function goes here
%   Detailed explanation goes here

hh=4/n;
xx=-2:hh:2;
yy=-2:hh:2;
zz=-2:hh:2;
X=zeros((n+1)*(n+1)*(n+1),3);
% for i=1:n+1
%     for j=1:n+1
%         for k=1:n+1
%             ind=k+(n+1)*(j-1)+(n+1)*(n+1)*(i-1);
%             X(ind,1)=xx(i);
%             X(ind,2)=yy(j);
%             X(ind,3)=zz(k);
%         end
%     end
% end

[YY,ZZ,XX] = meshgrid(yy,zz,xx);  % order this way to match for loop ordering
X(:,1) = reshape(XX,(n+1)^3,1);
X(:,2) = reshape(YY,(n+1)^3,1);
X(:,3) = reshape(ZZ,(n+1)^3,1);

dt = delaunayTriangulation(X);

ele=dt.ConnectivityList;
n_ele=length(ele);
node=dt.Points;
n_node=length(node);
end

