function [n_node,n_ele,node,ele ] = box_triangulation( n )
% triangulation on the unit cube
hh=1/n;
xx=0:hh:1;
yy=0:hh:1;
zz=0:hh:1;
X=zeros((n+1)*(n+1)*(n+1),3);

[YY,ZZ,XX] = meshgrid(yy,zz,xx);  % order this way to match for loop ordering
X(:,1) = reshape(XX,(n+1)^3,1);
X(:,2) = reshape(YY,(n+1)^3,1);
X(:,3) = reshape(ZZ,(n+1)^3,1);

dt = delaunayTriangulation(X);

ele=dt.ConnectivityList;
node=dt.Points;
n_ele=length(ele);
n_node=length(node);


end

