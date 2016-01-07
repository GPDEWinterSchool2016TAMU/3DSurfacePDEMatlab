function [  n_node,n_ele,node,ele] = triangulation_square( n )
h=1/n;
xx=0:h:1;
yy=0:h:1;
X=zeros((n+1)*(n+1),2);
[YY,XX]=meshgrid(yy,xx);
X(:,1)=reshape(XX,(n+1)^2,1);
X(:,2)=reshape(YY,(n+1)^2,1);
dt = delaunayTriangulation(X);

% Extract node and connectivity list of the square mesh.
ele=dt.ConnectivityList;
n_ele=length(ele);
node=dt.Points;
n_node=length(node);  % number of vertiecs on the surface.




end

