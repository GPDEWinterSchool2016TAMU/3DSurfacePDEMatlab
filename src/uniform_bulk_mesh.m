function [n_node,n_elem,node,elem ] = uniform_bulk_mesh( n,h ,return_nb)
% Generate 3D uniform mesh for the narrow band method 
%
% input:
%   h: constraints for narrow band
%   n: number of intervals for each direction of the cube.
%   return_nb: whether or not to extract narrow band and return instead of
%              full mesh
%   NOTE: usually fix n and let h=size_of_cube/n;
% output: 
%   n_node: number of vertices of the narrow band or full mesh.
%   n_ele: number of elements (cells) of the narrow band or full mesh.
%   node: vertex list of the narrow band or full mesh.
%   elem: connectivity list of the narror band or full mesh.
%
% We first generate a mesh for a box which contains the torus
% Here we use the function 'delaunayTriangulation' from MATLAB.
% Then select the cells intersects the narrow band
%
% D_h := {x \in U | I_h \Phi(x) < h}
%
% Here \Phi(x) is the distance function (Or level set function)
% and I_h is the lagrange interpolation on to continuous
% piecewise linear function space.
% To do this, we evaluate distance function on each vertices
% of the cell. We select this cell if one of the values is 
% between -h and h.
%
% Wenyu Lei
% Dec 28, 2015

if (nargin < 3)
    return_nb = 1;
end

% Generate the box mesh and get the vertex list and 
% and the connectivity list
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

elem=dt.ConnectivityList;
node=dt.Points;

if (return_nb)
    [n_node,n_elem,node,elem] = extract_narrow_band(node,elem,h);
else % return full mesh
    n_node = size(node,1);
    n_elem = size(elem, 1);
end

% Plot the narrow band mesh for part x>0. 
% Here the function 'plotmesh' is from ISO2MESH package.
plotmesh(node,elem,'x>0'); axis equal;

end