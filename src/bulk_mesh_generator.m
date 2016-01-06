function [n_node,n_elem,node,elem]=bulk_mesh_generator(h, return_nb)
% Generate 3D mesh for the narrow band method from ISO2MESH
%
% input:
%   h: constraints for the maximum volume of cells in the 
%   box mesh which contains the torus.
%   return_nb: whether or not to extract narrow band and return instead of
%              full mesh
% output: 
%   n_node: number of vertices of the narrow band or full mesh.
%   n_ele: number of elements (cells) of the narrow band or full mesh.
%   node: vertex list of the narrow band or full mesh.
%   elem: connectivity list of the narror band or full mesh.
%
% We first generate a mesh for a box which contains the torus
% Here we use the function 'meshabox' from ISO2MESH package.
% The input argument h indicates the maximum volume of cells
% in the box mesh (denote the box domain as U) is less than
% h^3. Then select the cells intersects the narrow band
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

% Generate the box mesh which contains the torus.The constraints 
% for the maximum volume is h^3.
% You need to modify the first two arguments when you have 
% different surface.

if (nargin < 2)
    return_nb = 1;
end

[node,~,elem]=meshabox([-2 -2 -1],[2 2 1],h*h*h);


if (return_nb)
    [n_node,n_elem,node,elem] = extract_narrow_band(node,elem);
else % return full mesh
    n_node = size(node,1);
    n_elem = size(elem, 1);
end


% Plot the narrow band mesh for part x>0. 
% Here the function 'plotmesh' is from ISO2MESH package.
plotmesh(node,elem,'x>0'); axis equal;

end