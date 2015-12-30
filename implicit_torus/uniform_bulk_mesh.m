function [n_nb_node,n_nb_ele,nb_node,nb_elem ] = uniform_bulk_mesh( n,h )
% Generate 3D uniform mesh for the narrow band method 
%
% input:
%   h: constraints for narrow band
%   n: number of intervals for each direction of the cube.
%   NOTE: usually fix n and let h=size_of_cube/n;
% output: 
%   n_nb_node: number of vertices of the narrow band.
%   n_nb_ele: number of elements (cells) of the narrow band.
%   nb_node: vertex list of the narrow band.
%   nb_elem: connectivity list of the narror band.
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

% Generate the box mesh and get the vertex list and 
% and the connectivity list
hh=4./n;
xx=-2:hh:2;
yy=-2:hh:2;
zz=-2:hh:2;
X=zeros((n+1)*(n+1)*(n+1),3);
for i=1:(n+1)
    for j=1:n+1
        for k=1:n+1
            ind=k+(n+1)*(j-1)+(n+1)*(n+1)*(i-1);
            X(ind,1)=xx(i);
            X(ind,2)=yy(j);
            X(ind,3)=zz(k);
        end
    end
end
dt = delaunayTriangulation(X);

elem=dt.ConnectivityList;
node=dt.Points;


n_ele = length (elem);
n_local_dofs = 4;

% Count the number of cells which intersects the narrow band
narrow_band_counter=0;
for cell = 1:n_ele
    for v= 1:n_local_dofs
        interp_dist=distfunc(node(elem(cell,v),:));
        if ((interp_dist)<(h) && (interp_dist)>(-h))
            narrow_band_counter=narrow_band_counter+1;
            break;
        end
    end
end

n_nb_ele = narrow_band_counter;
nb_elem=zeros(n_nb_ele,4);

% Get the connectivity list of the narrow band
narrow_band_counter=1;
for cell = 1:n_ele
    for v= 1:n_local_dofs
        interp_dist=distfunc(node(elem(cell,v),:));
        if ((interp_dist)<(h) && (interp_dist)>(-h))
            nb_elem(narrow_band_counter,:)=elem(cell,:);
            narrow_band_counter=narrow_band_counter+1;
            break;
        end
    end
end

% Generate the index set for the vertices of cells which
% intersect the narrow band
node_nb_flag =zeros(1,length(node));
for cell = 1:n_nb_ele
    for v = 1:4
        node_nb_flag(nb_elem(cell,v))=1;
    end
end
n_nb_node=length(nonzeros(node_nb_flag));

% Extract the relavant vertices
nb_node=zeros(n_nb_node,3);
count=1;
for i = 1:length(node)
    if(node_nb_flag(i))
        nb_node(count,:)=node(i,:);
        node_nb_flag(i)=count;
        count=count+1;
    end
end

% Update the connectivity list with extracted vertex list.
for cell = 1:n_nb_ele
    for v=1:4
        nb_elem(cell,v)=node_nb_flag(nb_elem(cell,v));
    end
end

% Plot the narrow band mesh for part x>0. 
% Here the function 'plotmesh' is from ISO2MESH package.
plotmesh(nb_node,nb_elem,'x>0'); axis equal;

end