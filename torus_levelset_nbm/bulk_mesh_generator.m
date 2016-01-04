function [n_nb_node,n_nb_ele,nb_node,nb_elem]=bulk_mesh_generator(h)
% Generate 3D mesh for the narrow band method from ISO2MESH
%
% input:
%   h: constraints for the maximum volume of cells in the 
%   box mesh which contains the torus.
% output: 
%   n_nb_node: number of vertices of the narrow band.
%   n_nb_ele: number of elements (cells) of the narrow band.
%   nb_node: vertex list of the narrow band.
%   nb_elem: connectivity list of the narror band.
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
[node,~,elem]=meshabox([-2 -2 -1],[2 2 1],h*h*h);
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