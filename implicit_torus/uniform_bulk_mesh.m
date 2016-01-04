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


n_ele = length (elem);
n_local_dofs = 4;

% % Count the number of cells which intersects the narrow band
% narrow_band_counter=0;
% for cell = 1:n_ele
%     if (sum(abs(distfunc(node(elem(cell,:),:))) < h) > 0) % check if any of the nodes are |dist|< h
%         narrow_band_counter=narrow_band_counter+1;
%     end
% end
% 
% n_nb_ele = narrow_band_counter;
% nb_elem=zeros(n_nb_ele,4);
% 
% % Get the connectivity list of the narrow band
% narrow_band_counter=0;
% for cell = 1:n_ele
%     if (sum(abs(distfunc(node(elem(cell,:),:))) < h) > 0)
%         narrow_band_counter=narrow_band_counter+1;
%         nb_elem(narrow_band_counter,:)=elem(cell,:);
%     end
% end


% count and store the narrow band connectivity list at the same time
% this takes more memory but 1/4 as much time as the above separated
% system.  In the end, we have the same memory required but we start the
% narrow band with as much memory as the full set.
%
% The bottleneck of uniform_bulk_mesh is the identifying which cells touch
% the band, mainly evaluating the distance function so we vectorize this
% which leaves a larger memory footprint again but speeds it up a lot
dist_vals = distfunc(node); % vectorize distance function evals
nb_elem=zeros(length(node),4);% start with full set, then cut down
narrow_band_counter = 0;
for cell = 1:n_ele
    if (sum( abs( dist_vals(elem(cell,:)) )<h ) >0) %identify if cell touches band (this is the bottleneck still)
        narrow_band_counter=narrow_band_counter+1;
        nb_elem(narrow_band_counter,:)=elem(cell,:);
    end
end
n_nb_ele = narrow_band_counter;
nb_elem = nb_elem(1:n_nb_ele,:); % cut it down to proper size


% Generate the index set for the vertices of cells which
% intersect the narrow band
node_nb_flag =zeros(1,length(node));
for cell = 1:n_nb_ele
    for v = 1:n_local_dofs
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
    for v=1:n_local_dofs
        nb_elem(cell,v)=node_nb_flag(nb_elem(cell,v));
    end
end

% Plot the narrow band mesh for part x>0. 
% Here the function 'plotmesh' is from ISO2MESH package.
plotmesh(nb_node,nb_elem,'x>0'); axis equal;

end