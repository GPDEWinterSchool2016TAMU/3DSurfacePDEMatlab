function [n_nb_node,n_nb_elem,nb_node,nb_elem] = extract_narrow_band(node,elem,h)
% [n_nb_node,n_nb_elem,nb_node,nb_elem] = extract_narrow_band(node,elem)
%
% Extracts the narrow band tetrahedral elements from bulk mesh and
% reorders them so that that is the new 

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
% n_nb_elem = narrow_band_counter;
% nb_elem=zeros(n_nb_elem,4);
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
nb_elem=zeros(size(node,1),4);% start with full set, then cut down
narrow_band_counter = 0;
for cell = 1:n_ele
    if (sum( abs( dist_vals(elem(cell,:)) )<h ) >0) %identify if cell touches band (this is the bottleneck still)
        narrow_band_counter=narrow_band_counter+1;
        nb_elem(narrow_band_counter,:)=elem(cell,:);
    end
end
n_nb_elem = narrow_band_counter;
nb_elem = nb_elem(1:n_nb_elem,:); % cut it down to proper size


% Generate the index set for the vertices of cells which
% intersect the narrow band
node_nb_flag =zeros(1,length(node));
for cell = 1:n_nb_elem
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
for cell = 1:n_nb_elem
    for v=1:n_local_dofs
        nb_elem(cell,v)=node_nb_flag(nb_elem(cell,v));
    end
end

