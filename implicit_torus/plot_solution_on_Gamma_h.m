function [] = plot_solution_on_Gamma_h(sol,nb_node,nb_ele,n_nb_ele, plt_edges_in_black)
% [] = plot_solution_on_Gamma_h(sol,nb_node,nb_ele,n_nb_ele, plt_edges_in_black)
%
%  This function takes the solution vector at the nodes and plots the
%  interpolated solution along the faces of Gamma_h.  It calls the
%  subdivide_Gamma_h function which extracts the faces of a cell T
%  intersected with Gamma_h and the interpolated solution values at those
%  vertices.
%
%  The only parameter is whether or not to plot the edges of the faces with
%  dark black lines.  
%
%  Spencer Patty
%  Dec 30, 2015

if (nargin < 5)
    plt_edges_in_black = 0;
end

figure;

for cell = 1:n_nb_ele
    
    % extract the info needed to determine if cell interesects Gamma_h and
    % while we are at it, extract the solution values too.
    cell_node_ind = nb_ele(cell,1:4);           % [1x4]
    vertices = nb_node(cell_node_ind, :);       % [4x3]
    dist_value = distfunc(vertices);            % [1x4]
    sol_val_at_vertices_T = sol(cell_node_ind); % [4x1]
    
    for subface = subdivide_Gamma_h(vertices,dist_value,sol_val_at_vertices_T')
        
        face_vertices=subface{1}{1}; % [3x3] extract the face points as [fv1; fv2; fv3] ([nx3])
        solution_at_v_face = subface{1}{2}; % [1x3] extract the solution values at face points
        
        x = face_vertices(:,1);  % [3x1]
        y = face_vertices(:,2);  % [3x1]
        z = face_vertices(:,3);  % [3x1]
        c = solution_at_v_face';  % [3x1]
        if (plt_edges_in_black)
            patch(x,y,z,c, 'EdgeColor','black');
        else
            patch(x,y,z,c, 'EdgeColor','none');
        end
    end
    
end

axis equal
xlabel('X'); ylabel('Y'); zlabel('Z');

end