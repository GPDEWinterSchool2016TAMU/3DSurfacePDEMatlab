function [error_L2_Gamma_h] = level_set_torus_function(n, h, quad_degree_bulk, quad_degree_face, plot_solution)

% mesh size

% Mesh generation
% [n_nb_node,n_nb_ele,nb_node,nb_ele]=bulk_mesh_generator(h); % unstructure
[n_nb_node,n_nb_ele,nb_node,nb_ele]=uniform_bulk_mesh(n,h); % structured

% Matrix and vector initialization

% We estimate a maximum of 15 interactions on each row for the uniform mesh
% case with an average of 13 interactions over all rows .
%
% For the tetgen case, a maximum is not calculable, but for h values tested in 
% [0.0625,1] we have on average 13 interactions. It is probably safe again
% to use bound 15 but you should be aware that running out of memory
% will slow you down alot.
A = spalloc(n_nb_node, n_nb_node, 15*n_nb_node);
% A = sparse([],[],[],n_nb_node,n_nb_node,15*n_nb_node);
F = zeros(n_nb_node,1);

% Assembling
for cell = 1:n_nb_ele
    
    cell_node_ind = nb_ele(cell,1:4);     % [1x4]
    vertices = nb_node(cell_node_ind, :); % [4x3]
    dist_value = distfunc(vertices);      % [1x4]

    [ lstiff,lrhs ] = local_assembling( vertices,dist_value,h,quad_degree_bulk );

    A(cell_node_ind,cell_node_ind) = A(cell_node_ind,cell_node_ind) + lstiff'; %[4x4]
    F(cell_node_ind) = F(cell_node_ind) + lrhs';  %[4x1]

end

% Apply back slash solver
sol=A\F;


% compute error of solution
error_L2_Gamma_h = compute_error_L2_Gamma_h(sol,nb_node,nb_ele,n_nb_ele,quad_degree_face);


% plotting solution on Gamma_h
if (plot_solution)
    plot_edges_in_black = 0;
    plot_solution_on_Gamma_h(sol,nb_node,nb_ele,n_nb_ele,plot_edges_in_black);
end

