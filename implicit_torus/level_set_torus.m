clear all; clc;
format long;

% mesh size
n=32; h=4./n;  

% quadrature degree of accuracy
quad_degree_bulk = 2; % for tetrahedra
quad_degree_face = 2; % for faces


% Mesh generation
fprintf('...Generating Mesh...\n');
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
fprintf('...Assembling Stiffness Matrix and RHS...\n');
for cell = 1:n_nb_ele
    
    cell_node_ind = nb_ele(cell,1:4);     % [1x4]
    vertices = nb_node(cell_node_ind, :); % [4x3]
    dist_value = distfunc(vertices);      % [1x4]

    [ lstiff,lrhs ] = local_assembling( vertices,dist_value,h,quad_degree_bulk );

    A(cell_node_ind,cell_node_ind) = A(cell_node_ind,cell_node_ind) + lstiff'; %[4x4]
    F(cell_node_ind) = F(cell_node_ind) + lrhs';  %[4x1]

end

% Apply back slash solver
fprintf('...Solving System...\n')
sol=A\F;

% Compute the L2 error i.e.
% \|u^e - u_h\|_{L^2(\Gamma_h)}
% Here \Gamma_h is the discrete surface defined by the level set function
% \Gamma_h := {x\in box mesh | I_h\Phi(x)=0}
% and u^e is extension of the solution.
% We compute this in a similar way as we did in the local assmebling.
% That is, we compute (u^e - u_h)^2 on the intersection between a cell and 
% \Gamma_h. The intersection may triangle or quadrilateral. So we cut
% the quadrilateral by two triangle and apply quadrature rule only on 
% triangles.
fprintf('...Computing L2(Gamma_h) Error...\n')
err_square = 0;
for cell = 1:n_nb_ele
    
    % extract the info needed to determine if cell interesects Gamma_h and
    % while we are at it, extract the solution values too.
    cell_node_ind = nb_ele(cell,1:4);           % [1x4]
    vertices = nb_node(cell_node_ind, :);       % [4x3]
    dist_value = distfunc(vertices);            % [1x4]
    sol_val_at_vertices_T = sol(cell_node_ind); % [4x1]
    
    for subface = subdivide_error(vertices,dist_value)
        % get quadrature info
        v_face=subface{1}; % extract the face points as [fv1; fv2; fv3] ([nx3])
        [nq,q,w]=facequad(v_face, quad_degree_face);
        exact_at_q = exact(q); % moved function_extension to inside exact() 
        for q_point = 1:nq
            shape_values_at_q_point = eval_basis_value(q(q_point,:),1:4,vertices); % [4x1]
            solution_at_q_point = shape_values_at_q_point'*sol_val_at_vertices_T; % [1x4][4x1]
            err_square=err_square + (exact_at_q(q_point)-solution_at_q_point)^2.*w(q_point);
        end
    end
end
% end
err=sqrt(err_square)