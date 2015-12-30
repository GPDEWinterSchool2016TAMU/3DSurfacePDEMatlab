clear all; clc;
format long;

% mesh size
n=32; h=4./n;  

% quadrature degree of accuracy
quad_degree_bulk = 2; % for tetrahedra
quad_degree_face = 2; % for faces


% Mesh generation
fprintf('Generating Mesh...\n');
%[n_nb_node,n_nb_ele,nb_node,nb_ele]=bulk_mesh_generator(h);
[n_nb_node,n_nb_ele,nb_node,nb_ele]=uniform_bulk_mesh(n,h);

% Matrix and vector initialization
A = spalloc(n_nb_node, n_nb_node, 9*n_nb_ele);
% A = sparse([],[],[],n_nb_node,n_nb_node,9*n_nb_ele);
F = zeros(n_nb_node,1);

% Assembling
fprintf('Assembling Stiffness Matrix and RHS...\n');
for cell = 1:n_nb_ele
    
    cell_ind = nb_ele(cell,1:4);     % 1x4
    vertices = nb_node(cell_ind, :); % 4x3
    dist_value = distfunc(vertices); % 1x4

    [ lstiff,lrhs ] = local_assembling( vertices,dist_value,h,quad_degree_bulk );

    A(cell_ind,cell_ind)=A(cell_ind,cell_ind)+lstiff';
    F(cell_ind)=F(cell_ind)+lrhs';

end

% Apply back slash solver
fprintf('Solving System...\n')
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
fprintf('Computing L2(Gamma_h) Error...\n')
err_square = 0;
for cell = 1:n_nb_ele
    % We first get the subdivision of \Gamma_h on this cell.    
    cell_ind = nb_ele(cell,1:4);     % 1x4
    vertices = nb_node(cell_ind, :); % 4x3
    dist_value = distfunc(vertices); % 1x4
    sol_val_at_vertices_T = sol(cell_ind); % (1x4)^T = 4x1
    
    for subface = subdivide_error(vertices',dist_value)
        % get quadrature info
        points_list=subface{1}'; % extract the face points
        [nq,q,w]=facequad(points_list, quad_degree_face);
        exact_at_q = exact(q); % moved function_extension to inside exact() 
        for q_point = 1:nq
            shape_value=zeros(1,4);
            for j = 1:4
                shape_value(j)=eval_basis_value(q(q_point,:),j,vertices);               
            end
            err_square=err_square...
                +(exact_at_q(q_point)-shape_value*sol_val_at_vertices_T)^2.*w(q_point);
        end
    end
end
% end
err=sqrt(err_square)