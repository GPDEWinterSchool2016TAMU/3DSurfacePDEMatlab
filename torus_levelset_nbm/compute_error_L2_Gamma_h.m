function [err] = compute_error_L2_Gamma_h(sol,nb_node,nb_ele,n_nb_ele,quad_degree_face)
% Compute the L2 error i.e. \|u^e - u_h\|_{L^2(\Gamma_h)}
% Here \Gamma_h is the discrete surface defined by the level set function
% \Gamma_h := {x\in box mesh | I_h d(\x)=0}
% and u^e is extension of the solution.
% We compute this in a similar way as we did in the local assembling.
% That is, we compute (u^e - u_h)^2 on the intersection between a cell and 
% \Gamma_h. The intersection may triangle or quadrilateral. So we cut
% the quadrilateral by two triangle and apply quadrature rule only on 
% triangles using the subdivid_error() function.
%
%

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

err=sqrt(err_square);

end