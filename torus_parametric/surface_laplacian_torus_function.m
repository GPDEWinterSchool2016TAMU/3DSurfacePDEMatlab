function [l2_error] = surface_laplacian_torus_function(n,plot_solutions)
% [l2_error] = surface_laplacian_torus_function(n,plot_solutions)
%
% Calls the appropriate mesh, assembly and error analysis functions and
% plots the solution and error if desired.
%
%  input:
%    n = number of discretizations in each direction of parametric space.
%    plot_solution = flag for plotting or not.
%
%


[ n_node,n_ele,pm_node,ele,global_ind,global_ind_inverse] = triangulation_surface( n );

% Initialization
A = sparse([],[],[],n_node,n_node,7*n_node);
MASS = sparse([],[],[],n_node,n_node,7*n_node);
rhs = zeros(n_node,1);

% Since we are going to compute the integral on the reference element,
% we need provide the quadrature rule for the reference triangle.
% Meanwhile, we can also provide infomation of shape functions on the 
% reference element, i.e. function values (hat_phi) and function gradients
% (hat_phix and hat_phiy) on each quadrature points.

% Quadrature on reference element
nq=4;
% quadrature weights [nqx1]
q_weights= [1./24,1./24,1./24,9./24]';
% quadrature points [nq x 2]
q_yhat = [0,1,0,1./3;... % x components  
          0,0,1,1./3]';  % y components
      
% shape value and shape gradient (x and y components) each are [nqx1]
[ hat_phistar_at_q,hat_phistarx_at_q,hat_phistary_at_q ] = FEEVAL( q_yhat,nq );

% Assembling
for cell = 1:n_ele
    % Get local stiffness matrix and local rhs
    cell_ind = ele(cell,1:3);     % [1x3]
    T_vertices = pm_node(cell_ind, :); % [3x2]
    cell_global_ind=global_ind(cell_ind);
    [ local_stiff,local_rhs ] ...
        = local_assembling( T_vertices,...
                            hat_phistar_at_q, hat_phistarx_at_q, hat_phistary_at_q,...
                            q_yhat,nq,q_weights,...
                            1,1,1); % a, beta, rhs_flag
                        
    % Copy local to global
    A(cell_global_ind,cell_global_ind) ...
           = A(cell_global_ind,cell_global_ind) + local_stiff; %[3x3]
    rhs(cell_global_ind) = rhs(cell_global_ind) + local_rhs;  %[3x1]
end

% Apply back slash solver
solution = A\rhs;


% L2 error computation
% Here we are going to compute L2 norm of  
% u_h - I_h u
% Here u_h is our numerical results, i.e. the linear combinatioin of 
% basis function with coefficients from the 'solution' vector we 
% just computed.
% I_h u is the largrange interpolation of the exact solution u.
% So I_h u is also the linear combination of the basis function 
% with the coefficients from the evaluation on all vertices in the mesh.
% (we call this coefficient vecter 'exact_sol')
% So the square of the l2 norm of u_h - I_h u should be
% (exact_sol-solution)'M(exact_sol-solution)
% Here M is the mass matrix.

% We first get the coefficients of I_h u
% This can be achieved by using the inverse mapping 
% of the global indexing (global_ind_inverse) to 
% extract unique nodes from 
% parametric node list
exact_sol = exact(pm_node(global_ind_inverse,:));

err_vec =exact_sol - solution;
%Assemble mass matrix
for cell=1 : n_ele
    % Local mass matrix
    cell_ind = ele(cell,1:3);     % [1x3]
    T_vertices = pm_node(cell_ind, :); % [3x2]
    cell_global_ind=global_ind(cell_ind);
    [local_mass,~] = ...
        local_assembling( T_vertices,...
                          hat_phistar_at_q,hat_phistarx_at_q,hat_phistary_at_q,...
                          q_yhat,nq,q_weights,...
                          0,1,0); %alpha, beta, rhs_flag
    % copy local to global
    MASS(cell_global_ind,cell_global_ind)...
                =MASS(cell_global_ind,cell_global_ind) + local_mass; %[3x3]
end

% print out the error
l2_error = sqrt(transpose(err_vec)*MASS*err_vec);



if (plot_solutions)
    % Visualization
    % Using the function 'patch' to visualize each triangle.
    % Colors are decided by the value on the vertices.

    Xnodes = parameterization(pm_node(global_ind_inverse,:));

    % change the element connectivity list to use the unique nodes of
    % global_ind instead of the repeated ones of ele and save as sele.
    sele=global_ind(ele);
    figure(1);
    axis([-2,2,-2,2,-2,2]); title('Solution'); colormap('default'); colorbar;
    for i=1:n_ele
        XX=Xnodes(sele(i,:),1);
        YY=Xnodes(sele(i,:),2);
        ZZ=Xnodes(sele(i,:),3);
        CC=solution(sele(i,:),1);
        patch(XX,YY,ZZ,CC,'EdgeColor','interp');
    end

    figure(2);
    axis([-2,2,-2,2,-2,2]); title('Error'); colormap('jet'); colorbar;
    for i=1:n_ele
        XX=Xnodes(sele(i,:),1);
        YY=Xnodes(sele(i,:),2);
        ZZ=Xnodes(sele(i,:),3);
        CC=err_vec(sele(i,:),1);
        patch(XX,YY,ZZ,CC,'EdgeColor','interp');
    end
end

end