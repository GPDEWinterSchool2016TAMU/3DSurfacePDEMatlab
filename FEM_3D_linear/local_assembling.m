function [ lstiff,lrhs ] = local_assembling( v, quad_degree, alpha,beta,rhs_flag )
% Local assembling routine (linear element)
% for a cell intersects with the narrow band with size h
%
% input:
%   v: vertex list (4 by 4 matrix)
%   quad_degree: degree of quadrature to apply on each tetrahedra
%   alpha: coefficient of the higher order term
%   beta: coefficient of the lower order term
%   rhs_flag: right hand side computation flag
% output:
%   lstiff: local stiffness matrix (4 by 4 for linear case)
%   lrhs: local right hand side (1 by 4 for linear case)
%
% Wenyu Lei
% Jan 7, 2016

if (nargin < 4)
    quad_degree = 2;
end

% Initialization
lstiff=zeros(4,4);
lrhs=zeros(1,4);


    
    % Get the quadarature info for each subdivision.
    [nq,q,w]=volquad(v, quad_degree);
    
    % preallocate rhs values at q
    % moved function_extension to inside rhs_eval()
    rhs_vals_at_q = rhs_eval(q); % [4x1]
        
    % Add the increment for local stiffness and right hand side on 
    % each quadrature point.
    for q_point = 1:nq
        % evaluate the basis functions and gradients at once per q_point
        [shape_values_at_q_point, shape_grads_at_q_point] ...
                 = eval_basis_value_grad(q(q_point,:),1:4,v); % [4x1], [4x3]
         
        % faster evaluation of local matrix contributions:
        % inner product the basis values and gradients all at once.
        shape_val_jk_matrix_at_q_point = shape_values_at_q_point*shape_values_at_q_point'; % [4x1][1x4]
        shape_grad_jk_matrix_at_q_point = shape_grads_at_q_point*shape_grads_at_q_point'; % [4x3][3x4]
        
        %
        % S = \int (gradphi_j . gradphi_k + phi_j*phi_k) dx
        %
        lstiff = lstiff + ( alpha*shape_grad_jk_matrix_at_q_point ...
                            + beta*shape_val_jk_matrix_at_q_point    ) ...
                          * w(q_point);

        %
        % F = \int f*phi_j dx
        %      
        if (rhs_flag)
            lrhs = lrhs + rhs_vals_at_q(q_point) ...
                        * shape_values_at_q_point' ... % [1x4]
                        * w(q_point);
        end
%         %  slower calculation of local matrix, but easier to understand
%         for j = 1:4
%             for k = 1:4
%                 %
%                 % S = \int (gradphi_j . gradphi_k + phi_j*phi_k) dx
%                 %
%                 lstiff(j,k)=lstiff(j,k) ...
%                      + ( shape_grads_at_q_point(j,:)*shape_grads_at_q_point(k,:)'...
%                          +shape_values_at_q_point(j)*shape_values_at_q_point(k)...
%                        ) ...
%                        * w(q_point);
%             end
%
%             %
%             % F = \int f*phi_j dx
%             %
%             lrhs(j)=lrhs(j) + rhs_vals_at_q(q_point)...
%                             * shape_values_at_q_point(j)...
%                             * w(q_point);
%         end

    end

end

