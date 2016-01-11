function [ local_stiff,local_rhs ] = local_assembling( v,hat_phi,...
                                                       grad_hat_phi_x1,grad_hat_phi_x2,...
                                                       q_yhat,nq,q_weights,...
                                                       alpha,beta,rhs_flag)

% [ local_stiff,local_rhs ] = local_assembling( v,hat_phi,grad_hat_phi_x1,grad_hat_phi_x2, q_yhat,nq,q_weights, alpha,beta,rhs_flag)
%
% construct local contributions from cell (defiend by vertices v to 
% stiffness matrix and rhs vector corresponding to
% strong form equation
%
%  -div( a grad u) + bu = f
%
%  with du/dn = 0 on boundary.
%

    %init
    local_stiff = zeros (3,3);
    local_rhs = zeros (3,1);
    
    %affine map info
    mat_B = [v(2,:)-v(1,:); v(3,:)-v(1,:)]';
    det_B = abs(det(mat_B));
    inv_B = [mat_B(2,2),-mat_B(1,2);-mat_B(2,1),mat_B(1,1)]/det_B;
    
    % precompute various things at the quadrature points.
    q_points = q_yhat*mat_B' + repmat(v(1,:),nq,1); % [nqx2] list of quadrature points in cell 
    rhs_vals = rhs_eval(q_points);                  % [nqx1] f(q_points)
                    
    
    %local stiff and rhs
    for q_index = 1:nq

        grad_phi_at_q_point=inv_B'*[grad_hat_phi_x1(q_index,:); grad_hat_phi_x2(q_index,:)]; % [2x3] 3 basis func gradients
  
        % exploiting matlab's speed with vector operations, we avoid doing
        % the for loops to construct the local matrix and rhs.  Below in
        % the commented section is the same code using for loops.
        grad_phi_ij_matrix = grad_phi_at_q_point'*grad_phi_at_q_point;  % = [grad phi_i . grad_phi_j ]_{ij}
        phi_ij_matrix = hat_phi(q_index,:)'*hat_phi(q_index,:);         % = [phi_i phi_j]_{ij}
        
        local_stiff=local_stiff + (alpha * grad_phi_ij_matrix' ...
                                   + beta * phi_ij_matrix' )...
                                * q_weights(q_index)...
                                * det_B;   
                            
        if(rhs_flag)
            local_rhs=local_rhs + rhs_vals(q_index)...    % f(q_point)
                                  *hat_phi(q_index,:)'... % basis on reference element at q_point [3x1]
                                  *q_weights(q_index)...  % quadrature weight
                                  *det_B;                 % reference area-element transform
        end
        
        
%         % The following shows the local assembling with the for loop.
%         for j = 1:3
%             for k = 1:3
%                 local_stiff(j,k)= local_stiff(j,k)...
%                     +(beta*hat_phi(q_index,k)*hat_phi(q_index,j)...
%                       +alpha*dot(grad_phi_at_q_point(:,k),grad_phi_at_q_point(:,j)))...
%                     *q_weights(q_index)*det_B;
%             end
%             if(rhs_flag)
%                 local_rhs(j) = local_rhs(j) +...
%                     (rhs_vals(q_index)*hat_phi(q_index,j))*q_weights(q_index)*det_B;
%             end
%         end
    end
    

end