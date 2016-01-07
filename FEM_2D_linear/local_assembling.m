function [ local_stiff,local_rhs ] = local_assembling( v,hat_phi,hat_phix,hat_phiy,hatx,haty,nq,q_weights,alpha,beta,rhs_flag )


    %init
    local_stiff = zeros (3,3);
    local_rhs = zeros (3,1);
    
    %affine map info
    mat_B = [v(2,:)-v(1,:); v(3,:)-v(1,:)]';
    det_B = abs(det(mat_B));
    inv_B = [mat_B(2,2),-mat_B(1,2);-mat_B(2,1),mat_B(1,1)]/det_B;
    
    % precompute various things at the quadpoints.
    q_points = [hatx; haty]'*mat_B' + repmat(v(1,:),nq,1);
    rhs_vals = rhs_eval(q_points);
    
   
    %local stiff and rhs
    for q_index = 1:nq

        grad_phi_star_at_q_point=inv_B'*[hat_phix(q_index,:); hat_phiy(q_index,:)]; %grad_phistar
        
        grad_phi_ij_matrix = grad_phi_star_at_q_point'*grad_phi_star_at_q_point;
        phi_ij_matrix = hat_phi(q_index,:)'*hat_phi(q_index,:);
        
        local_stiff=local_stiff + (alpha * grad_phi_ij_matrix...
                                + beta * phi_ij_matrix )...
                                * q_weights(q_index)...
                                * det_B;              
        if(rhs_flag)
            local_rhs=local_rhs + rhs_vals(q_index)...    % f(q_point)
                                  *hat_phi(q_index,:)'... % basis on reference element [3x1]
                                  *q_weights(q_index)...  % quadrature weight
                                  *det_B;                 % reference area-element transform
        end
        % The following shows the local assembling in the for loop. We also
        % use the original formula (which can be simplified).
%         for j = 1:3
%             for k = 1:3
%                 local_stiff(j,k)= local_stiff(j,k)...
%                     +(beta*hat_phi(q_index,k)*hat_phi(q_index,j)...
%                     +alpha*dot(...
%                     [hat_phix(q_index,k), hat_phiy(q_index,k)]*inv_B,...
%                     [hat_phix(q_index,j), hat_phiy(q_index,j)]*inv_B)...
%                     *q_weights(q_index)*det_B;
%             end
%             if(rhs_flag)
%             local_rhs(j) = local_rhs(j) +...
%                 (rhs_eval(q_point)*hat_phi(q_index,j))*q_weights(q_index)*det_B;
%             end
%         end
    end
    

end