function [ local_stiff,local_rhs ] = local_assembling( v,hat_phi,hat_phix,hat_phiy,hatx,haty,nq,q_weights,alpha,beta,rhs_flag )
% Local assembling routine on the cell \tao
%
% \int_{\tao} alpha*\nabla_\Gamma\phi_j\cdot\nabla_\Gamma\phi_k
%            +beta*\phi_j*\phi_k dx
% where \nabla_\Gamma is the surface gradient.
% The above formula could be simplified as (3.3.3) in the lecture notes.

% In our case, the local stiffness matrix is the case alpha=beta=1.
% We also need to compute the mass matrix for error computation, i.e.
% the case when alpha=0 and beta=1.
% In the input arguemnts, v1, v2, v3 are vertices of the cell, and
% hat_phi, hat_phix, hat_phiy, nq, weights are finite element and 
% quadrature information.
% The last argument indicates whether we compute the right hand side or
% not.

%   Wenyu Lei
%   Dec 30, 2015

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
    
    % 
    % X(yhat) = I_h(chi)(yhat) = chi_1 + s(chi_2 - chi_1) + t(chi_3 - chi_1)  [1x3]
    % 
    % grad_Ih_chi = grad X(yhat) = [chi_2-chi_1; chi_3- chi_1]'*inv(B)  [3x2]
    %
    % G = grad I_h chi * grad I_h chi  % constant matrix on each cell
    % Q = sqrt(det(Gint) )             % constant on each cell
    %
    x = parameterization(v);
    grad_Ih_chi = [ x(2,:)-x(1,:); x(3,:)-x(1,:)]'*inv_B;
    G=grad_Ih_chi'*grad_Ih_chi;
    Q=sqrt(det(G));

%     grad_chi = grad_pm(v);
%     G=grad_chi'*grad_chi;
%     Q=sqrt(det(G)); % chi area element
    
    %local stiff and rhs
    for q_index = 1:nq
%        q_point=AFFINE_MAPPING(v1,v2,v3,[hatx(q_index),haty(q_index)]);

        %
        % int_{chi(T)} (alpha*surfgrad_phi_i . surfgrad_phi_j + beta*phi_i*phi_j ) d\x
        %
        % note:  surfgrad_phi(x) = grad_phistar(y) * inv(G)(y) * grad_chi(y)
        %
        % G = grad_chi' grad_chi  ([2x3][3x2])  first fundamental tensor
        % Q = sqrt(det(G))                      chi map area-element transform
        %
        % so \int_{\chi(T)} gradphi_i(x) gradphi_j(x) d\x 
        %       = ... = \int_{T} grad_phistar_i(y)' * inv(G)(y) * grad_phistar_j(y) * Q dy
        %
        %   now mapping to reference element, we have
        %      grad_phistar(y) = inv(B)*grad_phistar_hat(yhat)
        %
        %       = \int_{hat(T)} (inv(B)*grad_phistarhat_i(yhat))' * inv(G)(y) 
        %                             * (inv(B)*grad_phistar_j(yhat)) * Q * det(B) dyhat
        %
        %  
        %
        
%         grad_chi = grad_pm(q_points(q_index,:));
%         G=grad_chi'*grad_chi;
%         Q=sqrt(det(G)); % chi area element
        grad_phi_star_at_q_point=inv_B'*[hat_phix(q_index,:); hat_phiy(q_index,:)]; %grad_phistar
        
        grad_phi_ij_matrix = grad_phi_star_at_q_point'*(G\grad_phi_star_at_q_point);
        phi_ij_matrix = hat_phi(q_index,:)'*hat_phi(q_index,:);
        
        local_stiff=local_stiff + (alpha * grad_phi_ij_matrix...
                                   + beta * phi_ij_matrix )...
                                  * q_weights(q_index)...
                                  * det_B...
                                  * Q;                 
        %
        % \int_{chi(T)}  f(x) phi_j(x) d\x
        %     = \int_{T}  f(x) phistar_j(y) Q(y) dy
        %         = \int_{hat(T)} f(x) phistar_hat_j(yhat) Q(y) det(B) dyhat
        %
        if(rhs_flag)
            local_rhs=local_rhs + rhs_vals(q_index)...    % f(q_point)
                                  *hat_phi(q_index,:)'... % basis on reference element [3x1]
                                  *q_weights(q_index)...  % quadrature weight
                                  *det_B...               % reference area-element transform
                                  *Q;                     % chi area-element transform
        end
        % The following shows the local assembling in the for loop. We also
        % use the original formula (which can be simplified).
%         for j = 1:3
%             for k = 1:3
%                 local_stiff(j,k)= local_stiff(j,k)...
%                     +(beta*hat_phi(q_index,k)*hat_phi(q_index,j)...
%                     +alpha*dot(...
%                     [hat_phix(q_index,k), hat_phiy(q_index,k)]*inv_B*invG*grad_chi',...
%                     [hat_phix(q_index,j), hat_phiy(q_index,j)]*inv_B*invG*grad_chi'))...
%                     *q_weights(q_index)*det_B*q;
%             end
%             if(rhs_flag)
%             local_rhs(j) = local_rhs(j) +...
%                 (rhs_eval(q_point)*hat_phi(q_index,j))*q_weights(q_index)*det_B*q;
%             end
%         end
    end
    

end