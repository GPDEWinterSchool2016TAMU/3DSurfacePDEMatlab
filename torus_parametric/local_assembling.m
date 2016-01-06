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
    
    %local stiff and rhs
    for q_index = 1:nq
%        q_point=AFFINE_MAPPING(v1,v2,v3,[hatx(q_index),haty(q_index)]);
        q_point=[hatx(q_index),haty(q_index)]*mat_B'+v(1,:);
        grad_chi=grad_pm(q_point);
        G=grad_chi'*grad_chi;
        q=sqrt(det(G));
        invBtgradphi_at_q=inv_B'*[hat_phix(q_index,:); hat_phiy(q_index,:)];
        local_stiff=local_stiff+ (alpha*invBtgradphi_at_q'*(G\invBtgradphi_at_q)...
                                  +beta*hat_phi(q_index,:)'*hat_phi(q_index,:))...
                                  *q_weights(q_index)*det_B*q;
        if(rhs_flag)
            local_rhs=local_rhs+rhs_eval(q_point)*hat_phi(q_index,:)'...
                       *q_weights(q_index)*det_B*q;
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