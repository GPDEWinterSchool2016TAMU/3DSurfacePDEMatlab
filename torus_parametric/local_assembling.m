function [ local_stiff,local_rhs ] = local_assembling( v1,v2,v3,hat_phi,hat_phix,hat_phiy,hatx,haty,nq,q_weights,alpha,beta )
% Local assembling routine on the cell \tao
%
% \int_{\tao} alpha*\nabla_\Gamma\phi_j\cdot\nabla_\Gamma\phi_k
%            +beta*\phi_j*\phi_k dx
% where \nabla_\Gamma is the surface gradient.
% In our case, the local stiffness is case alpha=beta=1.
% We also need to compute the mass matrix for error computation, i.e.
% the case when alpha=0 and beta=1.
% In the input arguemnts, v1, v2, v3 are vertices of the cell, and
% hat_phi, hat_phix, hat_phiy, nq, weights are finite element and 
% quadrature information.
%
%   Wenyu Lei
%   Dec 30, 2015

    %init
    local_stiff = zeros (3,3);
    local_rhs = zeros (3,1);
    
    %affine map info
    mat_B = [v2(1)-v1(1), v3(1)-v1(1); v2(2)-v1(2),v3(2)-v1(2)];
    det_B = abs(det(mat_B));
    inv_B = [mat_B(2,2),-mat_B(1,2);-mat_B(2,1),mat_B(1,1)]/det_B;
    
    %local stiff and rhs
    for q_index = 1:nq
        q_point=AFFINE_MAPPING(v1,v2,v3,[hatx(q_index),haty(q_index)]);
        g_pm=grad_pm(q_point);
        G=g_pm'*g_pm;
        q=sqrt(det(G));
        invG=[G(2,2),-G(1,2);-G(2,1),G(1,1)]/q^2;
        for j = 1:3
            for k = 1:3
                local_stiff(j,k)= local_stiff(j,k)...
                    +(beta*hat_phi(q_index,k)*hat_phi(q_index,j)...
                    +alpha*dot(...
                    [hat_phix(q_index,k), hat_phiy(q_index,k)]*inv_B*invG*g_pm',...
                    [hat_phix(q_index,j), hat_phiy(q_index,j)]*inv_B*invG*g_pm'))...
                    *q_weights(q_index)*det_B*q;
            end
            if(alpha~=0)
            local_rhs(j) = local_rhs(j) +...
                (rhs_eval(q_point)*hat_phi(q_index,j))*q_weights(q_index)*det_B*q;
            end
        end
    end
    

end