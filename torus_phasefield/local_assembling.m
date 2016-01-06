function [ local_stiff,local_rhs,me_sum ] = local_assembling( v,hat_phi,hat_phix,hat_phiy,hat_phiz,hatx,haty,hatz,nq,q_weights,alpha,beta,epsilon,rhs_flag )
%LOCAL_ASSEMBLING 
    
    %init
    local_stiff = zeros (4,4);
    local_rhs = zeros (4,1);
    

    
    %affine map info
    mat_B =[v(2,:)-v(1,:); v(3,:)-v(1,:); v(4,:)-v(1,:)]';
    det_B = abs(det(mat_B));
    inv_B = inv(mat_B);
    
    %local stiff and rhs
    me_sum=0;
    for q_index = 1:nq
        q_point=[hatx(q_index),haty(q_index),hatz(q_index)]*mat_B'+v(1,:); 
        psi=phase_field_function(q_point,epsilon);
        M_epsilon=1-psi^2;
        me_sum=me_sum+M_epsilon*q_weights(q_index);
        
        invBtgradphi_at_q = inv_B'*[hat_phix(q_index,:);
                                    hat_phiy(q_index,:);
                                    hat_phiz(q_index,:)];
        local_stiff=local_stiff...
                    +(alpha*invBtgradphi_at_q'*invBtgradphi_at_q...
                    +beta*hat_phi(q_index,:)'*hat_phi(q_index,:))...
                    *q_weights(q_index)*det_B*M_epsilon;
        if(rhs_flag)
            local_rhs=local_rhs+rhs_eval(q_point)*hat_phi(q_index,:)'...
                      *q_weights(q_index)*det_B*M_epsilon;
        end

    end
    

end