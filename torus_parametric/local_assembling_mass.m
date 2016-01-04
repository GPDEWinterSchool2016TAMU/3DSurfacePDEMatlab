function [ local_stiff] = local_assembling_mass( v1,v2,v3,hat_phi,nq,q_weights )
%LOCAL_ASSEMBLING Summary of this function goes here
%   Detailed explanation goes here

    local_stiff = zeros (3,3);
 
    mat_B = [v2(1)-v1(1), v3(1)-v1(1); v2(2)-v1(2),v3(2)-v1(2)];
    det_B = abs(det(mat_B));
%    inv_B = [mat_B(2,2),-mat_B(1,2);-mat_B(2,1),mat_B(1,1)]/det_B;
    for q_index = 1:nq
        for j = 1:3
            for k = 1:3
                local_stiff(j,k)= local_stiff(j,k)...
                    +(hat_phi(q_index,k)*hat_phi(q_index,j)...
                    )*q_weights(q_index)*det_B;
            end
        end
    end
    

end

