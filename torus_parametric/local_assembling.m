function [ local_stiff,local_rhs ] = local_assembling( y, hat_phistar,...
                                                       grad_hat_phistar_x1,grad_hat_phistar_x2,...
                                                       q_yhat,nq,q_weights,...
                                                       alpha,beta,rhs_flag )
% Local assembling routine on the cell X0(T) (with T a cell in parametric
% domain).  This is set to assemble stiffness and rhs for the weak form of 
%
%       -alpha surface laplacian u + beta u = f on the discrete Torus
%
% in particular, the stiffness matrix is
%
%  A_{ji} := \int_{X0(T)} alpha*\nabla_S\phi_j\cdot\nabla_S\phi_k
%                  + beta*\phi_j*\phi_k dx
%
%    = ... simplifications described below
%
%    = \int_{T}  (alpha*\nabla\phistar_j^T inv(G_Gamma) \nabla\phistar_i
%                  + beta*phistar_j*phistar_i  )
%                  * Q_Gamma dy
%
% where \nabla_S is the surface gradient on surface X0(T).
%
% The right hand side is defined on X0(T) so as to be consistent
% with the integral on chi(T) for ease of error analysis.  
% In particular,
%
%  F_j := \sum_{T} \int_{X0(T)} f(P0(x)) \phi_j(x) q(P0(x))/ Q_Gamma dx 
%       = \sum_{T} \int_{T}  fstar(y) phistar_j(y) qstar(y) dy  (the one implemented )
%       
%      (= \sum_{T} \int_{chi(T)} f(x) phi_j(invP0(x)) dx ) 
%
% See section 3.3 in the lecture notes for all of this above.
%
% In our case, the local stiffness matrix is the case alpha=beta=1. We also
% need to compute the mass matrix for error computation, i.e. the case when
% alpha=0 and beta=1.
%
%  Input:
%    y = [y1; y2; y3] [3x3] are vertices of the cell T
%    hat_phistar [nqx3]  = the 3 basis elements (columns) evaluated at nq points
%    grad_hat_phistar1 [nqx3] = the 3 basis x1 derivatives (columns) evaluated at nq points
%    grad_hat_phistar2 [nqx3] = the 3 basis x2 derivatives (columns) evaluated at nq points
%    nq = number of quadrature points
%    q_yhat = [nqx2] list of reference element quadrature points
%    q_weights = [nqx1] list of quadrature weights corresponding to q_yhat
%    alpha = surface laplacian coefficient
%    beta = mass coefficient
%    rhs_flag = compute right hand side or not
%  Output:
%    local_assembly [3x3] local contributions to stiffness matrix from cell T
%    local_rhs [3x1] local contributions to rhs from cell T
%
%  Reference Cell Notation:
%    Affine map from hat(T) to T is    F(yhat) = B*yhat + y1   
%                   with B = [y2-y1; y3-y1]' 
%
%  Basis functions:
%    hat_phistar defined on hat(T), 
%    phistar     defined on T     
%    phi         defined on X0(T)
%
%  Elements:
%    hat(T) = reference element in R^2
%      T    = parameter space element in R^2
%    chi(T) = curvilinear element on Torus in R^3
%    X0(T)  = linear interpolant of chi(T) with vertices on Torus in R^3
%
%  Drawings:
%
%  yhat3                y3                  x3       x2
%   *                   *           chi      *---------* 
%   |\        F         | \      --------->   \chi(T) /
%   | \   -------->     | T \                  \     /    a curvilinear element
%   |  \ hat(T)         *----*     \            \   /      with x in chi(T)
%   |   \              y1     y2    \            \ /    actually on the torus
%   *----* yhat2                     \            * x1
%  yhat1                              \           ^       with P0(xi) = xi
%      |                               \          |       
%      |                                \ X0      | P0       linear interpolant of  
%      |                                 \         \          element on torus   
%      |                                  \     x3  \     x2   X0(T) = I_h chi(T)
%      |                                   \      *-------*
%      |                X_{hat{T}}          \--->  \X0(T)/
%      ----------------------------------------->   \   /   
%                                                    \ /    
%                                                     * x1
%    with X_{hat{T}} (yhat) = X0( F(yhat) )
%
% Final note of warning:
%   In the notes, we never use the parametric space elements(T) and simply
%   map from (our defined here) X0(T) = K (in the notes) back to
%   hat(T) = hat(K) (in the notes).  Thus we introduce here the map X0 
%   to the parametric space and never use X_{hat{T}} map since all of 
%   our work is done on T and then computed on the reference element of T.
%
%
%   Wenyu Lei
%   Spencer Patty
%   Jan 7, 2016

    %init
    local_stiff = zeros (3,3);
    local_rhs = zeros (3,1);
    
    %affine map info
    mat_B = [y(2,:)-y(1,:); y(3,:)-y(1,:)]';
    det_B = abs(det(mat_B));
    inv_B = [mat_B(2,2),-mat_B(1,2);-mat_B(2,1),mat_B(1,1)]/det_B;
    
    % precompute various things
    q_points = q_yhat*mat_B' + repmat(y(1,:),nq,1);  % [nqx2] list of quadrature points in T
    rhs_vals = rhs_eval(q_points);                   % [nqx1] f(q_points)
    
    % Following Section 3.3.1 of notes:
    %
    % X_{hat{T}}(yhat) = chi_1 + yhat^{T} * [chi_2-chi1; chi3-chi1]  ( [1x2][2x3] = [1x3] )
    % X0(y) = chi_1 + (y-y1)^{T}*B^{-T} * [chi_2-chi1; chi3-chi1]    ( [1x2][2x2][2x3] = [1x3] )
    % 
    % grad_Ih_chi = grad_{y} X0(y) = [chi_2-chi_1; chi_3-chi_1]'*inv(B) [3x2][2x2] = [3x2]
    %
    % G_Gamma = grad I_h chi * grad I_h chi  % constant [2x2] matrix on each cell
    % Q_Gamma = sqrt(det(G_Gamma))           % constant value on each cell
    %
    chi = parameterization(y);
    grad_Ih_chi = [ chi(2,:)-chi(1,:); chi(3,:)-chi(1,:)]'*inv_B;
    G_Gamma=grad_Ih_chi'*grad_Ih_chi;  % Ih_chi metric tensor
    Q_Gamma=sqrt(det(G_Gamma));  % Ih_chi area element

    
    %local stiff and rhs
    for q_index = 1:nq

        % Simplifications for bilinear form:
        %
        % int_{X0(T)} (alpha*surfgrad_phi_i . surfgrad_phi_j + beta*phi_i*phi_j ) d\x
        %
        % note:  surfgrad_phi(x) = grad_phistar(y) * inv(G_Gamma)(y) * grad_X0(y)
        %
        % G_Gamma = grad_X0' grad_X0  ([2x3][3x2])  first fundamental tensor
        % Q_Gamma = sqrt(det(G_Gamma))              Ih chi map area-element transform
        %
        % so \int_{X0(T)} gradphi_i(x) gradphi_j(x) d\x 
        %       = ... = \int_{T} grad_phistar_i(y)' * inv(G_Gamma)(y) * grad_phistar_j(y) * Q dy
        %
        %   now mapping to reference element, we have
        %      grad_phistar(y) = inv(B)*grad_phistar_hat(yhat)
        %
        %       = \int_{hat(T)} (inv(B)*grad_phistarhat_i(yhat))' * inv(G_Gamma)(y) 
        %                             * (inv(B)*grad_phistar_j(yhat)) * Q_Gamma * det(B) dyhat
        %

        grad_phi_star_at_q_point=inv_B'*[grad_hat_phistar_x1(q_index,:); grad_hat_phistar_x2(q_index,:)];
        grad_phi_ij_matrix = grad_phi_star_at_q_point'*(G_Gamma\grad_phi_star_at_q_point);  % = [grad phi_i . grad_phi_j ]_{ij}
        phi_ij_matrix = hat_phistar(q_index,:)'*hat_phistar(q_index,:); % = [phi_i phi_j]_{ij}
        
        % Note: that we add a transpose so that the matrix is assembled
        % properly.  In our case, it is symmetric so it doesn't matter but
        % to be consistent, we note that for  u = sum_i  u_i phi_i and 
        % v = phi_j, the jth row of A[u_i] = F gets the contribution from
        % phi_j as test function.  ie, for a non symmetric problem
        %    
        %      -Delta u + b\cdot grad u + u = f
        %
        % we have the weak form:
        %
        %  sum_{i} u_i * \int grad phi_j . grad phi_i + (b.grad phi_i) phi_j + phi_i phi_j  d\x 
        %                     = \int f phi_j
        %
        % so that the matrix system A U = F has matrix
        %
        %  A_ji = \int grad phi_j . grad phi_i + (b.grad phi_i) * phi_j + phi_i * phi_j d\x
        %
        % Hence the transpose for consistency on the two matrix
        % contributions below.
        
        local_stiff=local_stiff + (alpha * grad_phi_ij_matrix'...% alpha [grad phi_i . grad_phi_j ]_{ji}
                                   + beta * phi_ij_matrix' )...  % beta [phi_i phi_j]_{ji}
                                  * q_weights(q_index)...        % quadrature weight
                                  * det_B...                     % reference area-element transform
                                  * Q_Gamma;                     % Ih_chi area-element transform
                              
        
        % Simplifications for RHS calculation:
        %
        % \int_{chi(T)}  f(x) phi_j(x) d\x
        %     = \int_{T}  fstar(y) phistar_j(y) q(y) dy
        %         = \int_{hat(T)} fstar_hat(yhat) phistar_hat_j(yhat) qhat(yhat) det(B) dyhat
        %
        if(rhs_flag)
            
            % It turns out that to be consistent for our error analysis, the
            % right hand side must still use the true q (chi area element)
            % instead of the Q_Gamma.  Analysis shows that 
            % |q - Q_Gamma|_{L^infty(T)} \leq ch |chi|_{W^{2,\infty}(T)}
            % so that using either one shouldn't change our convergence rates
            % since we get first order by default and then by a duality argument
            % gain an extra h.  However, using q on rhs makes the error 
            % analysis doable.
            grad_chi = grad_parameterization(q_points(q_index,:));
            G=grad_chi'*grad_chi; % chi metric tensor
            q=sqrt(det(G)); % chi area element
            
            local_rhs=local_rhs + rhs_vals(q_index)...        % f(q_point)
                                  *hat_phistar(q_index,:)'... % basis on reference element [3x1]
                                  *q_weights(q_index)...      % quadrature weight
                                  *det_B...                   % reference area-element transform
                                  *q;                         % chi area-element transform
        end
        
%         %
%         % The following shows the local assembling in the for loop. We also
%         % use the original formula (which can be simplified).
%         %
%         for j = 1:3
%             for i = 1:3
%                 local_stiff(j,i)= local_stiff(j,i)...
%                     +(beta*hat_phistar(q_index,i)*hat_phistar(q_index,j)...
%                     +alpha*dot(...
%                     [grad_hat_phistar1(q_index,i), grad_hat_phistar2(q_index,i)]*inv_B*invG*grad_chi',...
%                     [grad_hat_phistar1(q_index,j), grad_hat_phistar2(q_index,j)]*inv_B*invG*grad_chi'))...
%                     *q_weights(q_index)*det_B*q;
%             end
%             if(rhs_flag)
%             local_rhs(j) = local_rhs(j) +...
%                 (rhs_eval(q_point)*hat_phistar(q_index,j))*q_weights(q_index)*det_B*q;
%             end
%         end
    end
    

end