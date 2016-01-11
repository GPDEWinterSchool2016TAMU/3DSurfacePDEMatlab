% Standard finite element method on the unique square
% We approximate the following equation
% 2u-\Delta u = f in \Omega=[0,1]^2
%       du/dn = 0 on the boundary of \Omega
%
% where \Delta is the Laplace operator and du/dn is 
% the normal derivative on the boundary.
%
% Wenyu Lei
% Jan 7, 2016

clear all; close all; clc;
format long;
n=64;
% Generate the uniform mesh on the unit suqare
[ n_node,n_ele,node,ele] = triangulation_square( n );

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
[ hat_phi_at_q, hat_phix_at_q, hat_phiy_at_q ] = FEEVAL( q_yhat,nq );


% Assembling
for cell = 1:n_ele
    
    % Get local stiffness matrix and local rhs
    cell_ind = ele(cell,1:3);     % [1x3]
    vertices = node(cell_ind, :); % [3x2]  
    
    [ local_stiff,local_rhs ] ...
        = local_assembling(vertices,...
                           hat_phi_at_q, hat_phix_at_q, hat_phiy_at_q,...
                           q_yhat,nq,q_weights,...
                           1,2,1); % alpha, beta, rhs_flag
    
    % Copy local to global
    A(cell_ind,cell_ind) ...
        = A(cell_ind,cell_ind) + local_stiff;   %[3x3]
    rhs(cell_ind) = rhs(cell_ind) + local_rhs;  %[3x1]
end

% Apply back slash solver
solution = A\rhs;

% L2 error computation
exact_sol = exact(node);

err_vec =exact_sol - solution;
% Assemble mass matrix
for cell = 1:n_ele
    
    % Local mass matrix
    cell_ind = ele(cell,1:3);     % [1x3]
    vertices = node(cell_ind, :); % [3x2]
    [local_mass,~] = ...
        local_assembling( vertices,...
                          hat_phi_at_q, hat_phix_at_q, hat_phiy_at_q,...
                          q_yhat,nq,q_weights,...
                          0,1,0); % a, beta, rhs_flag
                      
    % copy local to global
    MASS(cell_ind,cell_ind)...
        =MASS(cell_ind,cell_ind) + local_mass; %[3x3]
end

% print out the error
l2_err = sqrt(transpose(err_vec)*MASS*err_vec)

% plot the solution
figure(1)
plot_solution(n,solution);
title('solution');
xlabel('X'); ylabel('Y');

% plot the error
figure(2)
plot_solution(n,err_vec);
title('error');
xlabel('X'); ylabel('Y');
