% We consider solving the following equation
% u-\Delta_\Gamma u = f on \Gamma
% in a phase field approach
% 
% The weak formulation is
% \int_D (1-\psi^2)(\nabla u\cdot\nabla v+uv-fv) dx=0
% for every v in H^1(D)
% Here \psi is a phase field function
% D is a box domian containing the torus

% Wenyu Lei
% Jan 6, 2016
clear all; close all; clc;
format long;
n=8; h=4./n; 
epsilon=sqrt(h);

%[n_node,n_ele,node,ele]=uniform_bulk_mesh(n,h); % unsturcture
[n_node,n_ele,node,ele]=uniform_box_mesh(n,h); % structed

% Initialization
A = sparse([],[],[],n_node,n_node,15*n_node);
rhs = zeros(n_node,1);

% Initialize shape value and gradient in the reference triangle
a = (5-sqrt(5))/20.;
q_weights= [1./24,1./24,1./24,1./24];
hatx=[a,1-3*a,a,a];
haty=[a,a,1-3*a,a];
hatz=[a,a,a,1-3*a];
nq=4;

[ hat_phi,hat_phix,hat_phiy,hat_phiz ] = FEEVAL( hatx,haty,hatz,nq );

% Assembling
% Here M_epsilon is defined as 
% M_epsilon = \int_D (1-\psi^2) dx
% This can be done by summing up local contribution 
% local_me = \int_\tao (1-\psi^2) dx
% We will use M_epsilon to compute the energy error.
M_epsilon=0;
for i=1 : n_ele
    % local stiff
    cell_node_ind = ele(i,1:4);
	v=node(cell_node_ind,:);
    
    [ local_stiff,local_rhs,local_me] = local_assembling( v,hat_phi,hat_phix,hat_phiy,hat_phiz,hatx,haty,hatz,nq,q_weights,1,1,epsilon,1);
    
    % copy local to global
    M_epsilon=M_epsilon+local_me;
    A(cell_node_ind,cell_node_ind) = A(cell_node_ind,cell_node_ind) + local_stiff'; %[4x4]
    rhs(cell_node_ind) = rhs(cell_node_ind) + local_rhs; %[4x1]
end

% back slash solver
solution = A\rhs;


exact_sol = exact(node);

% The energy norm is defined as
% \|v\|_\epsilon^2
% := (\int_D (1-\psi^2)(|\nabla v|^2+v^2)dx)/M_epsilon
% As we compute the l2-err in the parametric case, in FE space, 
% the error in the energy norm can be computed as the following:
err_energy = sqrt(transpose(exact_sol-solution)*A*(exact_sol-solution)/M_epsilon)

% We also report the l2 error on the Gamma_h
% (Recall the discrete the surface Gamma_h is defined by
% Gamma_h={x\in D : I_h\phi(x)==0}, where \phi is the level set function).
err_l2_gamma_h = compute_error_L2_Gamma_h(solution,node,ele,n_ele,2)

% visualize the nuerical results on Gamma_h
plot_solution_on_Gamma_h(solution,node,ele,n_ele,0);
