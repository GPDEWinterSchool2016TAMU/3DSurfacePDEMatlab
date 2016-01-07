% Implementation of finite element method on 3D
% Note we assembe the stiffness and mass matrices directly on the real space
%
% The equation we try to approximate is
% u-\Delta u = f in \Omega=[0,1]^3
%      du/dn = 0 on the boundary of \Omega
%
% Wenyu Lei
% Jan 7, 2016
clear all; close all; clc; 
format long;

% mesh size
n=16;   

% quadrature degree of accuracy
quad_degree_bulk = 2; % for tetrahedra

% Mesh generation
fprintf('...Generating Mesh...\n');
[n_node,n_ele,node,ele]=box_triangulation(n); % structured

% Matrix and vector initialization

% We estimate a maximum of 15 interactions on each row for the uniform mesh
% case with an average of 13 interactions over all rows .
%
% For the tetgen case, a maximum is not calculable, but for h values tested in 
% [0.0625,1] we have on average 13 interactions. It is probably safe again
% to use bound 15 but you should be aware that running out of memory
% will slow you down alot.
A = spalloc(n_node, n_node, 15*n_node);
MASS = spalloc(n_node, n_node, 15*n_node);
F = zeros(n_node,1);

% Assembling
fprintf('...Assembling Stiffness Matrix and RHS...\n');
for cell = 1:n_ele
    
    cell_node_ind = ele(cell,1:4);     % [1x4]
    vertices = node(cell_node_ind, :); % [4x3]

    [ lstiff,lrhs ] = local_assembling( vertices,quad_degree_bulk,1,1,1);

    A(cell_node_ind,cell_node_ind) = A(cell_node_ind,cell_node_ind) + lstiff'; %[4x4]
    F(cell_node_ind) = F(cell_node_ind) + lrhs';  %[4x1]

end

% Apply back slash solver
fprintf('...Solving System...\n')
sol=A\F;


% compute error of solution
% see the standard 2D code for more details
exact_sol = exact(node);
err_vec = exact_sol-sol;

% Assembling mass matrix
fprintf('...Assembling mass matrix...\n');
for cell = 1:n_ele
    
    cell_node_ind = ele(cell,1:4);     % [1x4]
    vertices = node(cell_node_ind, :); % [4x3]

    [ lstiff,~ ] = local_assembling( vertices,quad_degree_bulk,0,1,0);

    MASS(cell_node_ind,cell_node_ind) = MASS(cell_node_ind,cell_node_ind) + lstiff; %[4x4]
    F(cell_node_ind) = F(cell_node_ind) + lrhs';  %[4x1]

end

fprintf('...Computing L2 error...\n');
l2_err = sqrt(err_vec'*MASS*err_vec)

% plot the solution with tri-slice centered at point [0.25,0.25,0.25]
plot_solution( n,sol,[0.25,0.25,0.25] );
