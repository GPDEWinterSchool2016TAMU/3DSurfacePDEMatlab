% Finite element method on parametrized surfaces
% Consider the following equation
%
% u-\Delta_\Gamma u = f on \Gamma
%
% where \Gamma is a surface (torus in this code) in 3D and 
% \Delta_gamma is the Laplace-Beltrami operator
% We parameterize the torus by
% x_1 = (R+r\cos\theta)\cos\phi
% x_2 = (R+r\cos\theta)\sin\phi
% x_3 = r\sin\theta
% for \theta and \phi \in [-\pi,\pi].
%
% Wenyu Lei
% Dec 30, 2015

clear all; close all; clc;
format long;

% Get the triangulation on the parametric domain
n=64;
[ n_node,n_ele,node,ele,global_ind ] = triangulation_surface( n );

% Initialization
A = sparse([],[],[],n*n,n*n,9*n_ele);
MASS = sparse([],[],[],n*n,n*n,9*n_ele);
rhs = zeros(n*n,1);

% Since we are going to compute the integral on the reference element,
% we need provide the quadrature rule for the reference triangle.
% Meanwhile, we can also provide infomation of shape functions on the 
% reference element, i.e. function values (hat_phi) nad function gradients
% (hat_phix and hat_phiy) on each quadrature points.
q_weights= [1./24,1./24,1./24,9./24];
hatx=[0,1,0,1./3];
haty=[0,0,1,1./3];
nq=4;
[ hat_phi,hat_phix,hat_phiy ] = FEEVAL( hatx,haty,nq );

% Assembling
for i=1 : n_ele
    % Get local stiffness matrix and local rhs
    v1 = [node(ele(i,1),1), node(ele(i,1),2)];
    v2 = [node(ele(i,2),1), node(ele(i,2),2)];
    v3 = [node(ele(i,3),1), node(ele(i,3),2)];
    [ local_stiff,local_rhs ] = local_assembling( v1,v2,v3,hat_phi,hat_phix,hat_phiy,hatx,haty,nq,q_weights,1,1);
    % Copy local to global
    for j=1:3
        for k=1:3
            A(global_ind(ele(i,k)),global_ind(ele(i,j)))=A(global_ind(ele(i,k)),global_ind(ele(i,j)))+ local_stiff(j,k);
        end
        rhs(global_ind(ele(i,j)))=rhs(global_ind(ele(i,j)))+local_rhs(j);
    end
end

% Apply back slash solver
solution = A\rhs;


% L2 error computation
% Here we are going to compute L2 norm of  
% u_h - I_h u
% Here u_h is our numerical results, i.e. the linear combinatioin of 
% basis function with coefficients from the 'solution' vector we 
% just computed.
% I_h u is the largrange interpolation of the exact solution u.
% So I_h u is also the linear combination of the basis function 
% with the coefficients from the evaluation on all vertices in the mesh.
% (we call this coefficient vecter 'exact_sol')
% So the square of the l2 norm of u_h - I_h u should be
% (exact_sol-solution)'M(exact_sol-solution)
% Here M is the mass matrix.

% Get the coefficients of I_h u
exact_sol = zeros(n*n,1);
for i = 1:n
    for j=1:n
        exact_sol(j+(i-1)*n)=exact(node(j+(i-1)*(n+1),:));
    end
end
err_vec =exact_sol - solution;
%Assemble mass matrix
for i=1 : n_ele
    % Local mass matrix
    v1 = [node(ele(i,1),1), node(ele(i,1),2)];
    v2 = [node(ele(i,2),1), node(ele(i,2),2)];
    v3 = [node(ele(i,3),1), node(ele(i,3),2)];
    [local_mass,lrhs] = local_assembling( v1,v2,v3,hat_phi,hat_phix,hat_phiy,hatx,haty,nq,q_weights,0,1);
    % copy local to global
    for j=1:3
        for k=1:3
            MASS(global_ind(ele(i,k)),global_ind(ele(i,j)))=MASS(global_ind(ele(i,k)),global_ind(ele(i,j)))+ local_mass(j,k);
        end
    end
end

% print out the error
err = sqrt(transpose(err_vec)*MASS*err_vec)

% Visualization
% Using the function 'patch' to visualize each triangle.
% Colors are decided by the value on the vertices.
sv=zeros(n*n,3);
for i=1:n
    for j=1:n
        sv(j+(i-1)*n,:)=parameterization(node(j+(i-1)*(n+1),:));
    end
end
sele=zeros(n_ele,3);
for i=1:n_ele
    for j=1:3
        sele(i,j)=global_ind(ele(i,j));
    end
end

figure(1);
axis([-2,2,-2,2,-2,2]);
for i=1:n_ele
    XX=[sv(sele(i,1),1); sv(sele(i,2),1);sv(sele(i,3),1)];
    YY=[sv(sele(i,1),2); sv(sele(i,2),2);sv(sele(i,3),2)];
    ZZ=[sv(sele(i,1),3); sv(sele(i,2),3);sv(sele(i,3),3)];
    CC=[solution(sele(i,1));solution(sele(i,2));solution(sele(i,3))];
    patch(XX,YY,ZZ,CC,'EdgeColor','interp');
end