%clear all; close all; clc;
format long;
n=64;
[ n_node,n_ele,node,ele,global_ind ] = triangulation_surface( n );

A = sparse([],[],[],n*n,n*n,9*n_ele);
MASS = sparse([],[],[],n*n,n*n,9*n_ele);
rhs = zeros(n*n,1);
q_weights= [1./24,1./24,1./24,9./24];
hatx=[0,1,0,1./3];
haty=[0,0,1,1./3];
nq=4;
% q_weights=[1./2];
% hatx=[1./3];
% haty=[1./3];
% nq=1;
[ hat_phi,hat_phix,hat_phiy ] = FEEVAL( hatx,haty,nq );

% Assemble
ttime(1)=0;
ttime(2)=0;
for i=1 : n_ele
    % local stiff
    v1 = [node(ele(i,1),1), node(ele(i,1),2)];
    v2 = [node(ele(i,2),1), node(ele(i,2),2)];
    v3 = [node(ele(i,3),1), node(ele(i,3),2)];
    tic;
    [ local_stiff,local_rhs ] = local_assembling( v1,v2,v3,hat_phi,hat_phix,hat_phiy,hatx,haty,nq,q_weights,1,1);
    ttime(1)=ttime(1)+toc;
    % copy local to global
    tic;
    for j=1:3
        for k=1:3
            A(global_ind(ele(i,k)),global_ind(ele(i,j)))=A(global_ind(ele(i,k)),global_ind(ele(i,j)))+ local_stiff(j,k);
        end
        rhs(global_ind(ele(i,j)))=rhs(global_ind(ele(i,j)))+local_rhs(j);
    end
    ttime(2)=ttime(2)+toc;
end

tic
solution = A\rhs;
ttime(3)=toc;

% error
exact_sol = zeros(n*n,1);
for i = 1:n
    for j=1:n
        exact_sol(j+(i-1)*n)=exact(node(j+(i-1)*(n+1),:));
    end
end
err_vec =abs(exact_sol - solution);
%assemble mass
for i=1 : n_ele
    % local stiff
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

err = sqrt(transpose(err_vec)*MASS*err_vec)
% err_max = max(abs(err_vec))
ttime

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