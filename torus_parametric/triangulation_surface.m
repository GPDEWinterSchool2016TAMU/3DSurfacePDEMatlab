function [ n_node,n_ele,pm_node,ele,global_ind,global_ind_inverse ] = triangulation_surface( n )
% Generate the triangulation on the parametric domain.
% In our case, for a torus, this domain is a square [-\pi,\pi]*[-\pi,\pi].
% Note that, points on the parametric domain are periodic on both 
% vertical and horizontal directions. So we only have n*n vertieces 
% on the surface (not (n+1)*(n+1)). To get the index set on the surface,
% we build a mapping from the index set from the square to the index 
% set of the torus. For example, we should label vertices on the first 
% line x=-pi by 1, 2, 3, ... ,n ,1 (totally n+1 points on parametric 
% domain but only n points on the torus).
% We return this mapping as 'global_ind'
%
% Wenyu Lei
% Dec 30, 2015

% Triangulation on the square [-\pi,\pi]*[-\pi,\pi].
h=2*pi/n;
xx=-pi:h:pi;
yy=-pi:h:pi;
X=zeros((n+1)*(n+1),2);
% for i=1:(n+1)
%     for j=1:n+1
%         X(j+(n+1)*(i-1),1)=xx(i);
%         X(j+(n+1)*(i-1),2)=yy(j);
%     end
% end
[YY,XX]=meshgrid(yy,xx);
X(:,1)=reshape(XX,(n+1)^2,1);
X(:,2)=reshape(YY,(n+1)^2,1);
dt = delaunayTriangulation(X);

% Extract node and connectivity list of the square mesh.
ele=dt.ConnectivityList;
n_ele=length(ele);
pm_node=dt.Points;
n_node=n*n;  % number of vertiecs on the surface.

%
% identify the periodic nodes and create mapping to adjust the nodes
%
global_ind=zeros(n+1,n+1);
global_ind(1:n,1:n)=reshape(1:n^2,n,n);
global_ind(end,1:end-1)=global_ind(1,1:end-1);
global_ind(1:end-1,end)=1:n;
global_ind(end,end)=1;
global_ind=reshape(global_ind,1,(n+1)*(n+1));


global_ind_inverse=reshape(1:(n+1)^2,[n+1 n+1]);
global_ind_inverse=reshape(global_ind_inverse(1:n,1:n),[1 n*n]);
