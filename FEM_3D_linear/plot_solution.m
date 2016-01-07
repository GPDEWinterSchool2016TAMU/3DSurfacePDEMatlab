function [  ] = plot_solution( n,sol ,x)
% print the solution with tri-slice plot centered at the point x

hh=1/n;
xx=0:hh:1;
yy=0:hh:1;
zz=0:hh:1;

[XX,YY,ZZ] = meshgrid(xx,yy,zz); 

sol=reshape(sol,n+1,n+1,n+1);

slice(XX,YY,ZZ,sol,x(1),x(2),x(3)); 
end

