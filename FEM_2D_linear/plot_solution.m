function [  ] = plot_solution( n,solution )
% Plot the solution
h=1/n;
xx=0:h:1;
yy=0:h:1;
solution=reshape(solution,n+1,n+1);
mesh(xx,yy,solution);
end

