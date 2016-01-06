function [] = plot_tetrahedron(Tsub_vertices, T_vertices, color)
% This function takes a set of vertices of a tetrahedron Tsub which is a 
% subset of a tetrahedron T and plots the faces of Tsub with color and
% gives a wire frame of the original tetrahedron T to see it in
% perspective.  This tool is part of a test suite for the subdivide routine
% which takes a tetrahedron T with the distance function values on the 4
% vertices and returns a set of sub tetrahedra that span T\cap D_h  where
% D_h = { x | |d_h(x)| < h }.  In particular, we are assuming linear finite
% element approximations so that d_h = I_h d(x) is the linear interpolant
% of the distance function.
%
% Spencer Patty
% Dec 23, 2015

figure(1)
% indices of vertices for the four faces of tetrahedron
faces = [1 2 3; 1 3 4; 1 2 4; 2 3 4];

% plot the sub-tetrahedron (note that we need to transpose the
% vertices to give them in the form that patch expects them ie. a column
% vector of vertices size: [4x3] 
patch('Vertices', Tsub_vertices, 'Faces', faces, 'FaceColor',color, ...
      'EdgeColor','black')
%plot frame of original tetrahedron T
patch('Vertices', T_vertices, 'Faces', faces, 'FaceColor','none',...
      'EdgeColor','black', 'LineWidth',2)
  
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal
