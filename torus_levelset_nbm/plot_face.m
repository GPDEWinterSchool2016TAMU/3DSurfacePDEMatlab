function [] = plot_face(face_vertices, T_vertices, color)
% This function takes a set of vertices of a face on Gamma_h subset of the
% tetrahedron T_vertices and plots Gamma_h and the frame of T.
%
% This tool is part of a test suite for the subdivide_error routine
% which takes a tetrahedron T with the distance function values on the 4
% vertices and returns a set of sub faces that span T\cap Gamma_h.  In 
% particular, we are assuming linear finite element approximations so that 
% d_h = I_h d(x) is the linear interpolant of the distance function and
% Gamma_h = { x | d_h(x) = 0 }.
%
% Spencer Patty
% Dec 23, 2015

figure(1)

% indices of vertices for the four faces of tetrahedron
faces_T = [1 2 3; 1 3 4; 1 2 4; 2 3 4];
% indices of face on Gamma_h
faces_F = [1 2 3];

% plot the sub-tetrahedron (note that we need to transpose the
% vertices to give them in the form that patch expects them ie. a column
% vector of vertices size: [3x3] 
patch('Vertices', face_vertices, 'Faces', faces_F, 'FaceColor',color, ...
      'EdgeColor','black')
%plot frame of original tetrahedron T
patch('Vertices', T_vertices, 'Faces', faces_T, 'FaceColor','none',...
      'EdgeColor','black', 'LineWidth',2)
  
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal
