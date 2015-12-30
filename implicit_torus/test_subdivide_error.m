% test_subdivide_error.m
%
% This test file contains at least one test for every possible scenario of
% nodes inside on the gamma or outside the band D_h = {x | |d_h(x)| < h } 
% and the proper subface of the intersection T \cap D_h into 
% tetrahedra.  Here we assume that d_h(x) = I_h d(x) is the linear
% interpolant of the distance function d(x).
%
% It is best to go through the below scenarios section by section using 
% the "Run Section" command or the "Run and Advance" command.  
% 
% NOTE: This is not a script that will tell you anything if you merely
% run it top to bottom.
%
% The following scenarios are tested:
%
% (4p 0gamma 0m)  returns 0 face(s)
% (3p 1gamma 0m)  returns 0 face(s)
% (3p 0gamma 1m)  returns 1 face(s)
% (2p 2gamma 0m)  returns 0 face(s)
% (2p 1gamma 1m)  returns 1 face(s)
% (2p 0gamma 2m)  returns 2 face(s)
% (1p 3gamma 0m)  returns 1 face(s)
% (1p 2gamma 1m)  returns 1 face(s)
% (1p 1gamma 2m)  returns 1 face(s)
% (1p 0gamma 3m)  returns 1 face(s)
% (0p 3gamma 1m)  returns 1 face(s)
% (0p 2gamma 2m)  returns 0 face(s)
% (0p 1gamma 3m)  returns 0 face(s)
% (0p 0gamma 4m)  returns 0 face(s)
%
% and then a uniform random generator on d\in [-h,h] is also provided to
% test all sorts of cases and see the many possibilities of subfaces
%
% The test cases are all specified on the standard reference element
% tetrahedron called vertices here.
%
% Spencer Patty
% Dec 23, 2015

vertices = [[0,0,0]; [1,0,0]; [0,1,0]; [0,0,1]];  

%% test random values

h=1; % generate random numbers from uniform [-h,h]
d_at_vertices = -h + 2*h*rand(1,4)  
subface = subdivide_error(vertices, d_at_vertices);
colors = {'red','blue','green'};
i=1; 
clf;
for F=subface
    plot_face(F{1},vertices,colors{i});
    i=i+1;
end

%% test (4p 0gamma 0m)

d_at_vertices = [2, 3, 4, 5];
compare = {}; 
subface = subdivide_error(vertices, d_at_vertices);

%% test (3p 1gamma 0m)

d_at_vertices = [3, 4, 5, 0];
compare = {};
subface = subdivide_error(vertices, d_at_vertices);

%% test (3p 0gamma 1m)

d_at_vertices = [1, 1, -1, 1];
compare = {}; 
subface = subdivide_error(vertices, d_at_vertices);


%% test (2p 2gamma 0m)

d_at_vertices = [1, 1, 0, 0];
compare = {}; 
subface = subdivide_error(vertices, d_at_vertices);

%% test (2p 1gamma 1m)

d_at_vertices = [-1, 0, 1, 1];
compare = {[[0;0;0.5], [1;0;0], [0;0.5;0]]}; 
subface = subdivide_error(vertices, d_at_vertices);
compare{1}
clf
plot_face(subface{1},vertices,'red');

%% test (2p 0gamma 2m)

d_at_vertices = [1, -1, 1, -1];
t = sqrt(2)*abs(2-h)/(abs(2-h)+abs(0-h));
compare = {[[0;0;0.5], [0.5;0;0], [0.5;0.5;0]],...
           [[0;0;0.5], [0;0.5;0.5], [0.5;0.5;0]]}; 
subface = subdivide_error(vertices, d_at_vertices);
compare{1}
compare{2}

clf
plot_face(subface{1},vertices,'red');
plot_face(subface{2},vertices,'blue');

%% test (1p 3gamma 0m)

d_at_vertices = [1, 0, 0, 0];
compare = {[[1;0;0], [0;1;0], [0;0;1]]}; 
subface = subdivide_error(vertices, d_at_vertices);
compare{1}
clf
plot_face(subface{1},vertices,'red');

%% test (1p 2gamma 1m)

d_at_vertices = [1, 0, -1, 0];
compare = {[[1;0;0], [0;0.5;0], [0;0;1]]}; 
subface = subdivide_error(vertices, d_at_vertices);
compare{1}
clf
plot_face(subface{1},vertices,'red');

%% test (1p 1gamma 2m)

d_at_vertices = [-1, 0, 1, -1];
compare = {[[1;0;0], [0;0.5;0], [0;0.5;0.5]]}; 
subface = subdivide_error(vertices, d_at_vertices);
compare{1}
clf
plot_face(subface{1},vertices,'red');

%% test (1p 0gamma 3m)

d_at_vertices = [1, -2, -2, -2];
t = abs(1)/(abs(1)+abs(-2));
compare = {[[t;0;0], [0;t;0], [0;0;t]]}; 
subface = subdivide_error(vertices, d_at_vertices);
compare{1}
clf
plot_face(subface{1},vertices,'red');

%% test (0p 3gamma 1m)

d_at_vertices = [0, 0, -1, 0];
compare = {[[0;0;0], [1;0;0], [0;0;1]]}; 
subface = subdivide_error(vertices, d_at_vertices);
compare{1}
clf
plot_face(subface{1},vertices,'red');

%% test (0p 2gamma 2m)

d_at_vertices = [-1, 0, 0, -2];
compare = {}; 
subface = subdivide_error(vertices, d_at_vertices);

%% test (0p 1gamma 3m)

d_at_vertices = [-1, -1, 0, -1];
compare = {}; 
subface = subdivide_error(vertices, d_at_vertices);

%% test (0p 0gamma 4m)

d_at_vertices = [-1, -1, -2, -3];
compare = {}; 
subface = subdivide_error(vertices, d_at_vertices);


