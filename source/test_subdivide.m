% test_subdivide_suite.m
%
% Spencer Patty
% Dec 23, 2015
%
% This test file contains at least one test for every possible scenario of
% nodes inside on the edge or outside the band D_h = {x | |d_h(x)| < h } 
% and the proper subdivision of the intersection T \cap D_h into 
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
% (4out 0edge 0in)  returns 0 sub-tetrahedra
% (3out 1edge 0in)  returns 0 sub-tetrahedra
% (3out 0edge 1in)  returns 1 sub-tetrahedra
% (2out 2edge 0in)  returns 0 sub-tetrahedra
% (2out 1edge 1in)  returns 1 sub-tetrahedra
% (2out 0edge 2in)  returns 3 sub-tetrahedra
% (1out 3edge 0in)  returns 0 sub-tetrahedra
% (1out 2edge 1in)  returns 1 sub-tetrahedra
% (1out 1edge 2in)  returns 2 sub-tetrahedra
% (1out 0edge 3in)  returns 3 sub-tetrahedra
% (0out 3edge 1in)  returns 1 sub-tetrahedra
% (0out 2edge 2in)  returns 1 sub-tetrahedra
% (0out 1edge 3in)  returns 1 sub-tetrahedra
% (0out 0edge 4in)  returns 1 sub-tetrahedra
%
% and then a uniform random generator on d\in [0,2h] and [-2h,0] is also provided to
% test all sorts of cases and see the many possibilities of subdivisions
%
% The test cases are all specified on the standard reference element
% tetrahedron called vertices here.

vertices = [[0;0;0], [1;0;0], [0;1;0], [0;0;1]]; 
h = 1; % don't change from h=1


%% test random values of d positive

d_at_vertices = 0 + (2*h - 0)*rand(1,4)  % generate random numbers from uniform [0,2h]
subdivision = subdivide(vertices, d_at_vertices, h);
colors = {'red','blue','green'};
i=1; 
clf;
for T=subdivision
    plot_tetrahedron(T{1},vertices,colors{i});
    i=i+1;
end

%% test random values of d negative

d_at_vertices = 0 - (2*h - 0)*rand(1,4)  % generate random numbers from uniform [-2h, 0]
subdivision = subdivide(vertices, d_at_vertices, h);
colors = {'red','blue','green'};
i=1; 
clf;
for T=subdivision
    plot_tetrahedron(T{1},vertices,colors{i});
    i=i+1;
end


%% test (4out 0edge 0in)

d_at_vertices = [2, 3, 4, 5];
compare = {}; 
subdivision = subdivide(vertices, d_at_vertices, h);

%% test (3out 1edge 0in)

d_at_vertices = [3, 4, 5, h];
compare = {};
subdivision = subdivide(vertices, d_at_vertices, h);

%% test (3out 0edge 1in)  test around h

d_at_vertices = [2, 2, 0, 2];
t = sqrt(2)*abs(2-h)/(abs(2-h) + abs(0-h));
compare = {[[0;1;0], [0;1/2;0], [1-t;t;0], [0;t; 1-t]]  }; 
subdivision = subdivide(vertices, d_at_vertices, h);
compare{1}
subdivision{1}
clf
plot_tetrahedron(subdivision{1},vertices,'red');

%% test (3out 0edge 1in)  test around -h

d_at_vertices = [-2, -2, 0, -2];
t = sqrt(2)*abs(-2+h)/(abs(-2+h) + abs(0+h));
compare = {[[0;1;0], [0;2/3;0], [1-t;t;0], [0;t; 1-t]]}; 
subdivision = subdivide(vertices, d_at_vertices, h);
clf
plot_tetrahedron(subdivision{1},vertices,'red');


%% test (2out 2edge 0in)

d_at_vertices = [h, h, 2, 2];
compare = {}; 
subdivision = subdivide(vertices, d_at_vertices, h);

%% test (2out 1edge 1in)

d_at_vertices = [0, h, 2, 2];
t = sqrt(2)*abs(2-h)/(abs(2-h)+abs(0-h));
compare = {[[0;0;0], [1;0;0],   [0;0.5;0], [0;0;0.5]]}; 
subdivision = subdivide(vertices, d_at_vertices, h);
compare{1}
clf
plot_tetrahedron(subdivision{1},vertices,'red');

%% test (2out 0edge 2in)

d_at_vertices = [0, 0, 2, 2];
t = sqrt(2)*abs(2-h)/(abs(2-h)+abs(0-h));
compare = {[[0;0;0],   [0;0.5;0], [0;0;0.5], [t;1-t;0]],...
           [[0;0;0],   [1.0;0;0], [0;0;0.5], [t;1-t;0]],...
           [[0;0;0.5], [1;0;0],   [t;0;1-t], [t;1-t;0]]}; 
subdivision = subdivide(vertices, d_at_vertices, h);
compare{1}
compare{2}
compare{3}
clf
plot_tetrahedron(subdivision{1},vertices,'red');
plot_tetrahedron(subdivision{2},vertices,'blue');
plot_tetrahedron(subdivision{3},vertices,'green');

%% test (1out 3edge 0in)

d_at_vertices = [2, h, h, h];
compare = {}; 
subdivision = subdivide(vertices, d_at_vertices, h);

%% test (1out 2edge 1in)

d_at_vertices = [2, h, h, 0];
compare = {[[1;0;0], [0;1;0], [0;0;1], [0;0;0.5]]}; 
subdivision = subdivide(vertices, d_at_vertices, h);
compare{1}
clf
plot_tetrahedron(subdivision{1},vertices,'red');

%% test (1out 1edge 2in)

d_at_vertices = [2, h, 0, 0];
compare = {[[1;0;0], [0;1;0], [0;0.5;0], [0;0;0.5]],...
           [[1;0;0], [0;1;0], [0;0;1.0], [0;0;0.5]]}; 
subdivision = subdivide(vertices, d_at_vertices, h);
compare{1}
compare{2}
clf
plot_tetrahedron(subdivision{1},vertices,'red');
plot_tetrahedron(subdivision{2},vertices,'blue');

%% test (1out 0edge 3in)

d_at_vertices = [2, 0, 0, 0];
t = sqrt(2)*abs(2-h)/(abs(2-h)+abs(0-h));
compare = {[[1;0;0], [0;1;0],   [0;0;1],   [0;0;0.5]],...
           [[1;0;0], [0;1;0],   [0.5;0;0], [0;0;0.5]],...
           [[0;1;0], [0;0.5;0], [0.5;0;0], [0;0;0.5]]}; 
subdivision = subdivide(vertices, d_at_vertices, h);
compare{1}
compare{2}
compare{3}
clf
plot_tetrahedron(subdivision{1},vertices,'red');
plot_tetrahedron(subdivision{2},vertices,'blue');
plot_tetrahedron(subdivision{3},vertices,'green');

%% test (0out 3edge 1in)

d_at_vertices = [h, h, 0.75, h];
compare = {[[0;0;0], [1;0;0], [0;1;0], [0;0;1]]}; 
subdivision = subdivide(vertices, d_at_vertices, h);
compare{1}
clf
plot_tetrahedron(subdivision{1},vertices,'red');

%% test (0out 2edge 2in)

d_at_vertices = [0, h, 0.75, h];
compare = {[[0;0;0], [1;0;0], [0;1;0], [0;0;1]]}; 
subdivision = subdivide(vertices, d_at_vertices, h);
compare{1}
clf
plot_tetrahedron(subdivision{1},vertices,'red');

%% test (0out 1edge 3in)

d_at_vertices = [0, 0.5, 0.75, h];
compare = {[[0;0;0], [1;0;0], [0;1;0], [0;0;1]]}; 
subdivision = subdivide(vertices, d_at_vertices, h);
compare{1}
clf
plot_tetrahedron(subdivision{1},vertices,'red');

%% test (0out 0edge 4in)

d_at_vertices = [0, 0.5, 0.75, 0.25];
compare = {[[0;0;0], [1;0;0], [0;1;0], [0;0;1]]}; 
subdivision = subdivide(vertices, d_at_vertices, h);
compare{1}
clf
plot_tetrahedron(subdivision{1},vertices,'red');

