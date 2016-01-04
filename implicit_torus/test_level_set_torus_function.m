%
% Good for profiling the entire setup and solve.
%


quad_degree_bulk = 2;
quad_degree_face = 2;
plot_solution = 0;  %false
n = 32;
h = 4/n;
error_l2 = level_set_torus_function(n,h, quad_degree_bulk, quad_degree_face, plot_solution); 