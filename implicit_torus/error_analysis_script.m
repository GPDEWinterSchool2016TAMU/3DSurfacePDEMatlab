%
% This scropt will run a sequence of 
%
%
clear all;  close all;  clc;

numsteps = 4;  % num steps 
quad_degree_bulk = 2;
quad_degree_face = 2;
plot_solution = 0;  %false
error_l2 = zeros(1,numsteps);
h = zeros(1,numsteps);
for i = 1:numsteps
    n = 2^(3+i);  
    h(i) = 4/n;
%     error_l2(i) = rand;
    fprintf(' computing solution for h = %1.5e.\n', h(i));
    error_l2(i) = level_set_torus_function(n,h(i), quad_degree_bulk, quad_degree_face, plot_solution);    
end

% compute rates
rate_l2 = zeros(1,numsteps);
for i = 2:numsteps
    rate_l2(i) = log(error_l2(i)/error_l2(i-1))/log(h(i)/h(i-1));
end

% output table
fprintf('\n    h        error    rate \n');
for i = 1:numsteps
    fprintf('%1.3e  %1.3e  %1.2f\n', h(i), error_l2(i), rate_l2(i));
end




