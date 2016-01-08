% This script does the error analysis and creates the table of rates for
% you.  Specify the number of steps desired, it starts with n=8 and doubles
% from there.
%
%  Spencer Patty
%  Jan 7, 2016
%
clear all; close all; clc;

% number of steps to compute error for
num_steps = 6;

nvec = 2.^(3+(0:(num_steps-1)))';
error = zeros(num_steps,1);
for i = 1:num_steps
    fprintf('Step %d: Computing error for n = %d\n',i, nvec(i));
    error(i) = surface_laplacian_torus_function(nvec(i),0);
    fprintf('  Error for n = %d is %1.3e.\n', nvec(i), error(i));
end
    
% compute rate of convergence
rates = zeros(num_steps,1);
ind = 1:(num_steps-1);
rates(2:end) = log(error(ind)./error(ind+1))'/log(2);

% display error analysis table
fprintf('\n n    error   rate\n');
fprintf('%3d %1.3e %1.2f\n', [nvec, error, rates]')
