% test vectorization over basis function indices of eval_basis_grad()
%

vertices = [[0,0,0]; [1,0,0]; [0,1,0]; [0,0,1]];
ind = 1:4;

d = rand(1,4);  % pick random distance function values on T vertices

% construct a random point inside reference triangle element
x = zeros(1,3);
x(1) = rand;
x(2) = rand*(1-x(1));
x(3) = rand*(1-x(1)-x(2));

% evaluate the vector of indices
cell_basis_grad_at_x = eval_basis_grad(x, ind, vertices);

% evaluate the |I_h d(x)| on T
norm_Ih_d_at_x = grad_interp_levelset(cell_basis_grad_at_x, d);

% true |I_h d(x)| on T
norm_true_Ih_d_at_x = norm(d(1)*[-1;-1;-1] + d(2)*[1;0;0] + d(3)*[0;1;0] + d(4)*[0;0;1]);

% test for error
if (abs(norm_Ih_d_at_x - norm_true_Ih_d_at_x) > 1e-10)
    fprintf('Error, vectorization of eval_basis_grad is not working correctly.\n');
else
    fprintf('eval_basis_grad() index vectorization is working correctly!\n');
end