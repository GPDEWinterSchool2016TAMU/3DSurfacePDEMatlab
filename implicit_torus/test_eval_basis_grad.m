% test vectorization over basis function indices of eval_basis_grad()
%

vertices = [[0,0,0]; [1,0,0]; [0,1,0]; [0,0,1]];
ind = 1:4;

% construct a random point inside reference triangle element
x = zeros(1,3);
x(1) = rand;
x(2) = rand*(1-x(1));
x(3) = rand*(1-x(1)-x(2));

% we know single terms are working so compare vector to the single ones
cell_basis_true_grad_at_x = zeros(4,3);
cell_basis_true_grad_at_x(1,:) = eval_basis_grad(x,1, vertices);
cell_basis_true_grad_at_x(2,:) = eval_basis_grad(x,2, vertices);
cell_basis_true_grad_at_x(3,:) = eval_basis_grad(x,3, vertices);
cell_basis_true_grad_at_x(4,:) = eval_basis_grad(x,4, vertices);

% evaluate the vector of indices
cell_basis_grad_at_x = eval_basis_grad(x,ind, vertices);

if (norm(cell_basis_grad_at_x - cell_basis_true_grad_at_x) > 1e-10)
    fprintf('Error, vectorization of eval_basis_grad is not working correctly.\n');
else
    fprintf('eval_basis_grad() index vectorization is working correctly!\n');
end