function [ res ] = grad_interp_levelset( basis_grad_points,d)
% Return norm of the gradient of the linear interpolation of the level set
% function on a certain point.
% d a 1 by 4 vector which stores the distance function values for 4
% vertices of the cell.
% basis_grad_points stores 4 shape gradients on a given point.
% So the output res is the eveluation on this given point.
vec = zeros(1,3);
for i=1:4
    vec = vec+d(i)*basis_grad_points(i,:);
end

res = norm(vec);
end

