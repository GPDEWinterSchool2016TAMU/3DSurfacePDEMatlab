function [ norm_grad ] = grad_interp_levelset( basis_grad_points,d)
% Return norm of the gradient of the linear interpolation of the level set
% function on a certain point.
%
%  input:
%    basis_grad_points = [4x3] vector of 4 basis function gradients.  Note
%                              that the basis gradients are constant on the
%                              cell T so we do not need a point x in T in
%                              which to evaluate the end result.
%    d = [1x4] vector of distance function values at 4 vertices.
%  output:
%    norm_grad = [1x1] ell_2 norm of gradient of I_h d on cell T since
%                      basis gradients are constant on T.
%

% sum of d values times the corresponding basis gradient can be seen as a
% dot product or in this case a matrix multiplication
vec = d*basis_grad_points; %  [1x3] = [1x4][4x3] 

norm_grad = norm(vec);
end

