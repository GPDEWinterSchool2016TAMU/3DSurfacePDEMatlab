function [ grads ] = eval_basis_grad(x,i,v )
% Return ith shape gradient (linear element) on the point x.
% Note that the input argument v is a 4 by 3 matrix
% and jth row of v indicates the coordinates of jth vertex
% in this cell.
%
% input:
%   x = [1x3] real coordinate in R^3
%   i = [nx1] list of basis functions of which to evaluate basis elements
%             at x
%   v = [4x3] vector of 4 vertex coordinates in R^3
% output:
%   grads = [nx3] vector of basis function gradients i evaluated at x.
%

ind=[1,2,3,4,1,2,3];
vec1=v(ind(i+2),:)-v(ind(i+1),:); % [nx3]
vec2=v(ind(i+3),:)-v(ind(i+1),:); % [nx3]
vec3=v(ind(i),:)-v(ind(i+1),:);   % [nx3]

cp=cross(vec1,vec2); % [nx3]

% grads=cp/(vec3*cp');  %row vector  [nx3] / [nx3]*[3xn]
grads=cp./repmat(dot(vec3,cp,2),1,3);  %row vector  [nx3]./ [nx3]

end

