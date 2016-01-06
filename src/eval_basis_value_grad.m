function [vals, grads] = eval_basis_value_grad (x,i,v)
% Return ith shape value and gradient (linear element) on the point x.
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
%   vals  = [nx1] vector of basis function i evaluated at x.
%   grads = [nx3] vector of basis function gradients i evaluated at x.

ind=[1,2,3,4,1,2,3];
vec1=v(ind(i+2),:)-v(ind(i+1),:); % [nx3]
vec2=v(ind(i+3),:)-v(ind(i+1),:); % [nx3]
vec3=v(ind(i),:)-v(ind(i+1),:);   % [nx3]
vecx=repmat(x,length(i),1)-v(ind(i+1),:); % [nx3]

cp=cross(vec1,vec2,2); % [nx3]  2nd dimension is where cross is applied

% We dot product along second dimension to get out [nx1] shape instead of
% [1x3] shape.  Cross product is the same, especially since cross doesn't
% make sense for anything other than vectors of dimension 3.

% only do this dot product once
v3dotcp = dot(vec3,cp,2);

% construct solutions
vals = dot(vecx,cp,2)./v3dotcp;  % [nx1] ./ [nx1]
grads = cp./repmat(v3dotcp,1,3); % [nx3] ./ [nx3]

end