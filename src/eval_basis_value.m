function vals = eval_basis_value (x,i,v)
% Return ith shape value (linear element) on the point x.
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
%   vals = [nx1] vector of basis function i evaluated at x.
%

ind=[1,2,3,4,1,2,3];
vec1=v(ind(i+2),:)-v(ind(i+1),:); % [nx3]
vec2=v(ind(i+3),:)-v(ind(i+1),:); % [nx3]
vec3=v(ind(i),:)-v(ind(i+1),:);   % [nx3]
vecx=repmat(x,length(i),1)-v(ind(i+1),:); % [nx3]

cp=cross(vec1,vec2); % [nx3]

% dot product along second dimension to get out [nx1] shape instead of
% [1x3] shape.
vals = dot(vecx,cp,2)./dot(vec3,cp,2); % [nx1] ./ [nx1]

end