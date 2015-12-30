function [ res ] = eval_basis_grad(x,i,v )
% Return ith shape gradient (linear element) on the point x.
% Note that the input argument v is a 4 by 3 matrix
% and jth row of v indicates the coordinates of jth vertex
% in this cell.
ind=[1,2,3,4,1,2,3];
vec1=v(ind(i+2),:)-v(ind(i+1),:);
vec2=v(ind(i+3),:)-v(ind(i+1),:);
vec3=v(ind(i),:)-v(ind(i+1),:);


cp=cross(vec1,vec2);

res=cp/(vec3*cp');  %row vector

end

