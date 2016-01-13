function [ p ] = function_extension( x )
% The pull back function 
% which map the a point, x, in 3D space 
% back to the corresponding point on the surface.
% For a distance function d(x), this mapping is
% 
% p = x - d(x)*\nabla d(x)
% 
% where \nabla d is the gradient of the distant function
% (note that |\nabla d|=1)
%
%  input:  
%     x = [nx3] vector of n points in R^3
%  output:
%     p = [nx3] vector of n points on Torus

R=1.0; r=0.6;
value1 = sqrt(x(:,1).^2+x(:,2).^2)+1e-10;
value2 = value1-R;
value3 = sqrt(value2.^2+x(:,3).^2);
dist_value = value3 - r;
value4 = value2./(value1.*value3);

dist_grad=zeros(size(x));
dist_grad(:,1)=x(:,1).*value4;
dist_grad(:,2)=x(:,2).*value4;
dist_grad(:,3)=x(:,3)./value3;

n = length(dist_value);
p = x-spdiags(dist_value,0,n,n)*dist_grad;  % want nx3 matrix result so diagonalize
                                   % the column vector dist_value to
                                   % preserve proper sizes in multiplication.

end
