function [ res ] = function_extension( p )
% The pull back function 
% which map the a point in 3D space 
% back to the corresponding point a the surface.
% For a distance function d(x), this mapping is
% 
% res = p - d(p)*\nabla d(p)
% 
% where \nabla d is the gradient of the distant function
% (note that |\nabla d|=1)

R=1.0; r=0.6;
value1 = sqrt(p(1)^2+p(2)^2);
value2 = value1-R;
value3 = sqrt(value2^2+p(3)^2);
dist_value = value3 - r;
value4 = value2/(value1*value3);
dist_grad=zeros(1,3);
dist_grad(1)=p(1)*value4;
dist_grad(2)=p(2)*value4;
dist_grad(3)=p(3)/value3;

res = p-dist_value*dist_grad;
end

