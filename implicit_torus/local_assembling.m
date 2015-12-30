function [ lstiff,lrhs ] = local_assembling( v,d,h )
% Local assembling routine (linear element)
% for a cell intersects with the narrow band with size h
%
% input:
%   v: vertex list (4 by 4 matrix)
%   d: distance function evaluation on each vertex (1 by 4 vector)
%   h: size for narrow band (usually mesh size)
% output:
%   lstiff: local stiffness matrix (4 by 4 for linear case)
%   lrhs: local right hand side (1 by 4 for linear case)
%
% We want to compute: for each cell \tao
%
% lstiff(j,k)= int_{\tao\cap D_h}(\nabla\phi_k\cdot\nabla\phi_j
%                                 +\phi_k\phi_j)|\nabla I_h\Phi| dx
% and
%
% lrhs(j)= int_{\tao\cap D_h} f^e \phi_j|\nabla I_h\Phi| dx
% 
% where D_h is the narrow band, phi_j are shape functions,
% I_h\Phi is the interpolation of the level set function and
% f^e is the extension of the right hand side function.
% To do this, we apply quadrature rule for each subdivision of 
% the integral domain \tao\cap D_h.
%
% Wenyu Lei
% Dec 29, 2015

% Initialization
lstiff=zeros(4,4);
lrhs=zeros(1,4);

% Get the subdivision of the integral domain.
subd = subdivide(v',d,h);
n_subdivide= length(subd);

for i = 1:n_subdivide
    % Get the quadarature info for each subdivision.
    points_list=subd{i}';
    [nq,q,w]=volquad(points_list,2);
    
    % Add the increament for local stiffenss and right hand side on 
    % each quadrature point.
    for q_point = 1:nq
        shape_value=zeros(1,4);
        shape_grad=zeros(4,3);
        for j = 1:4
            shape_value(j)=eval_basis_value(q(q_point,:),j,v);
            shape_grad(j,:)=eval_basis_grad(q(q_point,:),j,v);
        end
        gls_interp=grad_interp_levelset(shape_grad,d);
        for j = 1:4
            for k = 1:4
                lstiff(j,k)=lstiff(j,k)...
                    +(shape_grad(j,:)*shape_grad(k,:)'+shape_value(j)*shape_value(k))...
                    *gls_interp*w(q_point);
            end
            lrhs(j)=lrhs(j)+rhs_eval(function_extension(q(q_point,:)))...
                *shape_value(j)*gls_interp*w(q_point);
        end
    end
end

end

