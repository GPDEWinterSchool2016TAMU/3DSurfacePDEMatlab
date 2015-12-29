clear all; clc;
format long;
n=32;h=4./n; 

% Mesh generation
%[n_nb_node,n_nb_ele,nb_node,nb_ele]=bulk_mesh_generator(h);
[n_nb_node,n_nb_ele,nb_node,nb_ele]=uniform_bulk_mesh(n,h);

% Matrix and vector initialization
A = sparse([],[],[],n_nb_node,n_nb_node,9*n_nb_ele);
F = zeros(n_nb_node,1);

% Assembling
for cell = 1:n_nb_ele
    vertices = zeros(4,3);
    dist_value = zeros(1,4);
    for v = 1:4
        vertices(v,:)=nb_node(nb_ele(cell,v),:);
        dist_value(v)=distfunc(vertices(v,:));
    end
    [ lstiff,lrhs ] = local_assembling( vertices,dist_value,h );
    for j = 1:4
        for k = 1:4
            A(nb_ele(cell,j),nb_ele(cell,k))=A(nb_ele(cell,j),nb_ele(cell,k))+lstiff(k,j);
        end
        F(nb_ele(cell,j))=F(nb_ele(cell,j))+lrhs(j);
    end
end

% Apply back slash solver
sol=A\F;

% Compute the L2 error i.e.
% \|u^e - u_h\|_{L^2(\Gamma_h)}
% Here \Gamma_h is the discrete surface defined by the level set function
% \Gamma_h := {x\in box mesh | I_h\Phi(x)=0}
% and u^e is extension of the solution.
% We compute this in a similar way as we did in the local assmebling.
% That is, we compute (u^e - u_h)^2 on the intersection between a cell and 
% \Gamma_h. The intersection may triangle or quadrilateral. So we cut
% the quadrilateral by two triangle and apply quadrature rule only on 
% triangles.
err_square = 0;
for cell = 1:n_nb_ele
    % We first get the subdivision of \Gamma_h on this cell.
    vertices = zeros(4,3);
    dist_value = zeros(1,4);
    ver_sol=zeros(1,4);
    for v = 1:4
        vertices(v,:)=nb_node(nb_ele(cell,v),:);
        dist_value(v)=distfunc(vertices(v,:));
        ver_sol(v)=sol(nb_ele(cell,v));
    end
    subd = subdivide_error(vertices',dist_value);
    n_subdivide= length(subd);
    % if we do have subdivisions, then we apply quadature rule on each
    % subdivision.
    if (n_subdivide ~= 0)
        for i = 1:n_subdivide
            % get quadrature info
            points_list=subd{i}';
            [nq,q,w]=facequad(points_list,2);
            for q_point = 1:nq
                shape_value=zeros(1,4);
                for j = 1:4
                    shape_value(j)=eval_basis_value(q(q_point,:),j,vertices);               
                end
                err_square=err_square...
                    +(exact(function_extension(q(q_point,:)))-shape_value*ver_sol')^2.*w(q_point);
            end
        end
    end
end
err=sqrt(err_square)