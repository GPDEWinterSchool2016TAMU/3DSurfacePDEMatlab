function [n_node,n_elem,node,elem] =bulk_mesh_generator_adaptive(h, return_nb)
% Generate 3D adaptive mesh for the narrow band method using Tetgen and
% Iso2mesh software.
%
% input:
%   h: constraints for narrow band,  length of an edge
%   return_nb: whether or not to extract narrow band and return instead of
%              full mesh
% output: 
%   n_node: number of vertices of the narrow band or full mesh.
%   n_ele: number of elements (cells) of the narrow band or full mesh.
%   node: vertex list of the narrow band or full mesh.
%   elem: connectivity list of the narror band or full mesh.
%
% We first generate a mesh for a box which contains the torus.
% Then select the cells intersects the narrow band
%
% D_h := {x \in U | I_h \Phi(x) < h}
%
% Here \Phi(x) is the distance function (Or level set function)
% and I_h is the lagrange interpolation on to continuous
% piecewise linear function space.
% To do this, we evaluate distance function on each vertices
% of the cell. We select this cell if one of the values is 
% between -h and h.
%
% Spencer Patty
% Jan 4, 2016

if (nargin < 2)
    return_nb = 1;  % true
end

%
% Construct background mesh which will define the adaptive shape of mesh.
% At each node of background mesh, we choose hmesh value which says how
% long the edges in the neighborhood around that point should be roughly.
% It is interpolated as the refinement continues by Tetgen engine.
%
n = min(floor(4/h),64);
h_tmp = 4/n;
[n_node,n_ele,node,elem ] = uniform_bulk_mesh( n, h_tmp ,0);

dist = abs(distfunc(node));
hmesh = (dist < 1.2*h)*h/2 + ( dist >= 1.2*h).*(h/2 + (dist-1.2*h)/5);

ISO2MESH_TEMP = '/Users/srobertp/software/geometric_pdes_matlab/tmp_iso2mesh_files';
ISO2MESH_SESSION = 'srp_';
preamble = [ISO2MESH_TEMP '/' ISO2MESH_SESSION 'post_vmesh'];

FID = fopen([preamble '.b.node'],'w');
fprintf(FID,'# Node count, 3 dim, no attribute, no boundary marker\n');
fprintf(FID,'%d 3 0 0\n', n_node);
fprintf(FID,'# Node index, node coordinates\n');
fprintf(FID,'%d %f %f %f\n', [(1:n_node)', node]' );
fclose(FID);

FID = fopen([preamble '.b.ele'],'w');
fprintf(FID,'# <# of tetrahedra>,  <nodes per tet. (4 or 10)>, <region attribute (0 or 1)>\n');
fprintf(FID,'%d 4 0\n', n_ele);
fprintf(FID,'#<tetrahedron #>, <node> <node> ... <node> [attribute]\n');
fprintf(FID,'%d %d %d %d %d\n', [(1:n_ele)' elem]' );
fclose(FID);

FID = fopen([preamble '.b.mtr'],'w');
fprintf(FID,'# <# of nodes> <size of metric (always 1)>\n');
fprintf(FID,'%d 1\n', n_node);
fprintf(FID,'#<value>\n');
fprintf(FID,'%f\n', hmesh' );
fclose(FID);


%
% Create initial mesh which will be refined using -r flag
%
n2 = 4;
h_tmp2 = 4/n2;
[n_node2,n_ele2,node2,elem2 ] = uniform_bulk_mesh( n2, h_tmp2 ,0);

FID = fopen([preamble '.node'],'w');
fprintf(FID,'# Node count, 3 dim, no attribute, no boundary marker\n');
fprintf(FID,'%d 3 0 0\n', n_node2);
fprintf(FID,'# Node index, node coordinates\n');
fprintf(FID,'%d %f %f %f\n', [(1:n_node2)', node2]' );
fclose(FID);

FID = fopen([preamble '.ele'],'w');
fprintf(FID,'# <# of tetrahedra>,  <nodes per tet. (4 or 10)>, <region attribute (0 or 1)>\n');
fprintf(FID,'%d 4 0\n', n_ele2);
fprintf(FID,'#<tetrahedron #>, <node> <node> ... <node> [attribute]\n');
fprintf(FID,'%d %d %d %d %d\n', [(1:n_ele2)' elem2]' );
fclose(FID);


%
% Construct actual adaptive mesh 
%

p0 = [-2,-2,-2];
p1 = [2,2,2];
maxradiustoedge = 1.414;
minangle = 20;
maxvol = h^3;
nodesize = 1;

% ISO2MESH_TETGENOPT = ['-A -q' num2str(maxradiustoedge) '/' num2str(minangle) ...
%                        '-m -V']; 

ISO2MESH_TETGENOPT = ['-r -q' num2str(maxradiustoedge) '/' num2str(minangle) '-m -V ' preamble];

[node,elem]=surf2mesh([],[],p0,p1,1,maxvol,[],[],nodesize);
elem=elem(:,1:4);
plotmesh(node,elem,'x>0'); axis equal;

%
% cut down to narrow band if desired
%
if (return_nb)
    [n_node,n_elem,node,elem] = extract_narrow_band(node,elem,h);
else % return full mesh
    n_node = size(node,1);
    n_elem = size(elem, 1);
end

% Plot the narrow band mesh for part x>0. 
% Here the function 'plotmesh' is from ISO2MESH package.
plotmesh(node,elem,'x>0'); axis equal;

end