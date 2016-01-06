clear all; clc; close all;

h = 0.03125;

n = min(floor(4/h),64);
h_tmp = 4/n;
[n_node,n_ele,node,elem ] = uniform_bulk_mesh( n, h_tmp ,0);

dist = abs(distfunc(node));
hmesh = (dist < 1.2*h)*h/2 + ( dist >= 1.2*h).*(h/2 + ((dist-1.2*h))/6);
clear dist;

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

clear n_node n_elem node elem;

% copyfile([preamble, '.b.node'] ,[preamble, '.node']);
% copyfile([preamble, '.b.ele'] ,[preamble, '.ele']);
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

clear n_node2 n_elem2 node2 elem2;

p0 = [-2,-2,-2];
p1 = [2,2,2];
maxradiustoedge = 1.414;
minangle = 10;
maxvol = h^3;
nodesize = 1;

% ISO2MESH_TETGENOPT = ['-A -q' num2str(maxradiustoedge) '/' num2str(minangle) ...
%                        '-m -V']; 

ISO2MESH_TETGENOPT = ['-r -q' num2str(maxradiustoedge) '/' num2str(minangle) ' -m -V ' preamble];

[node,elem,~]=surf2mesh([],[],p0,p1,1.0,maxvol,[],[],nodesize,1);
elem=elem(:,1:4);
plotmesh(node,elem,'x>0'); axis equal;


%
% cut down to narrow band
%
fprintf('...Extracting the Narrow Band...\n');
n_ele = length (elem);
n_local_dofs = 4;

% Count the number of cells which intersects the narrow band
narrow_band_counter=0;
for cell = 1:n_ele
    for v= 1:n_local_dofs
        interp_dist=distfunc(node(elem(cell,v),:));
        if ((interp_dist)<(h) && (interp_dist)>(-h))
            narrow_band_counter=narrow_band_counter+1;
            break;
        end
    end
end

n_nb_ele = narrow_band_counter;
nb_elem=zeros(n_nb_ele,4);

% Get the connectivity list of the narrow band
narrow_band_counter=1;
for cell = 1:n_ele
    for v= 1:n_local_dofs
        interp_dist=distfunc(node(elem(cell,v),:));
        if ((interp_dist)<(h) && (interp_dist)>(-h))
            nb_elem(narrow_band_counter,:)=elem(cell,:);
            narrow_band_counter=narrow_band_counter+1;
            break;
        end
    end
end

% Generate the index set for the vertices of cells which
% intersect the narrow band
node_nb_flag =zeros(1,length(node));
for cell = 1:n_nb_ele
    for v = 1:4
        node_nb_flag(nb_elem(cell,v))=1;
    end
end
n_nb_node=length(nonzeros(node_nb_flag));

% Extract the relavant vertices
nb_node=zeros(n_nb_node,3);
count=1;
for i = 1:length(node)
    if(node_nb_flag(i))
        nb_node(count,:)=node(i,:);
        node_nb_flag(i)=count;
        count=count+1;
    end
end

clear node;

% Update the connectivity list with extracted vertex list.
for cell = 1:n_nb_ele
    for v=1:4
        nb_elem(cell,v)=node_nb_flag(nb_elem(cell,v));
    end
end

clear elem

figure;
plotmesh(nb_node,nb_elem,'x>0'); axis equal;
