function [ out ] = subdivide(vertices, d_at_vertices, h )
%subdivide returns a list of sub-tetrahedra of T that are inside D_h
%   Subdivide takes in a tetrahedron element T denoted by its 4 vertices
%   and returns a list of subdivided tetrahedron elements that are 
%   inside the band D_h = { x | |d_h(x)| < h }.  Since we are assuming
%   d_h(x) = I_h d(x) is the linear interpolant of the distance function, 
%   it is sufficient to give the 4 values di = d_h(vi) to determine where 
%   the cuts need to be made.
%
%   vertices = [v1; v2; v3; v4]  with vi = [vi_x, vi_y, vi_z] row vectors
%
%   d_at_vertices = [d1, d2, d3, d4]  with di = d_h(vi) distance value
%
%   h > 0  width of band D_h = { x | |d_h(x)| < h }
%
%   out is a cell array of size 1xm where there are m subtetrahedra 
%   returned ( = {} if m=0 ).  Each element of the cell array is a set of
%   four vertices of the subtetrahedra, [sv1; sv2; sv3; sv4].  For example
%   out{i}(3,:) would be the 3rd vertex of the ith subtetrahedron.
%
% Systematically checks the following cases:
%
% (0in 0edge 4out)  returns 0 sub-tetrahedra
% (0in 1edge 3out)  returns 0 sub-tetrahedra
% (1in 0edge 3out)  returns 1 sub-tetrahedra
% (0in 2edge 2out)  returns 0 sub-tetrahedra
% (1in 1edge 2out)  returns 1 sub-tetrahedra
% (2in 0edge 2out)  returns 3 sub-tetrahedra
% (0in 3edge 1out)  returns 0 sub-tetrahedra
% (1in 2edge 1out)  returns 1 sub-tetrahedra
% (2in 1edge 1out)  returns 2 sub-tetrahedra
% (3in 0edge 1out)  returns 3 sub-tetrahedra
% (1in 3edge 0out)  returns 1 sub-tetrahedra
% (2in 2edge 0out)  returns 1 sub-tetrahedra
% (3in 1edge 0out)  returns 1 sub-tetrahedra
% (4in 0edge 0out)  returns 1 sub-tetrahedra
%
%  Some remarks must be made about the quality of mesh that we have here.
%  If the mesh is not sufficiently regular, it could be possible to end up
%  with vertices on the outside of both sides of the band D_h. We are 
%  assuming here that this does not happen and that our original mesh
%  was created with a sufficiently large minimal angle condition so that 
%  in conjunction with the area constraint |T| ~ h^3, we have tetrahedra
%  that have edges roughly of length h and not longer than 2h.
%
%  It is also important to note that even though the original mesh has
%  regularity properties, the subdivided tetrahedra do not.  D_h could
%  intersect T in many ways that leave just a sliver.  The resulting mass
%  or stiffness matrices will undoubtedly be poorly conditioned since they
%  will be integrated over these subdivided tetrahedra.  
%
%  Finally,
%  In a few cases we end up with a prism shape which needs to be subdivided
%  into 3 tetrahedra.  We use the following pattern of subdivision although
%  the final list of 4 vertices does not have a specific order. (We do not
%  fill in all the lines of each subtetrahedra due to graphical
%  difficulties :) but you get the idea)
%
%   1-----3      1            1            1-----3
%   |\   /|      |             \            \   /|
%   | \ / |  =>  |              \            \ / |
%   |  2  |      |               2            2  |
%   |  |  |      |               |               |
%   4  |  6      4-----6         |  6            6
%    \ | /        \   /          | /        
%     \|/          \ /           |/         
%      5            5            5          
%               [1 4 5 6]    [1 2 5 6]    [1 2 3 6]
%
% Spencer Patty
% Dec 23, 2015

%                                         
% reorder vertices so that they are listed ( in | edge | out )
[abs_d_at_xi, I] = sort(abs(d_at_vertices));
xi = vertices(I,:);  di = d_at_vertices(I);

% calculate which case we are in
num_in  = sum( (abs_d_at_xi < h) );
num_out = sum( (abs_d_at_xi > h) );
% num_edge = 4 - num_in - num_out;

% We proceed through all the possibilities of vertices (out edge in)
% and construct the correct number of subdivided cells.

if (num_in == 0)  %(0in 0edge 4out), (0in 1edge 3out), (0in 2edge 2out), (0in 3edge 1out)
    out = {};
    return;
end

if (num_out == 0) %(4in 0edge 0out), (3in 1edge 0out), (2in 2edge 0out), (1in 3edge 0out)
    out = cell(1,1);
    out(1) = {xi};
    return;
end

if (num_out == 3) %(1in 0edge 3out)
    x12 = find_xij(xi(1,:),xi(2,:), di(1), di(2), h); 
    x13 = find_xij(xi(1,:),xi(3,:), di(1), di(3), h);
    x14 = find_xij(xi(1,:),xi(4,:), di(1), di(4), h);
    out = cell(1,1);
    out(1) = {[xi(1,:); x12; x13; x14]};
    return;
end

if (num_out == 2)
    x13 = find_xij(xi(1,:),xi(3,:), di(1), di(3), h);
    x14 = find_xij(xi(1,:),xi(4,:), di(1), di(4), h);
    if (num_in == 1) %(1in 1edge 2out)
        out = cell(1,1);
        out(1) = {[xi(1,:); xi(2,:); x13; x14]};
        return;
    else % (2in 0edge 2out)
        x23 = find_xij(xi(2,:),xi(3,:), di(2), di(3), h);
        x24 = find_xij(xi(2,:),xi(4,:), di(2), di(4), h);
        out = cell(1,3);
        out(1) = {[xi(1,:); x13;     x14; x23]};
        out(2) = {[xi(1,:); xi(2,:); x14; x23]};
        out(3) = {[x24;     xi(2,:); x14; x23]};
        return;
    end
end

if (num_out == 1)
    if(num_in == 3) %(3in 0edge 1out)
        x14 = find_xij(xi(1,:),xi(4,:), di(1), di(4), h);
        x24 = find_xij(xi(2,:),xi(4,:), di(2), di(4), h);
        x34 = find_xij(xi(3,:),xi(4,:), di(3), di(4), h);
        out = cell(1,3);
        out(1) = {[xi(1,:); xi(2,:); xi(3,:);  x34]};
        out(2) = {[xi(1,:); xi(2,:); x14;      x34]};
        out(3) = {[x24;     xi(2,:); x14;      x34]};
        return;
    elseif(num_in == 2) %(2in 1edge 1out)
        x14 = find_xij(xi(1,:),xi(4,:), di(1), di(4), h);
        x24 = find_xij(xi(2,:),xi(4,:), di(2), di(4), h);
        out = cell(1,2);
        out(1) = {[xi(1,:); x24; xi(3,:); x14    ]};
        out(2) = {[xi(1,:); x24; xi(2,:); xi(3,:)]};
        return;
    else %(1in 2edge 1out)
        x14 = find_xij(xi(1,:),xi(4,:), di(1), di(4), h);
        out = cell(1,1);
        out(1) = {[xi(1,:); xi(2,:); xi(3,:);  x14]};
        return;
    end
end


end

