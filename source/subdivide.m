function [ out ] = subdivide(vertices, d_at_vertices, h )
%subdivide returns a list of sub-tetrahedra of T that are inside D_h
%   Subdivide takes in a tetrahedron element T denoted by its 4 vertices
%   and returns a list of subdivided tetrahedron elements that are 
%   inside the band D_h = { x | |d_h(x)| < h }.  Since we are assuming
%   d_h(x) = I_h d(x) is the linear interpolant of the distance function, 
%   it is sufficient to give the 4 values di = d_h(vi) to determine where 
%   the cuts need to be made.
%
%                                          [vi_x]
%   vertices = [v1, v2, v3, v4]  with vi = [vi_y] column vectors
%                                          [vi_z]
%
%   d_at_vertices = [d1, d2, d3, d4]  with di = d(vi) distance value
%
%   h > 0  width of band 
%
%
%
% (4out 0edge 0in)  returns 0 sub-tetrahedra
% (3out 1edge 0in)  returns 0 sub-tetrahedra
% (3out 0edge 1in)  returns 1 sub-tetrahedra
% (2out 2edge 0in)  returns 0 sub-tetrahedra
% (2out 1edge 1in)  returns 1 sub-tetrahedra
% (2out 0edge 2in)  returns 3 sub-tetrahedra
% (1out 3edge 0in)  returns 0 sub-tetrahedra
% (1out 2edge 1in)  returns 1 sub-tetrahedra
% (1out 1edge 2in)  returns 2 sub-tetrahedra
% (1out 0edge 3in)  returns 3 sub-tetrahedra
% (0out 3edge 1in)  returns 1 sub-tetrahedra
% (0out 2edge 2in)  returns 1 sub-tetrahedra
% (0out 1edge 3in)  returns 1 sub-tetrahedra
% (0out 0edge 4in)  returns 1 sub-tetrahedra
%
% Spencer Patty
% Dec 23, 2015

%                                         
% reorder vertices so that they are listed ( in | edge | out )
[abs_d_at_xi, I] = sort(abs(d_at_vertices));
xi = vertices(:,I);  di = d_at_vertices(I);


% calculate which case we are in
num_in  = sum( (abs_d_at_xi < h) );
num_out = sum( (abs_d_at_xi > h) );
% num_edge = 4 - num_in - num_out;

% We proceed through all the possibilities of vertices (out edge in)
% and construct the correct number of subdivided cells.  Then we store
% them in a cellArray called out.  It can be referenced by out{i} gives the
% array of vertices of the ith subdivided tetrahedron. For example,
% out{i}(:,3) would be the z coordinates of the ith sub-tetrahedron.

if (num_in == 0)  % covers (4out 0edge 0in) and (3out 1edge 0in)
    out = {};
    return;
end

if (num_out == 3) %(3out 0edge 1in)
    x12 = find_xij(xi(:,1),xi(:,2), di(1), di(2), h); 
    x13 = find_xij(xi(:,1),xi(:,3), di(1), di(3), h);
    x14 = find_xij(xi(:,1),xi(:,4), di(1), di(4), h);
    out = cell(1,1);
    out(1) = {[xi(:,1), x12, x13, x14]};
    return;
end

if (num_out == 2)
    x13 = find_xij(xi(:,1),xi(:,3), di(1), di(3), h);
    x14 = find_xij(xi(:,1),xi(:,4), di(1), di(4), h);
    if (num_in == 1) %(2out 1edge 1in)
        out = cell(1,1);
        out(1) = {[xi(:,1), xi(:,2), x13, x14]};
        return;
    else % (2out 0edge 2in)
        x23 = find_xij(xi(:,2),xi(:,3), di(2), di(3), h);
        x24 = find_xij(xi(:,2),xi(:,4), di(2), di(4), h);
        out = cell(1,3);
        out(1) = {[xi(:,1), x13,     x14, x23]};
        out(2) = {[xi(:,1), xi(:,2), x14, x23]};
        out(3) = {[x24,     xi(:,2), x14, x23]};
        return;
    end
end

if (num_out == 1)
    if(num_in == 3) %(1out 0edge 3in)
        x14 = find_xij(xi(:,1),xi(:,4), di(1), di(4), h);
        x24 = find_xij(xi(:,2),xi(:,4), di(2), di(4), h);
        x34 = find_xij(xi(:,3),xi(:,4), di(3), di(4), h);
        out = cell(1,3);
        out(1) = {[xi(:,1), xi(:,2), xi(:,3),  x34]};
        out(2) = {[xi(:,1), xi(:,2), x14,      x34]};
        out(3) = {[x24,     xi(:,2), x14,      x34]};
        return;
    elseif(num_in == 2) %(1out 1edge 2in)
        x14 = find_xij(xi(:,1),xi(:,4), di(1), di(4), h);
        x24 = find_xij(xi(:,2),xi(:,4), di(2), di(4), h);
        out = cell(1,2);
        out(1) = {[xi(:,1), x24, xi(:,3), x14    ]};
        out(2) = {[xi(:,1), x24, xi(:,2), xi(:,3)]};
        return;
    else %(1out 2edge 1in)
        x14 = find_xij(xi(:,1),xi(:,4), di(1), di(4), h);
        out = cell(1,1);
        out(1) = {[xi(:,1), xi(:,2), xi(:,3),  x14]};
        return;
    end
end

if (num_out == 0) %(0out 0edge 4in), (0out 1edge 3in), (0out 2edge 2in), (0out 3edge 1in)
    out = cell(1,1);
    out(1) = {xi};
    return;
end

end

