function [out] = subdivide_error(vertices, d_at_vertices)
%
%
% (4p 0gamma 0m)  returns 0 faces
% (3p 1gamma 0m)  returns 0 faces
% (3p 0gamma 1m)  returns 1 faces
% (2p 2gamma 0m)  returns 0 faces
% (2p 1gamma 1m)  returns 1 faces
% (2p 0gamma 2m)  returns 2 faces
% (1p 3gamma 0m)  returns 1 faces
% (1p 2gamma 1m)  returns 1 faces
% (1p 1gamma 2m)  returns 1 faces
% (1p 0gamma 3m)  returns 1 faces
% (0p 3gamma 1m)  returns 1 faces
% (0p 2gamma 2m)  returns 0 faces
% (0p 1gamma 3m)  returns 0 faces
% (0p 0gamma 4m)  returns 0 faces
%

% reorder vertices so that they are listed ( m | gamma | p )
[d_at_xi, I] = sort(d_at_vertices);
xi = vertices(:,I);  di = d_at_vertices(I);

% calculate which case we are m or p
num_m = sum( (d_at_xi < 0) );
num_p = sum( (d_at_xi > 0) );
% num_gamma = 4 - num_m - num_p;

out = {};


if (num_p == 4)  % (4p 0gamma 0m)
    out = {};
    return;
end

if (num_p == 3) 
    if (num_m == 1) % (3p 0gamma 1m)
        x12 = find_xij_gamma(xi(:,1),xi(:,2), di(1), di(2));
        x13 = find_xij_gamma(xi(:,1),xi(:,3), di(1), di(3));
        x14 = find_xij_gamma(xi(:,1),xi(:,4), di(1), di(4));
        out = cell(1,1);
        out(1) = {[x12, x13, x14]};
        return;
    else % (3p 1gamma 0m)
        return;
    end
end

if (num_p == 2)
    if (num_m == 0) % (2p 2gamma 0m)
        return;
    elseif (num_m == 1) % (2p 1gamma 1m)
        x13 = find_xij_gamma(xi(:,1),xi(:,3), di(1), di(3));
        x14 = find_xij_gamma(xi(:,1),xi(:,4), di(1), di(4));
        out = cell(1,1);
        out(1) = {[xi(:,2), x13, x14]};
        return;
    else % (2p 0gamma 2m)
        x13 = find_xij_gamma(xi(:,1),xi(:,3), di(1), di(3));
        x14 = find_xij_gamma(xi(:,1),xi(:,4), di(1), di(4));
        x23 = find_xij_gamma(xi(:,2),xi(:,3), di(2), di(3));
        x24 = find_xij_gamma(xi(:,2),xi(:,4), di(2), di(4));
        out = cell(1,2);
        out(1) = {[x13, x14, x23]};
        out(2) = {[x14, x23, x24]};
        return;
    end
end

if (num_p == 1)
    if(num_m == 3) % (1p 0gamma 3m)
        x14 = find_xij_gamma(xi(:,1),xi(:,4), di(1), di(4));
        x24 = find_xij_gamma(xi(:,2),xi(:,4), di(2), di(4));
        x34 = find_xij_gamma(xi(:,3),xi(:,4), di(3), di(4));
        out = cell(1,1);
        out(1) = {[x14, x24, x34]};
        return;
    elseif(num_m == 2) % (1p 1gamma 2m)
        x14 = find_xij_gamma(xi(:,1),xi(:,4), di(1), di(4));
        x24 = find_xij_gamma(xi(:,2),xi(:,4), di(2), di(4));
        out = cell(1,1);
        out(1) = {[xi(:,3), x14, x24 ]};
        return;
    else % (1p 2gamma 1m)
        x14 = find_xij_gamma(xi(:,1),xi(:,4), di(1), di(4));
        out = cell(1,1);
        out(1) = {[xi(:,2), xi(:,3), x14]};
        return;
    end
end

if (num_p == 0)
    if (num_m == 1) % (0p 3gamma 1m)
        out = cell(1,1);
        out(1) = {xi(:,2:4)};
        return;
    else % (0p 2gamma 2m), (0p 1gamma 3m), (0p 0gamma 4m) 
        return
    end
end

end

