function res = exact (chi)
% Exact solution.
res=sin(3*chi(:,2)).*cos(3*chi(:,1)+chi(:,2));
end