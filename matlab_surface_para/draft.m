
% Generate the mesh
n=8;
h=2*pi/n;
xx=-pi:h:pi;
yy=-pi:h:pi;
X=zeros((n+1)*(n+1),2);
for i=1:(n+1)
    for j=1:n+1
        X(j+(n+1)*(i-1),1)=xx(i);
        X(j+(n+1)*(i-1),2)=yy(j);
    end
end
dt = delaunayTriangulation(X);
triplot(dt);
