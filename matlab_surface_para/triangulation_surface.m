function [ n_node,n_ele,node,ele,global_ind ] = triangulation_surface( n )
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
%triplot(dt);
n_ele=dt.size(1);
ele=dt.ConnectivityList;
node=dt.Points;
n_node=length(node);

global_ind=zeros(1,(n+1)*(n+1));
for i=1:n
    for j=1:n
        global_ind(j+(i-1)*(n+1))=j+(i-1)*n;
    end
end

for i=1:n
    global_ind((n+1)*i)=(i-1)*(n)+1;
    global_ind(n*(n+1)+i)=i;
end

global_ind((n+1)*(n+1))=1;
end
