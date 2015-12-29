function [ y ] = get_parameterization( x )
% Return the parametric coordinates for 3D point x on the surface
r=0.6; R=1.0;
y=zeros(1,2);
if (x(1)<0 && x(2)>0)
    y(2)=pi+atan(x(2)/x(1));
else if (x(1)<0 && x(2)<0)
        y(2)=-pi+atan(x(2)/x(1));
    else
        y(2)=atan(x(2)/x(1));
    end
end
rad = sqrt(x(1)^2+x(2)^2);
if (rad<R && x(3)>0)
    y(1)=pi-asin(x(3)/r);
else if (rad<R && x(3)<0);
        y(1)=-pi-asin(x(3)/r);
    else
        y(1)=asin(x(3)/r);
    end
end

end

