function res = AFFINE_MAPPING( v1,v2,v3,x )

mat_B = [v2(1)-v1(1), v2(2)-v1(2);v3(1)-v1(1),v3(2)-v1(2)];
res =x*mat_B+v1;

end

