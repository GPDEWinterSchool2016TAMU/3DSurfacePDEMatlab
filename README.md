# 3DSurfacePDEMatlab
A finite element code for solving the Laplace-Beltrami operator on a 2D surface embedded in 3D.

There are two approaches we implement:  

1) parameterize the surface then solve normally on a surface mesh.

2) define the surface implicitly as the zero levelset of a signed distance function then solve an adjusted problem on a bulk mesh in a narrow band around the zero levelset.
