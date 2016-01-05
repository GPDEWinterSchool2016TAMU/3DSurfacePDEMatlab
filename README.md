# 3DSurfacePDEMatlab
A finite element code for solving the Laplace-Beltrami operator on a closed surface, Gamma, in 3D.

There are three approaches we implement:

1) parameterize the surface then solve normally on a surface mesh.

2) define the surface implicitly as the zero levelset of a signed distance function then solve an adjusted problem on a bulk mesh in a narrow band around the zero levelset. This is called the Narrow Banded Method (NBM) by Deckelnick, et al. in (2014, Unfitted finite element methods using bulk meshes for surface partial differential equations)

3) define the surface implicitly as a levelset of a phase field function, then solve the adjusted problem on a bulk mesh as described by Martin Burger in (2008, Finite element approximation of elliptic partial differential equations on implicit surfaces)
