
This section mainly shows the implementation of the standard finite element method in 3D. Note that we assemble the stiffness and mass matrices directly on the real space (without using reference cell and affine mapping). The main script is 

FEM_3D_linear.m

which calls the other functions to mesh, assemble, solve on the domain in 3D.

We approximate the solution of the problem

u-\Delta u = f in \Omega=[0,1]^3
     du/dn = 0 on the boundary of \Omega

where \Delta is the Laplace operator and du/dn is the normal derivative on the boundary.
Let the solution u to be cos(\pi x)cos(\pi y)cos(\pi z). So corresponding right hand side function f should be (3\pi^2+1)cos(\pi x)cos(\pi y)cos(\pi z).

Wenyu Lei
Spencer Patty
Jan 7, 2016