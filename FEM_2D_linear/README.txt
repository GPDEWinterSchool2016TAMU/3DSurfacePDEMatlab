
This section mainly shows the implementation of the standard finite element method in 2D using the reference cell and the affine mapping. The main script is 

FEM_2D_linear.m

which calls the other functions to mesh, assemble, solve on the domain in 2D.

We approximate the solution of the problem

2u-\Delta u = f in \Omega=[0,1]^2
      du/dn = 0 on the boundary of \Omega

where \Delta is the Laplace operator and du/dn is the normal derivative on the boundary.
Let the solution u to be cos(\pi x)cos(2\pi y). So corresponding right hand side function f should be (5\pi^2+2)cos(\pi x)cos(2\pi y).

About error analysis:
Let u_h to be the computation result, then u_h is the linear comblination of global shape functions (\phi), that is

u_h = \sum c_i \phi_i

We also denote I_h u to be the lagrange interpolation of the exact solution, so we have

I_h u = \sum d_i \phi with d_i to be the nodal values.

We report the L2 norm of u-I_h u since it is not hard to see that

\|u-I_h u\|_{L^2} = sqrt((c-d)'M(c-d))

where M is the mass matrix, i.e.

M_{i,j} = \int_\Omega \phi_j \phi_i dx.

The mass matrix can be assembled using the stiffness matrix assembly code.

Wenyu Lei
Spencer Patty
Jan 7, 2016