
The main script of this section is 

phase_field_torus.m

which calls the other function to mesh, assemble, compute energy error and visualize on the surface of torus.

The torus centered at origin can parametrized as 

y = [theta,phi] \in [-pi,pi]x[-pi,pi]
chi : y -> x\in Torus \subset R^3
where
  chi(y) = [ (R+r*cos(theta))*cos(phi), (R+r*cos(theta))*sin(phi), r*sin(theta) ]

Especially in the code, we let R = 1.0 and r = 0.6.

The goal is to solve the Laplace-Beltrami equation

- \Delta_{S} u + u = f on surface S

by the phase field approach.

Introduce the phase field function: for x in R^3

\psi(x) = p(d(x)/ \epsilon)+\epsilon q(x)+O(\epsilon^2)

Here d is he signed distance function and \epsilon>0 is small parameter. p is a monotone function such that p(-\infty)=-1 and p(\infty)=1.

I our case we take q = 0 and p(s) = 1-2/(1+exp(-s));

Let D be domain in R^3 containing the surface S. For example, D=[-2,2]x[-2,2]x[-1,1]. So the phase field weak formulation is given by

\int_D (1-\psi^2)(\nabla u\cdot\nabla v + uv -fv) dx =0 for all v in H^1(D).

Since this is standard 3D finite element implementation, we refer to 2D parametric implementation.

The corresponding energy norm is

\|v\|_\epsilon := \int_D (1-\psi^2)(|\nabla v|^2+v^2)dx /M_\epsilon

where M_\epsilon = \int_D (1-\psi^2) dx.


The error analysis is evaluated using 
l2_error = sqrt( error_at_nodes^T * A * error_at_nodes )/M_epsilon;
where A is the stiffness matrix and we can compute M_epsilon through local assembling.

At last, we note that the convergence rate is depended on \epsilon due to the definition of the phase field function.


Reference: 
Burger, M., 2009. Finite element approximation of elliptic partial differential equations on implicit surfaces. Computing and visualization in science, 12(3), pp.87-100.

Wenyu Lei
Spencer Patty
Jan 6, 2016
