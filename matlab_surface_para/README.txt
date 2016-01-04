
The main script of this section is

surface_laplacian_torus.m

which calls the other functions to mesh, assemble, solve and compute l2_error on surface of torus.  The torus is hard coded into the various functions using

R=1.0
r=0.6
y = [theta,phi] \in [-pi,pi]x[-pi,pi]
chi : y -> x\in Torus \subset R^3
where
  chi(y) = [ (R+r*cos(theta))*cos(phi), (R+r*cos(theta))*sin(phi), r*sin(theta) ]

Note that chi is our parameterization, so that this is in reality a 2D code in parameter space (theta, phi) and mapped to the surface in R^3.  

We solve the Laplace-Beltrami equation

- \Delta_{S} u + u = f

on surface S=torus with the surface laplacian \Delta_{S} defined using differential geometry terms as 

\nabla \chi(y)  = [dx/dtheta dy/dtheta dz/dtheta]   is [2x3] 
                  [dx/dphi    dy/dphi  dz/dphi  ]
                      
g(y) = \nabla chi(y) * \nabla\chi(y)^T  \in [2x2]   is the metric defined by our parameterization
q(y) = sqrt(|det(g(y))|)

we denote u(x) = u(chi(y)) = u^{*}(y)   where u^{*} is called the pull back of u.

Check these ? 
grad_{S} u(x) = \nabla_y u^{*}(y) * inv(g(y)) * \nabla_y \chi(y)    [1x2] [2x2] [2x3] 
div_{S} V(x) = 1/q(y) \divergence ( q(y) V^{*}(y) ) 
\Delta_{S} u(x) =  div_{S} grad_{S} u(x)  =  1/q(y) * \divergence ( q(y) * inv(g(y)) * \nabla_y u^{*}(y) )


The exact solution is 

u(theta,phi) = sin(3*phi)*cos(3*theta+phi).

The error analysis is evaluated using 
l2_error = sqrt( error_at_nodes^T * MASS * error_at_nodes );
which is a second order quadrature approximation to the error.


References:


Wenyu Lei
Spencer Patty
Jan 4, 2016
