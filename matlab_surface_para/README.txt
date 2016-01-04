
The main script of this section is

surface_laplacian_torus.m

which calls the other functions to mesh, assemble, solve and compute l2_error on surface of torus.  The torus is hard coded into the various functions using

R=1.0
r=0.6
theta,phi \in [-pi,pi]
chi = [(R+r*cos(theta))*cos(phi), (R+r*cos(theta))*sin(phi), r*sin(theta)  ] 
where chi : (theta,phi) -> Torus \subset R^3

Note that chi is our parameterization, so that this is in reality a 2D code in parameter space (theta, phi) and mapped to the surface in R^3.  

The exact solution is 

u(theta,phi) = sin(3*phi)*cos(3*theta+phi).

The error analysis is evaluated using 
l2_error = sqrt( error_at_nodes^T * MASS * error_at_nodes );
which is a second order quadrature approximation to the error.


References:


Wenyu Lei
Spencer Patty
Jan 4, 2016
