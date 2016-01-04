
The main script is 

level_set_torus_script.m

which steps through constructing the mesh, assembling the matrix system, solving and then computing the L2(Gamma_h) error.  Any modification or addition will likely start there.  This entire code is setup to solve a the laplace beltrami equation on the Torus in R^3

-\Delta_{\Gamma} + u = f.

The weak form of our problem where Gamma is defined implicitly using the distance function (described below) is:  Given D_h = \{ \x \in \R^3 | |I_h d(\x)| < h \} and a sequence of bulk triangulations {T_h}_{h>0} of tetrahedra containing D_h,  find u_h \in P^1(T_h \cap D_h ) such that

\int_{D_h} \grad u_h \cdot \grad v_h + u_h v_h ) |\grad I_h d(\x)| d\x = \int_{D_h} f(x) v_h |\grad I_h d(\x)| d\x for all v_h \in P^1(T_h\cap D_h).

as described in the the narrow banded method (NBM) of [1].

The torus shape is hard coded into many of the individual functions.  The torus we are using is parameterized by

R = 1.0;
r = 0.6;
T = \{ p\in R^3 |  p = [(R+r*cos(theta))*cos(phi), (R+r*cos(theta))*sin(phi), r*sin(theta)]  with  phi, theta \in [-pi, pi] \}

with the exact solution defined on T as 

u(p) = u(theta,phi) = sin(3*phi)*cos(3*theta + phi)  

and extended to a neighborhood of the torus being constant in the normal direction.  In other words, using the signed distance function d(\x) to the torus

d(\x) = sqrt( (sqrt( x^2+y^2)-R)^2 + z^2 ) - r

we can project \x in the neighborhood to the closest p on the torus along the normal '\grad d(\x)' (which is always orthogonal to the levelsets) via

p = \x - d(\x)\nabla d(\x).

so that u_ext( x) = u(p(\x)).  Likewise, the right hand side function is evaluated analytically using Maple as f = -\Delta_{\Gamma} u + u
and so we must extend it constant in the normal directions like the exact solution.


Most if not all the functions have been vectorized in some fashion or another to get as much speed out of the code without sacrificing too much readability.  There are a number of test_*.m files which constitute our test suite to verify that the many individual components are working properly.  Feel free to add other tests as you add functionality.  


References:

[1] Deckelnick, Klaus, Charles M. Elliott, and Thomas Ranner. "Unfitted finite element methods using bulk meshes for surface partial differential equations." SIAM Journal on Numerical Analysis 52.4 (2014): 2137-2162.


Wenyu Lei
Spencer Patty
Dec 30, 2015
