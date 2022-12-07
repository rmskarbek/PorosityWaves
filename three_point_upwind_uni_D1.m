      function [D] = three_point_upwind_uni_D1(z0,zL,n,v)
%...
%...  The MatMol Group (2009)
%...
%...  function three_point_upwind_uni_D1 returns the differentiation matrix 
%...  for computing the first derivative, xz, of a variable x over the spatial 
%...  domain z0 < z < zL from upwind three-point, first-order finite difference 
%...  approximations (this function replaces dss014)
%...
%...  argument list
%...
%...     z0      lower boundary value of z (input)
%...
%...     zL      upper boundary value of z (input)
%...
%...     n       number of grid points in the z domain including the
%...             boundary points (input)
%...
%...     v       fluid velocity (positive from left to right - only the sign is used) (input)
%...
%...  origin of the approximation
%...
%...  this function is an application of second-order directional
%...  differencing in the numerical method of lines. It is intended
%...  specifically for the analysis of convective systems modelled by
%...  first-order hyperbolic partial differential equations as dis-
%...  cussed in the two_point function. The coefficients of the finite
%...  difference approximations used herein are taken from Bickley, W.
%...  G., Formulae for numerical differentiation, The Mathematical
%...  Gazette, pp. 19-27, 1941, n = 2, m = 1, p = 0, 1, 2.
%...
%...  compute the spatial increment
        dz=(zL-z0)/(n-1);
        r2fdz=1/(2*dz);
%...
%...     (1)  finite difference approximation for positive v     
              if v > 0
%...
%...             sparse discretization matrix      
%...
%...             interior points      
		  e = ones(n,1);
                  D = spdiags([e -4*e 3*e], [-2:0], n, n);
%...
%...             boundary points      
                   D(1,1:3) = [-3 +4 -1];
                   D(2,1:3) = [-1  0 +1];
              end;
%...
%...     (2)  finite difference approximation for negative v
              if v < 0
%...
%...             sparse discretization matrix      
%...
%...             interior points      
		  e = ones(n,1);
                  D = spdiags([-3*e 4*e -e], [0:2], n, n);
 %...
%...             boundary points      
                   D(n-1,(n-2):n) = [-1 0 +1];
                   D(n,(n-2):n) = [+1 -4 +3];
              end;                
%...      
        D=r2fdz*D;
