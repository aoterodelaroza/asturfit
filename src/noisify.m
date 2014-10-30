% Copyright (C) 2010 Victor Lua~na and Alberto Otero-de-la-Roza
%
% This octave routine is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version. See <http://www.gnu.org/licenses/>.
%
% The routine distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.

function [xn,yn] = noisify (x, y, points=-0.1, xfac=0.1, yfac=0.1, LOG=0)
% function [xn,yn] = noisify (x, y {,points, xfac, yfac, LOG})
%
% noisify - add noise to the (x,y) data.
%
% Required input variables:
% x: vector with values of the independent variable.
% y: vector with values of the dependent variable.
%
% Optional input variables (all have default values):
% {points = -0.1}: points to which some noise will be added. A positive
%   integer represents the number of points to be modified. A negative number
%   represents (after changing the sign) the fraction of input points that
%   will be modified. The default, -0.1, means that ten percent of the input
%   data will be changed.
% {xfac = 0.1}: noise factor for the x. An input value, x(k), will be
%   converted into x(k) + x(k)*xfac*random, where the random number is
%   uniformly distributed in [0,1). Use xfac=0 to prevent changes to the
%   x values.
%   *Caution* no test will be done internally to prevent absurd values for
%   (xfac,yfac).
% {yfac = 0.1}: noise factor for the y.
% {LOG = 0}: print an internal information report.
%
% Output:
% [xn,yn]: vectors containing the [x,y] data with added noise.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: October 2010

   if (nargin < 2)
      print_usage ();
   elseif (length(x) != length(y))
      error('noisify: x and y must be vectors of the same length!')
   elseif (points==0 || (xfac==0 && yfac==0))
      error('noisify: nothing will be done!')
   endif

   ndata = length(x);

   if (points > 0)
      nnoise = min(ndata, points);
   else
      nnoise = min(ndata, fix(-points*ndata));
   endif

   xn = x; yn = y;
   for i = 1:nnoise
      i1 = fix(rand()*ndata+1);
      xn(i1) *= 1 + (2*rand()-1)*xfac;
      yn(i1) *= 1 + (2*rand()-1)*yfac;
   endfor

   if (LOG > 0)
      nd = find(xn!=x || yn!=y);
      nnd = length(nd);
      printf('\nnoisify: Added noise to %d points of %d\n', nnd, ndata);
      printf('--point-- ------xold----- ------xnew----- ------yold----- ------ynew-----\n');
      for i = 1:nnd
         i1 = nd(i);
         printf('%9d %15.6f %15.6f %15.6f %15.6f\n', i1, x(i1), xn(i1), y(i1), yn(i1));
      endfor
   endif

endfunction
