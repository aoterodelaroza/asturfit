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

function p = strainevalp (c, V0, V, strain='eulerian')
% function p = strainevalp (c, V0, V {, strain})
%
% strainevalp - Evaluation of the pressure (in GPa) derived from a polynomial
% function based on a given strain form.
% Use the more general straineval() for more properties.
%
% Required input variables:
% c: coefficients of the "energy versus f" polynomial. Usually this is
%    determined in a previous call to the strainfit() routine.
% V0: reference volume.
% V: vector containing the volumes at which the polynomial properties must
%    be calculated.
%
% Optional input variables (all have default values):
% {strain = 'eulerian'}: strain form.
%
% Output:
% p: vector with pressures (in GPa) at the V points.
%
% See also: strainfit.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: September 2010.
%

   hybohr3togpa = 2*14710.50498740275538944426;

   if (nargin < 3 | nargin > 4)
      print_usage ();
   endif

   # Determine the strain:
   if (strcmp(strain,'eulerian'))
      f = ((V./V0).^(-2/3)-1)/2;
      f2 = f + f + 1;
      f1v = f2.^(5/2) * (-1/(3*V0));
   elseif (strcmp(strain, 'natural'))
      f = log(V./V0)/3;
      f1v = 1./(V*3);
   elseif (strcmp(strain,'lagrangian'))
      f = ((V./V0).^(2/3)-1)/2;
      f1v = (f + f + 1).^(-1/2) / (3*V0);
   elseif (strcmp(strain,'infinitesimal'))
      f = -(V./V0).^(-1/3)+1;
      f1v = (-f + 1).^4 / (3*V0);
   elseif (strcmp(strain, 'quotient') | strcmp(strain, 'x1'))
      f = V./V0;
      f1v = ones(size(V)) * (1/V0);
   elseif (strcmp(strain, 'x3'))
      f = (V./V0).^(1/3);
      f1v = V.^(-2/3)/(3*V0^(1/3));
   elseif (strcmp(strain, 'xinv3'))
      f = (V./V0).^(-1/3);
      f1v = -V.^(-4/3) * (V0^(1/3)/3);
   elseif (strcmp(strain, 'V'))
      f = V;
      f1v = ones(size(V));
   else
      error('straineval: strain form requested is unknown!');
   endif

   # Evaluate E(f):
   c1 = polyder(c);
   p = - polyval(c1, f) .* f1v * hybohr3togpa;

endfunction
