% Copyright (C) 2011 Victor Lua~na and Alberto Otero-de-la-Roza
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

function p = pvstrainevalp (c, Vref, V, strain='eulerian')
% function p = pvstrainevalp (c, Vref, V {, strain})
%
% pvstrainevalp - Evaluation of the pressure from a polynomial p(V)
% function based on a given strain form.
%
% Required input variables:
% c: coefficients of the "energy versus f" polynomial. Usually this is
%    determined in a previous call to the strainfit() routine.
% Vref: reference volume.
% V: vector containing the volumes at which the polynomial properties must
%    be calculated.
%
% Optional input variables (all have default values):
% {strain = 'eulerian'}: strain form.
%
% Output:
% p: pressures evaluated at the V points.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: January 2011.
%

   if (nargin < 3 | nargin > 4)
      print_usage ();
   endif

   # Determine the strain:
   f = volume2strain(V, Vref, strain);

   # Evaluate E(f) and the derivatives of E versus f:
   p = polyval(c, f);
   
endfunction
