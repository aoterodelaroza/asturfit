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

function f = volume2strain(V, V0, strain='eulerian')
% function f = volume2strain(V, V0, strain='eulerian')
%
% volume2strain - convert a vector of volumes into a vector of strain values.
%
% Required input variables:
% V: vector containing the volumes to convert.
% V0: reference volume used in the definition of the strain.
%
% Optional input variables (all have default values):
% {strain = 'eulerian'}: strain form.
%
% Output:
% f: vector of strains.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: October 2010

   if (nargin < 2 || nargin > 3)
      print_usage ();
   endif

   # Determine the strain:
   if (strcmp(strain,'eulerian'))
      f = ((V./V0).^(-2/3)-1)/2;
   elseif (strcmp(strain, 'natural'))
      f = log(V./V0)/3;
   elseif (strcmp(strain,'lagrangian'))
      f = ((V./V0).^(2/3)-1)/2;
   elseif (strcmp(strain,'infinitesimal'))
      f = -(V./V0).^(-1/3)+1;
   elseif (strcmp(strain, 'quotient') || strcmp(strain, 'x1'))
      f = V./V0;
   elseif (strcmp(strain, 'x3'))
      f = (V./V0).^(1/3);
   elseif (strcmp(strain, 'xinv3'))
      f = (V./V0).^(-1/3);
   elseif (strcmp(strain, 'V'))
      f = V;
   else
      error('volume2strain: strain form requested is unknown!');
   endif

endfunction
