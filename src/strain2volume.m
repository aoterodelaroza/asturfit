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

function V = strain2volume(f, V0, strain='eulerian')
% function V = strain2volume(f, V0, strain='eulerian')
%
% strain2volume - convert a vector of strain into a vector of volume values.
%
% Required input variables:
% f: vector of strains.
% V0: reference volume used in the definition of the strain.
%
% Optional input variables (all have default values):
% {strain = 'eulerian'}: strain form.
%
% Output:
% V: vector containing the volumes.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: October 2010

   if (nargin < 2 | nargin > 3)
      print_usage ();
   endif

   # Determine the strain:
   if (strcmp(strain,'eulerian'))
      V = (f*2+1).^(-3/2) * V0;
   elseif (strcmp(strain, 'natural'))
      V = exp(f*3) * V0;
   elseif (strcmp(strain,'lagrangian'))
      V = (f*2+1).^(3/2) * V0;
   elseif (strcmp(strain,'infinitesimal'))
      V = (-f+1).^(-3) * V0;
   elseif (strcmp(strain, 'quotient') | strcmp(strain, 'x1'))
      V = f * V0;
   elseif (strcmp(strain, 'x3'))
      V = f.^3 * V0;
   elseif (strcmp(strain, 'xinv3'))
      V = f.^(-3) * V0;
   elseif (strcmp(strain, 'V'))
      V = f;
   else
      error('strain2volume: strain form requested is unknown!');
   endif

endfunction
