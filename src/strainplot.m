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

function strainplot (V,E,V0, c, strain='eulerian', fplot='plot.eps', type='v')
% function strainplot (V,E,V0,c {,strain,fplot,type})
%
% strainplot - Plot of the E(V) or E(f) data and the fitting polynomial.
%
% Required input variables:
% V: cell volume.
% E: cell energy.
% V0: volume reference. Usually the volume at the minimum of energy.
% c: polynomial coefficients.
%
% Optional input variables (all have default values):
% {strain = 'eulerian'}: strain form.
% {fplot = 'plot.eps'}: plot file name.
% {type = 'v'}: select between E(V) (defaulf) or E(f) (option 'f') plots.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: September 2010

   if (nargin < 4)
      print_usage ();
   endif

   vv = linspace(min(V), max(V), min(length(V)*4,101));
   ee = strainevalE(c, V0, vv, strain);

   if (strcmp(type,'v'))

      plot(V, E, 'o', vv, ee, '-');
      xlabel('V (bohr3)'); ylabel('E (Hy)');
      print(fplot, '-depsc');

   elseif (strcmp(type,'f'))

      f = volume2strain(V, V0, strain);
      ff = volume2strain(vv, V0, strain);

      plot(f, E, 'o', ff, ee, '-');
      xlabel('f'); ylabel('E (Hy)');
      print(fplot, '-depsc');

   endif

endfunction
