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

function E = strainstepevalE (c, V0, step, V, strain = 'eulerian')
%
% strainstepevalE - evaluate a strain polynomial augmented with Heaviside
%    step functions.
%
% Required input variables:
% c: coefficients of the strain polynomial.
% V0: reference volume for the strain.
% step: position and value of the Heaviside step functions.
%    The structure of this datum is:
%    for k = 1 : length(step)
%       step{k}.V
%       step{k}.E
%    endfor
% V: volume of the points to be evaluated.
%
% Optional input variables (all have default values):
% {strain = 'eulerian'}: strain form.
%
% Output:
% E: energy values.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: November 2010.
%
   if (nargin < 4 )
      print_usage ();
   endif
   E = strainevalE(c, V0, V, strain);
   for l = 1 : length(step)
      rk = find(V > step{l}.V);
      E(rk) += step{l}.E;
   endfor
endfunction
