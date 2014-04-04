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

function [s] = strainspinodal(c, V0, Vrange, strain='eulerian')
% function [s] = strainspinodal(c, V0, Vrange, strain='eulerian')
%
% strainspinodal - find the spinodal point of the strain polynomial in the
% volume range indicated. The spinodal point is characterized because the
% second derivative of the energy versus the volume becomes zero, and so
% does the bulk modulus.  In addition, the pressure has its most negative
% value.
%
% Required input variables:
% c: coefficients of the strain polynomial.
% V0: reference volume used in the definition of the strain.
% Vrange: vector containing the range of volumes in which the search must be
%         performed.
%
% Optional input variables (all have default values):
% {strain = 'eulerian'}: strain form.
%
% Output:
% s: data structure containing:
%    s.Vsp: spinodal volume.
%    s.psp: pressure at the spinodal point.
%    s.Bsp: bulk modulus at the spinodal point (it should be zero).
%    s.err: 0 if no error happened, !=0 otherwise.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: October 2010

   if (nargin < 3 | nargin > 4)
      print_usage ();
   endif

   LOG = 0;
   ra = Vrange;
   s.err = 0;
   MPTS = 101;
   for steps = 1 : 4
      rr = linspace(min(ra), max(ra), MPTS);
      prop = straineval(c, V0, rr, strain);
      [bmin,ibmin] = min(abs(prop.B));
      s.err = s.err | (length(find(prop.B <= 0)) < 1);
      ra = [rr(max(1,ibmin-1)), rr(min(MPTS,ibmin+1))];
   endfor
   s.Vsp = rr(ibmin);
   s.Esp = prop.E(ibmin);
   s.psp = prop.p(ibmin);
   s.Bsp = prop.B(ibmin);

endfunction
