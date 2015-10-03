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

function [c,s] = pvstrainfit (p,V,Vref, ndeg=4, strain='eulerian', LOG=1)
% function [c {,s}] = pvstrainfit (p,V,Vref {,ndeg,strain,LOG})
%
% pvstrainfit - Polynomial fitting of the pressure to one of several forms
% of crystal strain. Currently available strain forms are: eulerian, natural,
% lagrangian, infinitesimal, quotient, and x.
%
% Required input variables:
% p: pressure.
% V: cell volume.
% Vref: volume reference used in thhe strain.
%
% Optional input variables:
% {ndeg = 4}: degree of the polynomial.
% {strain = 'eulerian'}: strain form.
% {LOG = 1}: print internal information about the fitting.
%
% Minimum output:
% c: coefficients of the fitting polynomial.
%
% Additional (optional) output:
% s: data structure containing:
%    s.V0: recalculated volume at the minimum of energy.
%    s.B0[1:3]: bulk modulus and its pressure derivatives.
%    s.R2: determination coefficient (1 for a exact fitting).
%    s.S2: square sum of residuals.
%    s.Efit: vector containing the predicted energies.
%    s.p: vector containing the predicted pressure.
%    s.B: vector containing the predicted bulk modulus.
%    s.f: vector containing the strain values.
%    s.err: 0 if no error.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: January 2011

   if (nargin < 3)
      print_usage ();
   endif

   # Determine the strain:
   f = volume2strain(V, Vref, strain);

   # Polynomial fitting:
   c = polyfit(f, p, ndeg-1);

   if (LOG > 0 || nargout > 1)
      # Evaluate the fitting polynomial and get the determination coefficient:
      pfit = polyval(c, f);
      SSerr = sum((p-pfit).^2);
      SStot = sum((p-mean(p)).^2);
      R2 = 1 - SSerr / SStot;
      # Find the equilibrium point, if any, and get its properties:
      smin = pvstrainmin(c, Vref, [min(V), max(V)], strain);
      eval = pvstraineval(c, Vref, smin.Vmin, strain);

      if (LOG > 0)
         printf('\n\n');
         printf('***************\n');
         printf('* pvstrainfit * p(V) fitting to the %s polynomial strain\n', strain)
         printf('***************\n')
         printf('Reference volume (Vref):  %.6f\n', Vref);
         printf('Polynomial degree :     %d\n', ndeg);
         printf('Determin. coef. (R2):   %.12f\n', R2);
         printf('Eq. volume (Vmin,err):  %.6f %3d\n', smin.Vmin, smin.err);
         printf('Pressure at eq. (peq):  %.9g\n', eval.p);
         printf('Beq                  :  %.9g\n', eval.B);
         printf('B1peq                :  %.9g\n', eval.B1p);
         printf('B2peq                :  %.9g\n', eval.B2p);
         printf('B3peq                :  %.9g\n', eval.B3p);
         printf('All values in the input units\n');
      endif

      if (nargout > 1)
         s.V0 = smin.Vmin;
         s.p0 = eval.p;
         s.B0 = [eval.B, eval.B1p, eval.B2p, eval.B3p];
         s.R2 = R2;
         s.S2 = SSerr;
         s.err = smin.err;
      endif
   endif

endfunction
