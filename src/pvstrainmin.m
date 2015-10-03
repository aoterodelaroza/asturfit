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

function [s] = pvstrainmin(c, Vref, Vrange, strain='eulerian')
% function [s] = pvstrainmin(c, Vref, Vrange, strain='eulerian')
%
% pvstrainmin - find the equilibrium point of the strain polynomial in the
% volume range indicated.
%
% Required input variables:
% c: coefficients of the strain polynomial.
% Vref: reference volume used in the definition of the strain.
% Vrange: vector containing the range of volumes in which the minimum of
%         energy will be determined.
%
% Optional input variables (all have default values):
% {strain = 'eulerian'}: strain form.
%
% Output:
% s: data structure containing:
%    s.Vmin: volume at the minimum of energy.
%    s.fmin: strain at the minimum.
%    s.pmin: pressure at the minimum (it should be 0).
%    s.err:  error code. 0 means no problem.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: January 2011

   if (nargin < 3 || nargin > 4)
      print_usage ();
   endif
   frange = volume2strain(Vrange, Vref, strain);
   ngrid = 10;
   do
      ngrid *= 2;
      fgrid = linspace(min(frange), max(frange), ngrid+1);
      pgrid = polyval(c,fgrid);
      [pmin,imin] = min(abs(pgrid));
   until ((imin>1 && imin<length(pgrid)) || ngrid > 1000)
   if (abs(pgrid(imin)) <= 1e-5)
      # We have located the p=0 point by chance.
      s.err = 0;
      s.pmin = pmin;
      s.fmin = fgrid(imin);
      s.Vmin = strain2volume(s.fmin,Vref,strain);
      return
   elseif (imin<=1 || imin>=length(pgrid))
      # Problem 1: the equilibrium point, if any, is at the end of the range.
      s.err = 1;
      s.pmin = pmin;
      s.fmin = fgrid(imin);
      s.Vmin = strain2volume(s.fmin,Vref,strain);
      return
   elseif (pgrid(imin-1)*pgrid(imin+1) > 0)
      # Problem 2: root condition is not found
      s.err = 2;
      s.pmin = pmin;
      s.fmin = fgrid(imin);
      s.Vmin = strain2volume(s.fmin,Vref,strain);
      return
   endif
   # We start a bisection search for the root:
   finf = fgrid(imin-1); pinf = pgrid(imin-1);
   fsup = fgrid(imin+1); psup = pgrid(imin+1);
   fmid = finf + (fsup-finf)/2;
   while (abs(fmid-finf)>eps && abs(fmid-fsup)>eps)
      pmid = polyval(c,fmid);
      if (pinf*pmid < 0)
         fsup = fmid;
         psup = pmid;
      elseif (psup*pmid < 0)
         finf = fmid;
         pinf = pmid;
      elseif (abs(pmid) < 1e-20)
         # p==0 point located by chance.
         s.err = 0;
         s.pmin = pmid;
         s.fmin = fmid;
         s.Vmin = strain2volume(s.fmin,Vref,strain);
         return
      else
         # Problem 3: bisection iteration failed.
         s.err = 3;
         s.pmin = pmid;
         s.fmin = fmid;
         s.Vmin = strain2volume(s.fmin,Vref,strain);
         return
      endif
      fmid = finf + (fsup-finf)/2;
   endwhile
   # Bisection iteration converged:
   s.err = 0;
   s.fmin = fmid;
   s.pmin = polyval(c,fmid);
   s.Vmin = strain2volume(s.fmin,Vref,strain);
   return
endfunction
