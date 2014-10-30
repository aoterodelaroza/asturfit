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

function [s] = strainmin(c, V0, Vrange, strain='eulerian')
% function [s] = strainmin(c, V0, Vrange, strain='eulerian')
%
% strainmin - find the minimum of the strain polynomial in the volume
% range indicated.
%
% Required input variables:
% c: coefficients of the strain polynomial.
% V0: reference volume used in the definition of the strain.
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
%    s.Emin: energy at the minimum.
%    s.err:  error code. 0 means no problem.
%
%
% Uses: strainevalE().
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: October 2010

   if (nargin < 3 || nargin > 4)
      print_usage ();
   endif

   LOG = 0;
   frange = volume2strain(Vrange, V0, strain);

   %
   % Get the roots of the polynomial derivative:
   % Choose only those that are real and produce a positive second derivative.
   %
   c1 = polyder(c);
   c2 = polyder(c1);
   rr = roots(c1);
   ipos = find(abs(imag(rr)) <= 1e-15 && polyval(c2,rr) > 0);
   if (length(ipos) < 1)
      %
      % There are no left roots. An error in the roots() routine?
      %
      if (LOG)
         printf('polymin algorithmic error (1):\n');
         printf('coefs. of polynomial and derivative:\n');
         printf('-i----real(c)---imag(c)---real(c1)---imag(c1)---\n');
         for i = 1:length(c)-1
            printf('%3d %12.5e %12.5e %12.5e %12.5e\n',...
                   i, real(c(i)), imag(c(i)), real(c1(i)), imag(c1(i)));
         endfor
         printf('%3d %12.5e %12.5e\n', length(c), real(c(end)), imag(c(end)));
         printf('-iroot---real---imag---deriv2---\n');
         for i = 1:length(rr)
            printf('%6d %13.6e %13.6e %13.6e\n',...
                   i, real(rr(i)), imag(rr(i)), polyval(c2,rr(i)));
         endfor
      endif
      ###error('polymin: fitted polynomial has no minima!');
      s.err = 1;
   elseif (length(ipos) == 1)
      %
      % There is a single root. Check if it is inside range.
      %
      s.fmin = rr(ipos);
      s.err = s.fmin < min(frange)*0.98 || s.fmin > max(frange)*1.02;
   else
      %
      % There are several candidates. Evaluate the polynomial in a grid
      % of points within the range. Get the position of the minimum in the
      % grid and choost the root closest to that point.
      % Check again if it is inside range.
      %
      ff = linspace(min(frange), max(frange), 51);
      ee = polyval(c,ff);
      [eemin,ieemin] = min(ee);
      ffmin = ff(ieemin);
      rr = rr(ipos);
      rx = abs(rr .- ffmin).^2;
      [rxmin,irxmin] = min(rx);
      s.fmin = rr(irxmin);
      s.err = s.fmin < min(frange)*0.98 || s.fmin > max(frange)*1.02;
   endif

   if (s.err != 0)
      %
      % The above method has failed. Try Newton-Raphson.
      %
      if (exist('ieemin','var') == 0)
         % Evaluate the grid points if that has not been done before:
         ff = linspace(min(frange), max(frange), 51);
         ee = polyval(c,ff);
         [eemin,ieemin] = min(ee);
         ffmin = ff(ieemin);
      endif
      CONV = 1e-8; MAXIT = 30;
      fnew = ffmin; delta = 1e30; it = 0;
      do
         fold = fnew; dold = delta; it++;
         fnew = fold - polyval(c1,fold) / polyval(c2,fold);
         delta = abs(fnew-fold);
      until (delta < CONV || it > MAXIT)
      s.fmin = fnew;
      s.err = delta<CONV || s.fmin<min(frange)*0.98 || s.fmin>max(frange)*1.02;
   endif

   s.Vmin = strain2volume(s.fmin, V0, strain);
   s.Emin = polyval(c, s.fmin);

endfunction
