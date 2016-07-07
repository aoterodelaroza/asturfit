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

function [c,s] = strainfit (V,E,V0, ndeg=4, strain='eulerian', LOG=0)
% function [c {,s}] = strainfit (V,E,V0 {,ndeg,strain,LOG})
%
% strainfit - Polynomial fitting of the energy to one of several forms
% of crystal strain. Currently available strain forms are: eulerian, natural,
% lagrangian, infinitesimal, quotient, and x.
%
% Required input variables:
% V: cell volume.
% E: cell energy.
% V0: volume reference. Usually the volume at the minimum of energy.
%
% Optional input variables (all have default values):
% {ndeg = 4}: degree of the polynomial.
% {strain = 'eulerian'}: strain form.
% {LOG = 0}: print internal information about the fitting.
%
% Minimum output:
% c: coefficients of the fitting polynomial.
%
% Additional (optional) output:
% s: data structure containing:
%    s.V0: recalculated volume at the minimum of energy.
%    s.E0: energy at the minimum.
%    s.B0[1:3]: bulk modulus and its pressure derivatives.
%    s.R2: determination coefficient (1 for a exact fitting).
%    s.S2: square sum of residuals.
%    s.Efit: vector contained the predicted energies.
%    s.p: vector contained the predicted pressure.
%    s.B: vector contained the predicted bulk modulus.
%    s.f: vector containing the strain values.
%    s.err: 0 if no error.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: September 2010

   if (nargin < 3)
      print_usage ();
   endif

   # Determine the strain:
   if (strcmp(strain,'eulerian'))
      f = ((V./V0).^(-2/3)-1)/2;
      f0 = 0;
      if (LOG > 0 || nargout > 1)
         fpv = -(f + f + 1).^(5/2) / (3*V0);
         fppv = (f + f + 1).^4 * (5/(3*V0)^2);
      endif
   elseif (strcmp(strain, 'natural'))
      f = log(V./V0)/3;
      f0 = 0;
      if (LOG > 0 || nargout > 1)
         fpv = 1./(V.*3);
         fppv = -1./(V.**2 .*3);
      endif
   elseif (strcmp(strain,'lagrangian'))
      f = ((V./V0).^(2/3)-1)/2;
      f0 = 0;
      if (LOG > 0 || nargout > 1)
         fpv = (f + f + 1).^(-1/2) / (3*V0);
         fppv = -(f + f + 1).^(-2) / (3*V0)^2;
      endif
   elseif (strcmp(strain,'infinitesimal'))
      f = -(V./V0).^(-1/3)+1;
      f0 = 0;
      if (LOG > 0 || nargout > 1)
         fpv = (-f + 1).^4 / (3*V0);
         fppv = -(-f + 1).^7 * (4/(3*V0)^2);
      endif
   elseif (strcmp(strain, 'quotient') || strcmp(strain, 'x1'))
      f = V./V0;
      f0 = 1;
      if (LOG > 0 || nargout > 1)
         fpv = ones(size(V)) * (1/V0);
         fppv = zeros(size(V));
      endif
   elseif (strcmp(strain, 'x3'))
      f = (V./V0).^(1/3);
      f0 = 1;
      if (LOG > 0 || nargout > 1)
         fpv = 1./(f.^2*(3*V0));
         fppv = -f.^(-5) * (2/(3*V0)^2);
      endif
   elseif (strcmp(strain, 'xinv3'))
      f = (V./V0).^(-1/3);
      f0 = 1;
      if (LOG > 0 || nargout > 1)
         fpv = -f.^4 * (1/(3*V0));
         fppv = f.^7 * (4/(3*V0)^2);
      endif
   elseif (strcmp(strain, 'V'))
      f = V;
      if (LOG > 0 || nargout > 1)
         fpv = ones(size(V));
         fppv = zeros(size(V));
      endif
   else
      error('strainfit: strain form requested is unknown!');
   endif

   c = polyfit(f, E, ndeg);

   if (LOG > 0 || nargout > 1)
      # Evaluate the fitting polynomial and get the determination coefficient:
      Efit = polyval(c, f);
      SSerr = sum((E-Efit).^2);
      SStot = sum((E-mean(E)).^2);
      R2 = 1 - SSerr / SStot;
      # Find the polynomial minimum and its properties:
      # 1. evaluate the polynomial in a fine grid and get the best grid point
      # 2. get the roots of the polynomial derivative
      # 3. select the real roots that give a positive second derivative
      # 4. get the candidate from (3) that is closest to the best grid point
      ff = linspace(min(f), max(f), 51);
      ee = polyval(c, ff);
      [emin,imin] = min(ee);
      ffmin = ff(imin);
      c1 = polyder(c);
      c2 = polyder(c1);
      rr = roots(c1);
      ipos = find(abs(imag(rr)) <= 1e-15 & polyval(c2,rr) > 0);
      if (length(ipos) < 1)
         printf('strainfit algorithmic error (1):\n');
         printf('strain & ndeg: %s %d\n', strain, ndeg);
         printf('coefs. of polynomial and derivative:\n');
         printf('-icoef----real(c)---imag(c)---real(c1)---imag(c1)---\n');
         for i = 1:length(c)-1
            printf('%6d %22.15e %22.15e %22.15e %22.15e\n', i, real(c(i)), imag(c(i)), real(c1(i)), imag(c1(i)));
         endfor
         printf('%6d %22.15e %22.15e\n', i, real(c(end)), imag(c(end)));
         printf('-iroot---real---imag---deriv2---\n');
         for i = 1:length(rr)
             printf('%6d %13.6e %13.6e %13.6e\n', i, real(rr(i)), imag(rr(i)), polyval(c2,rr(i)));
         endfor
         ###error('strainfit: fitted polynomial has no minima!');
         s.err = 1;
         return
      elseif (length(ipos) == 1)
         fmin = rr(ipos);
      else
         rr = rr(ipos);
         rx = abs(rr .- ffmin).**2;
         [rxmin,irxmin] = min(rx);
         fmin = rr(irxmin);
      endif
      Emin = polyval(c, fmin);
      E1min = polyval(c1, fmin);
      c2 = polyder(c1); E2min = polyval(c2, fmin);
      c3 = polyder(c2); E3min = polyval(c3, fmin);
      if (strcmp(strain, 'eulerian'))
         Vmin = V0 * (1+2*fmin)^(-3/2);
         fpvmin = -(2*fmin+1)^(5/2)/(3*V0);
         fppvmin = (2*fmin+1)^4 * (5/(3*V0)**2);
         fpppvmin = -(2*fmin+1)^(11/2) * (40/(3*V0)**3);
      elseif (strcmp(strain, 'natural'))
         Vmin = V0 * exp(3*fmin);
         fpvmin = 1/(3*Vmin);
         fppvmin = -1/(3*Vmin**2);
         fpppvmin = 2/(3*Vmin**3);
      elseif (strcmp(strain, 'lagrangian'))
         Vmin = V0 * (1+2*fmin)^(3/2);
         fpvmin = (2*fmin+1)^(-1/2)/(3*V0);
         fppvmin = -(2*fmin+1)^(-2)/(3*V0)**2;
         fpppvmin = (2*fmin+1)^(-7/2)*(4/(3*V0)**3);
      elseif (strcmp(strain, 'infinitesimal'))
         Vmin = V0 * (1-fmin)^(-3);
         fpvmin = (1-fmin)^4/(3*V0);
         fppvmin = -(1-fmin)^7 * (4/(3*V0)**2);
         fpppvmin = (1-fmin)^10 * (28/(3*V0)**3);
      elseif (strcmp(strain, 'quotient') || strcmp(strain, 'x1'))
         Vmin = V0 * fmin;
         fpvmin = 1/V0;
         fppvmin = 0;
         fpppvmin = 0;
      elseif (strcmp(strain, 'x3'))
         Vmin = V0 * fmin^3;
         ss = -fmin^(-3)/(3*V0);
         fpvmin = fmin^(-2)/(3*V0);
         fppvmin = 2*fpvmin*ss;
         fpppvmin = 5*fppvmin*ss;
      elseif (strcmp(strain, 'xinv3'))
         Vmin = V0 / fmin^3;
         ss = -fmin^3/(3*V0);
         fpvmin = - fmin^4/(3*V0);
         fppvmin = 4*fpvmin*ss;
         fpppvmin = 7*fppvmin*ss;
      elseif (strcmp(strain, 'V'))
         Vmin = fmin;
         fpvmin = 1;
         fppvmin = 0;
         fpppvmin = 0;
      endif
      Bmin = Vmin * (E2min * fpvmin^2 + E1min * fppvmin);
      Bfmin = E3min*Vmin*fpvmin^2 + E2min*(fpvmin+3*fppvmin*Vmin)...
            + E1min*(fppvmin+Vmin*fpppvmin)/fpvmin;
      pfmin = -(E2min*fpvmin-E1min*fppvmin/fpvmin);
      Bpmin = Bfmin/pfmin;
   endif

   if (LOG > 0)
      printf('\n\n*************\n')
      printf('* strainfit * Polynomial fitting to the %s strain\n', strain)
      printf('*************\n')
      printf('Reference volume (V0):  %.6f\n', V0);
      printf('Polynomial degree :     %d\n', ndeg);
      printf('Determin. coef. (R2):   %.12f\n', R2);
      printf('Vol. at min. (Vmin):    %.6f\n', Vmin);
      printf('Energy at min. (Emin):  %.9g\n', Emin);
      printf('B at min.      (Bmin):  %.9g\n', Bmin);
      printf('Bp at min.    (Bpmin):  %.9g\n', Bpmin);
      printf('All values in the input units\n');
   endif

   if (nargout > 1)
      s.V0 = Vmin;
      s.E0 = Emin;
      s.B0 = [Bmin, Bpmin];
      s.R2 = R2;
      s.S2 = SSerr;
      s.Efit = Efit;
      s.f = f;
      s.err = 0;
   endif

endfunction
