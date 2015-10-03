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

function [c, stepout, Ecorr] = strainjumpfit (V, E, V0, stepin, ndeg=14 \
         , strain='eulerian', conv=1e-4, MAXSTEP=20, LOG=0)
% function [c, stepout {,Ecorr}] = strainjumpfit (V, E, V0, stepin {,ndeg=14 \
% , strain='eulerian', conv=1e-4, MAXSTEP=20, LOG=0})
%
% strainjumpfit - least squares fitting of a strain polynomial plus a set
% of step functions. Notice that a good guess of the step functions is very
% important, so calling this routine should be the second part after calling
% guessjump() or a similar task.
%
% Required input variables:
% (V,E): column vectors containing the volume and energy of the E(V) points.
% V0: reference volume for the strain.
% stepin: position and value of the Heaviside step functions. The value
%    will be optimized in this routine, but the number and position will
%    be maintained. The structure of this datum is:
%    for k = 1 : length(stepin)
%       stepin{k}.V
%       stepin{k}.E
%    endfor
%
% Optional input variables (all have default values):
% {ndeg = 14}: degree for the strain polynomial. If negative, an average of
%     strain poynomials up to degree abs(ndeg) will be used.
% {strain = 'eulerian'}: strain form.
%
% Output:
% c: coefficients of the fitting polynomial.
% stepout: optimized step functions. The structure is identical to that of
%    stepin, so the output step o f a run can be used as the input step for
%    a second round.
%
% Additional (optional) output:
% Ecorr: energy values with the known jumps taken out.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: November 2010.
%

   if (nargin < 4 )
      print_usage ();
   elseif (ndeg==0)
      error('strainjumpfit: ndeg cannot be zero!');
   endif

   # Initialization:
   m = length(stepin);
   for l = 1:m
      stepout{l}.V = stepin{l}.V;
      stepout{l}.E = stepin{l}.E;
   endfor
   if (LOG > 0)
      plot(V,E,'o');
      vv = linspace(min(V), max(V), 201);
      printf('DBG(strainjumpfit): Fitting a function with jumps');
      printf('\nIteration: %d\n', 0);
      printf('-it- -----stepV----- -----stepE-----\n');
      for l = 1 : m
         printf('%4d %15.6f %15.6f\n', l, stepout{l}.V, stepout{l}.E);
      endfor
   endif
   A = zeros(m,m); b = zeros(m,1); maxdif = 10*conv; nstp = 0;
   while (++nstp <= MAXSTEP && maxdif > conv)
      #
      # Remove the known steps:
      Erest1 = E;
      for l = 1 : m
         rk = find(V > stepout{l}.V);
         Erest1(rk) -= stepout{l}.E;
      endfor
      #
      # Polynomial fitting:
      if (ndeg > 0)
         c = strainfit(V, Erest1, V0, ndeg, strain, 0);
      elseif (ndeg < 0)
         c = avgstrainfit(V, Erest1, V0, -ndeg, :, strain, 0);
      endif
      Epol = strainevalE(c, V0, V, strain);
      Erest2 = E - Epol;
      #
      # Fit of the Heaviside step functions:
      for s = 1:m
         for r = 1:m
            rk = find(V > max(stepout{r}.V,stepout{s}.V));
            A(r,s) = length(rk);
         endfor
         rk = find(V > stepout{s}.V);
         b(s) = sum(Erest2(rk));
      endfor
      delta = A\b;
      maxdif = 0;
      for l = 1:m
         maxdif = max(maxdif, abs(stepout{l}.E-delta(l)));
         stepout{l}.E = delta(l);
      endfor
      if (LOG > 0)
         plot(vv,strainstepevalE(c,V0,stepout,vv,strain),'-r', V,E,'ob');
         printf('\nIteration: %d  maxdif: %.6e\n', nstp, maxdif);
         printf('-it- -----stepV----- -----stepE-----\n');
         for l = 1 : m
            printf('%4d %15.6f %15.6f\n', l, stepout{l}.V, stepout{l}.E);
         endfor
      endif
   endwhile

   if (nargout > 2)
      #
      # Remove the known steps:
      Ecorr = E;
      for l = 1 : m
         rk = find(V > stepout{l}.V);
         Ecorr(rk) -= stepout{l}.E;
      endfor
   endif
   
endfunction
