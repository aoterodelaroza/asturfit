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

function [cavg,savg] = pvavgstrainfit (p, V, Vref, nmax=16, MODE=1, strain='eulerian', LOG=1)
% function [cavg {,savg}] = pvavgstrainfit (p, V, Vref {,nmax,MODE,strain,LOG})
%
% pvavgstrainfit - Fit to an average strain polynomials.
%
% Required input variables:
% p: pressure.
% V: cell volume.
% Vref: volume reference. Usually the volume at the minimum of energy.
%
% Optional input variables (all have default values):
% {nmax = -1}: maximum degree of the polynomial.
% {MODE = 1}: polynomial weighting scheme.
% {strain = 'eulerian'}: strain form.
% {LOG = 0}: print internal information about the fitting if LOG>0.
%
% Minimum output:
% cavg: coefficients of the average fitting polynomial.
%
% Additional (optional) output:
% savg: data structure containing:
%    savg.eqmean[V, B, B1p, B2p, B3p]: mean value of the equilibrium
%       properties. If the average of polynomials has no suitable minimum
%       within the input volume range the "eq" properties will correspond
%       to the reference point.
%    savg.eqstd[V, B, B1p, B2p, B3p]: eq. props. standard deviation.
%    savg.R2: determination coefficient (1 for a exact fitting).
%    savg.pmean, savg.pstd: mean and standard error of the pressure at
%       all the input volumes.
%    savg.Bmean, savg.Bstd: mean and standard error of the bulk modulus at
%       all the input volumes.
%    savg.B1pmean, savg.B1pstd, savg.B2pmean, savg.B2pstd, savg.B3pmean,
%       savg.B3pstd, mean and standard error of the derivatives of the
%       bulk modulus.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: January 2011
%

   if (nargin < 3)
      print_usage ();
   elseif (length(V) != length(p))
      error('pvavgstrainfit: p and V must be vectors of the same length!')
   elseif (length(V) < 7)
      error('pvavgstrainfit: dataset must have at least 7 points!')
   elseif (MODE != 1 & MODE != 2)
      error('pvavgstrainfit: weighting mode must be 1 or 2!')
   endif

   ndata = length(V);
   Vrange = [min(V), max(V)];

   # Determine the maximum degree of the polynomials:
   if (nmax < 0)
      MaxDegree = min(ndata-5, fix(ndata/2));
   else
      MaxDegree = min(nmax,min(ndata-5, fix(ndata/2)));
   endif

   # Some statistics of the averaging proccess:
   Morder = npol = 0;
   for n = 2 : MaxDegree
      [c,sf] = pvstrainfit(p, V, Vref, n, strain, 0);
      [sm] = pvstrainmin(c, Vref, Vrange, strain);
      npol++;
      pol{npol}.c = c;
      s2(npol) = pol{npol}.SSerr = sf.S2;
      pol{npol}.order = n;
      pol{npol}.data = ndata;
      pol{npol}.smin = sm;
      Morder = max(n,Morder);
   endfor

   # Get the polynomial weights:
   SSmin = min(s2);
   Q = 0;
   if (MODE == 1)
      for k = 1 : npol
         w = (pol{k}.SSerr/SSmin) * (pol{k}.order / pol{k}.data);
         pol{k}.w = exp(-w);
         Q += pol{k}.w;
      endfor
   elseif (MODE == 2)
      for k = 1 : npol
         w = (pol{k}.SSerr/SSmin) * (pol{k}.order / pol{k}.data);
         pol{k}.w = exp(-w*w);
         Q += pol{k}.w;
      endfor
   endif
   ##ww = zeros(1,1:npol);
   for k = 1 : npol
      ww(k) = pol{k}.w = pol{k}.w / Q;
   endfor

   # Form the average polynomial:
   cavg = zeros(1,Morder);
   for k = 1 : npol
      n = pol{k}.order;
      cavg(Morder-n+1:Morder) += pol{k}.c .* pol{k}.w;
   endfor

   # Extra output and optional LOG record:
   # If the average of polynomials has a minimum, analize the equilibrium
   # geometry.
   # Otherwise analyze the reference volume.
   if (LOG > 0 | nargout > 1)
      [avgmin] = pvstrainmin(cavg, Vref, Vrange, strain);
      if (avgmin.err == 0)
         # Analyze and report the equilibrium geometry:
         if (LOG > 0)
            printf('\n\npvavgpolyfit: AVERAGE OF %s STRAIN POLYNOMIALS\n', strain);
            printf('Volume reference (Vref): %.6f\n', Vref);
            printf('Range of degrees: 2--%d\n', MaxDegree);
            printf('Number of polynomials: %d\n', npol);
            printf('\nProperties at the eq point of each polynomial:\n');
            printf('--i- npol data npar --SSerr-- ---w---- ----Vmin--- ----Bmin--- ---B1min--- ---B2min---\n');
         endif
         [srt,isrt] = sort(ww,'descend');
         pmean = pstd = zeros(5,1);
         for k = 1 : npol
            k1 = isrt(k);
            smin = pol{k1}.smin;
            [prop] = pvstraineval(pol{k1}.c, Vref, smin.Vmin, strain);
            if (LOG > 0 & k <= 25)
               # In the conversion of the bulk modulus and derivatives we
               # assume the units: volume (bohr^3), pressure (GPa).
               # hybohr3togpa = 2*14710.50498740275538944426;
               hybohr3togpa = 1;
               printf('%4d %4d %4d %4d', k, k1, pol{k1}.data, pol{k1}.order);
               printf(' %9.2e %8.6f', pol{k1}.SSerr, pol{k1}.w);
               printf(' %11.6f', smin.Vmin);
               printf(' %11.6f', prop.B * hybohr3togpa);
               printf(' %11.6f', prop.B1p);
               printf(' %11.9f', prop.B2p / hybohr3togpa);
               printf(' %11.9f', prop.B3p / hybohr3togpa^2);
               printf('\n');
            endif
            pr = [smin.Vmin; prop.B; prop.B1p; prop.B2p; prop.B3p];
            pmean = pmean .+ pr * pol{k1}.w;
            pstd = pstd .+ (pr.^2) * pol{k1}.w;
         endfor
         pstd = sqrt(pstd - pmean.^2);
         if (LOG > 0)
            printf('\nAverage properties (weighted polynomials):\n');
            printf('------ ---volume-- --B-(GPa)-- ----B1p---- B2p-(1/GPa) B3p--(1/GPa^2)\n');
            printf('-mean- %11.6f %11.6f %11.6f %11.6f %14.9f\n', pmean(1)...
                  , pmean(2)*hybohr3togpa, pmean(3), pmean(4)/hybohr3togpa...
                  , pmean(5)/hybohr3togpa^2);
            printf('stdvev %11.6f %11.6f %11.6f %11.6f %14.9f\n', pstd(1)...
                  , pstd(2)*hybohr3togpa, pstd(3)...
                  , pstd(4)/hybohr3togpa, pstd(5)/hybohr3togpa^2);
         endif
      else
         # Analyze and report the reference geometry:
         if (LOG > 0)
            printf('\n\npvavgpolyfit: AVERAGE OF %s STRAIN POLYNOMIALS\n', strain);
            printf('Volume reference (Vref): %.6f\n', Vref);
            printf('Range of degrees: 2--%d\n', MaxDegree);
            printf('Number of polynomials: %d\n', npol);
            printf('The average polynomial has no eq point in the data volume range.\n');
            printf('\nProperties at the reference volume: %.9f\n', Vref);
            printf('--i- npol data npar --SSerr-- ---w---- ----Bref--- ---B1ref--- ---B2ref---\n');
         endif
         [srt,isrt] = sort(ww,'descend');
         pmean = pstd = zeros(5,1);
         for k = 1 : npol
            k1 = isrt(k);
            [prop] = pvstraineval(pol{k1}.c, Vref, Vref, strain);
            if (LOG > 0 & k <= 25)
               # In the conversion of the bulk modulus and derivatives we
               # assume the units: volume (bohr^3), energy (Hy).
               # hybohr3togpa = 2*14710.50498740275538944426;
               hybohr3togpa = 1;
               printf('%4d %4d %4d %4d', k, k1, pol{k1}.data, pol{k1}.order);
               printf(' %9.2e %8.6f', pol{k1}.SSerr, pol{k1}.w);
               printf(' %11.6f', prop.B * hybohr3togpa);
               printf(' %11.6f', prop.B1p);
               printf(' %11.9f', prop.B2p / hybohr3togpa);
               printf(' %11.9f', prop.B3p / hybohr3togpa^2);
               printf('\n');
            endif
            pr = [Vref; prop.B; prop.B1p; prop.B2p; prop.B3p];
            pmean = pmean .+ pr * pol{k1}.w;
            pstd = pstd .+ (pr.^2) * pol{k1}.w;
         endfor
         pstd = sqrt(pstd - pmean.^2);
         if (LOG > 0)
            printf('\nAverage properties at the ref. volume: %.9f\n', Vref);
            printf('------ --B-(GPa)-- ----B1p---- B2p-(1/GPa) B3p--(1/GPa^2)\n');
            printf('-mean- %11.6f %11.6f %11.6f %14.9f\n'...
                  , pmean(2)*hybohr3togpa, pmean(3), pmean(4)/hybohr3togpa...
                  , pmean(5)/hybohr3togpa^2);
            printf('stdvev %11.6f %11.6f %11.6f %14.9f\n'...
                  , pstd(2)*hybohr3togpa, pstd(3)...
                  , pstd(4)/hybohr3togpa, pstd(5)/hybohr3togpa^2);
         endif
      endif

      pfit = pvstrainevalp(cavg, Vref, V, strain);
      SSerr = sum((p-pfit).^2);
      SStot = sum((p-mean(p)).^2);
      R2 = 1 - SSerr / SStot;
      if (nargout > 1)
         savg.eqmean = pmean;
         savg.eqstd = pstd;
         savg.R2 = R2;
         savg.pfit = pfit;
      endif
   endif

endfunction
