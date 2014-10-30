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

function [cavg,savg] = avgstrainfit (V, E, V0, nmax=16, MODE=1, strain='eulerian', LOG=0)
% function [cavg {,savg}] = avgstrainfit (V, E, V0 {,nmax,MODE,strain,LOG})
%
% avgstrainfit - Fit to an average strain polynomials.
%
% Required input variables:
% V: cell volume.
% E: cell energy.
% V0: volume reference. Usually the volume at the minimum of energy.
%
% Optional input variables (all have default values):
% {nmax = -1}: maximum degree of the polynomial.
% {MODE = 1}: polynomial weighting scheme.
% {strain = 'eulerian'}: strain form.
% {LOG = 0}: print internal information about the fitting if LOG>0.
%            Some special modes:
%            LOG > 1 produces an statistical analysis of the spinodal point,
%                    if it is present in the average polynomial.
%            LOG > 2 produces an statistical analysis of p, B, B1p, etc
%
% Minimum output:
% cavg: coefficients of the average fitting polynomial.
%
% Additional (optional) output:
% savg: data structure containing:
%    savg.eqmean[V, E, B, B1p, B2p, B3p]: mean value of the equilibrium
%       properties. If the average of polynomials has no suitable minimum
%       within the input volume range the "eq" properties will correspond
%       to the reference point.
%    savg.eqstd[V, E, B, B1p, B2p, B3p]: eq. props. standard deviation.
%    savg.R2: determination coefficient (1 for a exact fitting).
%    savg.Efit: vector containing the predicted energies.
%    savg.Estd: vector with the standard error for the energies.
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
% Created: September 2010
%

   if (nargin < 3)
      print_usage ();
   elseif (length(V) != length(E))
      error('avgstrainfit: V and E must be vectors of the same length!')
   elseif (length(V) < 7)
      error('avgstrainfit: dataset must have at least 7 points!')
   elseif (MODE != 1 && MODE != 2)
      error('avgstrainfit: weighting mode must be 1 or 2!')
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
      [c] = strainfit(V, E, V0, n, strain, 0);
      [sm] = strainmin(c, V0, Vrange, strain);
      Efit = strainevalE(c, V0, V, strain);
      SSerr = sum((E-Efit).**2);
      npol++;
      pol{npol}.c = c;
      s2(npol) = pol{npol}.SSerr = SSerr;
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
   cavg = zeros(1,Morder+1);
   for k = 1 : npol
      n = pol{k}.order;
      cavg(Morder-n+1:Morder+1) += pol{k}.c .* pol{k}.w;
   endfor

   # Extra output and optional LOG record:
   # If the average of polynomials has a minimum, analize the equilibrium
   # geometry.
   # Otherwise analyze the reference volume.
   if (LOG > 0 || nargout > 1)
      [avgmin] = strainmin(cavg, V0, Vrange, strain);
      if (avgmin.err == 0)
         # Analyze and report the equilibrium geometry:
         if (LOG > 0)
            printf('\n\navgpolyfit: AVERAGE OF %s STRAIN POLYNOMIALS\n', strain);
            printf('Volume reference (V0): %.6f\n', V0);
            printf('Range of degrees: 2--%d\n', MaxDegree);
            printf('Number of polynomials: %d\n', npol);
            printf('\nProperties at the minimum of each polynomial:\n');
            printf('--i- npol data npar --SSerr-- ---w---- ----Vmin--- ----Emin--- ----Bmin--- ---B1min--- ---B2min---\n');
         endif
         [srt,isrt] = sort(ww,'descend');
         pmean = pstd = zeros(6,1);
         for k = 1 : npol
            k1 = isrt(k);
            ###[smin] = strainmin(pol{k1}.c, V0, Vrange, strain);
            smin = pol{k1}.smin;
            [prop] = straineval(pol{k1}.c, V0, smin.Vmin, strain);
            if (LOG > 0 && k <= 25)
               # In the conversion of the bulk modulus and derivatives we
               # assume the units: volume (bohr^3), energy (Hy).
               hybohr3togpa = 2*14710.50498740275538944426;
               printf('%4d %4d %4d %4d', k, k1, pol{k1}.data, pol{k1}.order);
               printf(' %9.2e %8.6f', pol{k1}.SSerr, pol{k1}.w);
               printf(' %11.6f', smin.Vmin);
               printf(' %11.6f', prop.E);
               printf(' %11.6f', prop.B * hybohr3togpa);
               printf(' %11.6f', prop.B1p);
               printf(' %11.9f', prop.B2p / hybohr3togpa);
               printf(' %11.9f', prop.B3p / hybohr3togpa^2);
               printf('\n');
            endif
            pr = [smin.Vmin; prop.E; prop.B; prop.B1p; prop.B2p; prop.B3p];
            pmean = pmean .+ pr * pol{k1}.w;
            pstd = pstd .+ (pr.^2) * pol{k1}.w;
         endfor
         pstd = sqrt(pstd - pmean.^2);
         if (LOG > 0)
            printf('\nAverage properties (weighted polynomials):\n');
            printf('------ ---volume-- ---energy-- --B-(GPa)-- ----B1p---- B2p-(1/GPa) B3p--(1/GPa^2)\n');
            printf('-mean- %11.6f %11.6f %11.6f %11.6f %11.6f %14.9f\n', pmean(1)...
                  , pmean(2), pmean(3)*hybohr3togpa, pmean(4), pmean(5)/hybohr3togpa...
                  , pmean(6)/hybohr3togpa^2);
            printf('stdvev %11.6f %11.6f %11.6f %11.6f %11.6f %14.9f\n', pstd(1)...
                  , pstd(2), pstd(3)*hybohr3togpa, pstd(4)...
                  , pstd(5)/hybohr3togpa, pstd(6)/hybohr3togpa^2);
         endif
      else
         # Analyze and report the reference geometry:
         if (LOG > 0)
            printf('\n\navgpolyfit: AVERAGE OF %s STRAIN POLYNOMIALS\n', strain);
            printf('Volume reference (V0): %.6f\n', V0);
            printf('Range of degrees: 2--%d\n', MaxDegree);
            printf('Number of polynomials: %d\n', npol);
            printf('The average polynomial has no minimum in the data volume range.\n');
            printf('\nProperties at the reference volume: %.9f\n', V0);
            printf('--i- npol data npar --SSerr-- ---w---- ----Eref--- ----Bref--- ---B1ref--- ---B2ref---\n');
         endif
         [srt,isrt] = sort(ww,'descend');
         pmean = pstd = zeros(6,1);
         for k = 1 : npol
            k1 = isrt(k);
            [prop] = straineval(pol{k1}.c, V0, V0, strain);
            if (LOG > 0 && k <= 25)
               # In the conversion of the bulk modulus and derivatives we
               # assume the units: volume (bohr^3), energy (Hy).
               hybohr3togpa = 2*14710.50498740275538944426;
               printf('%4d %4d %4d %4d', k, k1, pol{k1}.data, pol{k1}.order);
               printf(' %9.2e %8.6f', pol{k1}.SSerr, pol{k1}.w);
               printf(' %11.6f', prop.E);
               printf(' %11.6f', prop.B * hybohr3togpa);
               printf(' %11.6f', prop.B1p);
               printf(' %11.9f', prop.B2p / hybohr3togpa);
               printf(' %11.9f', prop.B3p / hybohr3togpa^2);
               printf('\n');
            endif
            pr = [V0; prop.E; prop.B; prop.B1p; prop.B2p; prop.B3p];
            pmean = pmean .+ pr * pol{k1}.w;
            pstd = pstd .+ (pr.^2) * pol{k1}.w;
         endfor
         pstd = sqrt(pstd - pmean.^2);
         if (LOG > 0)
            printf('\nAverage properties at the ref. volume: %.9f\n', V0);
            printf('------ ---energy-- --B-(GPa)-- ----B1p---- B2p-(1/GPa) B3p--(1/GPa^2)\n');
            printf('-mean- %11.6f %11.6f %11.6f %11.6f %14.9f\n', pmean(2)...
                  , pmean(3)*hybohr3togpa, pmean(4), pmean(5)/hybohr3togpa...
                  , pmean(6)/hybohr3togpa^2);
            printf('stdvev %11.6f %11.6f %11.6f %11.6f %14.9f\n', pstd(2)...
                  , pstd(3)*hybohr3togpa, pstd(4)...
                  , pstd(5)/hybohr3togpa, pstd(6)/hybohr3togpa^2);
         endif
      endif

      Efit = strainevalE(cavg, V0, V, strain);
      SSerr = sum((E-Efit).^2);
      SStot = sum((E-mean(E)).^2);
      R2 = 1 - SSerr / SStot;
      if (nargout > 1)
         savg.eqmean = pmean;
         savg.eqstd = pstd;
         ###[prop] = straineval(cavg, V0, avgmin.Vmin, strain);
         ###savg.V0 = avgmin.Vmin;
         ###savg.E0 = prop.E;
         ###savg.B0 = [prop.B, prop.B1p, prop.B2p, prop.B3p];
         savg.R2 = R2;
         savg.Efit = Efit;
      endif
   endif

   if (LOG > 1)
      # Get the spinodal properties:
      nreject3 = Qnew = 0;
      pmean3 = pstd3 = zeros(2,1);
      [spavg] = strainspinodal(cavg,V0,Vrange,strain);
      if (spavg.err > 0)
         printf('\nSpinodal point missing in the average polynomial\n');
         printf('Check again your data, please!\n');
      else
         printf('\nSpinodal properties of the average polynomial:\n');
         printf('(Vsp,psp,Bsp): %.6f %.6f %.6f\n',...
                spavg.Vsp, spavg.psp*hybohr3togpa, spavg.Bsp*hybohr3togpa);
         printf('\nSpinodal properties of the main polynomials:\n');
         printf('--i- npol data npar --SSerr-- ---w---- ----Vsp---- ----psp---- ----Bsp----\n');
         for k = 1 : npol
            k1 = isrt(k);
            [sk] = strainspinodal(pol{k1}.c,V0,Vrange,strain);
            if (sk.err != 0)
               nreject3++;
            else
               Qnew += pol{k1}.w;
               if (k <= 25)
                  printf('%4d %4d %4d %4d', k, k1, pol{k1}.data, pol{k1}.order);
                  printf(' %9.2e %8.6f', pol{k1}.SSerr, pol{k1}.w);
                  printf(' %11.6f', sk.Vsp);
                  printf(' %11.6f', sk.psp*hybohr3togpa);
                  printf(' %11.6f', sk.Bsp*hybohr3togpa);
                  printf('\n');
               endif
               pr = [sk.Vsp; sk.psp];
               pmean3 = pmean3 .+ pr * pol{k1}.w;
               pstd3 = pstd3 .+ (pr.^2) * pol{k1}.w;
            endif
         endfor
         pmean3 *= 1/Qnew; pstd3 *= 1/Qnew; 
         pstd3 = sqrt(pstd3 - pmean3.^2);
         printf('Renormalization? Qnew: %.6e\n', Qnew);
         printf('nreject3: %d\n', nreject3);
         printf('avg-Vsp (mean,std): %.6f %.6f\n', pmean3(1), pstd3(1));
         printf('avg-psp (mean,std): %.6f %.6f\n',...
                pmean3(2)*hybohr3togpa, pstd3(2)*hybohr3togpa);
      endif
   endif

   if (nargout > 1)
      # Putting error bars to the pressure, bulk modulus, etc at all volumes
      savg.Emean = savg.Estd = zeros(size(V));
      savg.pmean = savg.pstd = zeros(size(V));
      savg.Bmean = savg.Bstd = zeros(size(V));
      savg.B1pmean = savg.B1pstd = zeros(size(V));
      savg.B2pmean = savg.B2pstd = zeros(size(V));
      savg.B3pmean = savg.B3pstd = zeros(size(V));
      for k = 1 : npol
         k1 = isrt(k);
         [prop] = straineval(pol{k1}.c, V0, V, strain);
         savg.Emean = savg.Emean .+ prop.E * pol{k1}.w;
         savg.Estd  = savg.Estd  .+ (prop.E.^2) * pol{k1}.w;
         savg.pmean = savg.pmean .+ prop.p * pol{k1}.w;
         savg.pstd  = savg.pstd  .+ (prop.p.^2) * pol{k1}.w;
         savg.Bmean = savg.Bmean .+ prop.B * pol{k1}.w;
         savg.Bstd  = savg.Bstd  .+ (prop.B.^2) * pol{k1}.w;
         savg.B1pmean = savg.B1pmean .+ prop.B1p * pol{k1}.w;
         savg.B1pstd  = savg.B1pstd  .+ (prop.B1p.^2) * pol{k1}.w;
         savg.B2pmean = savg.B2pmean .+ prop.B2p * pol{k1}.w;
         savg.B2pstd  = savg.B2pstd  .+ (prop.B2p.^2) * pol{k1}.w;
         savg.B3pmean = savg.B3pmean .+ prop.B3p * pol{k1}.w;
         savg.B3pstd  = savg.B3pstd  .+ (prop.B3p.^2) * pol{k1}.w;
      endfor
      savg.Estd = sqrt(savg.Estd - savg.Emean.^2);
      savg.pstd = sqrt(savg.pstd - savg.pmean.^2);
      savg.Bstd = sqrt(savg.Bstd - savg.Bmean.^2);
      savg.B1pstd = sqrt(savg.B1pstd - savg.B1pmean.^2);
      savg.B2pstd = sqrt(savg.B2pstd - savg.B2pmean.^2);
      savg.B3pstd = sqrt(savg.B3pstd - savg.B3pmean.^2);
      if (LOG > 2)
         f = hybohr3togpa;
         errorbar(savg.pmean*f, savg.Bmean*f, savg.pstd*f, savg.Bstd*f, '~>');
         grid('on'); xlabel('pressure (GPa)'); ylabel('bulk modulus (GPa)');
      endif
   endif

endfunction
