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

function [cb,sb] = strainbootstrap(V, E, V0, ndeg=12, nsample=100, strain='eulerian', LOG=0)
% function [cb{,sb}] = strainbootstrap(V, E, V0 {,ndeg, nsample, strain, LOG})
%
% strainbootstrap - bootstrap statistical method applied to the strain
% polynomials. A number of random samples are created from the initial {V,E}
% data. A collection of strain polynomials of different degree are fitted
% to each sample, and the equilibrium properties are averaged over the
% samples and polynomials. This is a form of Monte carlo statistics.
%
% Required input variables:
% V: cell volume.
% E: cell energy.
% V0: volume reference. Usually the volume at the minimum of energy.
%
% Optional input variables (all have default values):
% {ndeg = 12}: maximum degree of the polynomial.
%              Polynomials up to this degree wil be used and averaged.
% {nsample = 100}: number of random samples to create and average.
% {strain = 'eulerian'}: strain form.
% {LOG = 0}: print internal information about the fitting.
%          0 ....... no print
%          1 ....... print basic information
%          2 ....... time the run
%
% Minimum output:
% cb: coefficients of the average fitting polynomial.
%
% Additional (optional) output:
% sb: data structure containing:
%    sb.mean(1:5): mean values of equilibrium {V,E,B,B',B'',B'''}.
%    sb.std(1:5): std dev of equilibrium {V,E,B,B',B'',B'''} values.
%    sb.outliers: list of the points declared as outliers.
%    sb.Efit: fitted values for the energy.
%    sb.R2: determination coefficient (1 for a exact fitting).
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: October 2010

   if (nargin < 3)
      print_usage ();
   elseif (length(V) != length(E))
      error('strainbootstrap: V and E must be vectors of the same length!')
   elseif (length(V) < 7)
      error('strainbootstrap: dataset must have at least 7 points!')
   elseif (abs(ndeg) <= 1)
      error('strainbootstrap: polynomials must be of second or higher degree!')
   endif

   if (LOG>1)
      tic
   endif

   n = length(V);
   Vrange = [min(V), max(V)];
   nreject = 0;
   npol = 0;
   s2 = zeros(nsample,1);
   Morder = 0;
   for k = 1 : nsample
      rr = rand(n,1);
      r = find(rr >= 0.5);
      nr = length(r);
      for nd = 2 : ndeg
         if (nr < abs(nd)+5)
            # Too few points for a meaningful fit
            nreject++;
         else
            [c] = strainfit(V(r), E(r), V0, nd, strain, 0);
            Morder = max(Morder,length(c)-1);
            [sm] = strainmin(c, V0, Vrange, strain);
            Efit = strainevalE(c, V0, V(r), strain);
            SSerr = sum((E(r)-Efit).**2);
            npol++;
            pol{npol}.c = c;
            s2(npol) = pol{npol}.SSerr = SSerr;
            pol{npol}.deg  = nd;
            pol{npol}.data = nr;
            pol{npol}.V0 = sm.Vmin;
            pol{npol}.smin = sm;
         endif
      endfor
   endfor

   # Get the polynomial weights:
   SSmin = min(s2);
   Q = 0;
   for k = 1 : npol
      w = (pol{k}.SSerr/SSmin) * (pol{k}.deg / pol{k}.data);
      pol{k}.w = exp(-w);
      Q += pol{k}.w;
   endfor
   ##ww = zeros(1,1:npol);
   for k = 1 : npol
      ww(k) = pol{k}.w = pol{k}.w / Q;
   endfor

   # Form the average polynomial:
   cb = zeros(1,Morder+1);
   for k = 1 : npol
      n = length(pol{k}.c) - 1;
      cb(Morder-n+1:Morder+1) += pol{k}.c * pol{k}.w;
   endfor

   # Extra output and optional LOG record:
   # If the average of polynomials has a minimum, analize the equilibrium
   # geometry.
   # Otherwise analyze the reference volume.
   if (LOG || nargout > 1)
      [avgmin] = strainmin(cb, V0, Vrange, strain);
      if (avgmin.err == 0)
         if (LOG>0)
            printf('\n\nsb2: BOOTSTRAP ANALYSIS OF %s STRAIN POLYNOMIALS\n', strain);
            printf('Volume reference (V0): %.6f\n', V0);
            printf('Fitting mode: average of polynomials up to %d degree\n', ndeg);
            printf('Samples tried, Pol. rejected/accepted:: %d %d %d\n', nsample, nreject, npol);
            printf('\n Properties at the minimum of each polynomial:\n');
            printf('--i- --npol-- data npar --SSerr-- ---w---- ----Vmin--- ----Emin--- ----Bmin--- ---B1min--- ---B2min---\n');
         endif
         [srt,isrt] = sort(ww,'descend');
         pmean = pstd = pmean2 = pstd2 = zeros(6,1);
         for k = 1 : npol
            k1 = isrt(k);
            smin = pol{k1}.smin;
            Vmin = smin.Vmin;
            [prop] = straineval(pol{k1}.c, V0, Vmin, strain);
            if (LOG>0 && k <= 25)
               # In the conversion of the bulk modulus and derivatives we
               # assume the units: volume (bohr^3), energy (Hy).
               hybohr3togpa = 2*14710.50498740275538944426;
               printf('%4d %8d %4d %4d', k, k1, pol{k1}.data, pol{k1}.deg);
               printf(' %9.2e %8.6f', pol{k1}.SSerr, pol{k1}.w);
               printf(' %11.6f', Vmin);
               printf(' %11.6f', prop.E);
               printf(' %11.6f', prop.B * hybohr3togpa);
               printf(' %11.6f', prop.B1p);
               printf(' %11.9f', prop.B2p / hybohr3togpa);
               printf(' %11.9f', prop.B3p / hybohr3togpa^2);
               printf('\n');
            endif
            pr = [Vmin; prop.E; prop.B; prop.B1p; prop.B2p; prop.B3p];
            pmean = pmean .+ pr * pol{k1}.w;
            pstd = pstd .+ (pr.^2) * pol{k1}.w;
            pmean2 = pmean2 .+ pr * (1/npol);
            pstd2 = pstd2 .+ (pr.^2) * (1/npol);
         endfor
         pstd = sqrt(pstd - pmean.^2);
         pstd2 = sqrt(pstd2 - pmean2.^2);
         if (LOG>0)
            printf('\nAverage properties (weighted polynomials):\n');
            printf('------ ---volume-- ---energy-- --B-(GPa)-- ----B1p---- B2p-(1/GPa) B3p--(1/GPa^2)\n');
            printf('-mean- %11.6f %11.6f %11.6f %11.6f %11.6f %14.9f\n'...
                  , pmean(1), pmean(2), pmean(3)*hybohr3togpa, pmean(4)...
                  , pmean(5)/hybohr3togpa, pmean(6)/hybohr3togpa^2);
            printf('stdvev %11.6f %11.6f %11.6f %11.6f %11.6f %14.9f\n'...
                  , pstd(1), pstd(2), pstd(3)*hybohr3togpa, pstd(4)...
                  , pstd(5)/hybohr3togpa, pstd(6)/hybohr3togpa^2);
            printf('\nAverage properties (equally weighted polynomials):\n');
            printf('------ ---volume-- ---energy-- --B-(GPa)-- ----B1p---- B2p-(1/GPa) B3p--(1/GPa^2)\n');
            printf('-mean- %11.6f %11.6f %11.6f %11.6f %11.6f %14.9f\n'...
                  , pmean2(1), pmean2(2), pmean2(3)*hybohr3togpa, pmean2(4)...
                  , pmean2(5)/hybohr3togpa, pmean2(6)/hybohr3togpa^2);
            printf('stdvev %11.6f %11.6f %11.6f %11.6f %11.6f %14.9f\n'...
                  , pstd2(1), pstd2(2), pstd2(3)*hybohr3togpa, pstd2(4)...
                  , pstd2(5)/hybohr3togpa, pstd2(6)/hybohr3togpa^2);
         endif
      else
         # Analyze and report the reference geometry:
         if (LOG>0)
            printf('\n\nsb2: BOOTSTRAP ANALYSIS OF %s STRAIN POLYNOMIALS\n', strain);
            printf('Volume reference (V0): %.6f\n', V0);
            printf('Fitting mode: average of polynomials up to %d degree\n', ndeg);
            printf('Samples tried/rejected/accepted:: %d %d %d\n', nsample, nreject, npol);
            printf('The average polynomial has no minimum in the data volume range.\n');
            printf('\nProperties at the reference volume: %.9f\n', V0);
            printf('--i- --npol-- data npar --SSerr-- ---w---- ----Eref--- ----Bref--- ---B1ref--- ---B2ref---\n');
         endif
         [srt,isrt] = sort(ww,'descend');
         pmean = pstd = pmean2 = pstd2 = zeros(6,1);
         for k = 1 : npol
            k1 = isrt(k);
            [prop] = straineval(pol{k1}.c, V0, V0, strain);
            if (LOG>0 && k <= 25)
               # In the conversion of the bulk modulus and derivatives we
               # assume the units: volume (bohr^3), energy (Hy).
               hybohr3togpa = 2*14710.50498740275538944426;
               printf('%4d %8d %4d %4d', k, k1, pol{k1}.data, pol{k1}.deg);
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
            pmean2 = pmean2 .+ pr * (1/npol);
            pstd2 = pstd2 .+ (pr.^2) * (1/npol);
         endfor
         pstd = sqrt(pstd - pmean.^2);
         pstd2 = sqrt(pstd2 - pmean2.^2);
         if (LOG>0)
            printf('\nAverage properties (weighted polynomials):\n');
            printf('\nAverage properties at the ref. volume: %.9f\n', V0);
            printf('------ ---energy-- --B-(GPa)-- ----B1p---- B2p-(1/GPa) B3p--(1/GPa^2)\n');
            printf('-mean- %11.6f %11.6f %11.6f %11.6f %14.9f\n'...
                  , pmean(2), pmean(3)*hybohr3togpa, pmean(4)...
                  , pmean(5)/hybohr3togpa, pmean(6)/hybohr3togpa^2);
            printf('stdvev %11.6f %11.6f %11.6f %11.6f %14.9f\n'...
                  , pstd(2), pstd(3)*hybohr3togpa, pstd(4)...
                  , pstd(5)/hybohr3togpa, pstd(6)/hybohr3togpa^2);
            printf('\nAverage properties (equally weighted polynomials):\n');
            printf('------ ---energy-- --B-(GPa)-- ----B1p---- B2p-(1/GPa) B3p--(1/GPa^2)\n');
            printf('-mean- %11.6f %11.6f %11.6f %11.6f %11.6f %14.9f\n'...
                  , pmean2(1), pmean2(2), pmean2(3)*hybohr3togpa, pmean2(4)...
                  , pmean2(5)/hybohr3togpa, pmean2(6)/hybohr3togpa^2);
            printf('stdvev %11.6f %11.6f %11.6f %11.6f %11.6f %14.9f\n'...
                  , pstd2(1), pstd2(2), pstd2(3)*hybohr3togpa, pstd2(4)...
                  , pstd2(5)/hybohr3togpa, pstd2(6)/hybohr3togpa^2);
         endif
      endif

      Efit = strainevalE(cb, V0, V, strain);
      SSerr = sum((E-Efit).^2);
      SStot = sum((E-mean(E)).^2);
      R2 = 1 - SSerr / SStot;
      if (nargout > 1)
         sb.mean = pmean;
         sb.std = pstd;
         sb.R2 = R2;
         sb.Efit = Efit;
      endif

      if (LOG>0)
         printf('\n\nProperties of the average polynomial:\n');
         printf('Vmin, Emin: %.6f %.6f\n', Vmin, prop.E);
         printf('SSerr, SStot, R2, 1-R2: %.9e %.9e %.12f %.2e\n'...
               , SSerr, SStot, R2, 1-R2);
         nd = zeros(npol,1);
         for i = 1:npol
             nd(i) = pol{i}.data;
         endfor
         printf('Data (min/max/mode/mean/std): %d %d %d %.4f %.4f\n'...
               , min(nd), max(nd), mode(nd), mean(nd), std(nd));
      endif

      # To detect outliers we check the residuals for each point, i.e.
      # the difference between the actual energy and the value predicted
      # by the average polynomial. A robust measurement of what is normal
      # for a residual is provided by the median (but not the mean, which
      # too sensitive to outliear with large residuals):
      if (nargout > 1 || LOG > 0)
         lambda = 10.0;
         Edif = abs(Efit - E);
         Ecentral = median(Edif);
         Edisp = median(abs(Edif-Ecentral));
         iout = find(Edif > Ecentral+lambda*Edisp);
         sb.outliers = iout;
         if (LOG>0)
            printf('\nDetect outlier points (lambda=%.6f):\n', lambda);
            printf('Ecentral (median): %.6e\n', Ecentral);
            printf('Edisp    (median): %.6e\n', Edisp);
            printf('-npt- out ---volume-- ---energy-- ---error---\n');
            for i = 1 : length(V)
               printf('%5d ', i);
               printf('%3d ', (Edif(i) > Ecentral+lambda*Edisp));
               printf('%11.6f %11.6f ', V(i), E(i));
               printf('%11.4e', Edif(i));
               printf('\n');
            endfor
            printf('\nOutliers (%d):', length(iout));
            for i = 1 : length(iout)
               printf(' %d,', iout(i));
            endfor
            printf('\n');
         endif
      endif

   endif

endfunction
