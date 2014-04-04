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

function [steps,Ecorr] = guessjumps3(V,E, deltaE=0, LOG=0)
% function [steps{,Ecorr}] = guessjumps3(V,E {, deltaE=0, LOG=0})
%
% guessjumps3 - Detect jumps in the assumed continuous curve inherent
% to the (V,E) input dataset.
%
% Required input variables:
% (V,E): vectors containing points of a E(V) curve.
%
% Optional input variables (all have default values):
% {deltaE = 0}: neglect the steps smaller than this threshold. The default
%    is accepting all jumps, no matter how small.
% {LOG = 0}: print internal information about the fitting if (LOG>0).
%
% Minimum output:
% steps: cell array structure containing the detected jumps.
%    for k = 1 : nsteps
%    steps{k}.V ..... position of the jump.
%    steps{k}.E ..... estimated value of the step.
%
% Additional (optional) output:
% Ecorr: vector containing the corrected values for E.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: October 2010

   n = length(E);

   # linear interpolation
   E1 = zeros(size(E));
   r = 2:n;   E1(1) = interp1(V(r), E(r), V(1), 'linear', 'extrap');
   r = 1:n-1; E1(n) = interp1(V(r), E(r), V(n), 'linear', 'extrap');
   for i = 2 : n-1
      r = [1:i-1, i+1:n];
      E1(i) = interp1(V(r), E(r), V(i), 'linear');
   endfor

   # look for points that differ too much from the interpolation
   # (discard points for which the error is less than the average plus
   # one standard deviation):
   dE = E - E1;
   adE = abs(dE);
   mdE = abs(mean(dE));
   sdE = abs(std(dE));
   ierr = find(adE > mdE+sdE);

   if (LOG > 1)
      # slope
      sl = zeros(size(E));
      sl(1) = (E(2)-E(1))/(V(2)-V(1));
      sl(n) = (E(n)-E(n-1))/(V(n)-V(n-1));
      for i = 2 : n-1
         slp = (E(i)-E(i-1))/(V(i)-V(i-1));
         sln = (E(i+1)-E(i))/(V(i+1)-V(i));
         sl(i) = (slp+sln)/2;
      endfor
      plot(V,adE./abs(sl));
      xlabel('V (bohr^3)'); ylabel('Desv/slope');
      fileplot = 'gj3-step1.eps';
      print(fileplot, '-FHelvetica:28', '-depsc');
      printf('Plot file produced: %s\n', fileplot);
   endif

   # look for pairs of consecutive points
   ni = length(ierr);
   ip = find((ierr(2:ni)-ierr(1:ni-1)) == 1);
   ierr2 = ierr(ip);
   
   # determine position and step of the detected discontinuities
   jumps = 0;
   for i = 1 : length(ierr2)
      i1 = ierr2(i); i2 = i1 + 1;
      Vj = (V(i2)+V(i1))/2;
      left = find(V < Vj);
      right = find(V > Vj);
      if (length(left)>1 & length(right)>1)
         Ej = interp1(V(right),E(right),Vj,'linear','extrap') \
            - interp1(V(left),E(left),Vj,'linear','extrap');
         if (abs(Ej) >= deltaE)
            steps{++jumps}.V = Vj;
            Ejump(jumps) = steps{jumps}.E = Ej;
         endif
      endif
   endfor
   if (LOG > 0)
      [Esrt,isrt] = sort(Ejump, 'descend');
      printf('DBG(gj3):\n');
      printf('Interpolation error (mean,stddev): %.3e %.3e\n', mdE, sdE);
      printf('Jumps detected: %d\n', jumps);
      printf('       -jump- ----Vjump--- -----Ejump----\n');
      for j = 1 : jumps
         j1 = isrt(j);
         printf('%6d %6d %12.6f %14.6e\n', j, j1, steps{j1}.V, steps{j1}.E);
      endfor
      printf('\n');
      if (LOG > 1)
         plot(V,adE,';desv;');
         xlabel('V (bohr^3)'); ylabel('|E-E_{inter}| (Hy)');
         fileplot = 'gj3-step2.eps';
         print(fileplot, '-FHelvetica:28', '-depsc');
         printf('Plot file produced: %s\n', fileplot);
      endif
   endif
   if (nargout > 1)
      Ecorr = E;
      for j = 1 : jumps
         rango = find(V > steps{j}.V);
         Ecorr(rango) -= steps{j}.E;
      endfor
      if (LOG > 0)
         plot(V,E,'o', V,Ecorr,'-r');
         xlabel('V (bohr^3)'); ylabel('E (Hy)'); grid('on');
         fileplot = 'gj3-step3.eps';
         print(fileplot, '-FHelvetica:28', '-depsc');
         printf('Plot file produced: %s\n', fileplot);
      endif
   endif
endfunction
