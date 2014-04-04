#! /usr/bin/octave -q
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

% Compare the properties of bcc Li to Rb.

global nelectrons

addpath("../src/");

# elements
element = {"Na", "K", "Rb"};

hybohr3togpa = 2*14710.50498740275538944426;
fac = [1,1,hybohr3togpa,1,1/hybohr3togpa,1/hybohr3togpa^2];
 
# Table: equilibrium properties of the elements
printf("El");
printf(" -Eceil-");
printf(" --V0-(bohr^3)---");
printf(" --E0-(hartree)--");
printf(" ----B0-(GPa)----");
printf(" ------B1p0------");
printf(" --B2p0-(1/GPa)--");
printf("\n");
for i = 1 : length(element)
   el = tolower(element{i});
   filein = sprintf("w2k-lda-%s.dat", el);
   [v{i},e{i}] = readEVdata(filein,1);
   vref{i} = median(v{i});
   eceil{i} = ceil(max(e{i}));
   e{i} -= eceil{i};
   [cf{i},sf{i}] = avgstrainfit(v{i},e{i},vref{i},6,:,:,0);
   printf("%-2s", element{i});
   printf(" %7d", eceil{i});
   for k = 1 : 5
      printf(" %-16s", errformatf(sf{i}.eqmean(k)*fac(k), sf{i}.eqstd(k)*fac(k)));
   endfor
   printf("\n");
   rk = find(sf{i}.pmean > 0);
   vpos{i} = v{i}(rk);
endfor

# Figures will not appear on screen:
set(gcf(),"visible","off");

# Figure: E(V) for all the elements.
# We will make sure that all elements have the same number of points.
V = E = [];
for i = 1 : length(element)
   vv = linspace(min(v{i}), max(v{i}), 201)';
   ee = strainevalE(cf{i}, vref{i}, vv);
   V = [V, vv/sf{i}.eqmean(1)];
   E = [E, (ee-sf{i}.eqmean(2))/(sf{i}.eqmean(1)*sf{i}.eqmean(3))];
   key{i} = sprintf("-;%s;", element{i});
endfor
fileplot = "w2k-lda-ev.eps";
plot(V,E,key);
xlabel('V/V_0'); ylabel('(E-E_0)/(B_0V_0)'); grid('on');
print(fileplot, '-depsc');
printf("E(V/V_0) plot created in: %s\n", fileplot);

# Figures: B(p), B(V), B1p(p), B*B2p(p)
V = P = B = B1 = BB2 = [];
npts = 201; rk = 1:npts;
for i = 1 : length(element)
   vv = linspace(min(vpos{i}), max(vpos{i}), npts)';
   vv = vv(rk);
   s = straineval(cf{i}, vref{i}, vv);
   V = [V, vv/sf{i}.eqmean(1)];
   P = [P, s.p*hybohr3togpa];
   B = [B, s.B*hybohr3togpa];
   B1 = [B1, s.B1p];
   BB2 = [BB2, s.B.*s.B2p];
   key{i} = sprintf("-;%s;", element{i});
endfor

fileplot = "w2k-lda-bp.eps";
plot(P,B,key);
xlabel('p_{static} (GPa)'); ylabel('B (GPa)'); grid('on');
print(fileplot, '-depsc');
printf("B(p) plot created in: %s\n", fileplot);

fileplot = "w2k-lda-bv.eps";
plot(V,B,key);
xlabel('V/V_0'); ylabel('B (GPa)'); grid('on');
print(fileplot, '-depsc');
printf("B(V/V_0) plot created in: %s\n", fileplot);

fileplot = "w2k-lda-b1p.eps";
plot(P,B1,key);
xlabel('p_{static} (GPa)'); ylabel("B^'"); grid('on');
print(fileplot, '-depsc');
printf("B1p(p) plot created in: %s\n", fileplot);

fileplot = "w2k-lda-b2pbp.eps";
plot(P,BB2,key);
xlabel('p_{static} (GPa)'); ylabel('B^{''}B'); grid('on');
print(fileplot, '-depsc');
printf("B2pB(p) plot created in: %s\n", fileplot);
