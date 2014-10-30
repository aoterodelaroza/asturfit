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

% Test07: transition between two phases under high pressure

addpath("../src/");
global nelectrons
hybohr3togpa = 2*14710.50498740275538944426;

# read the data of the two phases
if (1)
   [vb1,eb1] = readEVdata('mgo-pbe-b1.dat',1);
   vrefb1 = median(vb1);
   [vb2,eb2] = readEVdata('mgo-pbe-b2.dat',1);
   vrefb2 = median(vb2);
else
   [vb1,eb1] = readEVdata('mgo-lda-b1.dat',1);
   vrefb1 = median(vb1);
   [vb2,eb2] = readEVdata('mgo-lda-b2.dat',1);
   vrefb2 = median(vb2);
endif

# fit both phases
[cb1,sb1] = avgstrainfit(vb1,eb1,vrefb1,10,:,:,1);
[cb2,sb2] = avgstrainfit(vb2,eb2,vrefb2,10,:,:,1);

# Figures will not appear on screen:
set(gcf(),"visible","off");

# plot the E(V) curves:
fileplot = 'test07-ev.eps';
plot(vb1,eb1,'or;;',vb1,sb1.Efit,'-r;B1 phase;'...
    ,vb2,eb2,'og;;',vb2,sb2.Efit,'-g;B2 phase;')
xlabel('V (bohr^3)'); ylabel('E (hartree)'); grid('on');
print(fileplot, '-FHelvetica:24', '-depsc');
printf('See the E(V) curve in the file: %s\n', fileplot);

# get the enthalpies and plot H(p) for both phases
hb1 = sb1.Efit + sb1.pmean .* vb1;
hb2 = sb2.Efit + sb2.pmean .* vb2;
pb1 = sb1.pmean*hybohr3togpa;
pb2 = sb2.pmean*hybohr3togpa;
fileplot = 'test07-hp.eps';
plot(pb1,hb1,'-r;B1 phase;', pb2,hb2,'-g;B2 phase;');
xlabel('p (GPa)'); ylabel('H (hartree)'); grid('on');
print(fileplot, '-FHelvetica:24', '-depsc');
printf('See the H(p) curve in the file: %s\n', fileplot);

# It is more effective to print the enthalpy difference between both
# phases, but we must be sure that both are defined in the same grid
# of pressures. We use interpolation in the {pb1,hb1} and {pb2,hb2} sets.
pp = linspace(max([min(pb1),min(pb2),0]), min([max(pb1),max(pb2)]), 101)';
hhb1 = interp1(pb1, hb1, pp);
hhb2 = interp1(pb2, hb2, pp);
fileplot = 'test07-dhp.eps';
plot(pp, hhb2-hhb1, '-b;B2 phase;');
xlabel('p (GPa)'); ylabel('H_{B2}-H_{B1} (hartree)'); grid('on');
print(fileplot, '-FHelvetica:24', '-depsc');
printf('See the DifH(p) curve in the file: %s\n', fileplot);

# Reactivate the viewing on screen of the plots:
set(gcf(),"visible","on");
