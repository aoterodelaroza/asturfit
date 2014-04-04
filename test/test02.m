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

function status = test02 (filein, rootname)
%
% TEST02: Standard fit to a set of data files.
%
% Run as:   test02.m file1.dat file2.dat ...

global nelectrons
hybohr3togpa = 2*14710.50498740275538944426;
status = 0;

addpath("../src/");

printf("+++++++++++++++++++++++++++++++\n");
printf("test02: filein and rootname --> %s, %s\n", filein, rootname);
printf("+++++++++++++++++++++++++++++++\n");

[v,e] = readEVdata(filein,1);
vref = median(v);

# Strain fitting. Decide up to wich degree do the average.
n = 4;
newerr = ones(4,1)*1e10;
rk = [1 3 4 5];
strain = "eulerian";
printf("Fit to polynomials of %s strain\n", strain);
printf("Select the best number of polynomials to average:\n");
printf("Checking err(X) = (std(X)/mean(X))^2 for:\n");
printf(" --n-- --volume-- ----B----- ----B1p--- ----B2p---\n");
do
   olderr = newerr;
   [cf,sf] = avgstrainfit(v,e,vref,n,:,strain,0);
   newerr = (sf.eqstd(rk)./sf.eqmean(rk)).^2;
   printf("%6d", n);
   printf("  %.3e", newerr);
   printf("\n");
   n++;
until (sum(newerr > 1.5*olderr)>1 | n > 18)
[cf,sf] = avgstrainfit(v,e,vref,--n,:,strain,1);

# Figures will not appear on screen:
set(gcf(),"visible","off");

# Figures: E(V), p(V), B(p), B(V), B1(p), B*B2(p)
fileplot = sprintf("%sev.eps", rootname);
plot(v,e,"ob", v,sf.Efit,"-r");
xlabel("V (bohr^3)"); ylabel("E (hartree)"); grid("on");
print(fileplot, "-FHelvetica:22", "-depsc");
printf("E(V) plot created in: %s\n", fileplot);

fileplot = sprintf("%spv.eps", rootname);
f = hybohr3togpa;
errorbar(v, sf.pmean*f, sf.pstd*f, "~");
xlabel("V (bohr^3)"); ylabel("p (GPa)"); grid("on");
print(fileplot, "-FHelvetica:22", "-depsc");
printf("p(V) plot created in: %s\n", fileplot);

fileplot = sprintf("%sbp.eps", rootname);
f = hybohr3togpa;
errorbar(sf.pmean*f, sf.Bmean*f, sf.pstd*f, sf.Bstd*f, "~>");
xlabel("p (GPa)"); ylabel("B (GPa)"); grid("on");
print(fileplot, "-FHelvetica:22", "-depsc");
printf("B(p) plot created in: %s\n", fileplot);

fileplot = sprintf("%sbv.eps", rootname);
plot(v,sf.Bmean*f,"-ob");
xlabel("V (bohr^3)"); ylabel("B (GPa)"); grid("on");
print(fileplot, "-FHelvetica:22", "-depsc");
printf("B(V) plot created in: %s\n", fileplot);

fileplot = sprintf("%sb1p.eps", rootname);
errorbar(sf.pmean*f, sf.B1pmean, sf.pstd*f, sf.B1pstd, "~>");
xlabel("p (GPa)"); ylabel("B^'"); grid("on");
print(fileplot, "-FHelvetica:22", "-depsc");
printf("B1(p) plot created in: %s\n", fileplot);

fileplot = sprintf("%sbb2p.eps", rootname);
errorbar(sf.pmean*f, sf.Bmean.*sf.B1pmean, sf.pstd*f, sf.Bstd.*sf.B1pstd, "~>");
xlabel("p (GPa)"); ylabel("B B^{''}"); grid("on");
print(fileplot, "-FHelvetica:22", "-depsc");
printf("B*B2p(p) plot created in: %s\n", fileplot);

# Reactivate the viewing on screen:
set(gcf(),"visible","on");

endfunction

## Use it as a script
if (!exist("argn"))
   if (nargin > 0)
      args = argv();
      for i = 1 : nargin
         test02(args{i}, sprintf("test02-c%02d-",i));
      endfor
   else
      printf("Use as: test02.m file(s)\n");
   endif
endif
