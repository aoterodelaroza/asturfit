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
function [data] = asturfit (filein)
% function [data] = asturfit (filein)
%      _____            __                 _____.__  __   
%     /  _  \   _______/  |_ __ __________/ ____\__|/  |_ 
%    /  /_\  \ /  ___/\   __\  |  \_  __ \   __\|  \   __\
%   /    |    \\___ \  |  | |  |  /|  | \/|  |  |  ||  |  
%   \____|__  /____  > |__| |____/ |__|   |__|  |__||__|  
%           \/     \/                                     
%
% asturfit - standard fitting task.
%
% Required input variables:
% filein: name of the input data file.
%
% Output (if any):
% data:
%

global nelectrons
hybohr3togpa = 2*14710.50498740275538944426;

if (~exist("avgstrainfit", "file"))
   error("AsturFit is not correctly installed! Check ReadMe!");
endif

tic;

s = {'       d8888          888                      .d888 d8b 888   ';
     '      d88888          888                     d88P"  Y8P 888   ';
     '     d88P888          888                     888        888   ';
     '    d88P 888 .d8888b  888888 888  888 888d888 888888 888 888888';
     '   d88P  888 88K      888    888  888 888P"   888    888 888   ';
     '  d88P   888 "Y8888b. 888    888  888 888     888    888 888   ';
     ' d8888888888      X88 Y88b.  Y88b 888 888     888    888 Y88b. ';
     'd88P     888  88888P"  "Y888  dY88888 888     888    888  "Y888';
     '                                                               ';
     '          Victor Lua~na and Alberto Otero-de-la-Roza           ';
     '                                                               '};
for i = 1 : length(s)
   printf("***  %s  ***\n", s{i});
endfor

n = max(strfind(filein,'.'));
if (!isempty(n))
  rootname = strtrunc(filein,n-1);
else
  rootname = filein
endif
printf("filein ----> %s\n", filein);
printf("rootname --> %s\n", rootname);

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
i = 0;
do
   olderr = newerr;
   [cf,sf] = avgstrainfit(v,e,vref,n,:,strain,0);
   degree(++i) = n; degerrs(i,:) = sf.eqstd(rk);
   newerr = (sf.eqstd(rk)./sf.eqmean(rk)).^2;
   printf("%6d", n);
   printf("  %.3e", newerr);
   printf("\n");
   n++;
until (sum(newerr > 1.5*olderr)>1 || n > 14)
[m,im] = min(degerrs);
nbest = degree(max(im));
printf("\nBest degree found for averaging: %d\n", nbest);
[cf,sf] = avgstrainfit(v,e,vref,nbest,:,strain,1);
newerr = (sf.eqstd(rk)./sf.eqmean(rk)).^2;

# Is this a noisy dataset?
worse = 1 / max(newerr);

# Create a file with the important data:
fileout = sprintf("%s-all.dat", rootname);
lu = fopen(fileout, "w+t");
fprintf(lu, "# volume bohr^3\n");
fprintf(lu, "# energy hartree\n");
fprintf(lu, "# pressure GPa\n");
fprintf(lu, "# bulk_modulus GPa\n");
fprintf(lu, "# z 1\n");
fprintf(lu, "# nelectrons %d\n", nelectrons);
prop = {"volume", "energy", "sigma(E)", "pressure", "sigma(p)", "B1p", "sigma(B1p)", "B2p", "sigma(B2p)"};
for nc = 1 : length(prop)
   fprintf(lu, "# column %02d: %s\n", nc, prop{nc});
endfor
for i = 1 : length(v)
   fprintf(lu, "%15.9f", v(i));
   fprintf(lu, "%18.9f", sf.Emean(i));
   fprintf(lu, "%12.9f", sf.Estd(i));
   fprintf(lu, "%15.6f", sf.pmean(i)*hybohr3togpa);
   fprintf(lu, "%9.6f",  sf.pstd(i)*hybohr3togpa);
   fprintf(lu, "%15.6f", sf.Bmean(i)*hybohr3togpa);
   fprintf(lu, "%9.6f",  sf.Bstd(i)*hybohr3togpa);
   fprintf(lu, "%15.6f", sf.B1pmean(i));
   fprintf(lu, "%9.6f",  sf.B1pstd(i));
   fprintf(lu, "%15.6f", sf.B2pmean(i)/hybohr3togpa);
   fprintf(lu, "%9.6f",  sf.B2pstd(i)/hybohr3togpa);
   fprintf(lu, "\n");
end
fclose(lu);
fprintf("File written with fitted data: %s\n", fileout);

# Figures will not appear on screen:
set(gcf(),"visible","off");

# Figures: E(V), p(V), B(p), B(V), B1(p), B*B2(p)
fileplot = sprintf("%s-ev.eps", rootname);
plot(v,e,'ob', v,sf.Efit,'-r');
xlabel('V (bohr^3)'); ylabel('E (hartree)'); grid('on');
print(fileplot, '-FHelvetica:22', '-depsc');
printf("E(V) plot created in: %s\n", fileplot);

fileplot = sprintf("%s-pv.eps", rootname);
f = hybohr3togpa;
errorbar(v, sf.pmean*f, sf.pstd*f, '~');
xlabel('V (bohr^3)'); ylabel('p (GPa)'); grid('on');
print(fileplot, '-FHelvetica:22', '-depsc');
printf("p(V) plot created in: %s\n", fileplot);

fileplot = sprintf("%s-bp.eps", rootname);
f = hybohr3togpa;
errorbar(sf.pmean*f, sf.Bmean*f, sf.pstd*f, sf.Bstd*f, '~>');
xlabel('p (GPa)'); ylabel('B (GPa)'); grid('on');
print(fileplot, '-FHelvetica:22', '-depsc');
printf("B(p) plot created in: %s\n", fileplot);

fileplot = sprintf("%s-bv.eps", rootname);
plot(v,sf.Bmean*f,'-ob');
xlabel('V (bohr^3)'); ylabel('B (GPa)'); grid('on');
print(fileplot, '-FHelvetica:22', '-depsc');
printf("B(V) plot created in: %s\n", fileplot);

fileplot = sprintf("%s-b1p.eps", rootname);
rk = find(sf.pmean > 0);
errorbar(sf.pmean(rk)*f, sf.B1pmean(rk), sf.pstd(rk)*f, sf.B1pstd(rk), '~>');
###errorbar(sf.pmean*f, sf.B1pmean, sf.pstd*f, sf.B1pstd, '~>');
xlabel('p (GPa)'); ylabel("B^'"); grid('on');
print(fileplot, '-FHelvetica:22', '-depsc');
printf("B1(p) plot created in: %s\n", fileplot);

fileplot = sprintf("%s-bb2p.eps", rootname);
errorbar(sf.pmean*f, sf.Bmean.*sf.B1pmean, sf.pstd*f, sf.Bstd.*sf.B1pstd, '~>');
xlabel('p (GPa)'); ylabel("B B^{''}"); grid('on');
print(fileplot, '-FHelvetica:22', '-depsc');
printf("B*B2p(p) plot created in: %s\n", fileplot);

# Reactivate the viewing on screen:
set(gcf(),"visible","on");

toc

endfunction

## Use it as a script
if (!exist("argn"))
   if (nargin > 0)
      args = argv();
      for i = 1 : nargin
         asturfit(args{i});
      endfor
   else
      printf('Use as: asturfit file(s)\n');
   endif
endif
