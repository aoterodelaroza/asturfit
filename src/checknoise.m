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

function status = checknoise (filein)
% function status = checknoise (filein)
%
% checknoise - check the noise of a suspect datafile.
%
% Required input variables:
% filein: name of the input data file.
%
% Output (if any):
% status: code of the type of problem found in the dataset.
%     0 ---> data appears to be smooth.
%    10 ---> the error bars indicate some noise.
%    20 ---> very large error bars.
%    +1 ---> outliers have been found.
%    +2 ---> jumps have been found.
%    These codes are added out to form the final status. A status of
%    23, for instance, would mean very large error bars, with outliers
%    and jumps detected.
%

global nelectrons
hybohr3togpa = 2*14710.50498740275538944426;
status = 0;

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
     '+-------------------------------------------------------------+';
     '|     checknoise: check the noise of a suspect datafile       |';
     '+-------------------------------------------------------------+'};
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

# Strain fitting. Check the errorbars:
rk = [1 3 4 5];
strain = "eulerian";
[cf,sf] = avgstrainfit(v,e,vref,10,:,strain,1);
printf('The first symptom of a noisy datafile is a large set of error bars\n');
printf('Check carefully the errors from the previous fit and the E(V) plot\n');

# Figures will not appear on screen:
set(gcf(),"visible","off");

# Figures: E(V).
fileplot = sprintf("%s-ev.eps", rootname);
plot(v,e,'ob', v,sf.Efit,'-r');
xlabel('V (bohr^3)'); ylabel('E (hartree)'); grid('on');
print(fileplot, '-FHelvetica:22', '-depsc');
printf("E(V) plot created in: %s\n", fileplot);

# Is this a noisy dataset?
fiterr = abs(sf.eqstd./sf.eqmean);
worse = max(fiterr);
if (worse < 0.01)
   status = 0;
   printf("Datafile appears very smooth. Check the E(V) plot, anyway.\n");
   return
elseif (worse > 1.0)
   status += 20;
   printf("Noise is a serious problem in your data!\n");
   sample = 10000;
else
   status += 10;
   printf("Some level of noise exist. Check carefully the next results.\n");
   sample = 1000;
endif

# Try first a bootstrap fit to detect possible outliers:
printf("\nA bootstrap fit will be tried to detect possible outliers:\n");
printf("Bootstrap sample size: %d\n", sample);
[cb,sb] = strainbootstrap(v, e, vref, :, sample, :, 1);
bad = sb.outliers;
good = 1:length(v); good(bad) = [];
fileplot = sprintf("%s-evn.eps", rootname);
plot(v,sb.Efit,'-r;BS fit;', v(good),e(good),'ob;good data;' \
    , v(bad),e(bad),'om;bad data;', 'markersize', 12);
xlabel('V (bohr^3)'); ylabel('E (hartree)'); grid('on');
print(fileplot, '-FHelvetica:24', '-depsc');
printf('See the E(V) curve in the file: %s\n', fileplot);
if (length(bad) > 0)
   status += 1;
endif

# Remove the detected outliers:
vg = v(good); eg = e(good);
printf("\nRemove the outliers and keep %d points.\n", length(vg));

# Try detecting jumps in the denoised dataset:
precision = figures(eg);
[stj,ej] = guessjumps3(vg,eg,20*precision,1);
if (length(stj) > 0)
   status += 2;
endif

# Reactivate the viewing on screen:
set(gcf(),"visible","on");

toc

endfunction

## Use it as a script
if (!exist("argn"))
   if (nargin > 0)
      args = argv();
      for i = 1 : nargin
         checknoise(args{i});
      endfor
   else
      printf('Use as: checknoise file(s)\n');
   endif
endif
