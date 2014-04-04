#! /usr/bin/octave3.2 -q
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

function status = allfits(filein)
%
% allfits: Perform a complete set of fittings to a set of data files.
%
% Run as:   allfits.m file1.dat file2.dat ...

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
     '|  Task: fitting (almost) all the EOS types to the E(V) data  |';
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

nltypes = {'bm3', 'bm4', 'bm5', 'pt3', 'pt4', 'pt5', 'murn', 'ap2'};

for i = 1 : length(nltypes)
   nlres{i} = nlf(v,e,nltypes{i},[],1,1);
endfor

sttypes = {'eulerian', 'natural', 'lagrangian', 'infinitesimal'};

for i = 1 : length(sttypes)
   # Strain fitting. Decide up to wich degree do the average.
   n = 4;
   newerr = ones(4,1)*1e10;
   rk = [1 3 4 5];
   strain = sttypes{i};
   printf("\n\n");
   printf("+++++++++++++++++++++++++++++++++++++++++\n");
   printf("Fit to polynomials of %s strain\n", strain);
   printf("+++++++++++++++++++++++++++++++++++++++++\n");
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
   until (sum(newerr > 1.2*olderr)>1 | n > 18)
   [cf,sf] = avgstrainfit(v,e,vref,--n,:,strain,1);
endfor

toc

endfunction

## Use it as a script
if (!exist("argn"))
   if (nargin > 0)
      args = argv();
      for i = 1 : nargin
         allfits(args{i});
      endfor
   else
      printf('Use as: allfits file(s)\n');
   endif
endif
