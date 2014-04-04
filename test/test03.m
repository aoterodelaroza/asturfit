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

function status = test03 (filein, rootname)
%
% TEST03: Perform a complete set of fittings to a set of data files.
%
% Run as:   test03.m file1.dat file2.dat ...

global nelectrons
hybohr3togpa = 2*14710.50498740275538944426;
status = 0;

addpath("../src/");

printf("+++++++++++++++++++++++++++++++\n");
printf("test03: filein and rootname --> %s, %s\n", filein, rootname);
printf("+++++++++++++++++++++++++++++++\n");

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

endfunction

## Use it as a script
###if (!exist("__nargout__"))
###   if (__nargin__ > 0)
###      args = argv();
###      for i = 1 : length(args)
###         test03(args{i}, sprintf("test03-c%02d-",i));
###      endfor
###   else
###      print_usage();
###   endif
###endif
if (!exist("argn"))
   if (nargin > 0)
      args = argv();
      for i = 1 : nargin
         test03(args{i}, sprintf("test03-c%02d-",i));
      endfor
   else
      printf('Use as: test03.m file(s)\n');
   endif
endif
