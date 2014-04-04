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

% Test09: Detecting jumps in a dataset

addpath("../src/");

# read in a smooth data set
[v,e] = readEVdata('mgo-sety1.dat',1);
vref = median(v);

# now we add artificially some jumps
# (keep the original data for further comparison)
vmin = min(v); vmax = max(v); vrange = vmax-vmin;
vj = v; ej = e;
jumps.v = [vmin+vrange/3, vmin+2*vrange/3];
jumps.e = [0.010, -0.011];
for i = 1 : length(jumps.v)
   printf('Adding a jump of %.6f for V>%.6f\n', jumps.e(i), jumps.v(i));
   rk = find(vj > jumps.v(i));
   ej(rk) += jumps.e(i);
endfor

# let's do the automatic detection
[st,ecorr] = guessjumps3(vj,ej,:,1);
fileplot = 'test09-ev.eps';
print(fileplot, '-FHelvetica:24', '-depsc');
printf('See the E(V) curve in the file: %s\n', fileplot);

# compare the original and the corrected energies
SSerr = sum((e-ecorr).^2);
SStot = sum((e-mean(e)).^2);
R2 = 1 - SSerr / SStot;
printf('\n');
printf('Global error of the corrected data:\n');
printf('Determination coef. (R2, 1 for a exact correction): %.9e\n', R2);
printf('1-R2: %.2e\n', SSerr / SStot);
