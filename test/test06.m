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

% Test06: working with subsets of the initial data.

addpath("../src/");
global nelectrons
hybohr3togpa = 2*14710.50498740275538944426;

# read the initial dataset, former by 129 points
[v,e] = readEVdata('mgo-sety1.dat',1);
vref = median(v);

# select one every n points, for n being [1, 2, 4, 8]:
# determine the equilibrium properties for each set
f = [1, 1, hybohr3togpa, 1, 1/hybohr3togpa, 1/hybohr3togpa^2];
printf("-n- data --V-(bohr^3)-- --E-(hartree)- ---B-(GPa)----");
printf("  ----B'----  B''(1/GPa)-\n");
for n = [1, 2, 4, 8]
   vset = v(1:n:end); eset = e(1:n:end); nset = length(vset);
   [cf,sf] = avgstrainfit(vset,eset,vref,10,:,:,0);
   printf("%3d %4d ", n, nset);
   printf(" %-14s", errformatf(sf.eqmean(1), sf.eqstd(1)));
   printf(" %-14s", errformatf(sf.eqmean(2), sf.eqstd(2)));
   printf(" %-14s", errformatf(sf.eqmean(3)*f(3), sf.eqstd(3)*f(3)));
   printf(" %-10s", errformatf(sf.eqmean(4)*f(4), sf.eqstd(4)*f(4)));
   printf(" %-10s", errformatf(sf.eqmean(5)*f(5), sf.eqstd(5)*f(5)));
   printf("\n");
endfor
