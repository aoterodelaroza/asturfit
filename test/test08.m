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

% Test08: Bootstrap fit to a data set with noise (disperse outlier points)

addpath("../src/");
global nelectrons
hybohr3togpa = 2*14710.50498740275538944426;

# read in a smooth data set
[v,e] = readEVdata('mgo-sety1.dat',1);
vref = median(v);

# fit to the smooth set
[c,s] = avgstrainfit(v,e,vref,10,:,:,1);

# now we add artificially some noise, only to the energy
noutliers = 10;
[vn,en] = noisify(v, e, noutliers, 0, 0.004, 1);

# lets use the bootstrap sampling and see if all outliers are detected
sample = 1000;
chance = (1-(1-0.5^noutliers)^sample);
printf("\n");
printf("Number of outliers (artificially added): %d\n", noutliers);
printf("Bootstrap sample size: %d\n", sample);
printf("Chance of finding all outliers is %.6e\n", chance);
[cb,sb] = strainbootstrap(vn, en, vref, :, sample, :, 1);

bad = sb.outliers;
good = 1:length(vn); good(bad) = [];
fileplot = 'test08-ev.eps';
plot(vn,sb.Efit,'-r;BS fit;', vn(good),en(good),'ob;good data;'...
    , vn(bad),en(bad),'om;bad data;', 'markersize', 12);
xlabel('V (bohr^3)'); ylabel('E (hartree)'); grid('on');
print(fileplot, '-FHelvetica:24', '-depsc');
printf('See the E(V) curve in the file: %s\n', fileplot);
