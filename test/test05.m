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

% Test: gluing together two data files and using an analytical EOS to
% extrapolate to smaller volumes.

global nelectrons

addpath("../src/");

# read the central dataset
[v1,e1] = readEVdata('mgo-sety1.dat',1);
vref = median(v1);

# determine the equilibrium properties
[cf,sf] = avgstrainfit(v1,e1,vref,:,:,:,1);
plot(v1,e1,'ob', v1,sf.Efit,'-r');
xlabel('V (bohr^3)'); ylabel('E (hartree)'); grid('on');
print('test05a.eps', '-FHelvetica:38', '-depsc');

# read the additional datasets for very small and very large volumes
[v2,e2] = readEVdata('mgo-sety1-compress.dat',1);
[v3,e3] = readEVdata('mgo-sety1-expand.dat',1);
v = [v2;v1;v3]; e = [e2;e1;e3];
# volumes are already sorted, otherwise:
# [v,isort] = sort(v); e = e(isort);

# evaluate (do not fit) Vinet and AP2 EOS with the equilibrium
# parameters on the whole volume range
pin = sf.eqmean([2,1,3,4])
rvinet = nlf(v, e, 'vinet', pin, 0, 1);
rap2 = nlf(v, e, 'ap2', pin, 0, 1);

plot(v,e,'ob;Data;', v,rvinet.Efit,'-r;Vinet;', v,rap2.Efit, '-g;AP2;');
xlabel('V (bohr^3)'); ylabel('E (hartree)'); grid('on');
print('test05b.eps', '-FHelvetica:38', '-depsc');
