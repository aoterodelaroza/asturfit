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

function [p,v,t] = skimPVT(p,v,t,mode='square',a1=0,a2=0)
% function [p,V,T] = skimPVT(p,V,T,{mode,a1,a2})
%
% skimPVT - discard p-v-t points based on some low-p high-t criterion.
%
% The p, v and t arrays must be 1D and have the same size. There are
% two discarding modes: 'square' (default), that discards points that
% are both low pressure (below a1 GPa) and high temperature (above a2
% K). a1 and a2 are required parameters for the 'square' mode.
%
% The 'debye' mode estimates the Debye temperature from the p(V)
% isotherm at temperature a1 (default min(t)) assuming poisson ratio
% 0.25. Then, for each pressure, points above the corresponding Debye
% temperature times a factor a2 are discarded. a2 defaults to 1.5.
% 
% Using the 'debye' mode requires two global variables: natoms (the
% number of atoms associated to the volume v) and mm (the molecular
% mass of the formula).

  %% Physical constants and conversion factors
  pckbau = 3.166815d-6; %% Boltzmann constant (Hy/K)
  pcamu = 1.660538782d-24; %% kg / amu
  pcme = 9.10938215d-28;   %% mass of an electron (kg)
  amu2au = pcamu/pcme;     %% amu to atomic units conversion 
  rybohr3togpa = 14710.50498740275538944426; %% ry/bohr^3 to GPa
  hybohr3togpa = rybohr3togpa * 2; %% atomic units to GPa

  if (strcmp(mode,'square'))
    %% Cut off a square of the input data in the region of high T and low p.
    %% The a1 is the pressure beyond which, no point is discarded.
    %% For the low-pressure region, temperatures above a2 are discarded.
    if (a1 == 0 || a2 == 0) 
      error("Linear mode requires 2 parameters (pmax, tmin)")
    endif
    pmax = a1;
    tmin = a2;
    kgood = find(p>pmax | t<tmin);
    v = v(kgood); p = p(kgood); t = t(kgood);
  elseif(strcmp(mode,'debye'))
    %% Use a scaled Debye temperature to discard a region of low p and high T
    %% The a1 is the reference temperature (Default: min(t))
    %% The a2 is the Debye temperature scaling factor for the cut (Def. 1.5)

    %% Global natoms and molecular mass to calculate Debye temperature
    global natoms
    global mm

    %% Set defaults for Debye cut.
    if (a1 == 0)
      tref = min(t);
    else
      tref = a1;
    endif
    if (a2 == 0)
      tdfac = 1.5;
    else
      tdfac = a2;
    endif

    %% Check that the user has set the variables
    if (!exist("natoms","var") || natoms < 1)
      error("skimPVT: requires global natoms")
    elseif (!exist("mm","var") || mm < 0.0001)
      error("skimPVT: requires global mm")
    endif

    %% Fit the cold isotherm using eulerian strain and Vref at the min. p
    iref = find(t == tref);
    if (length(iref) < 7)
      error('skimPVT: too few cold isotherm points in the dataset!')
    endif
    [pcold,icold] = sort(p(iref));
    vcold = v(iref)(icold);
    pcold = p(iref)(icold);
    vref = vcold(1);
    [cavg,savg] = pvavgstrainfit(pcold,vcold,vref,:,:,'eulerian',0);
    res = pvstraineval(cavg,vref,vcold,'eulerian');

    %% Calculate the Debye temperature using the bulk modulus
    td = (6*pi^2*natoms*vcold.^2).^(1/3) ./ pckbau .* 0.859949 .* sqrt(max(res.E2v/(mm*amu2au)/hybohr3togpa,0));

    %%for i = 1:length(td)
    %%  printf("%.15f %.15f\n",p(i),tdfac*td(i));
    %%endfor 
    %% Discard points above the tdfac * Debye temperature.
    kgood = false(size(v));
    for i = 1:length(pcold)
      kgood(find(p == pcold(i) & t<tdfac*td(i))) = true;
    endfor
    v = v(kgood); p = p(kgood); t = t(kgood);
  else
    error("Unknown mode: %s",mode)
  endif
  

endfunction
