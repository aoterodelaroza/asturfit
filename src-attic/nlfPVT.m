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

function res = nlfPVT(p, V, T, mode="none", pin=[], fit=1, LOG=1)
% function res = nlfPVT(p, V, T {, mode="none", pin=[], fit=1, LOG=1})
%
% nlfPVT (non-linear fitting) - fitting to a traditional p(V,T) EOS.
%
% Required input variables:
% (p,V,T): column vectors containing the points of the p(V,T) surface.
%
% Optional input variables (all have default values):
% {mode = "none"}: select the EOS. The default value will only provide an
%    error message. Implemented modes:
%    * Birch-Murnaghan family: 'bm3', 'bm4', 'bm5';
%    * Poirier-Tarantola family: 'pt3', 'pt4', 'pt5';
%    * Murnaghan EOS: 'murn';
%    * Holzapfel AP2: 'ap2';
% {pin = []}: column vector with the start value for the EOS parameters. If
%   no values are entered, the routine will try to determine some starting
%   point. If mode='ap2', pin(1) *must* contain the value for Z, even
%   though this is not a true parameter and will not be used for the
%   fitting.
% {fit = 1}: in default mode the analytical EOS will be fitted to the
%   (V,E) data. In the fit=0 mode, the EOS will be evaluated with the
%   provided pin parameters (required in this mode) at the volumes
%   contained in V.
% {LOG = 1}: print internal information about the fitting if (LOG>0).
%
% Output:
%   res.pout .... column vector with the best estimate of the parameters.
%   res.konv .... 1 if the fitting converged and 0 if not.
%   res.Efit .... column vector with the values for E predicted by the fit.
%   res.R2 ...... determination coefficient (1 for a exact fitting).
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: November 2010

   global nelectrons

   mode = tolower(mode);
   if (nargin < 3 | strcmp(mode,"none"))
      print_usage();
   elseif (length(p) != length(V) | length(p) != length(T))
      error('nlf: p, V, and T must be vectors of the same length!')
   elseif (exist("leasqr") ~= 2)
      error('nlf: leasqr, from octaveforge/optim, MUST be in your path!');
   endif

   rybohr3togpa = 14710.50498740275538944426;
   hybohr3togpa = 14710.50498740275538944426 * 2;

   ###if (length(pin) < 1)
   ###   Tmin = min(T); imin = find(T == Tmin);
   ###   [pcold,icold] = sort(p(imin));
   ###   vcold = V(imin)(icold);
   ###   vref = max(vcold);
   ###   [cavg,savg] = pvavgstrainfit(pcold,vcold,vref,:,:,:,0);
   ###endif
   if (strcmp(mode,'bm3'))
      Tmin = min(T); imin = find(T == Tmin);
      [pcold,icold] = sort(p(imin));
      vcold = V(imin)(icold);
      vref = max(vcold);
      if (length(pin) < 3)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            [cavg,savg] = pvavgstrainfit(pcold,vcold,vref,:,:,:,0);
            # pin : [V, B, Bprime]
            pin = savg.eqmean([1,2,3]);
         endif
      endif
      xvar = vcold; yvar = pcold;
      name = "BM3: 3rd order Birch-Murnaghan EOS";
      par_name = {'V0 (bohr^3)', 'B0 (GPa)', 'B1p'};
      #par_fact = [1, hybohr3togpa, 1];
      par_fact = [1, 1, 1];
      fun = @bm3_p;
   elseif (strcmp(mode,'bm3-mgd1'))
      Tmin = min(T); imin = find(T == Tmin);
      [pcold,icold] = sort(p(imin));
      vcold = V(imin)(icold);
      vref = max(vcold);
      if (length(pin) < 4)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            [cavg,savg] = pvavgstrainfit(pcold,vcold,vref,:,:,:,0);
            # pin: [V_0, B_0, B_0', Theta0, gamma0]
            #        1    2    3      4       5
            pin = savg.eqmean([1,2,3]);
            pin(4) = 500;
            pin(5) = 1.5;
         endif
      endif
      xvar = [V, T]; yvar = p;
      name = "BM3-Mie-Gruneisen-Debye-1 EOS";
      par_name = {'V0 (bohr^3)', 'B0 (GPa)', 'B1p', 'Theta0 (K)', 'gamma0'};
      par_fact = [1, 1, 1, 1, 1];
      fun = @bm3_mgd1_p;
   elseif (strcmp(mode,'bm3-mgd2'))
      Tmin = min(T); imin = find(T == Tmin);
      [pcold,icold] = sort(p(imin));
      vcold = V(imin)(icold);
      vref = max(vcold);
      if (length(pin) < 4)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            [cavg,savg] = pvavgstrainfit(pcold,vcold,vref,:,:,:,0);
            # pin: [V_0, B_0, B_0', Theta0, gamma0, q0]
            #        1    2    3      4       5     6
            pin = savg.eqmean([1,2,3]);
            pin(4) = 500;
            pin(5) = 1.5;
            pin(6) = 1.4;
         endif
      endif
      xvar = [V, T]; yvar = p;
      name = "BM3-Mie-Gruneisen-Debye-2 EOS";
      par_name = {'V0 (bohr^3)', 'B0 (GPa)', 'B1p', \
                  'Theta0 (K)', 'gamma0', 'q0'};
      par_fact = [1, 1, 1, 1, 1, 1];
      fun = @bm3_mgd2_p;
   elseif (strcmp(mode,'bm3-mgd3p'))
      Tmin = min(T); imin = find(T == Tmin);
      [pcold,icold] = sort(p(imin));
      vcold = V(imin)(icold);
      vref = max(vcold);
      if (length(pin) < 4)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            [cavg,savg] = pvavgstrainfit(pcold,vcold,vref,:,:,:,0);
            # pin: [V_0, B_0, B_0', Theta0, gamma0, a, b]
            #        1    2    3      4       5     6  7
            pin = savg.eqmean([1,2,3]);
            pin(4) = 500;
            pin(5) = 1.5;
            pin(6) = 0.8;
            pin(7) = 1.1;
         endif
      endif
      xvar = [V, T]; yvar = p;
      name = "BM3-Mie-Gruneisen-Debye-3' EOS";
      par_name = {'V0 (bohr^3)', 'B0 (GPa)', 'B1p', \
                  'Theta0 (K)', 'gamma0', 'a0', 'b0'};
      par_fact = [1, 1, 1, 1, 1, 1, 1];
      fun = @bm3_mgd3p_p;
   elseif (strcmp(mode,'bm4-mgd3p'))
      Tmin = min(T); imin = find(T == Tmin);
      [pcold,icold] = sort(p(imin));
      vcold = V(imin)(icold);
      vref = max(vcold);
      if (length(pin) < 4)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            [cavg,savg] = pvavgstrainfit(pcold,vcold,vref,:,:,:,0);
            # pin: [V_0, B_0, B_0', B_0'', Theta0, gamma0, a, b]
            #        1    2    3     4       5       6     7  8
            pin = savg.eqmean([1,2,3,4]);
            pin(5) = 500;
            pin(6) = 1.5;
            pin(7) = 0.8;
            pin(8) = 1.1;
         endif
      endif
      xvar = [V, T]; yvar = p;
      name = "BM4-Mie-Gruneisen-Debye-3' EOS";
      par_name = {'V0 (bohr^3)', 'B0 (GPa)', 'B1p', 'B2p (1/GPa)', \
                  'Theta0 (K)', 'gamma0', 'a0', 'b0'};
      par_fact = [1, 1, 1, 1, 1, 1, 1, 1];
      fun = @bm4_mgd3p_p;
   elseif (strcmp(mode,'bm3-mg1'))
      Tmin = min(T); imin = find(T == Tmin);
      [pcold,icold] = sort(p(imin));
      vcold = V(imin)(icold);
      vref = max(vcold);
      if (length(pin) < 4)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            [cavg,savg] = pvavgstrainfit(pcold,vcold,vref,:,:,:,0);
            # pin: [V_0, B_0, B_0', Theta0, gamma0, q0, a0, m]
            #        1    2    3      4       5      6   7  8
            pin = savg.eqmean([1,2,3]);
            pin(4) = 300;
            pin(5) = 4;
            pin(6) = 1.5;
            pin(7) = 1e-5;
            pin(8) = 3;
         endif
      endif
      xvar = [V, T]; yvar = p;
      name = "BM3-Mie-Gruneisen EOS";
      par_name = {'V0 (bohr^3)', 'B0 (GPa)', 'B1p', 'Theta0 (K)', \
                  'gamma0', 'q0', 'a0 (1/K)', 'm'};
      par_fact = [1, 1, 1, 1, 1, 1, 1, 1];
      fun = @bm3_mg1_p;
   else
      error('nlfPVT: EOS mode unknown!')
   endif

   if (fit==0)
      res.pfit = fun(xvar,pin);
      res.pout = pin;
      res.konv = 0;
      res.R2 = 0;
      if (LOG > 0)
         printf('\n(nlf) Evaluating %s\n', name);
         printf('Parameters (%d)\n', length(pin));
         for i = 1 : length(pin)
            printf('%-16s', par_name{i});
            printf(' %16.6f', pin(i)*par_fact(i));
            printf('\n');
         endfor
      endif
      if (size(yvar) == size(res.pfit))
         SSerr = sum((yvar-res.pfit).^2);
         SStot = sum((yvar-mean(yvar)).^2);
         res.R2 = 1 - SSerr / SStot;
         if (LOG>0)
            printf('SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
                  , SSerr, SStot, res.R2, 1-res.R2);
         endif
      endif
   else
      [f1, p1, kvg1, iter1, corp1, covp1, covr1, stdresid1, Z1, r21] \
      = leasqr(xvar, yvar, pin, fun, 1e-9, 30);
      SSerr = sum((yvar-f1).^2);
      SStot = sum((yvar-mean(yvar)).^2);
      R2 = 1 - SSerr / SStot;
      res.pout = p1;
      res.konv = kvg1;
      res.pfit = f1;
      res.R2 = R2;
      if (LOG > 0)
         printf('\nNon-linear fitting: %s\n', name);
         printf('Parameters (%d) start / converged\n', length(pin));
         for i = 1 : length(pin)
            printf('%-16s', par_name{i});
            printf(' %16.6f', pin(i)*par_fact(i));
            printf(' %16.6f', p1(i)*par_fact(i));
            printf('\n');
         endfor
         printf('Convergence (1=yes)? %d\n', kvg1);
         printf('Iterations: %d\n', iter1);
         printf('SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
               , SSerr, SStot, R2, 1-R2);
         printf('Correlation matrix of the parameters:\n');
         for i = 1 : rows(corp1)
            for j = 1 : i
               printf(' %12.6f', corp1(i,j));
            endfor
            printf('\n');
         endfor
      endif
   endif


endfunction

# Birch-Murnaghan, third order. p(V) expression.
# parameters: [V_0, B_0, B_0']
#               1    2    3
function pressure = bm3_p(V,par)
  xinv = par(1) ./ V;
  xinv23 = xinv.^(2/3)-1;
  xdif = xinv.^(7/3) - xinv.^(5/3);
  pressure = (xinv23 * (-3*(4-par(3))/4) + 1) .* xdif * (1.5*par(2));
endfunction

# BM3-Mie-Gruneisen model
# parameters: [V_0, B_0, B_0', Theta0, gamma0, q0, a0, m]
#               1    2    3      4       5      6   7  8
# Assumed units: V(bohr^3), T(K), p(GPa)
# Assumed reference: T0 = 300 K
function pressure = bm3_mg1_p(x,par)
  global natoms;
  T0 = 300;
  V = x(:,1); T = x(:,2);
  xx = V ./ par(1); xinv = par(1) ./ V;
  xinv23 = xinv.^(2/3)-1;
  xdif = xinv.^(7/3) - xinv.^(5/3);
  pr1 = (xinv23 * (-3*(4-par(3))/4) + 1) .* xdif * (1.5*par(2));
  tht = exp((-xx.^par(6)+1)*(par(5)/par(6)))*par(4);
  gamma = xx.^par(6) * par(5);
  R = 8.314472; bohr = 0.52917720859e-10; N_A = 6.02214179e23;
  n3R = xinv * (3*natoms*R/bohr^3/N_A/1e12/par(1));
  pr2 = (tht./(exp(tht./T)-1) .- tht./(exp(tht./T0)-1)) .* n3R .* gamma;
  pr3 = xx.^par(8) .* (T.*T-T0*T0) .* n3R * (par(7)*par(8)/2);
  pressure = pr1 + pr2 + pr3;
endfunction

# BM3-Mie-Gruneisen-Debye-1 model
# parameters: [V_0, B_0, B_0', Theta0, gamma0]
#               1    2    3      4       5
# Assumed units: V(bohr^3), T(K), p(GPa)
# Assumed reference: T0 = 300 K
function pressure = bm3_mgd1_p(x,par)
  global natoms;
  T0 = 300;
  V = x(:,1); T = x(:,2);
  x = V ./ par(1); xinv = par(1) ./ V;
  xinv23 = xinv.^(2/3)-1;
  xdif = xinv.^(7/3) - xinv.^(5/3);
  pr1 = (xinv23 * (-3*(4-par(3))/4) + 1) .* xdif * (1.5*par(2));
  R = 8.314472; bohr = 0.52917720859e-10; N_A = 6.02214179e23;
  gamma = par(5);
  theta = (par(1)./V).^gamma * par(4); 
  n3R = (T./V) * (gamma * 3*natoms*R/bohr^3/N_A/1e9);
  pr2 = (debye_3(theta./T) .- debye_3(theta/T0)) .* n3R;
  pressure = pr1 + pr2;
endfunction

# BM3-Mie-Gruneisen-Debye-2 model
# parameters: [V_0, B_0, B_0', Theta0, gamma0, q0]
#               1    2    3      4       5      6
# Assumed units: V(bohr^3), T(K), p(GPa)
# Assumed reference: T0 = 300 K
function pressure = bm3_mgd2_p(x,par)
  global natoms;
  T0 = 300;
  V = x(:,1); T = x(:,2);
  x = V ./ par(1); xinv = par(1) ./ V;
  xinv23 = xinv.^(2/3)-1;
  xdif = xinv.^(7/3) - xinv.^(5/3);
  pr1 = (xinv23 * (-3*(4-par(3))/4) + 1) .* xdif * (1.5*par(2));
  R = 8.314472; bohr = 0.52917720859e-10; N_A = 6.02214179e23;
  gamma = (V./par(1)).^par(6) * par(5);
  theta = exp(-(gamma.-par(5))/par(6)) * par(4);
  n3R = (T./V) .* gamma * (3*natoms*R/bohr^3/N_A/1e9);
  pr2 = (debye_3(theta./T) .- debye_3(theta/T0)) .* n3R;
  pressure = pr1 + pr2;
endfunction

# BM3-Mie-Gruneisen-Debye-3' model
# parameters: [V_0, B_0, B_0', Theta0, gamma0, a, b]
#               1    2    3      4       5     6  7
# Assumed units: V(bohr^3), T(K), p(GPa)
# Assumed reference: T0 = 300 K
function pressure = bm3_mgd3p_p(x,par)
  global natoms;
  T0 = 300;
  V = x(:,1); T = x(:,2);
  x = V ./ par(1); xinv = par(1) ./ V;
  xinv23 = xinv.^(2/3)-1;
  xdif = xinv.^(7/3) - xinv.^(5/3);
  pr1 = (xinv23 * (-3*(4-par(3))/4) + 1) .* xdif * (1.5*par(2));
  R = 8.314472; bohr = 0.52917720859e-10; N_A = 6.02214179e23;
  gamma = (((V./par(1)).^par(7)-1) * par(6) + 1) * par(5);
  theta = exp(-(gamma.-par(5))/par(7)) .* (V./par(1)).^((par(6)-1)*par(5)) * par(4);
  n3R = (T./V) .* gamma * (3*natoms*R/bohr^3/N_A/1e9);
  pr2 = (debye_3(theta./T) .- debye_3(theta/T0)) .* n3R;
  pressure = pr1 + pr2;
endfunction

# BM4-Mie-Gruneisen-Debye-3' model
# parameters: [V_0, B_0, B_0', B_0'', Theta0, gamma0, aa0, bb0]
#               1    2    3     4       5       6      7   8
# Assumed units: V(bohr^3), T(K), p(GPa)
# Assumed reference: T0 = 300 K
function pressure = bm4_mgd3p_p(x,par)
  global natoms;
  T0 = 300;
  v0 = par(1); b0 = par(2); b0p = par(3); b0pp = par(4);
  theta0 = par(5); gamma0 = par(6); aa0 = par(7); bb0 = par(8);
  V = x(:,1); T = x(:,2);
  f = ((v0./V).^(2/3)-1)/2;
  c2 = 4.5*b0*v0; c3 = c2*(b0p-4);
  c4 = (3/8)*b0*v0 * (9*(b0pp*b0+b0p^2) - 63*b0p + 143);
  pr1 = (f*2+1).^2.5 .* f .* (f.*(f*4*c4+3*c3)+2*c2) / (3*v0);
  R = 8.314472; bohr = 0.52917720859e-10; N_A = 6.02214179e23;
  gamma = (((V./v0).^bb0-1) * aa0 + 1) * gamma0;
  theta = exp(-(gamma.-gamma0)/bb0) .* (V./v0).^((aa0-1)*gamma0) * theta0;
  n3R = (T./V) .* gamma * (3*natoms*R/bohr^3/N_A/1e9);
  pr2 = (debye_3(theta./T) .- debye_3(theta./T0)) .* n3R;
  pressure = pr1 + pr2;
endfunction
