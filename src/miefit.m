% Copyright (C) 2011 Victor Lua~na and Alberto Otero-de-la-Roza
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

function res = miefit(p, V, T, Vref, Tref, mode1="none", mode2="none", pin1=[], pin2=[], fit=1, LOG=1)
% function res = miefit(p, V, T, Vref, Tref
%          {,mode1="none", mode2="none", pin1=[], pin2=[], fit=1, LOG=1})
%
% miefit (Mie-Gruneisen models fit) - fitting and evaluation of a
% Mie-Gruneisen p(V,T) EOS. The p(V,T) is made by adding up two
% different parts: 
%
%      p(V,T) = p_T0(V) + p_th(V,T)
%
% The p_T0(V) or "cold" isotherm EOS can be represented by many forms, like
% the Birch-Murnaghan, or Poirier-Tarantola, for instance. The
% "thermal" part is usually described in terms of a modified Debye form.
% This routine is designed to let the user combine many different forms for
% both EOS components.
%
% Required input variables:
% (p,V,T): column vectors containing the points of the p(V,T) surface.
% Vref: reference volume used in many cold EOS forms.
% Tref: reference temperature used in the thermal EOS part.
%       Notice that a subset of the (p,V,T) data (the points such that
%       T==Tref) will be used to fit the cold EOS.
%
% Required global variables:
%
% natoms: the number of atoms in the cell consistent with the volume v.
%         Required by Tange thermal pressure expressions.
% nelectrons: the number of electrons, for Holzapfel's AP2 equation of
%             state.
%
% Optional input variables (all have default values):
% {mode1 = "none"}: select the cold EOS. The default value only provides an
%    error message. Implemented modes:
%    * "bm" or "eulerian": average of Birch-Murnaghan polynomials.
%    * "bm3", "bm4", ...: Birch-Murnaghan family with a fixed order.
%    * "pt" or "natural": average of Poirier-Tarantola polynomials.
%    * "pt3", "pt4", ...: Poirier-Tarantola family with a fixed order.
%    Notice that no input parameters are required to fit the linear
%    forms (bm, pt).
% {pin1 = []}: column vector with the values (evaluation) or starting
%   values (fitting) of the cold EOS parameters.
%   In the case of fitting (fit = 1), if no values are entered, the routine will
%   try to determine some reasonable guess. If evaluating (fit = 0),
%   pin1 is mandatory.
% {mode2 = "none"}: select the thermal EOS. The default value only provides an
%    error message. Implemented modes:
%    * "tange1", "tange2", "tange3": hierarchical Debye models by Tange et
%      al [J. Geophys. Res. 114 (2009) B03208]. The parameters are:
%      * "tange1": [theta0, gamma0].
%      * "tange2": [theta0, gamma0, q0].
%      * "tange3": [theta0, gamma0, a0, b0].
%      * "jackson4": [a0, a1, a2, a3].
% {pin2 = []}: column vector with the start value for the cold EOS parameters.
%   If no values are entered, the routine will try to determine some starting
%   point.
% {fit = 1}: in default mode the analytical EOS will be fitted to the
%   (V,E) data. In the fit=0 mode, the pressure will be evaluated with the
%   provided pin parameters (required in this mode) at the volumes
%   contained in V and the temperatures in T.
% {LOG = 1}: print internal information about the fitting if (LOG>0).
%
% Output:
%   res.pout1 ... cold isotherm fitting parameters (equiv. to pin1)
%   res.pout2 ... thermal pressure fitting parameters (equiv. to pin2)
%   res.konv .... 1 if the fitting converged and 0 if not.
%   res.pfit .... column vector with the values of p predicted by the fit.
%   res.R2 ...... determination coefficient (1 for a exact fit).
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: January 2011

   %% Global variables
   global nelectrons
   global natoms
   global V0
   global T0

   %% Startup
   mode1 = tolower(mode1);
   mode2 = tolower(mode2);
   V0 = Vref; T0 = Tref;
   if (nargin < 5 | strcmp(mode1,"none") | strcmp(mode2,"none"))
      print_usage();
   elseif (length(p) != length(V) | length(p) != length(T))
      error('miefit: p, V, and T must be vectors of the same length!');
   elseif (exist("leasqr") ~= 2)
      error('miefit: leasqr, from octaveforge/optim, MUST be in your path (pth fit)!');
   endif

   %% Conversion factors
   rybohr3togpa = 14710.50498740275538944426;
   hybohr3togpa = 14710.50498740275538944426 * 2;

   if (fit == 0)
      %% Evaluate the cold isotherm for the provided pin1 and volume
      if (length(pin1) < 1)
         error('miefit: In fit==0 mode you *must* give the pin parameters!')
      endif
      if (strcmp(mode1,'eulerian') | regexp(mode1,'bm[0-9]*'))
         pfit1 = pvstrainevalp(pin1, Vref, V, 'eulerian');
      elseif (strcmp(mode1,'natural') | regexp(mode1,'pt[0-9]*'))
         pfit1 = pvstrainevalp(pin1, Vref, V, 'natural');
      else
         error('miefit: cold EOS mode unknown!')
      endif
      res.pout1 = pin1;
      pth = p - pfit1;
   else
      iref = find(T == Tref);
      if (length(iref) < 7)
         error('miefit: too few cold isotherm points in the dataset!')
      endif
      [pcold,icold] = sort(p(iref));
      vcold = V(iref)(icold);
      if (strcmp(mode1,'bm') | strcmp(mode1,'eulerian'))
         [cavg,savg] = pvavgstrainfit(pcold,vcold,Vref,10,:,'eulerian',LOG);
         res.pout1 = cavg;
         pfit1 = pvstrainevalp(cavg, Vref, V, 'eulerian');
         pth = p - pfit1;
      elseif (regexp(mode1,'bm[0-9]+'))
         norder = str2num(regexp(mode1, '[0-9]+','match', 'once'));
         if (norder < 2)
            error('miefit: BM order < 2!')
         endif
         [cf,sf] = pvstrainfit(pcold,vcold,Vref,norder,'eulerian',LOG);
         res.pout1 = cf;
         pfit1 = pvstrainevalp(cf, Vref, V, 'eulerian');
         pth = p - pfit1;
      elseif (strcmp(mode1,'pt') | strcmp(mode1,'natural'))
         [cavg,savg] = pvavgstrainfit(pcold,vcold,Vref,10,:,'natural',LOG);
         res.pout1 = cavg;
         pfit1 = pvstrainevalp(cavg, Vref, V, 'natural');
         pth = p - pfit1;
      elseif (regexp(mode1,'pt[0-9]+'))
         norder = str2num(regexp(mode1, '[0-9]+','match', 'once'));
         if (norder < 2)
            error('miefit: PT order < 2!')
         endif
         [cf,sf] = pvstrainfit(pcold,vcold,Vref,norder,'natural',LOG);
         res.pout1 = cf;
         pfit1 = pvstrainevalp(cf, Vref, V, 'natural');
         pth = p - pfit1;
      else
         error('miefit: cold EOS unknown!')
      endif
   endif

   if (strcmp(mode2,'tange1'))
      if (length(pin2) < 2)
         # pin2: [Theta0, gamma0]
         #          1       2
         pin2 = [500, 1.5];
      endif
      xvar = [V, T]; yvar = pth;
      name2 = "(Tange etal) Mie-Gruneisen-Debye-1 EOS";
      par_name2 = {'Theta0 (K)', 'gamma0'};
      par_fact2 = [1, 1];
      fun2 = @tange1_p;
   elseif (strcmp(mode2,'tange2'))
      if (length(pin2) < 3)
         # pin2: [Theta0, gamma0, q0]
         #          1       2      3
         pin2 = [500, 1.5, 1.1];
      endif
      xvar = [V, T]; yvar = pth;
      name2 = "(Tange etal) Mie-Gruneisen-Debye-2 EOS";
      par_name2 = {'Theta0 (K)', 'gamma0', 'q0'};
      par_fact2 = [1, 1, 1];
      fun2 = @tange2_p;
   elseif (strcmp(mode2,'tange3'))
      if (length(pin2) < 4)
         # pin2: [Theta0, gamma0, aa0, bb0]
         #          1       2      3    4
         pin2 = [500, 1.5, 0.8, 1.2];
      endif
      xvar = [V, T]; yvar = pth;
      name2 = "(Tange etal) Mie-Gruneisen-Debye-3 EOS";
      par_name2 = {'Theta0 (K)', 'gamma0', 'aa0', 'bb0'};
      par_fact2 = [1, 1, 1, 1];
      fun2 = @tange3_p;
   elseif (strcmp(mode2,'tange4'))
      if (length(pin2) < 3)
         # pin2: [Theta0, gamma0, bb0]
         #          1       2      3
         pin2 = [500, 1.5, 1.2];
      endif
      xvar = [V, T]; yvar = pth;
      name2 = "(Tange etal) Mie-Gruneisen-Debye-3 EOS";
      par_name2 = {'Theta0 (K)', 'gamma0', 'bb0'};
      par_fact2 = [1, 1, 1];
      fun2 = @tange4_p;
   elseif (strcmp(mode2,'jackson4'))
      if (length(pin2) < 4)
         # pin2: [a0, a1, a2, a3]
         #        1   2   3   4
         pin2 = [1.1, 1.2, 1.3, 1.4];
      endif
      xvar = [V, T]; yvar = pth;
      name2 = "Jackson-Rigden-4 EOS";
      par_name2 = {'a0', 'a1 (1/K)', 'a2 (1/K^2)', 'a3 (1/K^3)'};
      par_fact2 = [1, 1, 1, 1];
      fun2 = @jackson4_p;
   else
      error('miefit: thermal EOS mode unknown!')
   endif

   if (fit==0)
      pfit2 = fun2(xvar,pin2);
      res.pfit = pfit1 + pfit2;
      res.pout2 = pin2;
      res.konv = 0;
      res.R2 = 0;
      if (LOG > 0)
         printf('\n(miefit) Evaluating %s\n', name2);
         printf('Parameters (%d)\n', length(pin2));
         for i = 1 : length(pin2)
            printf('%-16s', par_name2{i});
            printf(' %16.6f', pin2(i)*par_fact2(i));
            printf('\n');
         endfor
      endif
      if (size(yvar) == size(pfit2))
         SSerr = sum((yvar-pfit2).^2);
         SStot = sum((yvar-mean(yvar)).^2);
         res.R2 = 1 - SSerr / SStot;
         if (LOG>0)
            printf('(Thermal) SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
                  , SSerr, SStot, res.R2, 1-res.R2);
         endif
         SSerr = sum((p-res.pfit).^2);
         SStot = sum((p-mean(p)).^2);
         res.R2 = 1 - SSerr / SStot;
         if (LOG>0)
            printf('(.Total.) SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
                  , SSerr, SStot, res.R2, 1-res.R2);
         endif
      endif
   else
      [f1, p1, kvg1, iter1, corp1, covp1, covr1, stdresid1, Z1, r21] \
      = leasqr(xvar, yvar, pin2, fun2, 1e-9, 30);
      SSerr2 = sum((yvar-f1).^2);
      SStot2 = sum((yvar-mean(yvar)).^2);
      R22 = 1 - SSerr2 / SStot2;
      res.pout2 = p1;
      res.konv = kvg1;
      res.pfit = pfit1 + f1;
      SSerr = sum((p-res.pfit).^2);
      SStot = sum((p-mean(p)).^2);
      R2 = 1 - SSerr / SStot;
      res.R2 = R2;
      if (LOG > 0)
         printf('\nNon-linear fitting: %s\n', name2);
         printf('Parameters (%d) start / converged\n', length(pin2));
         for i = 1 : length(pin2)
            printf('%-16s', par_name2{i});
            printf(' %16.6f', pin2(i)*par_fact2(i));
            printf(' %16.6f', p1(i)*par_fact2(i));
            printf('\n');
         endfor
         printf('Convergence (1=yes)? %d\n', kvg1);
         printf('Iterations: %d\n', iter1);
         printf('(Thermal) SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
               , SSerr2, SStot2, R22, 1-R22);
         printf('(.Total.) SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
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

# (Tange etal) Mie-Gruneisen-Debye-1 model
# parameters: [Theta0, gamma0]
#                1       2
# Assumed units: V(bohr^3), T(K), p(GPa)
function p_th = tange1_p(x,par)
  global natoms;
  global T0
  global V0
  V = x(:,1); T = x(:,2);
  R = 8.314472; bohr = 0.52917720859e-10; N_A = 6.02214179e23;
  theta0 = par(1);
  gamma = gamma0 = par(2);
  theta = (V0./V).^gamma * theta0; 
  n3R = (T./V) * (gamma * 3*natoms*R/bohr^3/N_A/1e9);
  p_th = (debye_3(theta./T) .- debye_3(theta/T0)) .* n3R;
endfunction

# (Tange etal) Mie-Gruneisen-Debye-2 model
# parameters: [Theta0, gamma0, q0]
#                1       2      3
# Assumed units: V(bohr^3), T(K), p(GPa)
function p_th = tange2_p(x,par)
  global natoms;
  global T0
  global V0
  V = x(:,1); T = x(:,2);
  R = 8.314472; bohr = 0.52917720859e-10; N_A = 6.02214179e23;
  theta0 = par(1); gamma0 = par(2); q0 = par(3);
  gamma = (V./V0).^q0 * gamma0;
  theta = exp(-(gamma.-gamma0)/q0) * theta0;
  n3R = (T./V) .* gamma * (3*natoms*R/bohr^3/N_A/1e9);
  p_th = (debye_3(theta./T) .- debye_3(theta/T0)) .* n3R;
endfunction


# (Tange etal) Mie-Gruneisen-Debye-3 model
# parameters: [Theta0, gamma0, aa0, bb0]
#                1       2      3    4
# Assumed units: V(bohr^3), T(K), p(GPa)
function p_th = tange3_p(x,par)
  global natoms;
  global T0
  global V0
  V = x(:,1); T = x(:,2);
  theta0 = par(1); gamma0 = par(2); aa0 = par(3); bb0 = par(4);
  R = 8.314472; bohr = 0.52917720859e-10; N_A = 6.02214179e23;
  gamma = ((((V./V0).^bb0-1) * aa0) + 1) * gamma0;
  theta = exp(-(gamma.-gamma0)/bb0) .* (V./V0).^((aa0-1)*gamma0) * theta0;
  n3R = (T./V) .* gamma * (3*natoms*R/bohr^3/N_A/1e9);
  p_th = (debye_3(theta./T) .- debye_3(theta/T0)) .* n3R;
endfunction

# (Tange etal) Mie-Gruneisen-Debye-4 model, same as tange3 with a=1 fixed.
# parameters: [Theta0, gamma0, bb0]
#                1       2     3
# Assumed units: V(bohr^3), T(K), p(GPa)
function p_th = tange4_p(x,par)
  global natoms;
  global T0
  global V0
  V = x(:,1); T = x(:,2);
  theta0 = par(1); gamma0 = par(2); aa0 = 1; bb0 = par(3);
  R = 8.314472; bohr = 0.52917720859e-10; N_A = 6.02214179e23;
  gamma = ((((V./V0).^bb0-1) * aa0) + 1) * gamma0;
  theta = exp(-(gamma.-gamma0)/bb0) .* (V./V0).^((aa0-1)*gamma0) * theta0;
  n3R = (T./V) .* gamma * (3*natoms*R/bohr^3/N_A/1e9);
  p_th = (debye_3(theta./T) .- debye_3(theta/T0)) .* n3R;
endfunction

# (Jackson-Rigden) Jackson-Ridge-4 model
# parameters: [a0, a1, a2, a3]
#              1   2   3   4
# Assumed units: V(bohr^3), T(K), p(GPa)
function p_th = jackson4_p(x,par)
  global natoms;
  global T0
  global V0
  V = x(:,1); tt = (x(:,2)-T0)/T0;
  #R = 8.314472; bohr = 0.52917720859e-10; N_A = 6.02214179e23;
  p_th = log(V./V0).*tt*par(1) .+ tt.*(tt.*(tt*par(3)+par(2))+par(1));
endfunction
