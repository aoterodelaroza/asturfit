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

function res = nlf(V, E, mode="none", pin=[], fit=1, LOG=1)
% function res = nlf(V, E {, mode="none", pin=[], fit=1, LOG=1})
%
% nlf (non-linear fitting) - fitting to a traditional E(V) EOS.
%
% Required input variables:
% (V,E): column vectors containing the points of the E(V) curve.
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
   if (nargin < 2 || strcmp(mode,"none"))
      print_usage();
   elseif (length(V) != length(E))
      error('nlf: V and E must be vectors of the same length!')
   elseif (exist("leasqr") ~= 2)
      error('nlf: leasqr, from octaveforge/optim, MUST be in your path!');
   endif

   rybohr3togpa = 14710.50498740275538944426;
   hybohr3togpa = 14710.50498740275538944426 * 2;

   if (length(pin) < 1)
      [Emin,Imin] = min(E); Vmin = V(Imin);
      [cavg,savg] = avgstrainfit(V,E,Vmin,:,:,:,:);
   endif
   if (strcmp(mode,'bm3'))
      if (length(pin) < 4)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            pin = savg.eqmean([2,1,3,4]);
            ###pin = [savg.E0; savg.V0; savg.B0(1); savg.B0(2)];
         endif
      endif
      name = "BM3: 3rd order Birch-Murnaghan EOS";
      par_name = {'E0 (Hy)', 'V0 (bohr^3)', 'B0 (GPa)', 'B1p'};
      par_fact = [1, 1, hybohr3togpa, 1];
      fun = @bm3_e;
   elseif (strcmp(mode,'bm4'))
      if (length(pin) < 5)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            pin = savg.eqmean([2,1,3,4,5]);
            ###pin = [savg.E0; savg.V0; savg.B0(1); savg.B0(2); savg.B0(3)];
         endif
      endif
      name = "BM4: 4th order Birch-Murnaghan EOS";
      par_name = {'E0 (Hy)', 'V0 (bohr^3)', 'B0 (GPa)', 'B1p', 'B2p (1/GPa)'};
      par_fact = [1, 1, hybohr3togpa, 1, 1/hybohr3togpa];
      fun = @bm4_e;
   elseif (strcmp(mode,'bm5'))
      if (length(pin) < 6)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            pin = savg.eqmean([2,1,3,4,5,6]);
            ###pin = [savg.E0; savg.V0; savg.B0(1); savg.B0(2); savg.B0(3); savg.B0(4)];
         endif
      endif
      name = "BM5: 5th order Birch-Murnaghan EOS";
      par_name = {'E0 (Hy)', 'V0 (bohr^3)', 'B0 (GPa)', 'B1p'...
               , 'B2p (1/GPa)', 'B3p (1/GPa^2)'};
      par_fact = [1, 1, hybohr3togpa, 1, 1/hybohr3togpa, 1/hybohr3togpa^2];
      fun = @bm5_e;
   elseif (strcmp(mode,'pt3'))
      if (length(pin) < 4)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            pin = savg.eqmean([2,1,3,4]);
            ###pin = [savg.E0; savg.V0; savg.B0(1); savg.B0(2)];
         endif
      endif
      name = "PT3: 3rd order Poirier-Tarantola EOS";
      par_name = {'E0 (Hy)', 'V0 (bohr^3)', 'B0 (GPa)', 'B1p'};
      par_fact = [1, 1, hybohr3togpa, 1];
      fun = @pt3_e;
   elseif (strcmp(mode,'pt4'))
      if (length(pin) < 5)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            pin = savg.eqmean([2,1,3,4,5]);
            ###pin = [savg.E0; savg.V0; savg.B0(1); savg.B0(2); savg.B0(3)];
         endif
      endif
      name = "PT4: 4th order Poirier-Tarantola EOS";
      par_name = {'E0 (Hy)', 'V0 (bohr^3)', 'B0 (GPa)', 'B1p', 'B2p (1/GPa)'};
      par_fact = [1, 1, hybohr3togpa, 1, 1/hybohr3togpa];
      fun = @pt4_e;
   elseif (strcmp(mode,'pt5'))
      if (length(pin) < 6)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            pin = savg.eqmean([2,1,3,4,5,6]);
            ###pin = [savg.E0; savg.V0; savg.B0(1); savg.B0(2); savg.B0(3); savg.B0(4)];
         endif
      endif
      name = "PT5: 5th order Poirier-Tarantola EOS";
      par_name = {'E0 (Hy)', 'V0 (bohr^3)', 'B0 (GPa)', 'B1p'...
               , 'B2p (1/GPa)', 'B3p (1/GPa^2)'};
      par_fact = [1, 1, hybohr3togpa, 1, 1/hybohr3togpa, 1/hybohr3togpa^2];
      fun = @pt5_e;
   elseif (strcmp(mode,'murn'))
      if (length(pin) < 4)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            pin = savg.eqmean([2,1,3,4]);
            ###pin = [savg.E0; savg.V0; savg.B0(1); savg.B0(2)];
         endif
      endif
      name = "Murn: Murnaghan EOS";
      par_name = {'E0 (Hy)', 'V0 (bohr^3)', 'B0 (GPa)', 'B1p'};
      par_fact = [1, 1, hybohr3togpa, 1];
      fun = @murn_e;
   elseif (strcmp(mode,'vinet'))
      if (length(pin) < 4)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            pin = savg.eqmean([2,1,3,4]);
            ###pin = [savg.E0; savg.V0; savg.B0(1); savg.B0(2)];
         endif
      endif
      name = "Vinet: Vinet EOS";
      par_name = {'E0 (Hy)', 'V0 (bohr^3)', 'B0 (GPa)', 'B1p'};
      par_fact = [1, 1, hybohr3togpa, 1];
      fun = @vinet_e;
   elseif (strcmp(mode,'antons'))
      if (length(pin) < 4)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            Esup = strainevalE(cavg,Vmin,max(V),:);
            pin = [Esup; savg.eqmean(1); savg.eqmean(3); savg.eqmean(4)];
            ###pin = [Esup; savg.V0; savg.B0(1); savg.B0(2)];
         endif
      endif
      name = "antons: Anton-Schmidt EOS";
      par_name = {'Einf (Hy)', 'V0 (bohr^3)', 'B0 (GPa)', 'B1p'};
      par_fact = [1, 1, hybohr3togpa, 1];
      fun = @antons_e;
   elseif (strcmp(mode,'ap2'))
      if (exist("gamma_inc") == 0)
         error('nlf: gamma_inc, from octaveforge/gsl, MUST be in your path!');
      elseif (!exist("nelectrons","var") || nelectrons < 1)
         error('nlf: AP2 EOS *requires* a non null global nelectrons!');
      elseif (length(pin) < 4)
         if (fit == 0)
            error('nlf: all parameters must be given if fit=0!');
         else
            pin = savg.eqmean([2,1,3,4]);
            ###pin = [savg.E0; savg.V0; savg.B0(1); savg.B0(2)];
         endif
      endif
      name = "ap2: Holzapfel AP2 EOS";
      par_name = {'E0 (Hy)', 'V0 (bohr^3)', 'B0 (GPa)', 'B1p'};
      par_fact = [1, 1, hybohr3togpa, 1];
      fun = @EHolzapfel_AP2;
   else
      error('nlf: EOS mode unknown!')
   endif

   if (fit==0)
      res.Efit = fun(V,pin);
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
      if (size(E) == size(V))
         SSerr = sum((E-res.Efit).^2);
         SStot = sum((E-mean(E)).^2);
         res.R2 = 1 - SSerr / SStot;
         if (LOG>0)
            printf('SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'...
                  , SSerr, SStot, res.R2, 1-res.R2);
         endif
      endif
   else
      [f1, p1, kvg1, iter1, corp1, covp1, covr1, stdresid1, Z1, r21] ...
      = leasqr(V, E, pin, fun);
      SSerr = sum((E-f1).^2);
      SStot = sum((E-mean(E)).^2);
      R2 = 1 - SSerr / SStot;
      res.pout = p1;
      res.konv = kvg1;
      res.Efit = f1;
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
         printf('SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'...
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

# Birch-Murnaghan, third order. E(V) expression.
# From ASE: (Hebbache2004) Ab initio study of high-pressure behavior of a low
#     compressibility metal and a hard material:u2003Osmium and diamond,
#     Phys. Rev. B 70, 224107 (2004)
# p(1) F_0
# p(2) V_0
# p(3) B_0
# p(4) B_0'
function y = bm3_e(x,p)
  xred = p(2) ./ x;
  xred23 = xred.^(2/3)-1;
  y = p(1) + 9*p(2)*p(3)/16 * ( xred23.^3 * p(4) +...
                                xred23.^2 .* (6-4*xred.^(2/3)));
endfunction

# Birch-Murnaghan, fourth order. E(V) expression.
# p(1) F_0
# p(2) V_0
# p(3) B_0
# p(4) B_0'
# p(5) B_0''
function y = bm4_e(x,p)
  v0=p(2);
  e0=p(1);
  b0=p(3);
  b0p=p(4);
  b0pp=p(5);
  t1=(v0./x).^(1.d0/3.d0);
  t2=t1.^2;
  t3=t2-1.d0;
  t4=t3.^2/4.d0;
  t5=b0p.^2;
  y=3.d0/8.d0*b0*v0*t4.*(9.d0*t4*b0*b0pp+9.d0*t4.*t5-63.d0*t4*b0p+143.d0*t4+...
                         6.d0*b0p*t3-24.d0*t2+36.d0)+e0;
endfunction

# Birch-Murnaghan, 5th order.
# p(1) F_0
# p(2) V_0
# p(3) B_0
# p(4) B_0'
# p(5) B_0''
# p(6) B_0'''
function y0 = bm5_e(x,p)
   E0 = p(1);
   V0 = p(2); f = ((V0./x).^(2/3) - 1)/2; ff1 = 2*f+1;
   B0 = p(3); c2 = 9*B0*V0/2;
   Bp0 = p(4); c3 = c2*(Bp0-4);
   Bpp0 = p(5); c4 = (3/8)*B0*V0*(9*(Bpp0*B0+Bp0*Bp0)-63*Bp0+143);
   Bppp0 = p(6);
   c5 = (432*c2*c3*c4+576*c2^2*c4-243*c3^3-648*c2*c3^2-1350*c2^2*c3...
      -2520*c2^3)/(180*c2*c2) + Bppp0*c2^3/(45*V0^2);
   y0 = f.*f.*(f.*(f.*(f.*c5+c4)+c3)+c2)+E0;
endfunction

# Poirier-Tarantola natural strain EOS, third order
# Poirier J-P and Tarantola A, Phys. Earth Planet Int. 109, p1 (1998)
# p(1) = F_0
# p(2) = V_0
# p(3) = B_0
# p(4) = B_0'
function y = pt3_e(x,p)
  xr = -log(x/p(2));
  y = p(1) + p(2)*p(3)*xr.^2/6 .* (3+xr*(p(4)-2));
endfunction

# Poirier-Tarantola natural strain EOS, fourth order
# Poirier J-P and Tarantola A, Phys. Earth Planet Int. 109, p1 (1998)
# p(1) = F_0
# p(2) = V_0
# p(3) = B_0
# p(4) = B_0'
# p(5) = B_0''
function y = pt4_e(x,p)
  v0=p(2);
  e0=p(1);
  b0=p(3);
  b0p=p(4);
  b0pp=p(5);
  t1=b0*v0;
  t2=log(v0./x);
  t3=t2.^2;
  t4=t3.^2;
  t5=b0^2;
  t6=b0p^2;
  t7=t3.*t2;
  y=t1.*t4/8.d0+t5*v0.*t4*b0pp/24.d0-t1.*t4*b0p/8.d0+t1.*t4.*t6/24.d0+...
      t1.*t7*b0p/6.d0-t1.*t7/3.d0+t1.*t3/2.d0+e0;
endfunction

# Poirier-Tarantola natural strain EOS, fifth order
# Poirier J-P and Tarantola A, Phys. Earth Planet Int. 109, p1 (1998)
# p(1) = F_0
# p(2) = V_0
# p(3) = B_0
# p(4) = B_0'
# p(5) = B_0''
# p(6) = B_0'''
function y = pt5_e(x,p)
  v0=p(2); f=log(x/v0)/3;
  e0=p(1);
  b0=p(3); c=4.5*b0*v0;
  b0p=p(4); d=c*(2-b0p);
  b0pp=p(5); ee=(9*(d*d-c*d+c*c)*v0+2*c**3*b0pp)/(12*c*v0);
  b0ppp=p(6); g=((432*(d-c)*c*ee-243*d**3+486*(d-c)*c*d+324*c**3)*v0*v0-4*c**5*b0ppp)/(180*(c*v0)**2);
  y = f.**2 .* (f.*(f.*(f*g+ee)+d)+c) + e0;
endfunction

# Murnaghan EOS
# F.D. Murnaghan, Proc. Natl. Acad. of Sci.: USA 30 (1944) 244-247
# * Murnaghan F D, Am. J. Math. 49, p235 (1937)
# * (Fu1983) C. -L. Fu and K. -M. Ho, First-principles calculation of
#   the equilibrium ground-state properties of transition metals:
#   Applications to Nb and Mo, Phys. Rev. B 28, 5480 (1983)
# p(1) F_0
# p(2) V_0
# p(3) B_0
# p(4) B_0'
function y = murn_e(x,p)
   xred = p(2) ./ x;
   y = p(1) + p(3)*x/p(4) .* (xred.^p(4) / (p(4)-1) + 1) ...
     - p(3)*p(2)/(p(4)-1);
endfunction

# Vinet EOS
# Vinet P et al., J. Phys.: Condens. Matter 1, p1941 (1989)
# p(1) = F_0
# p(2) = V_0
# p(3) = B_0
# p(4) = B_0'
function y = vinet_e(x,p)
  xr = (x ./ p(2)).^(1/3);
  y = p(1)-2*p(3)*p(2)*exp(-3/2*(p(4)-1)*(xr-1)).*(3*xr*p(4)-3*xr+5-3*p(4))/(p(4)-1)^2 +...
      4*p(3)*p(2)/(p(4)-1)^2;
endfunction

# Anton-Schmidt
#  B. Mayer, H. Anton, E. Botta, M. Methfessel,
#  J. Sticht, J. Harris and P. C. Schmidt, Ab-initio calculation of the
#  elastic constants and thermal expansion coefficients of Laves phases,
#  Intermetallics, 11, 23 (2003)
# p(1) = F_0
# p(2) = V_0
# p(3) = B_0
# p(4) = n
function y = antons_e(x,p)
  xr = (x / p(2));
  y = p(1) + p(3)*p(2)/(p(4)+1)*xr.^(p(4)+1).*(log(xr)-1/(p(4)+1));
endfunction

# Holzapfel AP2 EOS. Parameters: p = [E0, V0, B0, B0p]
# Important: use atomic units (hartree and bohr).
function E = EHolzapfel_AP2(V,p)
   global nelectrons
   eta = (V./p(2)).**(1/3);
   pFG = (3*pi^2)^(2/3)/5 * (nelectrons/p(2))^(5/3);
   c0 = -log(3*p(3)/pFG);
   c2 = 1.5 * (p(4)-3) - c0;
   z = eta*c0;
   ec0 = exp(c0);
   g2 = (gamma_inc(-2,z)-gamma_inc(-2,c0)) * (c0**2 * ec0);
   g1 = (gamma_inc(-1,z)-gamma_inc(-1,c0)) * (c0 * (c2-1) * ec0);
   g0 = (gamma_inc( 0,z)-gamma_inc( 0,c0)) * (-2 * c2 * ec0);
   gg = (exp(-z+c0)-1)*(c2/c0);
   E = (g2+g1+g0+gg) * (9*p(2)*p(3)) + p(1);
endfunction
