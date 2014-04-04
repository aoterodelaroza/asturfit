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

function [s] = straineval (c, V0, V, strain='eulerian')
% function [s] = straineval (c, V0, V {, strain})
%
% straineval - Evaluation of properties derived from a polynomial energy
% function based on a given strain form.
%
% Required input variables:
% c: coefficients of the "energy versus f" polynomial. Usually this is
%    determined in a previous call to the strainfit() routine.
% V0: reference volume.
% V: vector containing the volumes at which the polynomial properties must
%    be calculated.
%
% Optional input variables (all have default values):
% {strain = 'eulerian'}: strain form.
%
% Output:
% s: data structure containing:
%       s.E: vector with energies at the V points.
%       s.E1v: first derivative of E versus V.
%       s.E2v: second derivative of E versus V.
%       s.E3v: third derivative of E versus V.
%       s.E4v: fourth derivative of E versus V.
%       s.p: vector with pressures at the V points.
%       s.B: vector with bulk moduli at the V points.
%       s.B1p: first derivative of the bulk modulus versus p.
%       s.B2p: second derivative of the bulk modulus versus p.
%       s.B3p: third derivative of the bulk modulus versus p.
%
% See also: strainfit.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: September 2010.
%

   if (nargin < 3 | nargin > 4)
      print_usage ();
   endif

   # Determine the strain:
   if (strcmp(strain,'eulerian'))
      f = ((V./V0).^(-2/3)-1)/2;
      f2 = f + f + 1;
      ss = -f2.^(3/2)/(3*V0);
      f1v = f2.^(5/2) * (-1/(3*V0));
      f2v = f1v .* ss * 5;
      f3v = f2v .* ss * 8;
      f4v = f3v .* ss * 11;
      f5v = f4v .* ss * 14;
      f6v = f5v .* ss * 17;
      ###f2v = f2.^4 * (5/(3*V0)^2);
      ###f3v = f2.^(11/2) * (-40/(3*V0)^3);
      ###f4v = f2.^7 * (440/(3*V0)^4);
      ###f5v = f2.^(17/2) * (-6160/(3*V0)^5);
   elseif (strcmp(strain, 'natural'))
      f = log(V./V0)/3;
      ss = -1./V;
      f1v = 1./(V*3);
      f2v = f1v .* ss;
      f3v = f2v .* ss * 2;
      f4v = f3v .* ss * 3;
      f5v = f4v .* ss * 4;
      f6v = f5v .* ss * 5;
   elseif (strcmp(strain,'lagrangian'))
      f = ((V./V0).^(2/3)-1)/2;
      ss = -1./(V*3);
      f1v = (f + f + 1).^(-1/2) / (3*V0);
      f2v = ss .* f1v;
      f3v = ss .* f2v * 4;
      f4v = ss .* f3v * 7;
      f5v = ss .* f4v * 10;
      f6v = ss .* f5v * 13;
   elseif (strcmp(strain,'infinitesimal'))
      f = -(V./V0).^(-1/3)+1;
      ss = -1./(V*3);
      f1v = (-f + 1).^4 / (3*V0);
      f2v = f1v .* ss * 4;
      f3v = f2v .* ss * 7;
      f4v = f3v .* ss * 10;
      f5v = f4v .* ss * 13;
      f6v = f5v .* ss * 16;
   elseif (strcmp(strain, 'quotient') | strcmp(strain, 'x1'))
      f = V./V0;
      f1v = ones(size(V)) * (1/V0);
      f2v = zeros(size(V));
      f3v = zeros(size(V));
      f4v = zeros(size(V));
      f5v = zeros(size(V));
      f6v = zeros(size(V));
   elseif (strcmp(strain, 'x3'))
      f = (V./V0).^(1/3);
      ##ss = -1./(f.^3 * (3*V0));
      ##f1v = 1./(f.^2 * (3*V0));
      ss = -1./(V*3);
      f1v = V.^(-2/3)/(3*V0^(1/3));
      f2v = f1v .* ss * 2;
      f3v = f2v .* ss * 5;
      f4v = f3v .* ss * 8;
      f5v = f4v .* ss * 11;
      f6v = f5v .* ss * 14;
   elseif (strcmp(strain, 'xinv3'))
      f = (V./V0).^(-1/3);
      ss = -1./(V*3);
      f1v = -V.^(-4/3) * (V0^(1/3)/3);
      f2v = f1v .* ss * 4;
      f3v = f2v .* ss * 7;
      f4v = f3v .* ss * 10;
      f5v = f4v .* ss * 13;
      f6v = f5v .* ss * 16;
   elseif (strcmp(strain, 'V'))
      f = V;
      f1v = ones(size(V));
      f2v = zeros(size(V));
      f3v = zeros(size(V));
      f4v = zeros(size(V));
      f5v = zeros(size(V));
      f6v = zeros(size(V));
   else
      error('straineval: strain form requested is unknown!');
   endif

   # Evaluate E(f) and the derivatives of E versus f:
   s.E = polyval(c, f);
   c1 = polyderiv(c);  E1f = polyval(c1, f);
   c2 = polyderiv(c1); E2f = polyval(c2, f);
   c3 = polyderiv(c2); E3f = polyval(c3, f);
   c4 = polyderiv(c3); E4f = polyval(c4, f);
   c5 = polyderiv(c4); E5f = polyval(c5, f);

   # Get the derivatives of E versus V:
   s.E1v = E1f.*f1v;
   s.E2v = E2f.*f1v.^2 .+ E1f.*f2v;
   s.E3v = E3f.*f1v.^3 .+ E2f.*f1v.*f2v*3 .+ E1f.*f3v;
   s.E4v = E4f.*f1v.^4 .+ E3f.*f1v.^2.*f2v*6 .+ \
           E2f.*(f1v.*f3v*4.+f2v.^2*3) .+ E1f.*f4v;
   s.E5v = E5f.*f1v.^5 .+ E4f.*f1v.^3.*f2v*10 .+ \
           E3f.*(f1v.^2.*f3v*10.+f1v.*f2v.^2*15) .+ \
           E2f.*(f1v.*f4v*5.+f2v.*f3v*10) .+ E1f.*f5v;

   # Get the pressure, the bulk modulus, and its derivatives versus pressure:
   s.p = -s.E1v;
   s.B = V .* s.E2v;
   s.B1p = -(V .* s.E3v .+ s.E2v) ./ s.E2v;
   s.B2p = ((s.E4v.*s.E2v .- s.E3v.^2) .* V .+ s.E3v.*s.E2v) ./ s.E2v.^3;
   s.B3p = -((s.E2v.^2.*s.E5v.-s.E2v.*s.E3v.*s.E4v*4.+s.E3v.^3*3).*V \
         .+ s.E2v.^2.*s.E4v*2 .- s.E2v.*s.E3v.^2*3) ./ s.E2v.^5;
   
endfunction
