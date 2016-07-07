% Copyright (C) 2010 Victor Lua~na and Alberto Otero-de-la-Roza
%
% This octave routine is free software: you can redistribute it and/or
% modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version. See <http://www.gnu.org/licenses/>.
%
% The routine distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.

function s = errformatf(x, dx, nd=6)
% function s = errformatf(x, dx, nd=6)
%
% errformatf - compact x(dx) printing of a quantity and an estimated error.
% This is best explained with an example. If x=123.456789 and dx=0.0567
% the routine will produce 123.457(57), that should be interpreted as
% 123.457(+/-)0.057.
%
% Required input variables:
% x: quantity to print.
% dx: estimated error.
%
% Optional input variables (all have default values):
% {nd=6}: maximum number of decimals to use.
%
% Output:
% s: string containing the printed result.
%
% Examples:
%   errformatf(-123.4567890123456,0.99734)         -> -123.5(10)
%   errformatf(-123.4567890123456,0.099734)        -> -123.46(10)
%   errformatf(-123.4567890123456,0.09734)         -> -123.457(97)
%   errformatf(-123.4567890123456,0.00009734)      -> -123.456789(97)
%   errformatf(-123.4567890123456,0.000009734)     -> -123.456789(10)
%   errformatf(-123.4567890123456,0.0000009734)    -> -123.456789(1)
%   errformatf(-123.4567890123456,0.00000009734)   -> -123.456789(0)
%   errformatf(-123.4567890123456,0.00000009734,8) -> -123.45678901(10)
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: November 2010
%
   if (dx == 0)
      if (x == 0)
         s = '0(0)';
         return
      endif
      m = log10(figures(x));
      if (m >=0)
         s = sprintf('%d(0)', x);
      else
         s = sprintf('%.*f(0)', -m, x);
      endif
      return
   endif
   p = -fix(log10(abs(dx)));
   f=10^(p+2);
   tdx = round(dx*f)/f;
   f = 10^nd;
   sdx = sprintf('%d', round(tdx*f));
   ldx = length(sdx);
   if (ldx==1)
      s = sprintf('%.*f(%d)', nd, x, round(tdx*f));
   elseif (ldx<=nd+1)
      m = nd - ldx + 2;
      s = sprintf('%.*f(%d)', m, x, round(tdx*10^m));
   else
      s = sprintf('%d(%d)', x, round(tdx));
   endif
endfunction
