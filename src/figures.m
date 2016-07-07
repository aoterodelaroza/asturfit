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

function ur = figures(x)
%
% figures - determines the actual precision of the values contained in the
% vector x.
%
% Required input variables:
% x: data vector.
%
% Output:
% ur: resolution of x.
%
% Examples:
%    figures(1.1234000000)  -> 1e-4
%    figures(1.1234000001)  -> 1e-10
%    figures(-1.1234000001) -> 1e-10
%    figures(-11234000.001) -> 1e-3
%    figures(-11234000.001) -> 1e-3
%    figures(1.23e7)        -> 1e5
%    figures(1.23e-7)       -> 1e-9
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: November 2010.
%
   if (nargin != 1 )
      print_usage ();
   endif

   # number of decimals associated to the data type of x:
   ###p = -fix(log10(eps(x)));
   # number of decimals associated to the machine eps():
   peps = -fix(log10(eps()));

   # if the number of data is small analyze them all, otherwise check
   # a random percentage:
   n = length(x);
   if (n <= 10)
      rk = 1:n;
   else
      rk = find(rand(size(x)) >= 0.9);
   endif

   # print with p decimals and take out the number of terminal zeros:
   nz = 1000;
   for i = 1 : length(rk)
      ###p = -fix(log10(eps(x(rk(i)))));
      p = peps - fix(log10(x(rk(i))));
      s = sprintf('%.*f', p, x(rk(i)));
      fz = regexp(s,'[0\.]*$');
      nz = min(nz, length(s(fz:end))-(index(s(fz:end),'.')>0));
   endfor

   ur = 10^(-(p-nz));
endfunction
