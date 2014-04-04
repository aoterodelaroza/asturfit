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

% This is a test of the correct behavior of the readEVdata() routine.

global nelectrons

addpath("../src/");

# test list
tests = {"test01a.dat", # original source
	 "test01b.dat", # spaces and tabs in header and data
         "test01c.dat", # strange keywords in header
         "test01d.dat", # lines inside the body of the data
         "test01e.dat", # spurious fields in header and data
         "test01f.dat", # original source with DOS line terminators
	 "test01g.dat", # original source with MAC line terminators
         "test01h.dat",
	 };

# reference
load "test01_ref.dat";
angtobohr = 1.88972613288564;
vref = test01_ref(:,1) * angtobohr^3 / 4;
eref = test01_ref(:,2) / 2 / 4;
nelectronsref = 2;
clear test01_ref;

# do the tests
for i = 1:length(tests)
  [v,e] = readEVdata(tests{i},0);
  if (isequal(v,vref) && isequal(e,eref) && nelectrons == nelectronsref) 
    str = "pass";
  else
    str = "fail";
  endif
  printf("Test %d (%s): %s \n",i,tests{i},str);
endfor
