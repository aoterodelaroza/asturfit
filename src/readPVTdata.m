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

function [p, V, T] = readPVTdata (filename, LOG=1)
% function [p, V, T] = readPVTdata (filename, LOG=1)
%
% readPVTdata - read in the pressure, volume and temperature of a crystal phase.
%
% Required input variables:
% filename: name of the data file.
%
% Optional input variables (all have default values):
% {LOG = 1}: print information about the data read in if LOG>0.
%            LOG = 0  no output.
%            LOG = 1  number of points read in, volume and energy range.
%            LOG = 2  like 1 plus a complete list of the points read in.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: December 2010

  global nelectrons;
  global natoms;

  if (LOG>0)
    printf("readPVTdata: Reading input file -> %s\n", filename);
    tic;
  endif
  rybohr3togpa = 14710.50498740275538944426;
  hybohr3togpa = 14710.50498740275538944426 * 2;
  angtobohr = 1.88972613288564;
  bohrtoans = 0.52917720859;
  ev2hy = 0.0367493254;
  hy2ev = 27.21138386;   
  ry2ev = hy2ev/2;
  hytory = hy2ry = 2;
  rytohy = ry2hy = 1/hy2ry;

  ## Datafile can contain the units for some critical variables. Currently:
  ## * columns -----> column order of data. Default: "pvt".
  ## * pressure ----> default unit: GPa
  ## * volume ------> default unit: bohr^3
  ## * temperature -> default unit: K
  ## * z -----------> number of molecules in the unit cell (default 1)
  ## * nelectrons --> number of electrons in a molecule (required by some EOS)
  ##             This data is declared as global for the EOS that need it.
  ## * natoms ------> number of atoms in a molecule (required by the Mie-like
  ##             EOS.
  ## Those keywords will appear in comments in the head to the datafile.
  ## keywords
  Keyw={"columns","pressure","volume","temperature","nelectrons","natoms","z"};

  [fid,msg] = fopen(filename,"r");
  if (fid < 0 || ferror(fid)) 
    disp(msg)
    error("Could not find -- %s",filename);
  endif

  % *CAUTION:*
  % The number of electrons is required, at least, in the AP2 EOS.
  % Initializing it to 0 every time a new datafile is read in let the AP2
  % routine detect that this data is not set in and protest accordingly.
  % Keeping the last value read in is useful if several datasets for the
  % same system need to be read and worked with simultaneously. Only one
  % of the datafiles needs the data to be set.
  % Final decision: the danger of inheritance of a wrong 'nelectrons' from
  % a previous file is too large. It is better to initialize to a nonsense
  % value.
  %%%nelectrons = sprintf('%d',nelectrons);
  nelectrons = "0";
  natoms = "0";

  ## Read in data file:
  z = "1";
  columns = "pvt";
  nl = 1;
  datV = [];
  datP = [];
  datT = [];
  while (!feof(fid))
    line = fgetl(fid);

    ## headers and keywords
    [first,count] = sscanf(line," %1c",1);
    if (first == "\0" || first == " " || count == 0) 
      continue
    endif
    if (first == "#")
      # interpret the header as a sequence of lines:
      #   #  keyword  val
      # with keyword necessarily one of the keyw cell array.
      if (nl == 1) 
	# remove # chars and trailing blanks and tabs
	while (line(1) == "#" || line(1) == " " || line(1) == "\t")
	  line = line(2:end);
	endwhile    
	[first,sec] = sscanf(line,"%s %s","C");
	if (isa(first,"char") && isa(sec,"char"))
	  first = tolower(first); sec = tolower(sec);
          idx = ismember(Keyw,first);
	  if (any(idx) != 0)
	    eval(strcat(first,"='",sec,"' ;"));
	  endif
	endif
      endif
      continue
    endif

    ## read the line and advance
    [first,count] = sscanf(line,"%f %f %f",3);
    if (columns == "pvt")
       datP(nl) = first(1);
       datV(nl) = first(2);
       datT(nl) = first(3);
    elseif (columns == "ptv")
       datP(nl) = first(1);
       datV(nl) = first(3);
       datT(nl) = first(2);
    elseif (columns == "vpt")
       datP(nl) = first(2);
       datV(nl) = first(1);
       datT(nl) = first(3);
    elseif (columns == "vtp")
       datP(nl) = first(3);
       datV(nl) = first(1);
       datT(nl) = first(2);
    elseif (columns == "tpv")
       datP(nl) = first(2);
       datV(nl) = first(3);
       datT(nl) = first(1);
    elseif (columns == "tvp")
       datP(nl) = first(3);
       datV(nl) = first(2);
       datT(nl) = first(1);
    else
       error("readPVTdata: unknown order of data columns!");
    endif
    nl++;
  endwhile
  fclose(fid);

  ## after-read consistency checks
  if (nl < 3)
    error("too few points!")
  endif
  if (exist("volume"))
    if (LOG>0)
       printf("Input volume unit: %s\n", volume)
    endif
    if (strcmp(volume,"ang3") || strcmp(volume,"angstrom3") || strcmp(volume,"a3") ||\
        strcmp(volume,"ang^3") || strcmp(volume,"angstrom^3") || strcmp(volume,"a^3"))
      datV = datV * angtobohr^3;
    elseif (strcmp(volume,"bohr3") || strcmp(volume,"bohr^3"))
      datV = datV;
    else
      error('readPVTdata: input volume unit unknown!');
    endif
  else
    if (LOG>0)
       printf("Input volume unit: %s\n", 'bohr^3 assumed')
    endif
  endif
  if (exist("pressure"))
    if (LOG>0)
       printf("Input pressure unit: %s\n", pressure)
    endif
    if (strcmp(pressure,"gpa"))
      datP = datP;
    else
      error('readPVTdata: input pressure unit unknown!');
    endif
  else
    if (LOG>0)
       printf("Input pressure unit: %s\n", 'GPa assumed')
    endif
  endif
  if (exist("temperature"))
    if (LOG>0)
       printf("Input temperature unit: %s\n", temperature)
    endif
    if (strcmp(temperature,"k"))
      datT = datT;
    else
      error('readPVTdata: input temperature unit unknown!');
    endif
  else
    if (LOG>0)
       printf("Input temperature unit: %s\n", 'K assumed')
    endif
  endif
  datV = datV / str2num(z);
  datP = datP';
  datV = datV';
  datT = datT';
  nelectrons = str2num(nelectrons);
  natoms = str2num(natoms);

  ## ranges info
  if (LOG>0)
     printf("* Range of the PVT properties\n")
     printf("  p (GPa) : [%.5f--%.5f] dP = %.5f\n" \
           , min(datP), max(datP), max(datP)-min(datP));
     printf("  V (bohr^3) : [%.10f--%.10f] dV = %.10f\n" \
           , min(datV), max(datV), max(datV)-min(datV));
     printf("  T (K) : [%.4f--%.4f] dT = %.4f\n" \
           , min(datT), max(datT), max(datT)-min(datT));
     printf("  Number of data points: %d\n", length(datV));
     if (str2num(z) > 1)
        printf("  Number of molecules in a cell (Z): %d\n", str2num(z));
     endif
     if (nelectrons > 0)
        printf("  Number of electrons: %d\n", nelectrons);
     endif
     if (natoms > 0)
        printf("  Number of atoms: %d\n", natoms);
     endif
     if (LOG > 1)
        printf('--pt-- ----p-(GPa)----- ---V-(bohr3)---- -----T-(K)------\n');
        for i = 1 : length(datV)
           printf('%6d %16.9e %16.9e %16.9e\n', i, datP(i), datV(i), datT(i));
        endfor
     endif
     printf("\n");
     toc;
  endif

  p = datP; V = datV; T = datT;
  
endfunction
