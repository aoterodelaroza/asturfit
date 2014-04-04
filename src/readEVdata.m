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

function [V,E] = readEVdata (filename, LOG=1)
% function [V,E] = readEVdata (filename {, LOG=1})
%
% readEVdata - read in the volume and energy of a crystal phase.
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
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: September 2010

  global nelectrons

  if (LOG)
    printf("readEVdata: Reading input file -> %s\n", filename);
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
  ## * energy ------> default unit: hartree
  ## * volume ------> default unit: bohr^3
  ## * z -----------> number of molecules in the unit cell (default 1)
  ## * nelectrons --> number of electrons in a molecule (required by some EOS)
  ##             This data is declared as global for the EOS that need it.
  ## * x or length scale (it can be the volume too)
  ## Those keywords will appear in comments in the head to the datafile.
  ## keywords
  Keyw={"energy","volume","nelectrons","x","z"};

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

  ## input eps and ene from the elastic summary
  z = "1";
  nl = 1;
  eps = [];
  ene = [];
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
	  [idx, err] = cellidx(Keyw,first);
	  if (idx != 0)
	    eval(strcat(first,"='",sec,"' ;"));
	  endif
	endif
      endif
      continue
    endif

    ## read the line and advance
    [first,count] = sscanf(line,"%f %f",2);
    eps(nl) = first(1);
    ene(nl) = first(2);
    nl++;
  endwhile
  fclose(fid);

  ## after-read consistency checks
  if (nl < 3)
    error("too few points!")
  endif
  if (exist("x"))
    if (LOG>0)
       printf("Input length unit: %s\n", x)
    endif
    if (strcmp(x,"bohr") || strcmp(x,"au") || strcmp(x,"a.u."))
      eps = eps.^3;
    elseif (strcmp(x,"ang") || strcmp(x,"angstrom") || strcmp(x,"a"))
      eps = (eps * angtobohr).^3;
    elseif (strcmp(x,"ang3") || strcmp(x,"angstrom3") || strcmp(x,"a3") ||\
	    strcmp(x,"ang^3") || strcmp(x,"angstrom^3") || strcmp(x,"a^3"))
      eps = eps * angtobohr^3;
    else
      error('readEVdata: input length unit unknown!');
    endif
  endif
  if (exist("volume"))
    if (LOG>0)
       printf("Input volume unit: %s\n", volume)
    endif
    if (strcmp(volume,"ang3") || strcmp(volume,"angstrom3") || strcmp(volume,"a3") ||\
        strcmp(volume,"ang^3") || strcmp(volume,"angstrom^3") || strcmp(volume,"a^3"))
      eps = eps * angtobohr^3;
    elseif (strcmp(volume,"bohr3") || strcmp(volume,"bohr^3"))
      eps = eps;
    else
      error('readEVdata: input volume unit unknown!');
    endif
  else
    if (LOG>0)
       printf("Input volume unit: %s\n", 'bohr^3 assumed')
    endif
  endif
  if (exist("energy"))
    if (LOG>0)
       printf("Input energy unit: %s\n", energy)
    endif
    if (strcmp(energy,"hy") || strcmp(energy,"hartree") || strcmp(energy,"a.u."))
      ene = ene;
    elseif (strcmp(energy,"ev") || strcmp(energy,"eV"))
      ene = ene / hy2ev;
    elseif (strcmp(energy,"ry") || strcmp(energy,"rydberg"))
      ene = ene / hy2ry;
    else
      error('readEVdata: input energy unit unknown!');
    endif
  else
    if (LOG>0)
       printf("Input energy unit: %s\n", 'hartree assumed')
    endif
  endif
  ene = ene / str2num(z);
  eps = eps / str2num(z);
  eps = eps';
  ene = ene';
  nelectrons = str2num(nelectrons);

  ## ranges info
  if (LOG>0)
     xini = min(eps);
     xend = max(eps);
     xrange = xend - xini;
     erange = max(ene) - min(ene);
     printf("* Energy and volume ranges \n")
     printf("  V-range (bohr^3) : [%.10f , %.10f] dV = %.10f\n",xini,xend,xrange);
     printf("  E-range (hartree) : [%.10f , %.10f] dE = %.10f\n",min(ene),max(ene),erange);
     printf("  Number of data points: %d\n", length(eps));
     if (str2num(z) > 1)
        printf("  Number of molecules in a cell (Z): %d\n", str2num(z));
     endif
     if (nelectrons > 0)
        printf("  Number of electrons: %d\n", nelectrons);
     endif
     if (LOG > 1)
        printf('--pt-- ---V-(bohr3)---- -----E-(Hy)-----\n');
        for i = 1 : length(eps)
           printf('%6d %16.9e %16.9e\n', i, eps(i), ene(i));
        endfor
     endif
     printf("\n");
  endif

  V = eps; E = ene;
  
endfunction
