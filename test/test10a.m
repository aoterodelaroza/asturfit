#! /usr/bin/octave -q
#
# Fitting of the (p,V,T) data of MgO.
#

addpath("~/cvs/asturfit/src/");
addpath("/usr/share/octave/packages/3.2/optim-1.0.12/");

#### MgO pvt equation of state

[p,v,t] = readPVTdata('mgo_reduced.eos');

# Convert the T==0 data into 0.0001 K to avoid problems with the debye_3()
# routine.
t(abs(t)<0.0001) = 0.0001;

## Avoid the problematic (low p & high T)
global natoms
global mm
natoms = 2;
mm = 40.3044;
[p,v,t] = skimPVT(p,v,t,'debye',:,1.5);
## [p,v,t] = skimPVT(p,v,t,'square',10,1500);

## variability of the whole set of pressures.
SStot = sum((p-mean(p)).^2);

## Fit to a BM3-Tange3 model
tref = 0.0001;
vref = max(v(find(t==tref))); 
res = miefit(p, v, t, vref, tref, 'bm3', 'tange3', [], [], 1, 1);
res2 = miefit(p, v, t, vref, tref, 'bm3', 'tange3', res.pout1, res.pout2, 0, 0);
pdif = p - res2.pfit;
adif = abs(pdif);
printf("Pressure of max. deviation: %.5f\n",p(adif == max(adif)));
printf("Temperature of max. deviation: %.5f\n",t(adif == max(adif)));
printf("Volume of max. deviation: %.5f\n",v(adif == max(adif)));

# Get the parameters of the EOS
vmin = min(v(find(t==tref)));
res3 = pvstrainmin(res.pout1,vref,[vmin 2*vref],strain='eulerian');
vmin = res3.Vmin;
res3 = pvstraineval(res.pout1,vref,vmin,'eulerian');
printf("* PARAMETERS at zero pressure\n")
printf("p0   = %.15f\n",res3.p)
printf("V0   = %.15f\n",vmin)
printf("B0   = %.15f\n",res3.B)
printf("B0p  = %.15f\n",res3.B1p)
printf("B0pp = %.15f\n",res3.B2p)
printf("ThetaD = %.15f\n",res2.pout2(1))
printf("gamma0 = %.15f\n",res2.pout2(2))
printf("a = %.15f\n",res2.pout2(3))
printf("b = %.15f\n",res2.pout2(4))

# Compute deviations
SSerr = sum(pdif.^2);
R2 = 1 - SSerr / SStot;
printf('Total SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
      , SSerr, SStot, R2, 1-R2);
printf('Total Max,mean,sdev error: %.6e %.6e %.6e\n'
      , max(adif), mean(adif), std(adif))

#### Al pvt equation of state
[p,v,t] = readPVTdata('al.eos',0);

# Convert the T==0 data into 0.0001 K to avoid problems with the debye_3()
# routine.
tref = 0.0001;
t(t==0) += tref;

## Avoid the problematic (p<=10 GPa & T>=1500 K)
global natoms
global mm
natoms = 1;
mm = 26.981538;
[p,v,t] = skimPVT(p,v,t,'debye',:,1.5);

## Fit to a BM3-Tange3 model
SStot = sum((p-mean(p)).^2);
vref = max(v(find(t==tref))); 
res = miefit(p, v, t, vref, tref, 'bm3', 'tange3', [], [], 1, 1);
res2 = miefit(p, v, t, vref, tref, 'bm3', 'tange3', res.pout1, res.pout2, 0, 0);
pdif = p - res2.pfit;
adif = abs(pdif);
printf("Pressure of max. deviation: %.5f\n",p(adif == max(adif)));
printf("Temperature of max. deviation: %.5f\n",t(adif == max(adif)));
printf("Volume of max. deviation: %.5f\n",v(adif == max(adif)));

## Output parameters
vmin = min(v(find(t==tref)));
res3 = pvstrainmin(res.pout1,vref,[vmin 2*vref],strain='eulerian');
vmin = res3.Vmin;
res3 = pvstraineval(res.pout1,vref,vmin,'eulerian');
printf("* PARAMETERS at zero pressure\n")
printf("p0   = %.15f\n",res3.p)
printf("V0   = %.15f\n",vmin)
printf("B0   = %.15f\n",res3.B)
printf("B0p  = %.15f\n",res3.B1p)
printf("B0pp = %.15f\n",res3.B2p)
printf("ThetaD = %.15f\n",res2.pout2(1))
printf("gamma0 = %.15f\n",res2.pout2(2))
printf("a = %.15f\n",res2.pout2(3))
printf("b = %.15f\n",res2.pout2(4))

## Compute deviations
SSerr = sum(pdif.^2);
R2 = 1 - SSerr / SStot;
printf('Total SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
      , SSerr, SStot, R2, 1-R2);
printf('Total Max,mean,sdev error: %.6e %.6e %.6e\n'
      , max(adif), mean(adif), std(adif))

#### Diamond pvt equation of state
[p,v,t] = readPVTdata('c.eos',0);

# Convert the T==0 data into 0.0001 K to avoid problems with the debye_3()
# routine.
tref = 0.0001;
t(t==0) += tref;

## Avoid the problematic (p<=10 GPa & T>=1500 K)
global natoms
global mm
natoms = 2;
mm = 24.0214;
[p,v,t] = skimPVT(p,v,t,'debye',:,1.5);

## Fit to a BM3-Tange3 model
SStot = sum((p-mean(p)).^2);
vref = max(v(find(t==tref))); 
res = miefit(p, v, t, vref, tref, 'bm3', 'tange4', [], [], 1, 1);
res2 = miefit(p, v, t, vref, tref, 'bm3', 'tange4', res.pout1, res.pout2, 0, 0);
pdif = p - res2.pfit;
adif = abs(pdif);
printf("Pressure of max. deviation: %.5f\n",p(adif == max(adif)));
printf("Temperature of max. deviation: %.5f\n",t(adif == max(adif)));
printf("Volume of max. deviation: %.5f\n",v(adif == max(adif)));

# Get the parameters of the EOS
vmin = min(v(find(t==tref)));
res3 = pvstrainmin(res.pout1,vref,[vmin 2*vref],strain='eulerian');
vmin = res3.Vmin;
res3 = pvstraineval(res.pout1,vref,vmin,'eulerian');

## Output parameters
vmin = min(v(find(t==tref)));
res3 = pvstrainmin(res.pout1,vref,[vmin 2*vref],strain='eulerian');
vmin = res3.Vmin;
res3 = pvstraineval(res.pout1,vref,vmin,'eulerian');
printf("* PARAMETERS at zero pressure\n")
printf("p0   = %.15f\n",res3.p)
printf("V0   = %.15f\n",vmin)
printf("B0   = %.15f\n",res3.B)
printf("B0p  = %.15f\n",res3.B1p)
printf("B0pp = %.15f\n",res3.B2p)
printf("ThetaD = %.15f\n",res2.pout2(1))
printf("gamma0 = %.15f\n",res2.pout2(2))
printf("b = %.15f\n",res2.pout2(3))

## Compute deviations
SSerr = sum(pdif.^2);
R2 = 1 - SSerr / SStot;
printf('Total SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
      , SSerr, SStot, R2, 1-R2);
printf('Total Max,mean,sdev error: %.6e %.6e %.6e\n'
      , max(adif), mean(adif), std(adif))

