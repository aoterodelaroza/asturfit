#! /usr/bin/octave -q
#
# Fitting of the (p,V,T) data of MgO.
#

[p,v,t] = readPVTdata('mgo.eos');

# Convert the T==0 data into 0.0001 K to avoid problems with the debye_3()
# routine.
t(t==0) += 0.0001;

# Avoid the problematic (p<=10 GPa & T>=1500 K)
kgood = find(p>10 | t<1500);
pg = p(kgood); vg = v(kgood); tg = t(kgood); ng = length(vg);
printf("Removed %d problematic points of %d\n", length(v)-ng, ng);

# The dataset is too large for the fitting. Select a random subset.
ks = find(rand(size(vg))>0.90);
ps = pg(ks); vs = vg(ks); ts = tg(ks);
printf("Fit to a random set of %d out of %d\n", length(ks), ng);

# variability of the whole set of pressures.
SStot = sum((p-mean(p)).^2);
SStotg = sum((pg-mean(pg)).^2);

# Fit to a <BM>-Tange3 model.
vref = max(v(find(t==300))); tref = 300;
res = miefit(ps, vs, ts, vref, tref, 'bm', 'tange3', [], [], 1, 1);
res2 = miefit(p, v, t, vref, tref, 'bm', 'tange3', res.pout1, res.pout2, 0, 0);
pdif = p - res2.pfit;
adif = abs(pdif);
SSerr = sum(pdif.^2);
R2 = 1 - SSerr / SStot;
printf('Total SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
      , SSerr, SStot, R2, 1-R2);
printf('Total Max,mean,sdev error: %.6e %.6e %.6e\n'
      , max(adif), mean(adif), std(adif))
res2 = miefit(pg, vg, tg, vref, tref, 'bm', 'tange3', res.pout1, res.pout2, 0, 0);
pdif = pg - res2.pfit;
adif = abs(pdif);
SSerr = sum(pdif.^2);
R2 = 1 - SSerr / SStotg;
printf('Safe region SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
      , SSerr, SStot, R2, 1-R2);
printf('Safe region Max,mean,sdev error: %.6e %.6e %.6e\n'
      , max(adif), mean(adif), std(adif))

# Try a BM4-Tange3 model.
vref = max(v(find(t==300))); tref = 300;
res = miefit(ps, vs, ts, vref, tref, 'bm4', 'tange3', [], [], 1, 1);
res2 = miefit(p, v, t, vref, tref, 'bm4', 'tange3', res.pout1, res.pout2, 0, 0);
pdif = p - res2.pfit;
adif = abs(pdif);
SSerr = sum(pdif.^2);
R2 = 1 - SSerr / SStot;
printf('Total SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
      , SSerr, SStot, R2, 1-R2);
printf('Total Max,mean,sdev error: %.6e %.6e %.6e\n'
      , max(adif), mean(adif), std(adif))
res2 = miefit(pg, vg, tg, vref, tref, 'bm4', 'tange3', res.pout1, res.pout2, 0, 0);
pdif = pg - res2.pfit;
adif = abs(pdif);
SSerr = sum(pdif.^2);
R2 = 1 - SSerr / SStotg;
printf('Safe region SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
      , SSerr, SStot, R2, 1-R2);
printf('Safe region Max,mean,sdev error: %.6e %.6e %.6e\n'
      , max(adif), mean(adif), std(adif))

# Try a BM7-Tange3 model.
vref = max(v(find(t==300))); tref = 300;
res = miefit(ps, vs, ts, vref, tref, 'bm7', 'tange3', [], [], 1, 1);
res2 = miefit(p, v, t, vref, tref, 'bm7', 'tange3', res.pout1, res.pout2, 0, 0);
pdif = p - res2.pfit;
adif = abs(pdif);
SSerr = sum(pdif.^2);
R2 = 1 - SSerr / SStot;
printf('Total SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
      , SSerr, SStot, R2, 1-R2);
printf('Total Max,mean,sdev error: %.6e %.6e %.6e\n'
      , max(adif), mean(adif), std(adif))
res2 = miefit(pg, vg, tg, vref, tref, 'bm7', 'tange3', res.pout1, res.pout2, 0, 0);
pdif = pg - res2.pfit;
adif = abs(pdif);
SSerr = sum(pdif.^2);
R2 = 1 - SSerr / SStotg;
printf('Safe region SSerr, SStot, R2, 1-R2: %.6e %.6e %.12f %.2e\n'\
      , SSerr, SStot, R2, 1-R2);
printf('Safe region Max,mean,sdev error: %.6e %.6e %.6e\n'
      , max(adif), mean(adif), std(adif))
