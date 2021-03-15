function [fmap] = MakeScatFact90keV(s0,Scat,smap,eta,Atoms);

[UniqueAtoms,ia,ScatIndex] = unique(Atoms);

NUA=length(UniqueAtoms);

for i=1:NUA
fmap.(genvarname([UniqueAtoms{i}]))=interp1(s0,Scat.(genvarname([UniqueAtoms{i}])),smap,'spline')...
.*exp(1i*interp1(s0,eta.(genvarname([UniqueAtoms{i}])),smap,'spline'));
end