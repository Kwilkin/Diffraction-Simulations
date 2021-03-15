function [s0,Scat,ScatIndex] = LoadScattering90keV(Atoms);

[UniqueAtoms,ia,ScatIndex] = unique(Atoms);

NUA=length(UniqueAtoms);

for i=1:NUA
    load(['C:\Users\Kylej\Documents\MATLAB\Scattering\' UniqueAtoms{i} ' 90keV.mat'],'f','s');
    Scat.(genvarname([UniqueAtoms{i}]))=f;
end
s0=s;