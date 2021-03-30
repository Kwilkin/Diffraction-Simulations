function [s0,Scat,ScatIndex,eta] = LoadScattering90keV(Atoms);

energyKeV=90;
filename='C:\Users\Kyle\Documents\MATLAB\Scattering\form factors.xlsx';

[UniqueAtoms,ia,ScatIndex] = unique(Atoms);

NUA=length(UniqueAtoms);

for i=1:NUA
%     load(['C:\Users\Kyle\Documents\MATLAB\Scattering\' UniqueAtoms{i} ' 90keV.mat'],'f','s');
%     Scat.(genvarname([UniqueAtoms{i}]))=f;
    [s,Scat.(genvarname([UniqueAtoms{i}])),eta.(genvarname([UniqueAtoms{i}]))]=rdformfactor(filename,UniqueAtoms{i},energyKeV);
end
s0=s;