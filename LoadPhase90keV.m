function [eta,PhaseIndex] = LoadPhase90keV(Atoms);

[UniqueAtoms,~,PhaseIndex] = unique(Atoms);

NUA=length(UniqueAtoms);

for i=1:NUA
    load(['C:\Users\Kylej\Documents\MATLAB\Scattering\' UniqueAtoms{i} ' 90keV.mat'],'n');
    eta.(genvarname([UniqueAtoms{i}]))=n;
    
end