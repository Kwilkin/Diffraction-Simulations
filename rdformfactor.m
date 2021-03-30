function [s,sctamp,sctphs]=rdformfactor(filename,atom,eVenergy)
% This function read form factors, including scattering amplitude and
% phase, from a xls file; 
% eVenergy:electron energy with unit of KeV; available choices are 10, 40, 60, 90;
% sctamp,amplitude of form factor; sctphs, phase of form factor;
% For example,filename='D:\PHD\Research\programs\YW Xiong\Random Diffraction Simulation\form factors\form factors.xlsx';
% atom='H';[s,sctamp,sctphs]=rdformfactor(filename,atom,90);
% Written by Yanwei

sheet=atom;
A = xlsread(filename,sheet,'A3:I63');
s=A(:,1);

switch eVenergy
    case 10
        sctamp=A(:,2);
        sctphs=A(:,3);
    case 40
        sctamp=A(:,4);
        sctphs=A(:,5);
    case 60
        sctamp=A(:,6);
        sctphs=A(:,7);
    case 90
        sctamp=A(:,8);
        sctphs=A(:,9);
end

