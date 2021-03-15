% This program is for CF3I diffraction pattern calculation
% Calculates the diffraction pattern for a known angular distribution
% Inputs
% N is number of pixels
% AGD is population denisty (angular distribution x sin theta; should be same length as Ntheta)
% bb is the scattering factors
% Outputs
% I_mol is the molecular scattering
% I is the total scattering
% I_at is the atomic scattering
% smap gives the s-values at each pixel

function [I_mol,I,I_at,smap]=Diff_CF3I_Cylind_HCC_90keV_Theta_Quarter(N,theta)
if nargin < 2 || isempty(N)
    N=64;
end
Cs=load('/lustre/work/centurion/kylewilkin/Diffraction/C 90keV.mat');
Fs=load('/lustre/work/centurion/kylewilkin/Diffraction/F 90keV.mat');
Is=load('/lustre/work/centurion/kylewilkin/Diffraction/I 90keV.mat');

% Cs=load('G:\Users\kwilkin2\Documents\MATLAB\Scattering\Scattering 2\C 90keV.mat');
% Fs=load('G:\Users\kwilkin2\Documents\MATLAB\Scattering\Scattering 2\F 90keV.mat');
% Is=load('G:\Users\kwilkin2\Documents\MATLAB\Scattering\Scattering 2\I 90keV.mat');

%Ntheta=80;
Nphi=160;


I=zeros(N,N);I_at=I;

U=90e3; %Acceleration field in V
E=U*1.6022*1e-19;m=9.1094e-31;h=6.6261e-34;c=299792458;%Physical constants
lambda=h/sqrt(2*m*E)/sqrt(1+E/(2*m*c^2));% Electron wavelength
k=2*pi/lambda;
zd=0.55; % Distance to detector

dci=2.138e-10;dcf=1.330e-10;ang=108.07/180*pi;dff=sin(ang/2)*dcf*2;
ang2=pi-asin(dff/sqrt(3)/dcf);dfi=sqrt(dcf^2+dci^2-2*dcf*dci*cos(ang2));
tt=dff/sqrt(3);hh=sqrt(dcf^2-dff^2/3);


s0=Cs.s;
cl=0.04*3; % Size of phosphor screen
[x,y]=meshgrid(0:N-1,-N+1:0);x=x/N/2*cl;y=y/N/2*cl;[~,r]=cart2pol(x,y);
the=atan(r/zd);smap=2*k*sin(the/2)/1e10;
% fcmap=interp1(s0,sqrt(Cs.f),smap);ffmap=interp1(s0,sqrt(Fs.f),smap);fimap=interp1(s0,sqrt(Is.f),smap);
fcmap=interp1(s0,Cs.f,smap,'spline');ffmap=interp1(s0,Fs.f,smap,'spline');fimap=interp1(s0,Is.f,smap,'spline');

I_at=fcmap.^2+fimap.^2+ffmap.^2*3;

ecmap=interp1(s0,Cs.n,smap);efmap=interp1(s0,Fs.n,smap);
eimap=interp1(s0,Is.n,smap);
fcmap=fcmap.*exp(1i*ecmap);ffmap=ffmap.*exp(1i*efmap);fimap=fimap.*exp(1i*eimap);
%fimap=ffmap;

ccoord=[0,0,0];icoord=[dci,0,0];fcoord1=[-hh,tt,0];fcoord2=[-hh,tt,pi*2/3];
fcoord3=[-hh,tt,-pi*2/3]; % [Z,R,Phi]

phi=2*pi/Nphi:2*pi/Nphi:2*pi;
Distance=sqrt(x.^2+y.^2+zd^2);
dkx=k*(x./Distance+0.000);dky=k*(y./Distance+0.00);dkz=k*(zd./Distance-sqrt(1-0.025e-4));
dk=sqrt(dkx.^2+dky.^2+dkz.^2);anglefromh=pi/2;
RM=[1 0 0; 0 cos(anglefromh) -sin(anglefromh); 0 sin(anglefromh) cos(anglefromh)];

% Determine number of cylindrical harmonics required
Accuracy=1e-4;AccuracyM=1+Accuracy;
zz=numel(phi);  % Changed because I made phi sin dependent on the theta value
a1=RM*[sin(theta)*cos(phi(zz));sin(theta)*sin(phi(zz));cos(theta)];
ksi=dkx*a1(1)+dky*a1(2)+dkz*a1(3);
R=sqrt(dk.^2-ksi.^2);
% I_at=abs(fcmap).^2+abs(fimap).^2+abs(ffmap).^2*3;
m=35;
while AccuracyM > Accuracy
    m=m+1;Gm=0;
Gm=fcmap.*exp(1i*(ksi*ccoord(1)-m*ccoord(3)))*1i^m.*besselj(m,R*ccoord(2))+...
                fimap.*exp(1i*(ksi*icoord(1)-m*icoord(3)))*1i^m.*besselj(m,R*icoord(2))+...
                ffmap.*exp(1i*(ksi*fcoord1(1)-m*fcoord1(3)))*1i^m.*besselj(m,R*fcoord1(2))+...
                ffmap.*exp(1i*(ksi*fcoord2(1)-m*fcoord2(3)))*1i^m.*besselj(m,R*fcoord2(2))+...
                ffmap.*exp(1i*(ksi*fcoord3(1)-m*fcoord3(3)))*1i^m.*besselj(m,R*fcoord3(2));    
    PatternM=abs(Gm)./I_at;
    AccuracyM=max(PatternM(:));
end
disp(['Cylindrical harmonics orders needed for this calculation is ', num2str(m)]);
morder=m; 
tic;
    
    for zz=1:numel(phi)
        a1=RM*[sin(theta)*cos(phi(zz));sin(theta)*sin(phi(zz));cos(theta)];
        ksi=dkx*a1(1)+dky*a1(2)+dkz*a1(3);
        R=sqrt(dk.^2-ksi.^2);
        Gm=zeros(morder,N,N);Gmn=Gm;
        for m=0:morder-1
            %ksi=k*y/z;R=k*x/z;
            Gm(m+1,:,:)=fcmap.*exp(1i*(ksi*ccoord(1)-m*ccoord(3)))*1i^m.*besselj(m,R*ccoord(2))+...
                fimap.*exp(1i*(ksi*icoord(1)-m*icoord(3)))*1i^m.*besselj(m,R*icoord(2))+...
                ffmap.*exp(1i*(ksi*fcoord1(1)-m*fcoord1(3)))*1i^m.*besselj(m,R*fcoord1(2))+...
                ffmap.*exp(1i*(ksi*fcoord2(1)-m*fcoord2(3)))*1i^m.*besselj(m,R*fcoord2(2))+...
                ffmap.*exp(1i*(ksi*fcoord3(1)-m*fcoord3(3)))*1i^m.*besselj(m,R*fcoord3(2));
            
            Gmn(m+1,:,:)=fcmap.*exp(1i*(ksi*ccoord(1)+m*ccoord(3)))*1i^m.*besselj(m,R*ccoord(2))+...
                fimap.*exp(1i*(ksi*icoord(1)+m*icoord(3)))*1i^m.*besselj(m,R*icoord(2))+...
                ffmap.*exp(1i*(ksi*fcoord1(1)+m*fcoord1(3)))*1i^m.*besselj(m,R*fcoord1(2))+...
                ffmap.*exp(1i*(ksi*fcoord2(1)+m*fcoord2(3)))*1i^m.*besselj(m,R*fcoord2(2))+...
                ffmap.*exp(1i*(ksi*fcoord3(1)+m*fcoord3(3)))*1i^m.*besselj(m,R*fcoord3(2));
        end
        Gmn(1,:,:)=0;
        I=I+squeeze(sum(abs(Gm).^2,1)+sum(abs(Gmn).^2,1))/Nphi;
        if ~mod(zz,10)
            zz
            toc
        end
    end

I_mol=I-I_at;
