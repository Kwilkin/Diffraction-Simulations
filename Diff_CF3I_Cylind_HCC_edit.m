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

function [I_mol,I,I_at,smap]=Diff_CF3I_Cylind_HCC_edit(N,AGD,bb)
if nargin < 2 || isempty(N)
    N=64;
end
Ntheta=40;
Nphi=40;
% [~,AGD]=calc_angdis(cosdis,Ntheta);
%[~,bc]=textread('C:\Users\ASUN\Documents\MATLAB\Cross sections\scs c25.txt','%f %f');
%[~,bf]=textread('C:\Users\ASUN\Documents\MATLAB\Cross sections\scs f25.txt','%f %f');
%[a,bi]=textread('C:\Users\ASUN\Documents\MATLAB\Cross sections\scs i25.txt','%f %f');a=a/180*pi;
bc=bb.c25;bf=bb.f25;bi=bb.i25;a=bb.a*pi/180;
I=zeros(N,N);I_at=I;lambda=0.0766e-10;k=2*pi/lambda; zd=0.204;
dci=2.138e-10;dcf=1.330e-10;ang=108.07/180*pi;dff=sin(ang/2)*dcf*2;
ang2=pi-asin(dff/sqrt(3)/dcf);dfi=sqrt(dcf^2+dci^2-2*dcf*dci*cos(ang2));
tt=dff/sqrt(3);hh=sqrt(dcf^2-dff^2/3);
eta_C_40keV = [0.0541    0.0585    0.0713    0.0913    0.1165    0.1445    0.1730    0.2004    0.2256 0.2485    0.2689    0.2872    0.3037    0.3187    0.3326    0.3456    0.3577    0.3693 0.3802    0.3908    0.4009  ];
eta_C_10keV=[0.10201 0.11089 0.13649 0.17586 0.22496 0.2791 0.33405 0.38662 0.43505 0.47877 0.51788 0.55297 0.58472 0.61369 0.6404 0.66532 0.68876 0.71099 0.73218 0.75251 0.77206];
eta_F_40keV=[0.087053 0.090491 0.10064 0.11676 0.13781 0.16269 0.19031 0.2196 0.24963 0.2796 0.30888 0.33703 0.36382 0.38908 0.41278 0.43496 0.45571 0.47514 0.49337 0.51054 0.52675];
eta_F_10keV=[0.16411 0.17076 0.19016 0.22094 0.26122 0.30892 0.36194 0.41823 0.47594 0.53352 0.58974 0.64374 0.69511 0.74354 0.78898 0.83153 0.87138 0.90873 0.94383 0.9769 1.0082];
eta_I_40keV=[0.24909 0.27532 0.3506 0.46199 0.58889 0.71468 0.83374 0.94775 1.0598 1.1717 1.2833 1.3935 1.5009 1.6041 1.7021 1.7948 1.8821 1.9646 2.0427 2.1173 2.1888];
eta_I_10keV=[0.37139 0.41306 0.53506 0.72104 0.93788 1.155 1.3612 1.5594 1.7548 1.9497 2.1428 2.3314 2.5132 2.6864 2.8506 3.0062 3.1544 3.296 3.4323 3.5640 3.6917];
eta_C_25keV=(eta_C_40keV+eta_C_10keV)/2;
eta_F_25keV=(eta_F_40keV+eta_F_10keV)/2;
eta_I_25keV=(eta_I_40keV+eta_I_10keV)/2;
cl=7.4e-3*40/11*1.4;
[x,y]=meshgrid(-N/2+1:N/2,-N/2+1:N/2);x=x/N*cl;y=y/N*cl;[~,r]=cart2pol(x,y);
the=atan(r/zd);smap=2*k*sin(the/2)/1e10;
fcmap=interp1(a,sqrt(bc),the);ffmap=interp1(a,sqrt(bf),the);fimap=interp1(a,sqrt(bi),the);
ecmap=interp1(0:20,eta_C_25keV,smap);efmap=interp1(0:20,eta_F_25keV,smap);
eimap=interp1(0:20,eta_I_25keV,smap);
fcmap=fcmap.*exp(1i*ecmap);ffmap=ffmap.*exp(1i*efmap);fimap=fimap.*exp(1i*eimap);
%fimap=ffmap;
ccoord=[0-dci,0,0];icoord=[dci-dci,0,0];fcoord1=[-hh-dci,tt,0];fcoord2=[-hh-dci,tt,pi*2/3];
fcoord3=[-hh-dci,tt,-pi*2/3]; % [Z,R,Phi]
%theta=pi/Ntheta:pi/Ntheta:pi;
theta=linspace(pi/10-pi/10/2,pi-pi/10/2,10);
phi=2*pi/Nphi:2*pi/Nphi:2*pi;
Distance=sqrt(x.^2+y.^2+zd^2);
dkx=k*(x./Distance+0.000);dky=k*(y./Distance+0.00);dkz=k*(zd./Distance-sqrt(1-0.025e-4));
dk=sqrt(dkx.^2+dky.^2+dkz.^2);anglefromh=pi/2;
RM=[1 0 0; 0 cos(anglefromh) -sin(anglefromh); 0 sin(anglefromh) cos(anglefromh)];
% Determine number of cylindrical harmonics required
Accuracy=1e-4;AccuracyM=1+Accuracy;
z=numel(theta);
%phi=2*pi;
zz=numel(phi);  % Changed because I made phi sin dependent on the theta value
a1=RM*[sin(theta(z))*cos(phi(zz));sin(theta(z))*sin(phi(zz));cos(theta(z))];
ksi=dkx*a1(1)+dky*a1(2)+dkz*a1(3);
R=sqrt(dk.^2-ksi.^2);
I_at=abs(fcmap).^2+abs(fimap).^2+abs(ffmap).^2*3;
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
for z=1:numel(theta)
    tic;
    fprintf('\nTheta step %i ',z)
    %Nphi=round(Ntheta*2*sin(theta(z)));
    %phi=2*pi/Nphi:2*pi/Nphi:2*pi;
    Itheta=zeros(N,N);
    for zz=1:numel(phi)
        a1=RM*[sin(theta(z))*cos(phi(zz));sin(theta(z))*sin(phi(zz));cos(theta(z))];
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
        Itheta=Itheta+squeeze(sum(abs(Gm).^2,1)+sum(abs(Gmn).^2,1))/Nphi;
    end
    I=I+Itheta*AGD(z);	
    toc
end
I_mol=I-I_at;
