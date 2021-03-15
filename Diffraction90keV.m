function [I,I_mol,I_at,smap] = Diffraction90keV(path,center);


[Coor] = loadXYZ90keV(path);

[UniqueAtoms,~,~] = unique(Coor{4});

NA=length(Coor{1});
NUA=length(UniqueAtoms);

N=1024;

ds=0.0250; % inverse angstrom per pixel
if isempty(center)
    center=[round((N+1)/2),round((N+1)/2)];
end
% center = [581,539];

U=90e3; %Acceleration field in V
E=U*1.6022*1e-19;m=9.1094e-31;h=6.6261e-34;c=299792458;%Physical constants
lambda=h/sqrt(2*m*E)/sqrt(1+E/(2*m*c^2));% Electron wavelength
k=2*pi/lambda;

[s0,Scat,ScatIndex] = LoadScattering90keV(Coor{4});

[eta,~] = LoadPhase90keV(Coor{4});

[smap] = MakeSmap(N,ds,center);

% Set Variables
I=zeros(N,N);I_at=I;I_mol=0;

[fmap] = MakeScatFact90keV(s0,Scat,smap,eta,Coor{4});

for i=1:NA
    I_at=I_at+abs(fmap.(genvarname([UniqueAtoms{ScatIndex(i)}]))).^2;
end

for z1=1:NA
    for z2=1:NA
        if z1~=z2
            coij=[Coor{1}(z1) Coor{2}(z1) Coor{3}(z1)] - [Coor{1}(z2) Coor{2}(z2) Coor{3}(z2)];
            rij=sqrt(coij(1).^2+coij(2).^2+coij(3).^2);
            I_mol=I_mol+(fmap.(genvarname([UniqueAtoms{ScatIndex(z1)}])).*conj(fmap.(genvarname([UniqueAtoms{ScatIndex(z2)}])))).*sin(smap.*rij)./(smap.*rij);
        end
    end
end

I_mol=real(I_mol);
I=I_mol+I_at;
            
