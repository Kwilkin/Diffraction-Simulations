function [Coor] = loadXYZ90keV(path)
% Coor is a 1 x 5 cell {x y x type AtomicNumber}
% x, y, z, and AtomicNumber are 1 x N double vectors where N is the number
% of atoms you have in your molecule
% type is a 1 x N cell with each element a string and N is the number of
% atoms you have in your molecule
global FileName
if nargin < 1 || isempty(path)
    [FileName,PathName] = uigetfile('C:\Users\Kyle\Documents\MATLAB\Coordinates\*.xyz','Select the file with your coordinates');
    path=[PathName, FileName];
end

[Type,X,Y,Z]=textread(path,'%s %f %f %f','commentstyle','shell');


% %% Atom Types
% H=1;O=6;N=7;
% AtomicIndex=zeros(length(Type),1);
% for i=1:length(Type)
%     if strcmp(Type{i},'H')
%         AtomIndex(i)=1;
%     elseif strcmp(Type{i},'N')
%         AtomIndex(i)=7;
%     elseif strcmp(Type{i},'O')
%         AtomIndex(i)=8;
%     end
% end

% figure;scatter3(X,Y,Z,AtomIndex*10);title('Close to continue')
% uiwait

% Coor.X=X;
% Coor.Y=Y;
% Coor.Z=Z;
% Coor.Type=Type;
% Coor.AtomicIndex=AtomicIndex;
Type=Type(3:end);
X=X(3:end);
Y=Y(3:end);
Z=Z(3:end);
Coor={X Y Z Type};