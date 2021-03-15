function [smap] = MakeSmap(N,ds,center);

% Make smap
[x,y]=meshgrid((1:N)-center(2),(1:N)-center(1));[~,rr]=cart2pol(x,y);
smap=rr*ds;
smap(smap==0)=1e-10;