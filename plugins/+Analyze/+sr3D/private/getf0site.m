function [f0,PSFx0,PSFy0]=getf0site(locs,p)
p.ploton=0;
[zas,zn]=stackas2z(locs.PSFxpix,locs.PSFypix,locs.frame,locs.phot,p);
if isnan(zas)
    PSFx0=NaN;
    PSFy0=NaN;
else
    ind=find(locs.frame<=zas,1,'last');
    if isempty(ind)
        ind=1;
    end
PSFx0=locs.PSFxpix(ind);
PSFy0=locs.PSFypix(ind);
end
% drawnow
f0=zas;
end