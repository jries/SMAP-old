function calibrate3D(p)
% p.filelist
% p.outputfile
% p.dz
% p.modality
% p.zcorr
% p.ROIxy
% p.ROIz
% p.smoothxy
% p.smoothz
% p.gaussrange
% p.filter;
% p.zcorrframes
%p.gaussroi

%get bead positions
disp('get bead positions')
f=figure;
p.tabgroup=uitabgroup(f);
[beads,p]=images2beads_so(p);



p.ploton=false;
if contains(p.modality,'astig')
    %determine sx,sy
    disp('fit beads to get sx,sy')
    for k=1:length(beads)
        stackh=single(beads(k).stack.image);
        s=size(stackh); 
        d=round((s(1)-p.gaussroi)/2);
        stack=stackh(d+1:end-d,d+1:end-d,:);
        P=callYimingFitter(stack,1,100,4,0,1);
        beads(k).loc.PSFxpix=P(:,5);
        beads(k).loc.PSFypix=P(:,6);
        beads(k).loc.phot=P(:,3);
        beads(k).loc.bg=P(:,4);
        beads(k).f0=stackas2z_so(beads(k).loc.PSFxpix,beads(k).loc.PSFypix,beads(k).loc.frames,beads(k).loc.phot,p);
        ind=find(beads(k).loc.frames<=beads(k).f0,1,'last');
        if isnan(beads(k).f0)||isempty(ind)
            ind=1;
        end
        beads(k).psfx0=beads(k).loc.PSFxpix(ind);
        beads(k).psfy0=beads(k).loc.PSFypix(ind);
    end
    badind=isnan([beads(:).f0]);
    beads(badind)=[];
else
    f0g=round(size(beads(1).stack,3)/2);
    for k=1:length(beads)
        beads(k).f0=f0g;
    end
end
if contains(p.modality,'astig')
    [gausscal,indgoodc]=getgausscal_so(beads,p); 
else
    indgoodc=true(size(beads));
    gausscal=[];
end
[csplinecal,indgoods]=getstackcal_so(beads(indgoodc),p);
save(p.outputfile,'gausscal','csplinecal');
end



