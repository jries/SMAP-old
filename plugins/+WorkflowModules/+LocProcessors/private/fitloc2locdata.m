function locdat=fitloc2locdata(obj,locs,indin)
fieldsremove={'xerrpix','yerrpix','PSFxpix','PSFypix','xpix','ypix'};
keeptype={'frame','filenumber'};
fn=fieldnames(locs);
keepfields=setdiff(fn,fieldsremove);

for k=1:length(keepfields)
    if any(strcmp(keepfields{k},keeptype))
    locdat.(keepfields{k})=locs.(keepfields{k})(indin);
    else
        locdat.(keepfields{k})=single(locs.(keepfields{k})(indin));
    end
end
if isprop(obj,'fileinfo')
    pixelsize=obj.fileinfo.cam_pixelsize_um*1000;
    if myisfield(obj.fileinfo,'roi')&&~isempty(obj.fileinfo.roi)
    roi=obj.fileinfo.roi;
    else
        roi=zeros(2);
    end
else
    pixelsize=[1 1];
    roi=zeros(2);
end

locdat.xnm=(locs.xpix(indin)+roi(1))*pixelsize(1);
locdat.ynm=(locs.ypix(indin)+roi(2))*pixelsize(end);

if isfield(locs,'xerrpix')
locdat.xerr=locs.xerrpix(indin)*pixelsize(1);
else
    locdat.xerr=locdat.xnm*0+1;
end

if isfield(locs,'yerrpix')
locdat.yerr=locs.yerrpix(indin)*pixelsize(end);
else
    locdat.yerr=locdat.xerr;
end
if isfield(locs,'PSFxpix')
locdat.PSFxnm=locs.PSFxpix(indin)*pixelsize(1);
else
    locdat.PSFxnm=0*locdat.xnm+100;
end
if isfield(locs,'PSFypix')
    locdat.PSFynm=locs.PSFypix(indin)*pixelsize(end);
else
    locdat.PSFynm=locdat.PSFxnm;
end
locdat.locprecnm=sqrt((locdat.xerr.^2+locdat.yerr.^2)/2);
if isfield(locs,'zerr')
    locdat.locprecznm=locs.zerr(indin);
end


locdat.filenumber=uint8(0*locdat.xnm+obj.filenumber);
locdat.channel=0*locdat.xnm;
end