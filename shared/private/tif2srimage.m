function im=tif2srimage(file,p)
rangext=[p.sr_pos(1)-p.sr_size(1) p.sr_pos(1)+p.sr_size(1)];
rangeyt=[p.sr_pos(2)-p.sr_size(2) p.sr_pos(2)+p.sr_size(2)]; 

rangex=rangext-p.shiftxy_min;
rangey=rangeyt-p.shiftxy_max;
fs=p.ch_filelist.Value;
tnum=p.renderp1.Value;

%  file=lp.files.file;
 
 fileh=file(fs);
 
 if length(fileh.tif)<tnum||isempty(fileh.tif(tnum).image)
     im.image=zeros(p.sr_sizeRecPix(1),p.sr_sizeRecPix(2));
     im.rangex=rangex+p.shiftxy_min;
     im.rangey=rangey+p.shiftxy_max;
     return
 end
 
 if isfield(fileh.tif(tnum).info,'pixsize') 
    pixsize=fileh.tif(tnum).info.pixsize;
    roi=fileh.tif(tnum).info.roi;
 else
    pixsize=fileh.info.pixsize;
    roi=fileh.info.roi;
 end
 
srec=round(p.sr_sizeRecPix);
rangex=rangex+pixsize*1000/2;rangey=rangey+pixsize*1000/2;
rangexpix=rangex/1000/pixsize;rangeypix=rangey/1000/pixsize;
% rpixrx=[floor(rangexpix(1)) ceil(rangexpix(2))];rpixry=[floor(rangeypix(1)) ceil(rangeypix(2))];
rpixrx=[round(rangexpix(1)) round(rangexpix(2))];rpixry=[round(rangeypix(1)) round(rangeypix(2))];
position=[rpixrx(1),rpixrx(2),rpixry(1),rpixry(2)];

positionc=position;
positionc(1:2)=position(1:2)-roi(1);
positionc(3:4)=position(3:4)-roi(2);

coim=cutoutim(permute(double(fileh.tif(tnum).image),[2 1 3]),positionc);

magnification=pixsize/p.sr_pixrec*1000;
srcoim=imresize(coim,magnification,'nearest');

psr=[rangex(:)/p.sr_pixrec-position(1)*magnification+1; rangey(:)/p.sr_pixrec-position(3)*magnification+1];
psr(2)=psr(1)+srec(1)-1;
psr(4)=psr(3)+srec(2)-1;
srfinal=cutoutim(srcoim,psr);

im.image=permute(srfinal,[2,1,3]);
im.rangex=rangext;
im.rangey=rangeyt;