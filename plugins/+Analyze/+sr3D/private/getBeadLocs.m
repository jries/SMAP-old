function beadlocs=getBeadLocs(locs,p)
pos.x=locs.xnm;
pos.y=locs.ynm;
% dx=p.camPixSizeNm;
dx=p.cam_pixelsize_nm;
rangex=[min(pos.x) max(pos.x)];
rangey=[min(pos.y) max(pos.y)];
im=histrender(pos,rangex,rangey,dx,dx);

sigma=2;
h=fspecial('gaussian',3*sigma,sigma);
im1f=imfilter(im,h);
maxima=NMS2DBlockCcall(im1f,7);

minbeadsgauss=8;

co=minbeadsgauss/(sigma^2*pi*2);
% co=0.5;
indg=maxima(:,3)>co;
my=maxima(indg,1);
mx=maxima(indg,2);
mynm=my*dx+rangey(1);
mxnm=mx*dx+rangex(1);
    
    

imagesc(rangex,rangey,im1f);
hold on
plot(mxnm,mynm,'wo')
hold off

beadlocs.x=mxnm;
beadlocs.y=mynm;
end
