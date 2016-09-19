
function imout=drawerSMAP(him,p)
if nargin==0
    imout={'imaxtoggle','imax_min','sr_sizeRecPix','lut','colorfield_min','colorfield_max','gamma','normalizeFoV'};
    return
end

% si=round(p.sr_sizeRecPix([2 1]));

img=him.image;
si=size(img);
if isfield(p,'normalizeFoV')&&~isempty(p.normalizeFoV)&&p.normalizeFoV>0
    s=size(img)/2;
    x=round(s(1)-p.normalizeFoV:s(1)+p.normalizeFoV);
    y=round(s(2)-p.normalizeFoV:s(2)+p.normalizeFoV);
    imgnorm=img(x,y,:);
else
    imgnorm=img;
end

[imgn,norm]=normalizeImage(img,p.imaxtoggle,p.imax_min,imgnorm);
if p.gamma ~=1
    imgn=imgn.^p.gamma;
end

[iml,lut]=applyLut(imgn,p.lut.selection,p.colorfield_min,p.colorfield_max);
imout.image=iml;
if him.istiff==0
    imout.mask=double(makemask(imgn));
else
    him.numberOfLocs=0;
    imout.mask=zeros(si(1),si(2),3,'double');
end
imout.istiff=him.istiff;
imout.rangex=him.rangex;
imout.rangey=him.rangey;
imout.lut=lut;
imout.imax=norm;
imout.numberOfLocs=him.numberOfLocs;      
end


function [imout,norm]=normalizeImage(img,imaxtoggle,imax,imgnorm)
    if imaxtoggle %quantile
        if imax<0
            imax=1-10^imax;
        end
         norm=myquantilefast(imgnorm(:),imax,30/(1-imax));
        if norm==0
            norm=max(imgnorm(:));
        end
    else
        norm=imax;
    end
    if norm~=0
    imout=img/norm;
    imout(imout>1)=1;
    else
        imout=img;
    end
end

function [imo,lut]=applyLut(im,lutname,pmin,pmax)
s=size(im);
    lut=mymakelut(lutname);
if length(s)==2||s(3)==1
    im=im-pmin;
    im=im/(pmax-pmin);
    im(im<0)=0;
    im(im>1)=1;
    
%     im(im<pmin)=pmin;
%     im(im>pmax)=pmax;
%     im(1)=pmin;
%     im(2)=pmax;

    imo=(ind2rgb(uint8(im*2^8),lut));
%     whos imo im
else
%     if isa(im,'single')||isa(im,'double')
%         imo=uint8(im*2^8);
%     else
    imo=double(im);
%     end
end


end

function mask3=makemask(image)
maskfactor=3;
im2D=sum(image,3);
mx=max(im2D(:));
if mx~=0
im2D=im2D/mx;
mask=im2D*maskfactor;
mask(mask>1)=1;
else
mask=0*im2D;   
end
mask3=repmat(mask,1,1,3);
end
