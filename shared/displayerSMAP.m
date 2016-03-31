function [imout,sr_imagehandle]=displayerSMAP(layers,p)

if nargin==0
    %input parameters
    imout={'sr_layerson','sr_axes','sr_sizeRecPix','roihandle','sr_pixrec','rotationangle','sr_pos','sr_size'};
    return          
end

rim=0;

show=0;
tiffthere=0;
txtN='';

for k=1:length(layers)

    if p.sr_layerson(k)
        if isfield(layers(k),'images') && ~isempty(layers(k).images)
            if show==0
                si=size(layers(k).images.finalImages.image);
                 imall=zeros(si(1)+rim,si(2)+rim,3);
                 mask=zeros(si(1)+rim,si(2)+rim,3);
                 imtiff=zeros(si(1)+rim,si(2)+rim,3);
%                  imref=layers(k).images.finalImages;

                show=1;
            end
         fi=layers(k).images.finalImages;
         if fi.istiff==1
                imtiff=imtiff+fi.image;
                tiffthere=1;
         else
             imall=imall+fi.image;
             mask=mask+fi.mask;
%              txtN=[txtN 'N'  num2str(k) '=' shortnumber(fi.numberOfLocs) ', '];
         end

        end
    end
end
%make color bars
if ~show
    imout=[];
    sr_imagehandle=[];
    return
end


if tiffthere
    sm=size(mask);
    if min(sm(1:2))>4
        mask(1:4,:,:)=1;
        mask(:,1:4,:)=1;
        mask(end-4:end,:,:)=1;
        mask(:,end-4:end,:)=1;
    end
    mask(mask>1)=1;
    
    imfinal=mask.*imall+(1-mask).*imtiff;
else
    imfinal=imall;
end

compimage=imfinal;    
%rotate
if isfield(p,'rotationangle')&&~isempty(p.rotationangle)&&p.rotationangle~=0
    imfinal=imrotate(imfinal,p.rotationangle,'nearest','crop');
end

for k=1:min(4,length(layers))
    if p.sr_layerson(k)&&~isempty(layers(k).images)
        imfinal=addcolorbar(imfinal,layers(k).images.finalImages.lut,k);
        rangexplot=layers(k).images.finalImages.rangex;
         rangeyplot=layers(k).images.finalImages.rangey;
    end
end

    imfinal=addscalebar(imfinal,p.sr_pixrec(1));
    if isfield(p,'sr_axes')&&~isempty(p.sr_axes)&&ishandle(p.sr_axes)
        sr_imagehandle=image(rangexplot/1000,rangeyplot/1000,imfinal,'Parent',p.sr_axes,'Pickable','none','HitTest','off');
%                     plotovim=1;
        title(p.sr_axes,txtN)
        set(p.sr_axes,'Xlim',rangexplot/1000)
        set(p.sr_axes,'Ylim',rangeyplot/1000)
        set(p.sr_axes,'YDir','reverse')
        p.sr_axes.HitTest='on';
        p.sr_axes.PickableParts='all';
        axes(p.sr_axes) %bring to forground
        imout.handle=sr_imagehandle;
        drawnow limitrate nocallbacks
    else
        sr_imagehandle=[];
    end
    imout.image=imfinal;
    imout.composite=compimage;
    imout.rangex=rangexplot/1000;
    imout.rangey=rangeyplot/1000;

end


function imout=addcolorbar(imin,lut,layer)
s=size(imin);
if layer==1||layer==3
    l=s(1);
else
    l=s(2);
end
if layer<4
    ll=4;
else
    ll=6;
end
rim=10;
x=ceil((1:l-rim)/(l-rim)*255);
y=1:ll;
[X,~]=meshgrid(x,y);
lutim=ind2rgb(uint8(X),lut);
imout=imin;
switch layer
    case 3
        imout(rim/2+1:end-rim/2,1:4,:)=permute(lutim,[2 1 3]);
    case 2
%         imout(end-3:end,rim/2+1:end-rim/2,:)=lutim;
        imout(1:4,rim/2+1:end-rim/2,:)=lutim;
    case 1
        imout(rim/2+1:end-rim/2,end-3:end,:)=permute(lutim,[2 1 3]);
    case 4
%         imout(1:6,rim/2+1:end-rim/2,:)=lutim;
        imout(end-5:end,rim/2+1:end-rim/2,:)=lutim;
end

end

function imin=addscalebar(imin,pixrec,fac)
sim=size(imin);
if nargin<3
    fac=1;
end
   lennm=10^floor(log10(sim(2)*pixrec*0.7*fac));
   lenpix=round(lennm/pixrec);
   if lenpix<sim(2)-12&&sim(1)>4
   imin(end-4:end-1,end-11-lenpix:end-9,:)=0;
   imin(end-3:end-2,end-10-lenpix:end-10,:)=1;
   end
end
