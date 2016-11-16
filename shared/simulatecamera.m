function [img,simulpar]=simulatecamera(locs,p,frames)
%locs.x locs.y locs.z locs.frames locs.phot
%p.xrange p.yrange p.pixelsize p.bg p.psfmode p.calfile
%p.EMon

% p.xrange=[0 20000];p.yrange=[0 20000];p.pixelsize=100;
% p.xrange=[min(locs.x)-1000 max(locs.x)+1000];
% p.yrange=[min(locs.y)-1000 max(locs.y)+1000];
% p.EMon=0;
% p.offset=100; p.conversion=5; p.emgain=100;
% p.background=10;

p.sizex=round((p.xrange(2)-p.xrange(1))/p.pixelsize);
p.sizey=round((p.yrange(2)-p.yrange(1))/p.pixelsize);
im0=zeros(p.sizex,p.sizey,'single');
if isempty(locs)
    simulpar=p;
    img=[];
    return
end
%assume frames sorted, otherwise sort
if nargin<3
    frames=min(locs.frame):max(locs.frame);
end

numlocs=length(locs.x);
%2D in focus
locs.s=zeros(numlocs,1)+100;
img=zeros(p.sizex,p.sizey,length(frames),'single');
ind1=1;
for k=1:length(frames)
    while ind1<=numlocs&&locs.frame(ind1)<frames(k)
        ind1=ind1+1;
    end
    if locs.frame(ind1)>frames(k)
        range=[];
    else
        ind2=ind1;
        while ind2<=numlocs&&locs.frame(ind2)==frames(k)
            ind2=ind2+1;
        end    
        range=ind1:ind2-1;
        ind1=ind2;
    end
    
    if ~isempty(range)
        imh=psfimage(locs,range,p);
    else
        imh=im0;
    end
    imh2=addbg(imh,p.background);
    imh3=int2phot(imh2,p.EMon);   
    imh4=phot2adu(imh3,p);
    img(:,:,k)=imh4;
    
    if isfield(p,'plotaxis') &&isvalid(p.plotaxis)
        figure(88);
        imagesc([ imh4],'Parent',p.plotaxis)
        colorbar
        drawnow
    end
%     waitforbuttonpress;
    
    
end
simulpar=p;
end

function imo=phot2adu(imin,p)
if ~p.EMon
    emgain=1;
else
    emgain=p.emgain;
end
imo=imin*emgain/p.conversion+p.offset;
end

function img=psfimage(locs,range,p)
locsh.x=locs.x(range);
locsh.y=locs.y(range);
locsh.s=locs.s(range);
locsh.N=locs.phot(range);
[img,nlocs,Gc]=gaussrender(locsh,p.xrange, p.yrange, p.pixelsize, p.pixelsize);
end

function imh2=addbg(imh,bg)
imh2=imh+bg;
end

function imh3=int2phot(imh2,emon)
if emon
    em=2;
else
    em=1;
end
imh3=poissrnd(imh2/em)*em;
end