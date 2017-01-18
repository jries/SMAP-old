classdef CME3DDSpherefit_DC<interfaces.SEEvaluationProcessor
    properties
        savedevals
    end
    methods
        function obj=CME3DDSpherefit_DC(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function makeGui(obj)
            makeGui@interfaces.SEEvaluationProcessor(obj);
            obj.guihandles.saveimagesb.Callback={@saveimagesb_callback,obj};
        end
        function out=run(obj,p)
            out=runintern(obj,p);     
        end
        
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.layer.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|layer5|layer6');
pard.layer.position=[1,1];
pard.layer.Width=2;

pard.text1.object=struct('Style','text','String','slice thickness (nm)');
pard.text1.position=[2,1];
pard.text1.Width=2;
pard.slice.object=struct('Style','edit','String','50');
pard.slice.position=[2,3];
pard.slice.Width=2;

pard.maskon.object=struct('Style','checkbox','String','Use mask. Cutoff:','Value',1);
pard.maskon.position=[3,1];
pard.maskon.Width=2;
pard.maskcutoff.object=struct('Style','edit','String','1');
pard.maskcutoff.position=[3,3];
pard.maskcutoff.Width=2;

pard.textm.object=struct('Style','text','String','Filter mask (pix)');
pard.textm.position=[4,1];
pard.textm.Width=2;
pard.sigmamask.object=struct('Style','edit','String','1');
pard.sigmamask.position=[4,3];
pard.sigmamask.Width=2;

pard.refit.object=struct('Style','checkbox','String','refit sphere','Value',1);
pard.refit.position=[5,1];
pard.refit.Width=2;

pard.saveon.object=struct('Style','checkbox','String','save figure','Value',1);
pard.saveon.position=[6,1];
pard.saveon.Width=2;

pard.zrangecheck.object=struct('Style','checkbox','String','set z center to:','Value',1);
pard.zrangecheck.position=[7,1];
pard.zrangecheck.Width=2;
pard.zrange.object=struct('Style','edit','String','100 ');
pard.zrange.position=[7,3];
pard.zrange.Width=2;

pard.text8.object=struct('Style','text','String','sigma for rendering');
pard.text8.position=[8,1];
pard.text8.Width=2;
pard.renderfilter.object=struct('Style','edit','String','.8');
pard.renderfilter.position=[8,3];
pard.renderfilter.Width=2;


pard.saveimagescheck.object=struct('Style','checkbox','String','');
pard.saveimagescheck.position=[10,1];
pard.saveimagescheck.Width=.3;
pard.saveimagesb.object=struct('Style','pushbutton','String','save as:');
pard.saveimagesb.position=[10,1.3];
pard.saveimagesb.Width=1;
pard.saveimagesdir.object=struct('Style','edit','String','');
pard.saveimagesdir.position=[10,2.4];
pard.saveimagesdir.Width=2.6;
pard.plugininfo.type='ROI_Evaluate';
end

function saveimagesb_callback(a,b,obj)
p=obj.getGuiParameters.par;
f=p.saveimagesdir;
d=uigetdir(f);
if d
obj.guihandles.saveimagesdir.String=d;
end
end

function out=runintern(obj,p)
% obj.site.sePar.Settings
roisize=obj.site.sePar.Settings.siteroi/2;
locs=obj.getLocs({'xnm','ynm','znm'},'layer',1,'size',roisize*2);
locs2=obj.getLocs({'xnm','ynm','znm'},'layer',2,'size',roisize*2);

lenbar=20;
ranger=[0 roisize];
rangex=[-roisize roisize];

h=obj.setoutput('image');hp=h.Parent;delete(hp.Children)

medianz=median(locs.znm);
rangez=[medianz-roisize medianz+roisize];
if p.zrangecheck
    rangez=[p.zrange-roisize p.zrange+roisize];
end
pixelsx=5;
pixelsz=5;

x=locs.xnm-obj.site.pos(1);
y=locs.ynm-obj.site.pos(2);
z=locs.znm-medianz;

xr=locs.xnmrot-obj.site.pos(1);
yr=locs.ynmrot-obj.site.pos(2);

[im,mask]=findInMask(x,y,z,rangex,rangez,p);
if ~p.maskon
    im=true(length(x),1);
    mask=ones(size(mask));
end
 
out.numberOfLocsin3DMask=sum(im);
% m is unfiltered, n is filtered
xn=x(im);yn=y(im);zn=z(im);

%determine quantiles.
qx=myquantile(xn,[0.05,0.5,0.95]);
qy=myquantile(yn,[0.05,0.5,0.95]);
qz=myquantile(zn,[0.05,0.5,0.95]);

%positions corrected for median
dx=qx(2);dy=qy(2);
[fi,rn]=cart2pol(xn-dx,yn-dy);
[fi,r]=cart2pol(x-dx,y-dy);

qr=myquantile(rn,0.5);
qrangeall=0.05:0.05:0.95;
out.allqx=myquantile(xn,qrangeall);
out.allqy=myquantile(yn,qrangeall);
out.allqz=myquantile(zn,qrangeall);
out.allqr=myquantile(rn,qrangeall);
out.allrange=qrangeall;
out.qx=qx;out.qy=qy;out.qz=qz;out.qr=qr;

immask=makeprojections(xn-dx,yn-dy,zn,rn,p.slice,ranger,rangex,rangez,pixelsx,pixelsz,p.renderfilter,lenbar);
imout= makeprojections(x,y,z,r,p.slice*4,ranger,rangex,rangez,pixelsx,pixelsz,p.renderfilter,lenbar);
% plotimages1

%%spherefit
quantiles=[qx(2),qy(2),qz(2),qr];
out.spherefit=spherefit(xn,yn,zn,quantiles);
rSphere=out.spherefit(1);
fitpsphere=out.spherefit; %later rename, fitp is too generic

xw=fitpsphere(2)-dx; yw=fitpsphere(3)-dy;
xppi4=cos(pi/4)*xw+sin(pi/4)*yw; ppi4=cos(pi/4)*yw-sin(pi/4)*xw;
xmpi4=cos(-pi/4)*xw+sin(-pi/4)*yw; ympi4=cos(-pi/4)*yw-sin(-pi/4)*xw;

xc=xn-fitpsphere(2); yc=yn-fitpsphere(3); zc=zn-fitpsphere(4);

[tc,pc,rc]=cart2sph(zc,yc,xc);
out.allqtheta=myquantile(pc,qrangeall);

hsc=axes('Parent',hp);subplot(3,4,11,hsc);
hmap2=axes('Parent',hp);subplot(3,4,7,hmap2);
hmap3=axes('Parent',hp);subplot(3,4,3,hmap3);
hmaprot=axes('Parent',hp);subplot(3,4,8,hmaprot);

%fitcocerage
fitpsc=spherecoverage(tc,pc,rc,hsc);
out.fitcoverage_thetatop=min(max(-pi/2,fitpsc(2)),pi/2);
out.fitcoverage_thetabottom=min(max(-pi/2,fitpsc(3)),pi/2);

x2=locs2.xnm-obj.site.pos(1)-fitpsphere(2);
y2=locs2.ynm-obj.site.pos(2)-fitpsphere(3);
z2=locs2.znm-medianz-fitpsphere(4);
[tc2,pc2,rc2]=cart2sph(z2,y2,x2);

r2maxrel=2;
r2minrel=0;
badind=rc2<fitpsphere(1)*r2minrel|rc2>fitpsphere(1)*r2maxrel;
x2(badind)=[];y2(badind)=[];z2(badind)=[];

out.map2D=mapanalysis2D(xc,yc,zc,rSphere,hmap2,x2,y2,z2);
out.map3D=mapanalysis3D(xc,yc,zc,rSphere,hmap3,x2,y2,z2,0,[],1.1);

out.pos1=struct('x',xc,'y',yc,'z',zc,'r',rc,'p',pc,'t',tc);
out.pos2=struct('x',x2,'y',y2,'z',z2,'r',rc2,'p',pc2,'t',tc2);

%radial distribution
[tr,pr,rr]=cart2sph(xc,yc,zc);
[tr2,pr2,rr2]=cart2sph(x2,y2,z2);
hnor=axes('Parent',hp);subplot(3,4,4,hnor);
imn1p=renderpolar(rr,pr,rSphere);
imn2p=renderpolar(rr2,pr2,rSphere);
simn=size(imn1p);
imp=zeros(simn(1),simn(2),3);
imp(:,:,1)=imn1p/max(imn1p(:));
imp(:,:,2)=imn2p/max(imn2p(:));
imagesc([-2 2],[-2 2],imp,'Parent',hnor);

% correct rotation
axisp1=out.map3D.mainCentroid

[xr,yr,zr]=rotatecoord(xc,yc,zc,axisp1);
[xr2,yr2,zr2]=rotatecoord(x2,y2,z2,axisp1);
% xr=xc;yr=yc;zr=zc;xr2=x2;yr2=y2;zr2=z2;

rotateiterations=5;
isfull=out.map3D.coverageFraction<0.5;
 if isfull
     cutoff=2;
     cent=pi/2;
 else
     cutoff=0.3;
     cent=pi/2;
 end
    
% else
%     cent=pi/2;
% end

if isfull
    [xr2,yr2,zr2]=rotatecoord(xr2,yr2,zr2,[ pi/2 0]);
    [xr,yr,zr]=rotatecoord(xr,yr,zr,[pi/2 0]);
end

for k=1:rotateiterations
out.map3D=mapanalysis3D(xr,yr,zr,rSphere,hmaprot,xr2,yr2,zr2,cent,isfull,cutoff);
axisp3=out.map3D.mainCentroid
[xr,yr,zr]=rotatecoord(xr,yr,zr,axisp3);
[xr2,yr2,zr2]=rotatecoord(xr2,yr2,zr2,axisp3);
drawnow
% pause(1)

% waitforbuttonpress
end


% axisp=[0 0]
% tcr=tc-axisp(1);
% pcr=pc-axisp(2);
% 
% tcr2=tc2-axisp(1);
% pcr2=pc2-axisp(2);
% 
% tcx=pi;pcx=pi/2;
% pcr(pcr>pcx)=pcr(pcr>pcx)-2*pcx;
% pcr(pcr<-pcx)=pcr(pcr<-pcx)+2*pcx;
% tcr(tcr>tcx)=tcr(tcr>tcx)-2*tcx;
% tcr(tcr<-tcx)=tcr(tcr<-tcx)+2*tcx;
% 
% pcr2(pcr2>pcx)=pcr2(pcr2>pcx)-2*pcx;
% pcr2(pcr2<-pcx)=pcr2(pcr2<-pcx)+2*pcx;
% tcr2(tcr2>tcx)=tcr2(tcr2>tcx)-2*tcx;
% tcr2(tcr2<-tcx)=tcr2(tcr2<-tcx)+2*tcx;
% % tcr=tc+0.2;
% % pcr=pc+0.1;
% [zr,yr,xr]=sph2cart(tcr,pcr,rc);
% [zr2,yr2,xr2]=sph2cart(tcr2,pcr2,rc2);
out.map3D=mapanalysis3D(xr,yr,zr,rSphere,hmaprot,xr2,yr2,zr2,cent,isfull,cutoff);
axisp4=out.map3D.mainCentroid

[xr,yr,zr]=rotatecoord(xr,yr,zr,axisp4);
[xr2,yr2,zr2]=rotatecoord(xr2,yr2,zr2,axisp4);

if isfull
    [xr2,yr2,zr2]=rotatecoord(xr2,yr2,zr2,[0 pi]);
    [xr,yr,zr]=rotatecoord(xr,yr,zr,[0 pi]);
end


[tr,pr,rr]=cart2sph(xr,yr,zr);
[tr2,pr2,rr2]=cart2sph(xr2,yr2,zr2);
hnor2=axes('Parent',hp);subplot(3,4,12,hnor2);
imn1p=renderpolar(rr,pr,rSphere);
imn2p=renderpolar(rr2,pr2,rSphere);
simn=size(imn1p);
imp=zeros(simn(1),simn(2),3);
imp(:,:,1)=imn1p/max(imn1p(:));
imp(:,:,2)=imn2p/max(imn2p(:));
imagesc([-2 2],[-2 2],imp,'Parent',hnor2);




plotimages1
  



    

if p.saveimagescheck
    maind=p.saveimagesdir;
    lut=hot(256);
    for k=1:6
        thisd=[maind filesep 'view' num2str(k)];
        if ~exist(thisd,'dir')
            mkdir(thisd)
        end
        fo=[thisd filesep num2str(obj.site.indList) '_S' num2str(obj.site.ID) 'C' num2str(obj.site.info.cell) 'F' num2str(obj.site.info.filenumber) '.tif']
        imageo=uint8(immask(k).image*255);
        imrgb=ind2rgb(imageo,lut);
        imwrite(imrgb,fo)
    
    end
end

%filter from layer: or add with filter size and check box
    function plotimages1
        ima1=[imout(6).image,imout(6).image(:,1)*0+1, imout(2).image,imout(6).image(:,1)*0+1, imout(3).image];
        ima2=[imout(1).image,imout(1).image(:,1)*0+1, imout(4).image,imout(1).image(:,1)*0+1, imout(5).image];
        imall=[ima1;ima1(end,:)*0+1;ima2];
        %plot

        hs1=axes('Parent',hp);subplot(4,6,1,hs1);imagesc(rangex,rangez,immask(6).image,'Parent',hs1); colormap(hs1,hot);
        hs2=axes('Parent',hp);subplot(4,6,2,hs2);imagesc(rangex,rangez,immask(2).image,'Parent',hs2); colormap(hs2,hot);
        hs3=axes('Parent',hp);subplot(4,6,3,hs3);imagesc(rangex,rangez,immask(3).image,'Parent',hs3);colormap(hs3,hot);
        hs4=axes('Parent',hp);subplot(4,6,7,hs4);imagesc(rangex,rangex,immask(1).image,'Parent',hs4); colormap(hs4,hot);
        hs5=axes('Parent',hp);subplot(4,6,8,hs5);imagesc(rangex,rangez,immask(4).image,'Parent',hs5); colormap(hs5,hot);
        hs6=axes('Parent',hp);subplot(4,6,9,hs6);imagesc(rangex,rangez,immask(5).image,'Parent',hs6); colormap(hs6,hot);

        hs4.NextPlot='add';
        rectangle('Position',[qx(1)-dx,qy(1)-dy,qx(3)-qx(1),qy(3)-qy(1)],'EdgeColor',[1 1 1],'Parent',hs4)
        rectangle('Position',[-qr,-qr,2*qr,2*qr],'EdgeColor',[1 1 1],...
            'Curvature',[1,1],...
            'EdgeColor','g','Parent',hs4)
        plot(qx(2)-dx,qy(2)-dy,'c+','Parent',hs4)

        hs1.NextPlot='add';
        rectangle('Position',[-qr,qz(1),2*qr,qz(3)-qz(1)],'EdgeColor',[1 1 1],'Parent',hs1)
        plot(0,qz(2),'c+','Parent',hs1)

        colormap(hs1,hot)
        sim=size(imout(1).image);
        maskx=bwperim(imresize(squeeze((sum(mask,2))),[sim(1) sim(2)]))';
        masky=bwperim(imresize(squeeze((sum(mask,1))),[sim(1) sim(2)]))';
        maskz=bwperim(imresize(squeeze((sum(mask,3))),[sim(1) sim(2)]))';

        imxy=imout(1).image;imxy(maskz==1)=1;
        imxz=imout(2).image;imxz(maskx==1)=1;
        imyz=imout(3).image;imyz(masky==1)=1;

        ht1=axes('Parent',hp);subplot(4,6,13,ht1);imagesc(rangex,rangez,imout(6).image,'Parent',ht1); colormap(ht1,hot);
        ht2=axes('Parent',hp);subplot(4,6,14,ht2);imagesc(rangex,rangez,imxz,'Parent',ht2);colormap(ht2,hot);
        ht3=axes('Parent',hp);subplot(4,6,15,ht3);imagesc(rangex,rangez,imyz,'Parent',ht3); colormap(ht3,hot);
        ht4=axes('Parent',hp);subplot(4,6,19,ht4);imagesc(rangex,rangex,imxy,'Parent',ht4); colormap(ht4,hot);
        ht5=axes('Parent',hp);subplot(4,6,20,ht5);imagesc(rangex,rangez,imout(4).image,'Parent',ht5); colormap(ht5,hot);
        ht6=axes('Parent',hp);subplot(4,6,21,ht6);imagesc(rangex,rangez,imout(5).image,'Parent',ht6); colormap(ht6,hot);
 %%%%
        hs2.NextPlot='add';
        rectangle('Parent',hs2,'Position',[fitpsphere(2)-fitpsphere(1)-dx,fitpsphere(4)-fitpsphere(1),2*fitpsphere(1),2*fitpsphere(1)],'EdgeColor',[1 1 1], 'Curvature',[1,1])   
        plotwedge(fitpsphere(2)-dx,fitpsphere(4),fitpsphere(1),fitpsc(2),hs2,'c')
        plotwedge(fitpsphere(2)-dx,fitpsphere(4),fitpsphere(1),fitpsc(3),hs2,'c')
        hs3.NextPlot='add';
        rectangle('Parent',hs3,'Position',[fitpsphere(3)-fitpsphere(1)-dy,fitpsphere(4)-fitpsphere(1),2*fitpsphere(1),2*fitpsphere(1)],'EdgeColor',[1 1 1], 'Curvature',[1,1])
        plotwedge(fitpsphere(2)-dx,fitpsphere(4),fitpsphere(1),out.allqtheta(1),hs2,'m')
        plotwedge(fitpsphere(2)-dx,fitpsphere(4),fitpsphere(1),out.allqtheta(end),hs2,'m')    
        hs4.NextPlot='add';
        rectangle('Parent',hs4,'Position',[fitpsphere(2)-fitpsphere(1)-dx,fitpsphere(3)-fitpsphere(1)-dy,2*fitpsphere(1),2*fitpsphere(1)],'EdgeColor',[1 1 1], 'Curvature',[1,1])
        hs5.NextPlot='add';
        rectangle('Parent',hs5,'Position',[xppi4-fitpsphere(1),fitpsphere(4)-fitpsphere(1),2*fitpsphere(1),2*fitpsphere(1)],'EdgeColor',[1 1 1], 'Curvature',[1,1])
        hs6.NextPlot='add';
        rectangle('Parent',hs6,'Position',[xmpi4-fitpsphere(1),fitpsphere(4)-fitpsphere(1),2*fitpsphere(1),2*fitpsphere(1)],'EdgeColor',[1 1 1], 'Curvature',[1,1])
    end
end

function [im,mask]=findInMask(xn,yn,zn,rangex,rangez,p)
%3D reconstruction for mask - change these if masking doesn't work
%correctly
prxy=10;
prz=10;
% cutoff=1;
cutofffactor=p.maskcutoff;
% p.sigmamask=1; %to parameters????
dilationPixels=1; %integer, dilation. Two times along each direction
rangerz=rangez(1):prz:rangez(2);
% rangerx=rangex(1):prxy:rangex(2);
srz=length(rangerz);
srx=round((rangex(2)-rangex(1))/prxy);
outim=zeros(srx,srx,srz);

    
%     posr.x=0;posr.y=0;
%     posr.s=sigmarec*prxy;
%     [~,~,G]=gaussrender(posr,rangex, rangex, prxy, prxy);
    him=fspecial('gaussian',5,p.sigmamask);
      
for k=1:srz-1
    ig=zn>rangerz(k)&zn<rangerz(k+1);
    posr.x=xn(ig);posr.y=yn(ig);
%     posr.s=0*posr.x+sigmarec*prxy;
%     outim(:,:,k)=gaussrender(posr,rangex, rangex, prxy, prxy,[],[],G);
out1=histrender(posr,rangex, rangex, prxy, prxy);
 outim(:,:,k)=filter2(him,out1');
    
end
% cutoff=graythresh(outim)*2;
cutoff=quantile(outim(:),0.995)*cutofffactor;
imbw=0*outim;imbw(outim>cutoff)=1;
cc=bwconncomp(imbw);

lc=zeros(cc.NumObjects,1);
for k=1:cc.NumObjects
    lc(k)=length(cc.PixelIdxList{k});
end
[~,ix]=max(lc);
imbw2=0*outim;
imbw2(cc.PixelIdxList{ix})=1;
% imbw2=bwareaopen(imbw,3);
se = strel('disk',dilationPixels);
imbw3 = imdilate(imbw2>0,se);
imbw3 = imdilate(permute(imbw3,[2 3 1]),se);
imbw3 = imdilate(permute(imbw3,[2 3 1]),se);
imbw3=permute(imbw3,[2 3 1]);
% imbw3(imbw3>0)=1;
%  look3d(double(imbw3));

im=withinmask(imbw3,(xn-rangex(1))/prxy,(yn-rangex(1))/prxy,(zn-rangez(1))/prz);
mask=imbw3;
% sum(im)
end


function xslice=makeprojection(x,y,z,sliced,rangex,rangez,pixelsx,pixelsz,angle,renderfilter,lenbar)

xn=cos(angle)*x+sin(angle)*y;
yn=cos(angle)*y-sin(angle)*x;

confinex=[-sliced sliced];
ig=yn>confinex(1)& yn<confinex(2);
posp.x=xn(ig);
posp.y=z(ig);
% posp.s=posp.x*0+120;

xslice=histrender(posp,rangex, rangez, pixelsx, pixelsz);
if renderfilter>0
    h=fspecial('gaussian',3*ceil(renderfilter),renderfilter);
    xslice=imfilter(xslice,h);
end

% imagesc(rangex,rangez,xslice')!^!
xslice=xslice/max(xslice(:));

xslice(end-1,end-lenbar-1:end-1)=1;
end

function imout= makeprojections(xn,yn,zn,rn,sliced,ranger,rangex,rangez,pixelsx,pixelsz,renderfilter,lenbar)

posp.x=rn;
posp.y=zn;
posp.s=posp.x*0+120;

srim=histrender(posp,ranger, rangez, pixelsx, pixelsz);
[Y,X]=meshgrid(rangez(1)+pixelsz/2:pixelsz:rangez(2),ranger(1)+pixelsx/2:pixelsx:ranger(2));
% subplot(2,3,1)

imn=srim./(X'+pixelsx/4);
% imn=srim';
imn2=[imn(:,end) imn(:,end:-1:2) imn];
if renderfilter>0
    h=fspecial('gaussian',3*ceil(renderfilter),renderfilter);
    imn2=imfilter(imn2,h);
end
imn2=imn2/max(imn2(:));
imn2(end-1,end-lenbar-1:end-1)=1;
imout(6).image=imn2;
imout(6).angle=0;
imout(6).title='radialAverage_rz';


imout(1).image=makeprojection(xn,zn,yn,10,rangex,rangex,pixelsx,pixelsx,0,renderfilter,lenbar);
imout(1).angle=0;
imout(1).title='xy';

si=size(imout(1).image);
imout(1).image(si(1)/2,si(2)/2)=1;
imout(1).image(end-6:end-2,2:6)=1;


% subplot(2,3,2)
imout(2).image=makeprojection(xn,yn,zn,sliced,rangex,rangez,pixelsx,pixelsz,0,renderfilter,lenbar);
imout(2).angle=0;
imout(2).title='xz';
imout(2).image((end-5),2:12)=1;
% subplot(2,3,3)
imout(3).image=makeprojection(xn,yn,zn,sliced,rangex,rangez,pixelsx,pixelsz,pi/2,renderfilter,lenbar);
imout(3).angle=90;
imout(3).title='xz';
imout(3).image(end-12:end-2,5)=1;

% subplot(2,3,5)
imout(4).image=makeprojection(xn,yn,zn,sliced,rangex,rangez,pixelsx,pixelsz,pi/4,renderfilter,lenbar);
imout(4).angle=45;
imout(4).title='xz';
imout(4).image(end-12,2)=1;imout(4).image(end-11,3)=1;imout(4).image(end-10,4)=1;imout(4).image(end-9,5)=1;
% 
% subplot(2,3,6)
imout(5).image=makeprojection(xn,yn,zn,sliced,rangex,rangez,pixelsx,pixelsz,-pi/4,renderfilter,lenbar);
imout(5).angle=135;
imout(5).title='xz';
imout(5).image(end-12,6)=1;imout(5).image(end-11,5)=1;imout(5).image(end-10,4)=1;imout(5).image(end-9,3)=1;


end

%  
function    plotwedge(midx,midy,radius,angle,handle,appearence)
thetm=max(-pi/2,min(pi/2,angle));
plot([midx midx+radius*cos(thetm)],[midy midy+radius*sin(thetm)],appearence,'Parent',handle) 
plot([midx midx-radius*cos(thetm)],[midy midy+radius*sin(thetm)],appearence,'Parent',handle) 
end

function stat=mapanalysis2D(xc,yc,zc,rSphere,hmap,x2,y2,z2)
nump=256;
roisize=500; %half size (nm)
cutofffactor=1.5;
sigma=2.5;

rangex=[-roisize roisize];
rangey=rangex;
pixelsize=[2*roisize/nump,2*roisize/nump];
      
imh=myhist2(xc,yc,pixelsize(1),pixelsize(2),rangex,rangey);

imh2=myhist2(x2,y2,pixelsize(1),pixelsize(2),rangex,rangey);

hk=fspecial('gauss',20,sigma);
imhf=imfilter(sqrt(imh),hk,'replicate');
% imhf=imhf(nump/4+1:2.5*nump/2,:);
imh2f=imfilter(sqrt(imh2),hk,'replicate');

cutoffone=max(hk(:))*cutofffactor;
imbw=imhf>cutoffone;
imbw=imfill(imbw,'holes');
areapix=sum(imbw(:));
areanm=areapix*pixelsize(1)*pixelsize(2);
stat.areanm=areanm;

s=size(imhf);
imp=zeros(s(2),s(1),3);
imp(:,:,1)=imhf'/max(imhf(:));
imp(:,:,2)=imh2f'/max(imh2f(:));

if ~isempty(hmap)
    imagesc(rangex,rangey,(imp),'Parent',hmap);
    hmap.NextPlot='add';
%     colormap(hmap,'jet')
%     imagesc(rangex,rangey,(1-imbw')*max(imp(:)*.5),'AlphaData',0.25,'Parent',hmap);
%     plot(hmap,[-pi/2 pi/2],[0 0],'wo')
%     plot(hmap,cbx,cby,'wx')
%      plot(hmap,cdx,cdy,'w+')
    hmap.NextPlot='replace';
%     title(['coverage area (nm^2): ' num2str(stat.coverageArea,'%6.0f') ', fraction: ' num2str(stat.coverageFraction)],'Parent',hmap);
end
end

function stat=mapanalysis3D(xc,yc,zc,rSphere,hmap,x2,y2,z2,centerangle,isfull,cutofffactor)
if isempty(isfull)
    isfull=false;
end

[tc2,pc2,rc2]=cart2sph(zc,yc,xc);
tc2=tc2+centerangle;
[pc2,tc2]=rewindangles(pc2,tc2);

% tc2(tc2>pi)=tc2(tc2>pi)-2*pi;
% tc2(tc2<-pi)=tc2(tc2<-pi)+2*pi;
% pc2(pc2>pi)=tc2(tc2>pi)-2*pi;
% pc2(pc2<-pi)=tc2(tc2<-pi)+2*pi;

nump=128; %side length of reconstruction in pixels
sigma=3.5;
% cutofffactor=0.1;%1.7;

rangex=[-pi pi];
rangey=[-1 1];
pixelsize=[2*pi/nump,2/nump];
      
imh=myhist2(tc2,sin(pc2),pixelsize(1),pixelsize(2),rangex,rangey);
% sum(imh(:))
imh=[imh ;imh];
hk=fspecial('gauss',ceil(4*sigma),sigma);
imhf=imfilter(sqrt(imh),hk,'replicate');
imhf=imhf(nump/4+1:2.5*nump/2,:);
%  figure(88);imagesc(imh);waitforbuttonpress

[tc2,pc2,rc2]=cart2sph(z2,y2,x2);
tc2=tc2+centerangle;
[pc2,tc2]=rewindangles(pc2,tc2);
% tc2(tc2>2*pi)=tc2(tc2>2*pi)-2*pi;
imh=myhist2(tc2,sin(pc2),pixelsize(1),pixelsize(2),rangex,rangey);
imh=[imh ;imh];
hk=fspecial('gauss',ceil(4*sigma),sigma/2);
imhf2=imfilter(sqrt(imh),hk,'replicate');
imhf2=imhf2(nump/4+1:2.5*nump/2,:);

minco=myquantile(imhf(:),128/nump^2);

cutoffone=max(minco,max(hk(:))*cutofffactor);
imbw=imhf>cutoffone;

%         coveragepix=sum(imbw(:));
%         coverageFraction=coveragepix/numel(imbw);
statb=areastats(imbw);
stat.coverageFraction=statb.coverageFraction;
stat.coverageArea=stat.coverageFraction*4*pi*rSphere.^2;
stat.rSphere=rSphere;

cbx=statb.centroidx*pixelsize(1)-pi;
cby=statb.centroidy*pixelsize(2)-1;

statd=areastats(~imbw);
cdx=statd.centroidx*pixelsize(1)-pi;
cdy=statd.centroidy*pixelsize(2)-1;
% statb.main
% statd.main
if statd.main.istop == statb.main.istop
    disp('assignment top bottom not clear')
end

stat.topcoverage=~statd.main.istop;

stat.mainFraction=1-statd.main.areacombined/numel(imbw);
stat.mainArea=stat.mainFraction*4*pi*rSphere.^2;

%to outside
if ~isfull
cen=statd.main.Centroid;
else
cen=statb.main.Centroid;
%   cen(2)=cen(2)-nump/2;
end
    

t=cen(2)*pixelsize(1)-pi+pi/2-centerangle;
% cen(1)*pixelsize(2)-1
p=asin(cen(1)*pixelsize(2)-1);

stat.mainCentroid=[t p];
      
if ~isempty(hmap)
    perim=bwperim(imbw');
    [i,j]=ind2sub(size(perim),find(perim));
    s=size(imhf);
    imp=zeros(s(2),s(1),3);
    imp(:,:,1)=imhf'/max(imhf(:));
    imp(:,:,2)=imhf2'/max(imhf2(:));
    ind3=sub2ind(size(imp),i,j,ones(length(i),1)*3);
    
    imp(ind3)=1;
    imagesc(rangex,rangey,(imp),'Parent',hmap);
    hmap.NextPlot='add';
    colormap(hmap,'jet')
%     imagesc(rangex,rangey,(1-imbw')*max(imhf(:)*.5),'AlphaData',0.25,'Parent',hmap);
    plot(hmap,[-pi/2+centerangle pi/2+centerangle],[0 0],'wo')
    plot(hmap,cbx,cby,'wx')
    plot(hmap,cdx,cdy,'w+')
%      plot(hmap,cdx,cdy,'w+')
    hmap.NextPlot='replace';
    title(['coverage area (nm^2): ' num2str(stat.coverageArea,'%6.0f') ', fraction: ' num2str(stat.coverageFraction)],'Parent',hmap);
end
end

function stat=areastats(imbw)
main=struct('area',0,'areacombined',0,'istop',0);
stat=struct('coveragepix',0,'coverageFraction',0,'main',main,'centroidx',[],'centroidy',[],'area',[],'valid',false);

sim=size(imbw);
nump=sim(1);
bottompix=[nump/4,nump/2];
toppix=[3*nump/4,nump/2];
      
struc=regionprops(imbw,'Centroid','Extrema','Area');
if ~isempty(struc)
for k=1:length(struc)
    stat.centroidx(k)=struc(k).Centroid(2);stat.centroidy(k)=struc(k).Centroid(1);

    ex=struc(k).Extrema;
    stat.touchedge{k}=[min(ex(:,1))<1,max(ex(:,1))>nump,min(ex(:,2))<1,max(ex(:,2))>nump];
    stat.area(k)=struc(k).Area;
    stat.valid=true;
end
[~,indlargest]=max(stat.area);
dtop=sqrt((stat.centroidx(indlargest)-toppix(1)).^2+(stat.centroidy(indlargest)-toppix(2)).^2);
dbottom=sqrt((stat.centroidx(indlargest)-bottompix(1)).^2+(stat.centroidy(indlargest)-bottompix(2)).^2);
stat.main.istop=dtop<dbottom;
stat.coveragepix=sum(imbw(:));
stat.coverageFraction=stat.coveragepix/numel(imbw);
stat.main.area=stat.area(indlargest);
stat.main.Centroid=struc(indlargest).Centroid;
areacombined=stat.area(indlargest);
touchref=stat.touchedge{indlargest};
for k=1:length(stat.touchedge)
    if ~(k==indlargest)
        tochhere=stat.touchedge{k};
        if (touchref(1)&tochhere(1))||(touchref(2)&tochhere(2))||(touchref(3)&tochhere(4))||(touchref(4)&tochhere(3))
            areacombined=areacombined+stat.area(k);
        end  
    end
end

stat.main.areacombined=areacombined;
end
end


function fitp=spherefit(xn,yn,zn,quantiles)
    fh=@sphere_implicit;
    
    startptot=[quantiles(4),quantiles(1),quantiles(2),quantiles(3)+50];
    
    dz=[150  1000 -1000];
%     lb=[0 -2 -2 -2]*1000*0.25;
%     ub=[8 2 2 2]*1000*0.25;
    lb=[0    -200 -200 -2000];
    ub=[2000    200 200 2000];

    tic
    [fitp,r]=implicitfit(fh,startptot,xn,yn,zn,lb,ub,1);
    indf=1;
    for k=1:length(dz)
        startp=startptot; startp(1)=sqrt(startp(1)^2+dz(k)^2);startp(4)=startp(4)+dz(k);
        [fitph,rh]=implicitfit(fh,startp,xn,yn,zn,lb,ub,1);
        if rh<r
            indf=k;
        end
    end
    toc
    tic
    startp=startptot; startp(1)=sqrt(startp(1)^2+dz(indf)^2);startp(4)=startp(4)+dz(indf);
    [fitp,r]=implicitfit(fh,startp,xn,yn,zn,lb,ub,0);
    toc
%     fitp
end

function img=renderpolar(r,p,rSphere)
Nmax=10;
sx=0.05;
rmax=2;
px=0.03;
rnorm1=r/rSphere;
[pos1.x,pos1.y]=pol2cart(p,rnorm1);
pos1.N=1./cos(p);pos1.N(pos1.N>Nmax)=Nmax;
pos1.s=0*pos1.x+sx;
rangenx=[0 1]*rmax;rangeny=[-1 1]*rmax;
imn1=gaussrender(pos1,rangenx,rangeny, px,px);
img=horzcat(imn1(:,end:-1:1) ,imn1);    
end


function [xr,yr,zr]=rotatecoord(x,y,z,axisp)


[tc,pc,rc]=cart2sph(x,z,y);
 tcr=tc+axisp(2);
  [x,z,y]=sph2cart(tcr,pc,rc);

[tc,pc,rc]=cart2sph(z,y,x);
% pcr=pc-axisp(2);
 tcr=tc-axisp(1);
 
 [zr,yr,xr]=sph2cart(tcr,pc,rc);
%  [zm,ym,xm]=sph2cart(tcr,pc,rc);
%  
%  
%  [tcr,pcr,rc]=cart2sph(ym,xm,zm);
%  
%  
% tcr=tcr-axisp(2);
% % [pcr,tcr]=rewindangles(pcr,tcr);
% [yr,xr,zr]=sph2cart(tcr,pcr,rc);
end

function [pcr,tcr]=rewindangles(pcr,tcr)
tcx=pi;pcx=pi/2;
pcr(pcr>pcx)=pcr(pcr>pcx)-2*pcx;
pcr(pcr<-pcx)=pcr(pcr<-pcx)+2*pcx;
tcr(tcr>tcx)=tcr(tcr>tcx)-2*tcx;
tcr(tcr<-tcx)=tcr(tcr<-tcx)+2*tcx;
end