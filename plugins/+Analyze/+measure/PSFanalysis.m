classdef PSFanalysis<interfaces.DialogProcessor
%     PSFanalysis Analyzes the point-spread function from a z-stack of
%     images
    properties
        figure
    end
    methods
        function obj=PSFanalysis(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_layerson','linewidth_roi','sr_roihandle','znm_min','znm_max','numberOfLayers'};
            obj.history=false;
            obj.showresults=false;
        end
        
        function out=run(obj,p)  
            if isempty(obj.figure)||~isvalid(obj.figure)
               obj.figure=figure('Position',[1 1 800 800],'Units','Pixels');
            end
            p.figure=obj.figure;
            analyze_PSF(obj.locData,p) 
            out=[];
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end






function analyze_PSF(locData,p) %analyze PSF 
allzdist=3;
% allzstart=1;
DZ=p.dz/1000;
fileinfo=locData.files.file(1).info;
file=locData.files.file(1).name;

fd=[file(1:end-8) '/'];
allf=dir([fd '*.tif']);
if isempty(allf)||~exist([fd allf(1).name],'file')
    [~,fd]=uigetfile([fileparts(file) filesep '*.tif'],'select raw image');
    allf=dir([fd '*.tif']);
end
l=length(allf);
img=mytiffreadgeneral([fd allf(1).name],1:l);

[~,filen]=fileparts(file);
%
offset=min(img(:));
rp=p.sr_roihandle.getPosition;
croi=fileinfo.roi;
c=rp/fileinfo.pixsize;
% c=par.roicoords/par.pixelsize;
ccentx=round(c(1)+c(3)/2)-croi(1);
ccenty=round(c(2)+c(4)/2)-croi(2);
wins=round((p.roisize-1)/2);

xx=(-wins:wins)*fileinfo.pixsize*1000;
% xx=(-wins:wins)*pixtrue;
smallim=double(img(ccenty-wins:ccenty+wins,ccentx-wins:ccentx+wins,:)-offset);
ssm=size(smallim);

[~,ind]=max(smallim(:));
[mx,my,mz]=ind2sub(size(smallim),ind);

%take only centralpart
window=round(p.numframes/2);
minf=max(1,mz-window);
maxf=min(mz+window,ssm(3));

smallim=smallim(:,:,minf:maxf);

ssm=size(smallim);
[~,ind]=max(smallim(:));
[mx,my,mz]=ind2sub(size(smallim),ind);

positions=locData.getloc({'frame','xnm','ynm','PSFxnm','PSFynm'},'position','roi');

badind=(positions.frame<=minf)|(positions.frame>maxf);
positions.xnm(badind)=[];positions.ynm(badind)=[];
positions.frame(badind)=[];
positions.PSFxnm(badind)=[];positions.PSFynm(badind)=[];
positions.frame=positions.frame-minf+1;

% zz=((1:ssm(3))-mz)*DZ;

profilex=squeeze((sum(smallim(:,:,:),2)));
profiley=squeeze((sum(smallim(:,:,:),1)));
% profilez=squeeze((sum(smallim(:,:,:),3)));


dummypic=zeros(ssm(3));
%make one montage
combinedim=[squeeze(smallim(:,:,mz)) squeeze(smallim(:,my,:))
    squeeze(smallim(mx,:,:))' dummypic];
figure(p.figure)
subplot(3,3,1)
imagesc(combinedim)
axis off
title(filen,'Interpreter','none')

subplot(3,3,4)
nind0=mz-allzdist*4:allzdist:mz+allzdist*4;
nind0(nind0<1)=1;nind0(nind0>l)=l;
nind=ones(9,1);
nind(1:length(nind0))=nind0;
% size(smallim)

cuts=[smallim(:,:,nind(1)) smallim(:,:,nind(2)) smallim(:,:,nind(3))
    smallim(:,:,nind(4)) smallim(:,:,nind(5)) smallim(:,:,nind(6))
    smallim(:,:,nind(7)) smallim(:,:,nind(8)) smallim(:,:,nind(9))];
imagesc((cuts))
% imagesc(sqrt(cuts))
axis off
title(['dz = ' num2str(p.dz*allzdist)  ' nm'])

sp=size(profiley);
% mini=min(profiley(:));
maxi=max(profiley(:));
startp=[maxi, 0, 100,0];
fitpx=zeros(4,sp(2));fitpy=zeros(4,sp(2));
for k=1:sp(2)
    fitpx(:,k)=mygaussfit(xx,profilex(:,k)',startp);
    fitpy(:,k)=mygaussfit(xx,profiley(:,k)',startp);
end

mzo=mz;
[~,mzi]=max(fitpx(1,mz-3:mz+3)+fitpy(1,mz-3:mz+3));
[~,mzsig]=min(fitpx(3,mz-3:mz+3).^2+fitpy(3,mz-3:mz+3).^2);
mz=mzo+mzi-4;
mzsig=mzo+mzsig-4;

if mz<8
    mz=mzo;
end
indz=find(positions.frame>=mz,1,'first');
% mx=positions.xnm(indz);
% my=positions.ynm(indz);
zz=((1:ssm(3)))*DZ;

zzs=zz(mz-3:mz+3);
p1=fitpx(3,:);
p1s=p1(mz-3:mz+3);
p=polyfit(zzs,p1s,2);
m1=-p(2)/2/p(1);

p1=fitpy(3,:);
p1s=p1(mz-3:mz+3);
p=polyfit(zzs,p1s,2);
m2=-p(2)/2/p(1);

% zzss=zz(positions.frame);

m3=0;

subplot(3,3,8)
plot(zz,fitpx(3,:))
hold on
plot(zz,fitpy(3,:),'r')
plot(zz(positions.frame),positions.PSFxnm,'k')
if ~isempty(positions.PSFynm)
    plot(zz(positions.frame),positions.PSFynm,'k')
end
plot([zz(mz),zz(mz)],[min(fitpy(3,:)),max(fitpy(3,:))],'r')
plot([zz(mzsig),zz(mzsig)],[min(fitpy(3,:)),max(fitpy(3,:))],'g')
hold off
axis tight
ylim([0 0.5]*1000)
%legend('x','y','2Dfit','Location','NorthWest')
title(['\sigma_x = ' num2str(fitpx(3,mz),3) ', \sigma_y = ' num2str(fitpy(3,mz),3)...
    ', \sigma_f = ' num2str(positions.PSFxnm(indz),3)])

% minx=min(positions.xnm);miny=min(positions.ynm);
minx=positions.xnm(indz);miny=positions.ynm(indz);
x=(positions.xnm-minx);
y=(positions.ynm-miny);

subplot(3,3,5)
plot(positions.frame,x,positions.frame,y)
dx1=myquantilefast((x),[0.05 .95]);dy1=myquantilefast(abs(y),[0.05 .95]);
ylim(dy1);
subplot(3,3,6)

plot(x, y,'+')
inframes=(positions.frame>mz-5)&(positions.frame<mz+5);
inframes1=(positions.frame>mz-10)&(positions.frame<mz+10);
hold on
plot(x(inframes1), y(inframes1),'g+')
plot(x(inframes), y(inframes),'r+')
hold off

dxs=max(positions.xnm(inframes1))-min(positions.xnm(inframes1));
dys=max(positions.ynm(inframes1))-min(positions.ynm(inframes1));
dxa=max(positions.xnm)-min(positions.xnm);
dya=max(positions.ynm)-min(positions.ynm);

title(['d_{1000} = ' num2str(sqrt(dxs.^2+dys.^2),3),...
    ', d_{all} = ' num2str(sqrt(dxa.^2+dya.^2),4) ' nm'])
axis tight

r=[1 0 0]; k=[0 0 0]; b=[0 1 0];
col=[r;r;r;r;k;b;b;b;b];

subplot(3,3,2)
for k=1:9
plot(xx,profilex(:,mz-k+5),'Color',col(k,:))
hold on
end
plot(xx,profilex(:,mz),'k')
hold off
axis tight
set(gca,'YTickLabel','')
title('profile x')


subplot(3,3,3)
plot(xx,profilex(:,mz))
hold on
fitf=mygaussforfit(fitpx(:,mz),xx);
plot(xx,fitf,'r')
hold off
title(['x(z0), If= ' num2str(fitpx(1,mz))] )


subplot(3,3,7)
plot(zz,fitpx(1,:),zz,0*fitpx(4,:),zz,fitpy(1,:),zz,0*fitpy(4,:))
hold on
[~,inx]=max(profilex(:));
[mxx,mzz]=ind2sub(size(profilex),inx);

profz=profilex(mxx,:);
fitpz=mygaussfit(zz,profz,[max(profz),1,0.5,0]);
profzf=mygaussforfit(fitpz,zz);
plot(zz,profz,zz,profzf);
hold off

mf=max(profzf);
ylim([-mf*.2 mf*1.5])
title(['z: \sigma = ' num2str(abs(fitpz(3)))]);
xlabel('z')

plotpos=max(1,minf-20);

subplot(3,3,9)
imagesc(img(:,:,plotpos))
hold on
plot(ccentx,ccenty,'ow')
hold off
axis off

title([num2str(m1-1,3) ', ' num2str(m2-1,3)  ', ' num2str(m3-1,3)])

fileout=[fd(1:end-1) '_PSF' ];
export_fig(gcf,[fileout '.pdf'],'-pdf');
end


function pard=guidef

pard.text0.object=struct('String','Select ROI around single bead','Style','text');
pard.text0.position=[1,1];
pard.text0.Width=4;

pard.text1.object=struct('String','dz (nm):','Style','text');
pard.text1.position=[2,1]; 
pard.dz.object=struct('String','50','Style','edit');
pard.dz.position=[2,2];

pard.text2.object=struct('String','ROI size (pixel):','Style','text');
pard.text2.position=[3,1]; 
pard.roisize.object=struct('String','11','Style','edit');
pard.roisize.position=[3,2];

pard.text3.object=struct('String','number of frames:','Style','text');
pard.text3.position=[4,1]; 
pard.numframes.object=struct('String','32','Style','edit');
pard.numframes.position=[4,2];
pard.exportpdf.object=struct('String','Export results as PDF','Style','checkbox','Value',0);
pard.exportpdf.position=[5,1];
% 
pard.plugininfo.name='PSFanalysis';
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='PSFanalysis Analyzes the point-spread function from a z-stack of images';
end
