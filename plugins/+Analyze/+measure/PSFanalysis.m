classdef PSFanalysis<interfaces.DialogProcessor
    methods
        function obj=PSFanalysis(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_layerson','linewidth_roi','sr_roihandle','znm_min','znm_max','numberOfLayers'};
        end
        
        function out=run(obj,p)           
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
 allzstart=1;
 DZ=p.dz/1000;
%  DZ=0.1
% par.pixelsize=0.13;
% par.pixelsize=0.13;
% pixtrue=0.13;
% intmax=par.analv1;
% file=par.info.file;
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
wins=11;

xx=(-wins:wins)*fileinfo.pixsize;
% xx=(-wins:wins)*pixtrue;
smallim=double(img(ccenty-wins:ccenty+wins,ccentx-wins:ccentx+wins,:)-offset);
ssm=size(smallim);

[~,ind]=max(smallim(:));
[mx,my,mz]=ind2sub(size(smallim),ind);

%take only centralpart
window=16;
minf=max(1,mz-window);
maxf=min(mz+window,ssm(3));
% if minf==1
%     maxf=min(60,ssm(3));
% end
% if maxf==ssm(3)
%     minf=max(1,ssm(3)-60);
% end
smallim=smallim(:,:,minf:maxf);

ssm=size(smallim);
[~,ind]=max(smallim(:));
[mx,my,mz]=ind2sub(size(smallim),ind);

positions=locData.getloc({'frame','xnm','ynm','PSFxnm'},'position','roi');

badind=(positions.frame<=minf)|(positions.frame>maxf);
positions.xnm(badind)=[];positions.ynm(badind)=[];
positions.frame(badind)=[];positions.PSFxnm(badind)=[];
positions.frame=positions.frame-minf+1;



zz=((1:ssm(3))-mz)*DZ;

profilex=squeeze((sum(smallim(:,:,:),2)));
profiley=squeeze((sum(smallim(:,:,:),1)));
profilez=squeeze((sum(smallim(:,:,:),3)));
figure('Position',[1 1 800 800],'Units','Pixels')
% subplot(3,3,1)
% imagesc(xx,xx,squeeze(smallim(:,:,mz)))
% 
% subplot(3,3,2)
% imagesc(squeeze(smallim(:,my,:)))
% subplot(3,3,4)
% imagesc(squeeze(smallim(mx,:,:))')
dummypic=zeros(ssm(3));
%make one montage
combinedim=[squeeze(smallim(:,:,mz)) squeeze(smallim(:,my,:))
    squeeze(smallim(mx,:,:))' dummypic];
subplot(3,3,1)
imagesc(combinedim)
axis off
title(filen,'Interpreter','none')

subplot(3,3,5)
nind=allzstart:allzdist:ssm(3);
% size(smallim)

cuts=[smallim(:,:,nind(1)) smallim(:,:,nind(2)) smallim(:,:,nind(3))
    smallim(:,:,nind(4)) smallim(:,:,nind(5)) smallim(:,:,nind(6))
    smallim(:,:,nind(7)) smallim(:,:,nind(8)) smallim(:,:,nind(9))];
imagesc(sqrt(cuts))
axis off
title('sqrt, dz = 200 nm')

sp=size(profiley);
mini=min(profiley(:));
maxi=max(profiley(:));
startp=[maxi, 0, 1,0];
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
mx=positions.xnm(indz);
my=positions.ynm(indz);
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

zzss=zz(positions.frame);

% zzs=zzss(mz-3:mz+3);
% p1=double(positions(:,6))';
% p1s=p1(mz-3:mz+3);
% p=polyfit(zzs,p1s,2);
% m3=-p(2)/2/p(1)
m3=0;

subplot(3,3,8)
plot(zz,fitpx(3,:))
hold on
plot(zz,fitpy(3,:),'r')
plot(zz(positions.frame),positions.PSFxnm,'k')
plot([zz(mz),zz(mz)],[min(fitpy(3,:)),max(fitpy(3,:))],'r')
plot([zz(mzsig),zz(mzsig)],[min(fitpy(3,:)),max(fitpy(3,:))],'g')
hold off
axis tight
ylim([0.1 0.5])
%legend('x','y','2Dfit','Location','NorthWest')
title(['\sigma_x = ' num2str(fitpx(3,mz),3) ', \sigma_y = ' num2str(fitpy(3,mz),3)...
    ', \sigma_f = ' num2str(positions.PSFxnm(indz),3)])


subplot(3,3,6)

% posnm=positions;
% posnm(:,2:3)=posnm(:,2:3)*pixtrue*1000;
% posnm(:,2:3)=posnm(:,2:3)*par.pixelsize*1000;
minx=min(positions.xnm);miny=min(positions.ynm);

x=(positions.xnm-minx);
y=(positions.ynm-miny);

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
% dxa=max(posnm(:,2))-min(posnm(:,2));
% dya=max(posnm(:,3))-min(posnm(:,3));

title(['d_{1000} = ' num2str(sqrt(dxs.^2+dys.^2),3),...
    ', d_{all} = ' num2str(sqrt(dxa.^2+dya.^2),4) ' nm'])
axis tight

r=[1 0 0]; k=[0 0 0]; b=[0 1 0];
col=[r;r;r;r;k;b;b;b;b];

% 
% figure(33)
% plot(x(inframes),y(inframes),'+')
% 
% adsf
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

% subplot(3,3,4)
% 
% for k=1:9
% plot(xx,profiley(:,mz-k+5),'Color',col(k,:))
% hold on
% end
% plot(xx,profiley(:,mz),'k')
% hold off
% axis tight
% set(gca,'YTickLabel','')
% title('profile y')


subplot(3,3,3)
plot(xx,profilex(:,mz))
hold on
fitf=mygaussforfit(fitpx(:,mz),xx);
plot(xx,fitf,'r')
hold off
%set(gca,'YTickLabel','')
% axis tight
title(['x(z0), If= ' num2str(fitpx(1,mz))] )

% subplot(3,3,7)
% plot(xx,profiley(:,mz))
% hold on
% fitf=mygaussforfit(fitpy(:,mz),xx);
% plot(xx,fitf,'r')
% hold off
% set(gca,'YTickLabel','')
% axis tight
% title('profile y at best z')

subplot(3,3,7)
plot(zz,fitpx(1,:),zz,0*fitpx(4,:),zz,fitpy(1,:),zz,0*fitpy(4,:))
% title('I,BG vs z')
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




% subplot(3,3,4)
% 
% inbest=find(pos.original.fitpos(:,1)<mz+2+minf&pos.original.fitpos(:,1)>mz-2+minf);
% n=0:par.analv2:par.analv1;
% ints=pos.original.fitpos(inbest,4);
% hn=hist(ints,n);
% 
% 
% 
% plot(n,hn)
% hold on
% 
% [mh,mih]=max(hn);
% fitp=mygaussfit(n,hn,[mh,n(mih),n(mih)/10,0])
% plot(n,mygaussforfit(fitp,n))
% hold off
% axis tight
% s='%5.0f'
% title(['med=' num2str(median(ints),s) ', mean=' num2str(mean(ints),s) ',fit=' num2str(fitp(2),s)])
% % plot(pos.original.fitpos(inbest,4),pos.original.fitpos(inbest,5),'+');

% subplot(3,3,4)
% hold off
% plot(zz(positions(:,1)),positions(:,4),'k')
% [mh,mih]=max(positions(:,4));
% fitp=mygaussfit(zz(positions(:,1)),positions(:,4)',[mh,1,.5,0]);
% hold on
% plot(zz(positions(:,1)),mygaussforfit(fitp,zz(positions(:,1))))
% hold off
% title(fitp(3))


fileout=[fd(1:end-1) '_PSF' ];
export_fig(gcf,[fileout '.pdf'],'-pdf');
end


function pard=guidef
pard.text1.object=struct('String','dz (nm):','Style','text');
pard.text1.position=[1,1];
% 
pard.dz.object=struct('String','50','Style','edit');
pard.dz.position=[1,2];
% 
pard.plugininfo.name='PSFanalysis';
pard.plugininfo.type='ProcessorPlugin';
% pard.text2.object=struct('String','fitmodel:','Style','text');
% pard.text2.position=[2,1];
% 
% pard.fitmodel.object=struct('String','Gauss|Flat|Disk|Ring|Distance','Style','popupmenu');
% pard.fitmodel.position=[2,2];
% 
% pard.restrictsigma.object=struct('String','sigma=<locp>','Style','checkbox');
% pard.restrictsigma.position=[2,3];
% 
% pard.linelengthcheck.object=struct('String','length (nm)','Style','checkbox');
% pard.linelengthcheck.position=[3,1];
% 
% pard.linelength.object=struct('String','250','Style','edit');
% pard.linelength.position=[3,2];
end
