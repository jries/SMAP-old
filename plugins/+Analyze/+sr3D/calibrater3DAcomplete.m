classdef calibrater3DAcomplete<interfaces.DialogProcessor
    properties

    end
    methods
        function obj=calibrater3DAcomplete(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.showresults=true;
        end
        function intiGui(obj)
            beaddistribution_callback(0,0,obj)           
        end
        function out=run(obj,p)
            out=[];
            locsall=obj.locData.getloc({'frame','xnm','ynm','PSFxnm','PSFynm','filenumber','phot'},'position','roi','layer',1,'removeFilter','filenumber');
%             locsall=obj.locData.getloc('frame','xnm','ynm','PSFxnm','PSFynm','filenumber','phot');
            if isempty(locsall.PSFynm)
                error('no PSFy found')
            end
            locsall.PSFxpix=locsall.PSFxnm/p.cam_pixelsize_nm;
            locsall.PSFypix=locsall.PSFynm/p.cam_pixelsize_nm;
            locsall.zfnm=locsall.frame*p.dz/1000;
            
            initaxis(p.resultstabgroup,'found beads');
            numfiles=max(locsall.filenumber);
            maxd=p.cam_pixelsize_nm*2;
            locsall.beadnum=zeros(size(locsall.xnm));
            for k=1:numfiles
                indf=locsall.filenumber==k;
                beadlocs(k)=getBeadLocs(locsall.xnm(indf),locsall.ynm(indf),p);
                [beadnum,numlocs]=associatelocs(beadlocs(k).x,beadlocs(k).y,locsall.xnm(indf),locsall.ynm(indf),maxd);
                beadn=beadnum>0;
                locsall.beadnum(indf&beadn)=beadnum(beadn)+max(locsall.beadnum);
            end
            
            %make bead structure
            
%             [ztrue,I0]=getTrueZ(locsall,p);
             initaxis(p.resultstabgroup,'get true z')
            for k=max(locsall.beadnum):-1:1
                thisbead=(locsall.beadnum==k);
                bead(k).PSFxpix=double(locsall.PSFxpix(thisbead));
                bead(k).PSFypix=double(locsall.PSFypix(thisbead));
                bead(k).phot=double(locsall.phot(thisbead));
                bead(k).xpos=median(double(locsall.xnm(thisbead)));
                bead(k).ypos=median(double(locsall.ynm(thisbead)));
                bead(k).filenumber=double(locsall.filenumber(find(thisbead,1)));
                bead(k).numlocs=sum(thisbead);
                bead(k).frame=double(locsall.frame(thisbead));
                bead(k).zfnm=double(locsall.zfnm(thisbead));
                %derived properties
                bead(k).ztrue=stackas2z(bead(k).PSFxpix,bead(k).PSFypix,bead(k).zfnm,bead(k).phot,p.showresults);
                
                indz0=find(bead(k).zfnm>bead(k).ztrue,1);
                rangez0=indz0-15:indz0+15;rangez0(rangez0<1)=1;rangez0(rangez0>length(bead(k).zfnm))=length(bead(k).zfnm);
                bead(k).minSx=min(bead(k).PSFxpix(rangez0));
                bead(k).minSy=min(bead(k).PSFypix(rangez0));
%                 ztrue(k);
                bead(k).I0=max(bead(k).phot);
            end
            
            % if on glass: correct position based on average, not taking
            % into account position. This might create errors. change
            % later?
            if p.beaddistribution.Value==1 %glass
                zf1=robustMean([bead([bead(:).filenumber]==1).ztrue]);
                for k=2:max([bead(:).filenumber])
                        bi=find([bead(:).filenumber]==k);
                        dzh=robustMean([bead(bi).ztrue])-zf1;
                        for b=1:length(bi)
                            bead(b).ztrue=bead(b).ztrue-dzh;
                            bead(b).zfnm=bead(b).zfnm-dzh;
                        end
                end
            end
            
            %sort beads according to X,Y
            
            Xrange=[ p.Xmin:p.Xd:p.Xmax].*p.cam_pixelsize_nm;
            Yrange=[p.Ymin:p.Xd:p.Ymax].*p.cam_pixelsize_nm;
            Zrange=p.Zmin:p.Zd:p.Zmax;
            
            if length(Xrange)==1
                Xrange(2)=p.Xmax*p.cam_pixelsize_nm;
            end
            if length(Yrange)==1
                Yrange(2)=p.Ymax*p.cam_pixelsize_nm;
            end
            if length(Zrange)==1
                Zrange(2)=p.Zmax*p.cam_pixelsize_nm;
            end
            
            for X=1:length(Xrange)-1
                for Y=1:length(Yrange)-1
                    indh=[bead(:).xpos]>Xrange(X)&[bead(:).xpos]<Xrange(X+1)&[bead(:).ypos]>Yrange(Y)&[bead(:).ypos]<Yrange(Y+1);
                    if sum(indh)==0
                        continue
                    end
                    if p.beaddistribution.Value==1
                        Z=1;
                        SXY(X,Y,Z).spline=getcleanspline(bead(indh));
                    else %get gel curves
                        splinez=getcleangel(bead(indh),p);
                        for Z=1:length(splinez)
                            SXY(X,Y,Z).spline=splinez(Z);
                        end
                        SXY(X,Y,Z).bead=bead(indh);
                    end 
                end
            end
            sp=[SXY(:).spline];
            zt=0:0.1:2;figure(88);hold off;for k=1:numel(sp), plot(zt,sp(k).x(zt),zt,sp(k).y(zt));hold on;end

            %if in gel: make Z.dependent curves
            
            % calculate average curves 
            
            %remove some curves based on intensity, deviation from average
            %Sxy
            
            %recalculate average curves
            
            %calculate LUT
%             calculate sx.^2-sy.^2 
            
%             zglass=min(ztrue)
%             if p.check_glassposb
%                 zglass=p.cal_glasspos*p.dz/1000;
%             else
%                 zs=sort(ztrue(~isnan(ztrue)));
%                 zglass=mean(zs(1:3));
%             end
%             ztrue=ztrue-zglass;
%             title(['zglass = ' num2str(zglass)])
%             
%             B0=double(~p.B0);
%             
%             frame=(p.cal_fstart+p.cal_fstop)/2;
%             [sx,sy,zcorr]=getCalibrationCurve(locsall,ztrue,frame,p);
%             initaxis(p.resultstabgroup,'calibration curve fit');
% 
%             fitp=fitCalibrationCurve(sx./p.cam_pixelsize_nm,sy/p.cam_pixelsize_nm,zcorr,frame*p.dz/1000-zglass,B0,p.cal_zrange/2);          
%             obj.outforfit=real(fitp([2 4 5 6 7 8 1 3]));
%             drawnow
% %             obj.cal3D.zglass=zglass;
%             allframes=p.cal_fstart:p.cal_deltaf:p.cal_fstop;
%             
% 
%             initaxis(p.resultstabgroup,'sx2-sy2');
%             obj.outsx2sy2=fitsx2sy2(sx/p.cam_pixelsize_nm,sy/p.cam_pixelsize_nm,zcorr,frame*p.dz/1000-zglass,p.cal_zrange/2);
%             
% %             recgui.initaxis(p.resultstabgroup,'sx-sy');
% %             plot(zcorr,sx-sy,'.')
%         
%             for k=1:length(allframes)
%                     initaxis(p.resultstabgroup,['c' num2str(k)]);
%                 obj.cal3D(k).zglass=zglass;
%                 frame=allframes(k);
%                 [sx,sy,zcorr]=getCalibrationCurve(locsall,ztrue,frame,p);
% %                 recgui.initaxis(p.resultstabgroup,'calibration curve fit');
%                 fitp=fitCalibrationCurve(sx./p.cam_pixelsize_nm,sy/p.cam_pixelsize_nm,zcorr,frame*p.dz/1000-zglass,B0,p.cal_zrange/2);          
%                 obj.cal3D(k).fit=real(fitp([2 4 5 6 7 8 1 3]));
%                 obj.cal3D(k).midpoint=-fitp(9);
%                 obj.cal3D(k).sx2sy2=fitsx2sy2(sx./p.cam_pixelsize_nm,sy./p.cam_pixelsize_nm,zcorr,frame*p.dz/1000-zglass,p.cal_zrange/2);
%                 
%                 drawnow
%                 
%             end
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end
function s=getcleanspline(b)

%filter based on Sx, Sy and phot
minSx=[b(:).minSx];
minSy=[b(:).minSy];
[mSx,smSx]=robustMean(minSx);
[mSy,smSy]=robustMean(minSy);
phot=[b(:).I0];
[mp,sp]=robustMean(phot);
indgood=minSx<mSx+2*smSx & minSy<mSy+2*smSy & phot<mp+2*sp;

if sum(indgood)==0
    indgood=true(size(minSx));
end

for k=1:length(b)
    b(k).zfnm=b(k).zfnm-b(k).ztrue;
end

z=vertcat(b(indgood).zfnm);
splinex=getspline(vertcat(b(indgood).PSFxpix),z,1./vertcat(b(indgood).PSFxpix));
spliney=getspline(vertcat(b(indgood).PSFypix),z,1./vertcat(b(indgood).PSFypix));

%how good does it fit with spline?
figure(99)
hold off
for k=length(b):-1:1
    w=(b(k).phot).^2;
%     w=1;
    errh=(b(k).PSFxpix-splinex(b(k).zfnm)).^2.*w+(b(k).PSFypix-spliney(b(k).zfnm)).^2.*w;
    err(k)=sqrt(sum(errh)/sum(w));
    plot(b(k).zfnm,b(k).PSFxpix,'k.')
    hold on
    plot(b(k).zfnm,b(k).PSFypix,'k.')
%     plot(b(k).zfnm,w/max(w))
    
end

%  
[em,es]=robustMean(err);
indgood2=indgood&err<em+2*es;
if sum(indgood2)==0
    indgood2=indgood;
end
zg=vertcat(b(indgood2).zfnm);
sxg=vertcat(b(indgood2).PSFxpix);
syg=vertcat(b(indgood2).PSFypix);

splinex2=getspline(sxg,zg,1./(sxg-splinex(zg)).^2);
spliney2=getspline(syg,zg,1./(syg-spliney(zg)).^2);

zt=0:0.01:2;
for k=1:length(b)
    if indgood(k)
    plot(b(k).zfnm,b(k).PSFxpix,'.')
    
    plot(b(k).zfnm,b(k).PSFypix,'.')
    end
end

plot(zt,splinex2(zt),'k')
plot(zt,spliney2(zt),'k')
plot(zt,splinex(zt),'r')
plot(zt,spliney(zt),'r')
drawnow
s.x=splinex2;s.y=spliney2;
end
function spline=getspline(S,z,w)
spline=fit(z,S,'smoothingspline','Weights',w); %weigh smaller more
% 
% [zs,zind]=sort(z);Ss=S(zind);
% ds=Ss-spline(zs);
% spline2=fit(zs,Ss,'smoothingspline','Weights',1./ds.^2); %weigh smaller more
end

% 
% function savecalfile(a,b,obj)
% fn=[obj.locData.files.file(1).name(1:end-4) '_3DAcal.mat'];
% [f,p]=uiputfile(fn);
% if f
%     outforfit=obj.outforfit;
%     cal3D=obj.cal3D;
%     outsx2sy2=obj.outsx2sy2;
%     save([p f],'outforfit','cal3D','outsx2sy2')
% end
% end
% 
% function fitpsx=fitsx2sy2(sx,sy,zcorr,framez,zrange)
% 
% indf=abs(zcorr-framez)<zrange;
% if sum(indf)>5
% fitpsx=fit(sx(indf).^2-sy(indf).^2,zcorr(indf),'poly3','Robust','LAR');
% 
% % fitpsx=polyfit(zcorr,sx.^2-sy.^2,4);
% plot(sx.^2-sy.^2,zcorr,'r.')
% hold on
% plot(sx(indf).^2-sy(indf).^2,zcorr(indf),'b.')
% 
% sxsort=sort(sx.^2-sy.^2);
% zsort=feval(fitpsx,sxsort);
% 
% plot(sxsort,zsort,'k')
% % plot(zcorr,polyval(fitpsx,zcorr),'.')
% hold off
% xlim([-6 6])
% else
%     fitpsx=zeros(2,1);
% end
% end

% function fitp=fitCalibrationCurve(sx,sy,z,framez,B0,zrange,startp)
% if isempty(sx)
%     fitp=zeros(9,1);
% else
% if nargin<7
%     startp=[    0.3    1.0    1.0000  0   0        0         0  0.307   -framez];
% end
% indf=abs(z-framez)<zrange;
% fitp=lsqnonlin(@sbothfromsigmaerr,startp,[],[],[],[z(indf) z(indf)],[sx(indf) sy(indf)],B0);
% fitp=real(fitp);
% 
% subplot('Position',[0.05,0.65,.9,.3])
% plot(z,sx,'c.',z,sy,'m.')
% hold on;
% plot(z(indf),sx(indf),'b.',z(indf),sy(indf),'r.')
% sxf=sigmafromz(fitp([1 2 4 6 8 9]),z,B0);
% plot(z,sxf,'k.')
% fpy=fitp([1 3 5 7 8 9]);
% fpy(5)=-fpy(5);
% syf=sigmafromz(fpy,z,B0);
% plot(z,syf,'k.')
% hold off
% title(num2str(fitp(:)',2))
% ylim([0 max(max(sxf),max(syf))])
% 
% subplot('Position',[0.05,0.45,.9,.15])
% plot(z(indf),sx(indf)-sxf(indf),'b.')
% hold on
% plot(z(indf),0*sy(indf),'k.')
% hold off
% ylim([-1 1]*.7)
% 
% subplot('Position',[0.05,0.25,.9,.15])
% plot(z(indf),sy(indf)-syf(indf),'r.')
% hold on
% plot(z(indf),0*sy(indf),'k.')
% hold off
% ylim([-1 1]*.7)
% 
% subplot('Position',[0.05,0.05,.9,.15])
% end
% end
function [sx,sy,zcorr]=getCalibrationCurve(locs,ztrue,frame,p)
frameregion=p.cal_framerange/2;        

df=locs.frame-frame;

indf=abs(df)<=frameregion;
indbead=locs.beadnum>0;

sx=double(locs.PSFxnm(indf&indbead));
sy=double(locs.PSFynm(indf&indbead));
fn=locs.frame(indf&indbead);
beadnum=locs.beadnum(indf&indbead);
dff=df(indf&indbead);

zbead=ztrue(beadnum);
zcorr=zbead-dff*p.dz/1000;
% 
%  recgui.initaxis(p.resultstabgroup,'calibration curve')
% plot(zcorr,sx,'.',zcorr,sy,'.')

indg=~isnan(zcorr);
zcorr=zcorr(indg);sx=sx(indg);sy=sy(indg);
end

% function [z,I0]=getTrueZ(locs,p)
% mmax=max(locs.beadnum);
% z=zeros(mmax,1);
% I0=zeros(mmax,1);
%  initaxis(p.resultstabgroup,'get true z')
% 
% for k=1:mmax
%     ind=locs.beadnum==k;
% %     zf=locs.frame(ind)*p.dz/1000;
%     [zas,zn]=stackas2z(locs.PSFxpix(ind),locs.PSFypix(ind),locs.zfnm(ind),locs.phot(ind),p.showresults);
%     z(k)=zas;
%     I0(k)=max(locs.phot(ind));
% %     waitforbuttonpress
% end
%  w=warning('off');
%  warning(w);
% end

function calibrateAstig3D(locs,p)
global zt sxf syf
sxt=double(locs.PSFxnm)/p.cam_pixelsize_nm;syt=double(locs.PSFynm)/p.cam_pixelsize_nm;framet=double(locs.frame);
zt=framet*p.dz/1000;

B0=double(~p.B0);


recgui.initaxis(p.resultstabgroup,'select range');
plot(framet,sxt,'ro')
hold on
plot(framet,syt,'bo')
hold off
title('select range (two points)')
[indi,y]=ginput(2);

% indi=round(indum/p.dz*1000);

%rangef=round(max(1,min(indi))):round(min(max(indi),length(sxt)));
rangef=round(max(min(indi))):round(min(max(indi)));

range=find(framet>=rangef(1),1,'first'):find(framet<=rangef(end),1,'last');

sx=sxt(range);sy=syt(range);frame=framet(range);
z=frame*p.dz/1000;



[~,ix]=min(sx);
[~,iy]=min(sy);
% 
midpreal=(z(ix)+z(iy))/2;


indframez0=find(p.framez0<=framet,1,'first');
midp=zt(indframez0);

z=z-midp;zt=zt-midp;
% recgui.initaxis(p.resultstabgroup,'fitted points');

plot(zt,sxt,'r.')
hold on
plot(zt,syt,'b.')
plot(z,sx,'ro')
plot(z,sy,'bo')

mpf=midpreal-midp;
% parx= [d sx0 sy0 Ax Ay Bx By g mp]
startp=[    0.3    1.0    1.0000  0   0        0         0  0.307   -mpf];

% fitp=lsqcurvefit(@sbothfromsigma,startp,[z z],[sx sy])
fitp=lsqnonlin(@sbothfromsigmaerr,startp,[],[],[],[z z],[sx sy],B0);

% fitfunction=@(par,x)sbothfromsigmafit(par,x);
% fit2=fit([z; z],[sx; sy],fitfunction);

sxf=sigmafromz(fitp([1 2 4 6 8 9]),zt,B0);
plot(zt,sxf,'k')
fpy=fitp([1 3 5 7 8 9]);
fpy(5)=-fpy(5);
syf=sigmafromz(fpy,zt,B0);
plot(zt,syf,'k')
hold off

axis tight

% 
% figure(2)
% plot(sx,sy,'.')
% hold on
% plot(sxf,syf,'r')
% plot(syf,sxf,'c')
% hold off


%zpar=[sigma0x,Ax,Ay,Bx,By,gamma,d,sigma0y)
outforfit=real(fitp([2 4 5 6 7 8 1 3]));
button = questdlg('Is the fit good?','3D astigmatism bead calibration','Refit','Save','Cancel','Refit') ;
switch button
    case 'Refit'
        calibrateAstig3D(locs,p)
    case 'Save'
        fn=[p.fileout(1:end-4) '_3DAcal.mat']
        [f,p]=uiputfile(fn);
        if f
            save([p f],'outforfit')
        end
end
end
function s=sigmafromz(par,z,B0)
global gamma
par=real(par);
% parx= [d sx0 Ax Bx g mp]
s0=par(2);d=par(1);A=par(3);B=par(4)*B0;g=par(5);mp=par(6);

% s=s0*sqrt(1+(z-g+mp).^2/d^2);
s=s0*sqrt(1+(z-g+mp).^2/d^2+A*(z-g+mp).^3/d^3+B*(z-g+mp).^4/d^4);
s=real(s);
end


function s=sbothfromsigma(par,z,B0)
% parx= [d sx0 sy0 Ax Ay Bx By g mp]
par=real(par);
px=par([1 2 4 6 8 9]);
py=par([1 3 5 7 8 9]);
 py(5)=-py(5);
zh=z(:,1);
s=real([sigmafromz(px,zh,B0) sigmafromz(py,zh,B0)]);
end

function err=sbothfromsigmaerr(par,z,sx,B0)
sf=sbothfromsigma(par,z,B0);
err=sf-sx;

% stderr=std(err);
err=err./sqrt(abs(err));
end

function beaddistribution_callback(a,b,obj)
gel={'Zt','Zmin','Zd','Zmax','framerangecombinet','framerangecombine'};
glass={'ztoframe','ztoframet'};
p=obj.getSingleGuiParameter('beaddistribution');
switch p.Value 
    case 1%glass
        tgel='off';tglass='on';
    case 2 %gel
        tgel='on';tglass='off';
end
for k=1:length(gel)
    obj.guihandles.(gel{k}).Visible=tgel;
end
for k=1:length(glass)
    obj.guihandles.(glass{k}).Visible=tglass;
end
end

function lut_callback(a,b,obj)
end
function pard=guidef(obj)
pard.dzt.object=struct('Style','text','String','dz (nm)'); 
pard.dzt.position=[1,1];
pard.dzt.Width=.5;
pard.dz.object=struct('Style','edit','String','50'); 
pard.dz.position=[1,1.5];
pard.dz.Width=.5;

pard.beaddistribution.object=struct('String',{{'Beads on Glass','Beads in Gel'}},'Style','popupmenu','Callback',{{@beaddistribution_callback,obj}});
pard.beaddistribution.position=[2,1];
pard.beaddistribution.Width=1.5;

tp=2.6;tmin=3.1;td=3.4;tmax=3.7;
w=0.3;
wp=0.5;
pard.tt1.object=struct('String','grid','Style','text');
pard.tt1.position=[2,tp];
pard.tt1.Width=wp;
pard.mint.object=struct('String','min','Style','text');
pard.mint.position=[2,tmin];
pard.mint.Width=w;
pard.dxt.object=struct('String','delta','Style','text');
pard.dxt.position=[2,td];
pard.dxt.Width=w;
pard.maxt.object=struct('String','max','Style','text');
pard.maxt.position=[2,tmax];
pard.maxt.Width=w;

pard.Xt.object=struct('String','X (pix)','Style','text');
pard.Xt.position=[3,tp];
pard.Xt.Width=wp;

pard.Xmin.object=struct('String','0','Style','edit');
pard.Xmin.position=[3,tmin];
pard.Xmin.Width=w;
pard.Xd.object=struct('String','128','Style','edit');
pard.Xd.position=[3.5,td];
pard.Xd.Width=w;
pard.Xmax.object=struct('String','512','Style','edit');
pard.Xmax.position=[3,tmax];
pard.Xmax.Width=w;

pard.Yt.object=struct('String','Y (pix)','Style','text');
pard.Yt.position=[4,tp];
pard.Yt.Width=wp;

pard.Ymin.object=struct('String','0','Style','edit');
pard.Ymin.position=[4,tmin];
pard.Ymin.Width=w;
pard.Ymax.object=struct('String','512','Style','edit');
pard.Ymax.position=[4,tmax];
pard.Ymax.Width=w;

pard.Zt.object=struct('String','Z (nm)','Style','text');
pard.Zt.position=[5,tp];
pard.Zt.Width=wp;
pard.Zmin.object=struct('String','0','Style','edit');
pard.Zmin.position=[5,tmin];
pard.Zmin.Width=w;
pard.Zd.object=struct('String','1000','Style','edit');
pard.Zd.position=[5,td];
pard.Zd.Width=w;
pard.Zmax.object=struct('String','3000','Style','edit');
pard.Zmax.position=[5,tmax];
pard.Zmax.Width=w;

pard.lutb.object=struct('String','save Z(Sx,Sy,X,Y,Z)','Style','pushbutton','Callback',{{@lut_callback,obj}});
pard.lutb.position=[6,1];
pard.lutb.Width=1.5;

pard.St.object=struct('String','S (pix)','Style','text');
pard.St.position=[6,tp];
pard.St.Width=wp;
pard.Smin.object=struct('String','0','Style','edit');
pard.Smin.position=[6,tmin];
pard.Smin.Width=w;
pard.Sd.object=struct('String','.02','Style','edit');
pard.Sd.position=[6,td];
pard.Sd.Width=w;
pard.Smax.object=struct('String','3','Style','edit');
pard.Smax.position=[6,tmax];
pard.Smax.Width=w;

pard.framerangecombinet.object=struct('String','frames to combine','Style','text');
pard.framerangecombinet.position=[3,1];
pard.framerangecombinet.Width=1;
pard.framerangecombine.object=struct('String','10','Style','edit');
pard.framerangecombine.position=[3,2];
pard.framerangecombine.Width=.5;

pard.ztoframet.object=struct('String','set z to frame','Style','checkbox');
pard.ztoframet.position=[3,1];
pard.ztoframet.Width=1;
pard.ztoframe.object=struct('String','21','Style','edit');
pard.ztoframe.position=[3,2];
pard.ztoframe.Width=.5;

% pard.savebutton.object=struct('String','save','Style','pushbutton');
% pard.savebutton.position=[7,3];
pard.inputParameters={'cam_pixelsize_nm'};
pard.plugininfo.type='ProcessorPlugin';


end
