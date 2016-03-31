classdef calibrateAstigDeep<interfaces.DialogProcessor
    properties
        outforfit
        cal3D
        outsx2sy2
    end
    methods
        function obj=calibrateAstigDeep(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.inputParameters={'cam_pixelsize_nm'};
        end
        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
            set(obj.guihandles.savebutton,'Callback',{@savecalfile,obj})
        end
        function out=run(obj,p)
            out=[];
            locsall=obj.locData.getloc({'frame','xnm','ynm','PSFxnm','PSFynm','phot'},'position','roi');
%             locsall=obj.locData.getloc('frame','xnm','ynm','PSFxnm','PSFynm','filenumber','phot');
            if isempty(locsall.PSFynm)
                error('no PSFy found')
            end
            initaxis(p.resultstabgroup,'found beads');
            beadlocs=getBeadLocs(locsall,p);
            maxd=p.cam_pixelsize_nm*5;
            [locsall.beadnum,numlocs]=associatelocs(beadlocs.x,beadlocs.y,locsall.xnm,locsall.ynm,maxd);
            
            ztrue=getTrueZ(locsall,p);
%             zglass=min(ztrue)
            if p.check_glasspos
                zglass=p.cal_glasspos*p.dz/1000;
            else
                zs=sort(ztrue(~isnan(ztrue)));
                zglass=mean(zs(1:3));
            end
            ztrue=ztrue-zglass;
            title(['zglass = ' num2str(zglass)])
            
            B0=double(~p.B0);
            
            frame=(p.cal_fstart+p.cal_fstop)/2;
            [sx,sy,zcorr]=getCalibrationCurve(locsall,ztrue,frame,p);
            initaxis(p.resultstabgroup,'calibration curve fit');

            fitp=fitCalibrationCurve(sx./p.cam_pixelsize_nm,sy/p.cam_pixelsize_nm,zcorr,frame*p.dz/1000-zglass,B0,p.cal_zrange/2);          
            obj.outforfit=real(fitp([2 4 5 6 7 8 1 3]));
            drawnow
%             obj.cal3D.zglass=zglass;
            allframes=p.cal_fstart:p.cal_deltaf:p.cal_fstop;
            

            initaxis(p.resultstabgroup,'sx2-sy2');
            obj.outsx2sy2=fitsx2sy2(sx/p.cam_pixelsize_nm,sy/p.cam_pixelsize_nm,zcorr,frame*p.dz/1000-zglass,p.cal_zrange/2);
            
%             recgui.initaxis(p.resultstabgroup,'sx-sy');
%             plot(zcorr,sx-sy,'.')
        
            for k=1:length(allframes)
                    initaxis(p.resultstabgroup,['c' num2str(k)]);
                obj.cal3D(k).zglass=zglass;
                frame=allframes(k);
                [sx,sy,zcorr]=getCalibrationCurve(locsall,ztrue,frame,p);
%                 recgui.initaxis(p.resultstabgroup,'calibration curve fit');
                fitp=fitCalibrationCurve(sx./p.cam_pixelsize_nm,sy/p.cam_pixelsize_nm,zcorr,frame*p.dz/1000-zglass,B0,p.cal_zrange/2);          
                obj.cal3D(k).fit=real(fitp([2 4 5 6 7 8 1 3]));
                obj.cal3D(k).midpoint=-fitp(9);
                obj.cal3D(k).sx2sy2=fitsx2sy2(sx./p.cam_pixelsize_nm,sy./p.cam_pixelsize_nm,zcorr,frame*p.dz/1000-zglass,p.cal_zrange/2);
                
                drawnow
                
            end
            
        end
        function pard=pardef(obj)
            pard=pardef;
        end
    end
end

function savecalfile(a,b,obj)
fn=[obj.locData.files.file(1).name(1:end-4) '_3DAcal.mat'];
[f,p]=uiputfile(fn);
if f
    outforfit=obj.outforfit;
    cal3D=obj.cal3D;
    outsx2sy2=obj.outsx2sy2;
    save([p f],'outforfit','cal3D','outsx2sy2')
end
end

function fitpsx=fitsx2sy2(sx,sy,zcorr,framez,zrange)

indf=abs(zcorr-framez)<zrange;
if sum(indf)>5
fitpsx=fit(sx(indf).^2-sy(indf).^2,zcorr(indf),'poly3','Robust','LAR');

% fitpsx=polyfit(zcorr,sx.^2-sy.^2,4);
plot(sx.^2-sy.^2,zcorr,'r.')
hold on
plot(sx(indf).^2-sy(indf).^2,zcorr(indf),'b.')

sxsort=sort(sx.^2-sy.^2);
zsort=feval(fitpsx,sxsort);

plot(sxsort,zsort,'k')
% plot(zcorr,polyval(fitpsx,zcorr),'.')
hold off
xlim([-6 6])
else
    fitpsx=zeros(2,1);
end
end

function fitp=fitCalibrationCurve(sx,sy,z,framez,B0,zrange,startp)
if isempty(sx)
    fitp=zeros(9,1);
else
if nargin<7
    startp=[    0.3    1.0    1.0000  0   0        0         0  0.307   -framez];
end
indf=abs(z-framez)<zrange;
fitp=lsqnonlin(@sbothfromsigmaerr,startp,[],[],[],[z(indf) z(indf)],[sx(indf) sy(indf)],B0);
fitp=real(fitp);

subplot('Position',[0.05,0.65,.9,.3])
plot(z,sx,'c.',z,sy,'m.')
hold on;
plot(z(indf),sx(indf),'b.',z(indf),sy(indf),'r.')
sxf=sigmafromz(fitp([1 2 4 6 8 9]),z,B0);
plot(z,sxf,'k.')
fpy=fitp([1 3 5 7 8 9]);
fpy(5)=-fpy(5);
syf=sigmafromz(fpy,z,B0);
plot(z,syf,'k.')
hold off
title(num2str(fitp(:)',2))
ylim([0 max(max(sxf),max(syf))])

subplot('Position',[0.05,0.45,.9,.15])
plot(z(indf),sx(indf)-sxf(indf),'b.')
hold on
plot(z(indf),0*sy(indf),'k.')
hold off
ylim([-1 1]*.7)

subplot('Position',[0.05,0.25,.9,.15])
plot(z(indf),sy(indf)-syf(indf),'r.')
hold on
plot(z(indf),0*sy(indf),'k.')
hold off
ylim([-1 1]*.7)

subplot('Position',[0.05,0.05,.9,.15])
end
end
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

function z=getTrueZ(locs,p)
mmax=max(locs.beadnum);
z=zeros(mmax,1);
 initaxis(p.resultstabgroup,'get true z')

for k=1:mmax
    ind=locs.beadnum==k;
    zf=locs.frame(ind)*p.dz/1000;
    [zas,zn]=stackas2z(locs.PSFxnm(ind),locs.PSFynm(ind),zf,locs.phot(ind),p.showresults);
    z(k)=zas;
end
 w=warning('off');
 warning(w);
end

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

function pard=pardef
% pard.d3_color.object=struct('Style','checkbox','String','render in color');
% pard.d3_color.position=[2,1];

pard.text2.object=struct('String','calculate calibration for','Style','text');
pard.text2.position=[1,1];
pard.text2.Width=2;
pard.text4.object=struct('String','first frame','Style','text');
pard.text4.position=[2,1];
pard.text4.Width=1.5;
pard.cal_fstart.object=struct('Style','edit','String','30'); 
pard.cal_fstart.position=[2,2.5];
pard.cal_fstart.Width=.5;

pard.text5.object=struct('String','last frame','Style','text');
pard.text5.position=[3,1];
pard.text5.Width=1.5;
pard.cal_fstop.object=struct('Style','edit','String','30'); 
pard.cal_fstop.position=[3,2.5];
pard.cal_fstop.Width=.5;

pard.text6.object=struct('String','calculate for every F frames','Style','text');
pard.text6.position=[4,1];
pard.text6.Width=1.5;
pard.cal_deltaf.object=struct('Style','edit','String','5'); 
pard.cal_deltaf.position=[4,2.5];
pard.cal_deltaf.Width=.5;

pard.text9.object=struct('String','data info','Style','text');
pard.text9.position=[5,1];
pard.text9.Width=2;
pard.check_glasspos.object=struct('String','set frame focused on coverslip to: ','Style','checkbox');
pard.check_glasspos.position=[6,1];
pard.check_glasspos.Width=1.5;
pard.cal_glasspos.object=struct('Style','edit','String','5'); 
pard.cal_glasspos.position=[6,2.5];
pard.cal_glasspos.Width=.5;

pard.text3.object=struct('String','dz (nm)','Style','text');
pard.text3.position=[7,1];
pard.text3.Width=1.5;
pard.dz.object=struct('Style','edit','String','50'); 
pard.dz.position=[7,2.5];
pard.dz.Width=.5;

pard.text10.object=struct('String','fit parameters','Style','text');
pard.text10.position=[1,3];
pard.text10.Width=2;

pard.text11.object=struct('String','z range (um) for fit','Style','text');
pard.text11.position=[2,3];
pard.text11.Width=1.5;
pard.cal_zrange.object=struct('Style','edit','String','1'); 
pard.cal_zrange.position=[2,4.5];
pard.cal_zrange.Width=.5;

pard.text12.object=struct('String','range of frames to evaluate','Style','text');
pard.text12.position=[3,3];
pard.text12.Width=1.5;
pard.cal_framerange.object=struct('Style','edit','String','10'); 
pard.cal_framerange.position=[3,4.5];
pard.cal_framerange.Width=.5;

pard.B0.object=struct('String','set B = 0','Style','checkbox');
pard.B0.position=[4,3];

pard.savebutton.object=struct('String','save','Style','pushbutton');
pard.savebutton.position=[6,3];

% 
% pard.d3_color.object=struct('Style','checkbox','String','render in color');
% pard.d3_color.position=[5,1];
% 
% 
% pard.pixauto.object=struct('Style','checkbox','String','pixelsize in z auto','Value',1);
% pard.pixauto.position=[4,1];
% pard.pixrecset.object=struct('Style','edit','String','5'); 
% pard.pixrecset.position=[4,2];
% 
% pard.text4.object=struct('String','filter  \sigma in z (planes)','Style','text');
% pard.text4.position=[6,1];
% pard.filterz.object=struct('Style','edit','String','1'); 
% pard.filterz.position=[6,2];
% 
% 
% pard.preview.object=struct('Style','checkbox','String','preview');
% pard.preview.position=[2,4];
% 
% pard.N0_fit.object=struct('String','N0','Style','radiobutton');
% pard.N0_fit.position=[2,2];
% 
% pard.N0_v.object=struct('String','10','Style','edit');
% pard.N0_v.position=[2,3];
% pard.N0_v.isnumeric=1;
% 
% 
% pard.pmature_fit.object=struct('String','p mature','Style','radiobutton');
% pard.pmature_fit.position=[3,2];
% 
% pard.pmature_v.object=struct('String','.5','Style','edit');
% pard.pmature_v.position=[3,3];
% pard.pmature_v.isnumeric=1;
% 
% 


end
