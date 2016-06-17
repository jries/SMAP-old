classdef calibrateAstig<interfaces.DialogProcessor
    properties
        zold
    end
    methods
        function obj=calibrateAstig(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.inputParameters={'cam_pixelsize_nm'};
             obj.zold.changed=0;
             obj.showresults=true;

        end
        function out=run(obj,p)
            out=[];
            if isfield(obj.locData.loc,'gradient3Dellipticity')
                locs=obj.locData.getloc({'frame','gradient3Dellipticity'},'layer',1,'position','roi','removeFilter',{'PSFxnm','PSFynm'});
                p.fileout=obj.locData.files(1).file.name;
                out=calibrategradient3D(locs,p)  ;  
                obj.setPar('fit_gradient3Dellipticity',(out))
            elseif isfield(obj.locData.loc,'PSFynm')
                locs=obj.locData.getloc({'frame','PSFxnm','PSFynm'},'layer',1,'position','roi','removeFilter',{'PSFxnm','PSFynm'});

                p.fileout=obj.locData.files(1).file.name;
                calibrateAstig3D(locs,p)         
            else
                error('no 3D data found')
            end
            
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        

    end
    methods(Static)
    end
end


function ttxt=calibrategradient3D(locs,p)  
global zt 
framet=double(locs.frame);
zt=framet*p.dz/1000;
eps=locs.gradient3Dellipticity;
epsl=log(eps);
B0=double(~p.B0);
initaxis(p.resultstabgroup,'select range');
plot(framet,epsl,'ro')
% hold on
% plot(framet,syt,'bo')
% hold off
title('select range (two points)')
[indi,y]=ginput(2);

rangef=round(max(min(indi))):round(min(max(indi)));
range=find(framet>=rangef(1),1,'first'):find(framet<=rangef(end),1,'last');
epslr=epsl(range);
frame=framet(range);
z=frame*p.dz/1000;

indz0=find(epslr>0,1,'first');
midp=z(indz0);

% 
% midpreal=(z(ix)+z(iy))/2;
% 
% indframez0=find(p.framez0<=framet,1,'first');
% midp=zt(indframez0);
% 
z=z-midp;zt=zt-midp;

plot(zt,epsl,'r.')
hold on

plot(z,epslr,'ro')



% fitp=lsqnonlin(@sbothfromsigmaerr,startp,[],[],[],[z z],[sx sy],B0);

E=[ones(length(epslr),1) epslr];
b=E\z;
bplot(1)=-b(1)/b(2);
bplot(2)=1/b(2);
plot(zt,bplot(1)+bplot(2)*zt)
ttxt=sprintf([num2str(b(1)) '\t' num2str(b(2))]);
title(ttxt)
%zpar=[sigma0x,Ax,Ay,Bx,By,gamma,d,sigma0y)
outforfit=b;
button = questdlg('Is the fit good?','3D astigmatism bead calibration','Refit','Save','Cancel','Refit') ;
switch button
    case 'Refit'
        calibrategradient3D(locs,p)
    case 'Save'
        fn=[p.fileout(1:end-4) '_3DAGcal.mat'];
        [f,p]=uiputfile(fn);
        if f
            save([p f],'outforfit')
        end
end
end


function calibrateAstig3D(locs,p)
global zt sxf syf
sxt=double(locs.PSFxnm)/p.cam_pixelsize_nm;syt=double(locs.PSFynm)/p.cam_pixelsize_nm;framet=double(locs.frame);
% if p.reverse
%    p.dz=-p.dz;
% end

zt=framet*p.dz/1000;



B0=double(~p.B0);
initaxis(p.resultstabgroup,'select range');
plot(framet,sxt,'ro')
hold on
plot(framet,syt,'bo')
hold off
title('select range (two points)')
[indi,y]=ginput(2);

rangef=round(max(min(indi))):round(min(max(indi)));
range=find(framet>=rangef(1),1,'first'):find(framet<=rangef(end),1,'last');

sx=sxt(range);sy=syt(range);frame=framet(range);
z=frame*p.dz/1000;

[~,ix]=min(sx);
[~,iy]=min(sy);

midpreal=(z(ix)+z(iy))/2;

indframez0=find(p.framez0<=framet,1,'first');
midp=zt(indframez0);

z=z-midp;zt=zt-midp;

plot(zt,sxt,'r.')
hold on
plot(zt,syt,'b.')
plot(z,sx,'ro')
plot(z,sy,'bo')

mpf=midpreal-midp;
% parx= [d sx0 sy0 Ax Ay Bx By g mp]
startp=[    0.3    1.0    1.0000  0   0        0         0  0.307   -mpf];

fitp=lsqnonlin(@sbothfromsigmaerr,startp,[],[],[],[z z],[sx sy],B0);


sxf=sigmafromz(fitp([1 2 4 6 8 9]),zt,B0);
plot(zt,sxf,'k')
fpy=fitp([1 3 5 7 8 9]);
fpy(5)=-fpy(5);
syf=sigmafromz(fpy,zt,B0);
plot(zt,syf,'k')
hold off
axis tight
ylabel('PSFx,PSFy')

%zpar=[sigma0x,Ax,Ay,Bx,By,gamma,d,sigma0y)
outforfit=real(fitp([2 4 5 6 7 8 1 3]));
% initaxis(p.resultstabgroup,'sx^2-sy^2');
outsx2sy2=fitsx2sy2(sx,sy,z,1);
% outsx2sy2=obj.outsx2sy2;
button = questdlg('Is the fit good?','3D astigmatism bead calibration','Refit','Save','Cancel','Refit') ;
switch button
    case 'Refit'
        calibrateAstig3D(locs,p)
    case 'Save'
        fn=[p.fileout(1:end-4) '_3DAcal.mat'];
        [f,p]=uiputfile(fn);
        if f
            save([p f],'outforfit','outsx2sy2')
        end
end
end
function s=sigmafromz(par,z,B0)
% global gamma
% parx= [d sx0 Ax Bx g mp]
s0=par(2);d=par(1);A=par(3);B=par(4)*B0;g=par(5);mp=par(6);

s=s0*sqrt(1+(z-g+mp).^2/d^2+A*(z-g+mp).^3/d^3+B*(z-g+mp).^4/d^4);
end


function s=sbothfromsigma(par,z,B0)
% parx= [d sx0 sy0 Ax Ay Bx By g mp]
px=par([1 2 4 6 8 9]);
py=par([1 3 5 7 8 9]);
 py(5)=-py(5);
zh=z(:,1);
s=[sigmafromz(px,zh,B0) sigmafromz(py,zh,B0)];
end

function err=sbothfromsigmaerr(par,z,sx,B0)
sf=sbothfromsigma(par,z,B0);
err=sf-sx;

% stderr=std(err);
err=err./sqrt(abs(err));
end

function fitpsx=fitsx2sy2(sx,sy,z,zrange)

indf=abs(z)<zrange;
if sum(indf)>5
fitpsx=fit(sx(indf).^2-sy(indf).^2,z(indf),'poly8','Robust','LAR');
yyaxis right
% hold on
% % fitpsx=polyfit(zcorr,sx.^2-sy.^2,4);
% plot(z,sx.^2-sy.^2,'r.')
hold on
plot(z(indf),sx(indf).^2-sy(indf).^2,'.','Color',[0 0.5 0])

sxsort=sort(sx.^2-sy.^2);
zsort=feval(fitpsx,sxsort);

plot(zsort,sxsort,'k')
ylabel('sx^2-sy^2')
% plot(zcorr,polyval(fitpsx,zcorr),'.')
hold off
% xlim([-6 6])
else
    fitpsx=zeros(2,1);
end
end


function pard=guidef
% pard.d3_color.object=struct('Style','checkbox','String','render in color');
% pard.d3_color.position=[2,1];
% pard.mode3D.object=struct('String',{{'astigmatic PSFx/y for MLE','astigmatic gradient fit','bi-plane'}},'Style','popupmenu');
% pard.mode3D.position=[2,3];
% pard.mode3D.Width=3;

pard.text2.object=struct('String','calibrate 3D ','Style','text');
pard.text2.position=[1,1];
pard.text2.Width=1.5;
pard.text3.object=struct('String','dz (nm)','Style','text');
pard.text3.position=[3,1];
% 
pard.dz.object=struct('Style','edit','String','50'); 
pard.dz.position=[3,2.5];


pard.text4.object=struct('String','frame of Zero position','Style','text');
pard.text4.position=[4,1];
pard.text4.Width=1.5;

pard.framez0.object=struct('Style','edit','String','21'); 
pard.framez0.position=[4,2.5];

pard.B0.object=struct('String','set B = 0','Style','checkbox','Value',0);
pard.B0.position=[5,1];

% pard.reverse.object=struct('String','reverse z axis','Style','checkbox','Value',0);
% pard.reverse.position=[4,3.5];

pard.plugininfo.name='calibrate 3D';
pard.plugininfo.type='ProcessorPlugin';
end
