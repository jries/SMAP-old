classdef CheckCalibration<interfaces.DialogProcessor
    properties
        zold
    end
    methods
        function obj=CheckCalibration(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.inputParameters={'cam_pixelsize_nm'};
             obj.zold.changed=0;

        end

        function out=run(obj,p)
            locsROI=obj.locData.getloc({'frame','xnm','ynm','znm','phot'},'layer',1,'position','roi');
            locsall=obj.locData.getloc({'frame','xnm','ynm','znm','filenumber','phot'});
%             plotzvsframe(locsROI,p)
            beadlocs=getBeadLocs(locsROI,p);
            maxd=p.cam_pixelsize_nm*2;
            [locsall.beadnum,numlocs]=associatelocs(beadlocs.x,beadlocs.y,locsall.xnm,locsall.ynm,maxd);
            plotmanybeads(locsall,p)

            out=0;
        end
        
        function pard=pardef(obj)
            pard=pardef;
        end
        

    end
end

function plotmanybeads(locs,p)
mmax=max(locs.beadnum);
indf1=locs.filenumber==1;
indf2=locs.filenumber==2;
ax=initaxis(p.resultstabgroup,'zplots');
axis(ax)
for k=1:mmax
    indb=locs.beadnum==k;
    bf1=indb&indf1;
    bf2=indb&indf2;
    locsb1.frame=locs.frame(bf1);
    locsb1.znm=locs.znm(bf1);
    locsb1.phot=locs.phot(bf1);
    
    locsb2.frame=locs.frame(bf2);
    locsb2.znm=locs.znm(bf2);
    locsb2.phot=locs.phot(bf2);
    pf=plotzvsframe(locsb1,p);
    slope1(k)=pf(1);off1(k)=-pf(2)/pf(1);
    hold on
    pf=plotzvsframe(locsb2,p);
    slope2(k)=pf(1);off2(k)=-pf(2)/pf(1);
    hold off
    
    mx1(k)=mean(locs.xnm(bf1));
    my1(k)=mean(locs.ynm(bf1));
%     drawnow
%     waitforbuttonpress
end
initaxis(p.resultstabgroup,'slope vs z')
indb=abs(slope1)>1.3|abs(slope2)>2|abs(off1)>10000|abs(off2)>10000|abs(slope1)<.8;

plot(off1(~indb),slope1(~indb),'ro')
hold on
plot(off2(~indb),slope2(~indb),'bx')
hold off

initaxis(p.resultstabgroup,'slope in image')
scatter3(mx1(~indb),my1(~indb),slope1(~indb),[],slope1(~indb))
end


function pf=plotzvsframe(locs,p)

diffz=50;
refractiveIndexFactor=1.25;

z=locs.frame*diffz/refractiveIndexFactor;
plot(z,locs.znm,'o')
% 
% nd=7;
% dz=locs.znm;
% for k=1:nd
%     dz=diff(dz);
% end
 [~,ind]=max(locs.phot);
 
 hold on
plot(z(ind),locs.znm(ind),'kd')
 
 
 
%   [~,ind]=min(abs(dz))
nd=0;


dz1=diff(locs.znm);
i2=find((dz1(ind:end)-50)<-100,1,'first')+ind-2;
i1=find((dz1(1:ind-1)-50)<-100,1,'last')+1;

if isempty(i1)
    i1=1;
end
if isempty(i2)
    i2=length(locs.znm);
end
zlin=locs.znm(i1:i2);
flin=z(i1:i2);

pf=polyfit(flin,zlin,1);


plot(flin,zlin,'x')
plot(flin,polyval(pf,flin),'r')
hold off
title(['slope: ' num2str(pf(1))])
end

function pard=pardef
% pard.d3_color.object=struct('Style','checkbox','String','render in color');
% pard.d3_color.position=[2,1];

pard.text2.object=struct('String','check astigmatism calibration','Style','text');
pard.text2.position=[2,1];
% pard.text3.object=struct('String','dz (nm)','Style','text');
% pard.text3.position=[3,1];
% 
% pard.dz.object=struct('Style','edit','String','50'); 
% pard.dz.position=[3,2];
% 
% 
% pard.text4.object=struct('String','frame of Zero position','Style','text');
% pard.text4.position=[4,1];
% 
% pard.framez0.object=struct('Style','edit','String','21'); 
% pard.framez0.position=[4,2];
% % 
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
