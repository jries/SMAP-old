classdef NPCLabelingQuantify<interfaces.SEEvaluationProcessor
    properties
        savedevals
    end
    methods
        function obj=NPCLabelingQuantify(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function makeGui(obj)
            makeGui@interfaces.SEEvaluationProcessor(obj);
%             obj.guihandles.saveimagesb.Callback={@saveimagesb_callback,obj};
        end
        function out=run(obj,p)
            try
            out=runintern(obj,p);
            catch err
                err
                out=[];
            end
         
        end
        
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
% pard.layer.object=struct('Style','popupmenu','String','layer1|layer2');
% pard.layer.position=[1,1];
% pard.layer.Width=2;

% pard.dualchannel.object=struct('Style','checkbox','String','dual channel','Value',0);
% pard.dualchannel.position=[2,1];
% pard.dualchannel.Width=2;
% 
% 
% pard.sqrtfit.object=struct('Style','checkbox','String','fit sqrt(img)','Value',1);
% pard.sqrtfit.position=[3,1];
% pard.sqrtfit.Width=2;
% 
% 
% pard.t2.object=struct('Style','text','String','gaussfac for imfit (l1, l2)');
% pard.t2.position=[4,1];
% pard.t2.Width=3;
% 
% pard.gaussfac_imfit.object=struct('Style','edit','String','1, 1');
% pard.gaussfac_imfit.position=[4,4];
% pard.gaussfac_imfit.Width=1;
% 
% pard.fit_sigma.object=struct('Style','checkbox','String','Fit sigma of ring','Value',0);
% pard.fit_sigma.position=[5,1];
% pard.fit_sigma.Width=2;

pard.plugininfo.type='ROI_Evaluate';
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end


function out=runintern(obj,p)

locs=obj.getLocs({'xnm','ynm','xnm_gt','ynm_gt','locprecnm'},'layer',1,'size',p.se_siteroi);
xm=locs.xnm-obj.site.pos(1);
ym=locs.ynm-obj.site.pos(2);


step=2*pi/8;

[tha,rhoa]=cart2pol(xm,ym);
R=50;
dR=20;

minlp=step*R*.4;
inr=rhoa>R-dR&rhoa<R+dR;
inr=inr&locs.locprecnm<minlp;
th=tha(inr);rho=rhoa(inr);

[th_gt,rho_gt]=cart2pol(locs.xnm_gt(inr)-obj.site.pos(1),locs.ynm_gt(inr)-obj.site.pos(2));
locptheta=locs.locprecnm(inr)./rho;
% dth=mod(th,step);
mdt=cyclicaverage(th,step,1./locptheta);
if mdt>pi/16
    mdt=mdt-step;
end
[ths,inds]=sort(th);
ths(end+1)=ths(1)+2*pi;
dth=diff(ths);
lp=locptheta(inds);
lp(end+1)=lp(1);
dstep=sqrt(lp(1:end-1).^2+lp(2:end).^2); %uncertaintty
mindistance=1;
fullspace=ceil(dth-mindistance); fullspace( fullspace<0)=0;
numfound=8-sum(fullspace);
%when considering maximum gap, also look at next corners
%correlation: when uncertain and moves in one direction, distance to othr
%gets larger
% figure(188);histogram(th,-pi:pi/8:pi)
figure(188);
% histogram(dth/step,30)
subplot(2,1,1)
plot(dth,dstep,'+')
title((mean(mdt))/pi*180)
subplot(2,1,2)
plot(locs.xnm(inr),locs.ynm(inr),'x',locs.xnm_gt(inr),locs.ynm_gt(inr),'o')
axis equal
tg=sort(th_gt);tg(end+1)=tg(1)+2*pi;
numlocs=sum(diff(tg)>step*.7);
title(['corners: ' int2str(numlocs) ', found: ' int2str(numfound)])


end

