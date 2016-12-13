classdef CenterSites<interfaces.SEEvaluationProcessor
    properties
        
    end
    methods
        function obj=CenterSites(varargin)        
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
pard.centermode.object=struct('Style','popupmenu','String',{{'median','mean','mask','fitring'}});
pard.centermode.position=[1,1];
pard.centermode.Width=2;

pard.iterationst.object=struct('Style','text','String','iterations: ');
pard.iterationst.position=[2,1];
pard.iterationst.Width=1;

pard.iterations.object=struct('Style','edit','String',3);
pard.iterations.position=[2,2];
pard.iterations.Width=1;

pard.cxy.object=struct('Style','checkbox','String','correct x,y');
pard.cxy.position=[1,3];
pard.cxy.Width=2;

pard.cz.object=struct('Style','checkbox','String','correct z');
pard.cz.position=[2,3];
pard.cz.Width=2;

pard.plugininfo.type='ROI_Evaluate';
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end


function out=runintern(obj,p)
for k=1:p.iterations
locs=obj.getLocs({'xnm','ynm','znm','locprecznm','locprecnm'},'layer',1,'size',p.se_siteroi);
if obj.display
ax=obj.setoutput('images');
end
switch p.centermode.selection
    case 'median'
        dx=median(locs.xnm)-obj.site.pos(1);
        dy=median(locs.ynm)-obj.site.pos(2);
        if ~isempty(locs.znm)
        dz=median(locs.znm)-obj.site.pos(3);
        else
            dz=0;
        end
    case 'mean'
        dx=mean(locs.xnm)-obj.site.pos(1);
        dy=mean(locs.ynm)-obj.site.pos(2);
        if ~isempty(locs.znm)
        dz=mean(locs.znm)-obj.site.pos(3);
        else
            dz=0;
        end
    case 'mask'
        img=obj.site.image.image;
        ns=round(-p.se_sitefov/2:p.se_sitepixelsize:p.se_sitefov/2);
        [Xs,Ys]=meshgrid(ns,ns);
        mask=Xs.^2+Ys.^2<(p.se_siteroi/2)^2;
        
        pixrec=p.se_sitepixelsize;
        rim=20;
        diameterNPC=110;
        n=-1*diameterNPC:pixrec:1*diameterNPC;
        [X,Y]=meshgrid(n,n);

        h=exp(-abs(X.^2+Y.^2-(diameterNPC/2)^2)/4/rim^2);
        h=h/sum(h(:));
        imsq=sqrt(sum(img,3));
        imsq(~mask)=0;
        imf=filter2(h,imsq);
        [~,mind]=max(imf(:));
        [my,mx]=ind2sub(size(imf),mind);
        dx=ns(mx)/2;
        dy=ns(my)/2;   
        dz=0;
    case 'fitring'
        [dx,dy]=fitposring(locs.xnm,locs.ynm,R);

end 

obj.site.pos=obj.site.pos+[dx*p.cxy dy*p.cxy dz*p.cz];
end
if obj.display
imagesc(ax,[imsq/max(imsq(:)) imf/max(imf(:))])
title(ax,[dx dy])
end
out=[];
end

