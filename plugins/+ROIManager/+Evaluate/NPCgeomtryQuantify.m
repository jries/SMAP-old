classdef NPCgeomtryQuantify<interfaces.SEEvaluationProcessor
    properties
        savedevals
    end
    methods
        function obj=NPCgeomtryQuantify(varargin)        
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
pard.Rt.object=struct('Style','text','String','R (nm):');
pard.Rt.position=[1,1];
pard.Rt.Width=1;

pard.R.object=struct('Style','edit','String','50');
pard.R.position=[1,2];
pard.R.Width=1;

pard.dRt.object=struct('Style','text','String','dR:');
pard.dRt.position=[1,3];
pard.dRt.Width=1;

pard.dR.object=struct('Style','edit','String','20');
pard.dR.position=[1,4];
pard.dR.Width=1;

pard.plugininfo.type='ROI_Evaluate';
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
end


function out=runintern(obj,p)
locs=obj.getLocs({'xnm','ynm','znm','locprecznm','locprecnm'},'layer',1,'size',p.se_siteroi);
dz=5;
z=-400:dz:400;
hz=hist(locs.znm,z);
ac=myxcorr(hz,hz);
if obj.display
ax1=obj.setoutput('profile');
plot(ax1,z,hz)
ax2=obj.setoutput('correlation');
plot(ax2,z,ac)
end
out.ac=ac;
out.dz=dz;
out.hist=hz;
out.z=z;
end

