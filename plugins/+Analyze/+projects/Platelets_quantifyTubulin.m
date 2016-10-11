classdef Platelets_quantifyTubulin<interfaces.DialogProcessor
    % LINEPROFILE Calculates profiles along a linear ROI and fits it with a model of choice
    methods
        function obj=Platelets_quantifyTubulin(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_layerson','linewidth_roi','znm_min','znm_max','sr_pixrec'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)           
            results=make_lineprofiles(obj.locData,p);
            out.clipboard=results;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function addmt_callback(a,b,obj)
p=obj.getAllParameters;
locs=obj.locData.getloc({'xnmline','ynmline','znm'},'layer',1,'position','roi');
ax=initaxis(obj.resultstabgroup,'profile');
plot(ax,locs.xnmline,locs.znm,'.');
lls=p.loclist.String;
if isempty(lls{1})
    newline=1;
else
    newline=length(lls)+1;
end
lls{newline}=num2str(length(locs.xnmline));
obj.setPar(struct('loclist',struct('String'),lls));
end



function pard=guidef(obj)
pard.addmt.object=struct('String','add single MT (locs/um)','Style','pushbutton','Callback',{{@addmt_callback,obj}});
pard.addmt.position=[1,1];
pard.addmt.Width=1.5;

pard.loclist.object=struct('String',{{''}},'Style','listbox');
pard.loclist.position=[4,1];
pard.loclist.Width=1.5;
pard.loclist.Height=3;



pard.plugininfo.name='Platelets_quantifyTubulin';
% pard.plugininfo.description=sprintf('Calculates profiles along a linear ROI and fits it with a model of choice. \n Flat: step function convolved with Gaussian (=Erf). \n Disk: Projection of a homogeneously filled disk, convolved with Gaussian. \n Ring: Projection of a ring, convolved with Gaussian. \n Distance: Two Gaussians in a distance d.');
pard.plugininfo.type='ProcessorPlugin';
end