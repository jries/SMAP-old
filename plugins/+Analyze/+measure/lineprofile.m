classdef lineprofile<interfaces.DialogProcessor
    methods
        function obj=lineprofile(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_layerson','linewidth_roi','znm_min','znm_max'};
            obj.showresults=true;
%             obj.maxresultstabs=6;
        end
        
        function out=run(obj,p)           
            make_lineprofiles(obj.locData,p)
            out.clipboard={'test',1};
        end
        function pard=pardef(obj)
            pard=pardef;
        end
    end
end




function pard=pardef
pard.text1.object=struct('String','binwidth (nm):','Style','text');
pard.text1.position=[1,1];

pard.binwidth.object=struct('String','2','Style','edit');
pard.binwidth.position=[1,2];

pard.text2.object=struct('String','fitmodel:','Style','text');
pard.text2.position=[2,1];

pard.fitmodel.object=struct('String','Gauss|Flat|Disk|Ring|Distance','Style','popupmenu');
pard.fitmodel.position=[2,2];

pard.restrictsigma.object=struct('String','sigma=<locp>','Style','checkbox');
pard.restrictsigma.position=[2,3];

pard.linelengthcheck.object=struct('String','length (nm)','Style','checkbox');
pard.linelengthcheck.position=[3,1];

pard.linelength.object=struct('String','250','Style','edit');
pard.linelength.position=[3,2];

pard.plugininfo.name='lineprofile';
end