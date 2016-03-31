classdef clusterMIiSR<interfaces.DialogProcessor
    methods
        function obj=clusterMIiSR(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_layerson','linewidth_roi','sr_roihandle','znm_min','znm_max','numberOfLayers'};
        end
        
        function out=run(obj,p)           
           disp('not implemented')
            out=0;
        end
        function pard=pardef(obj)
            pard=pardef;
        end
    end
end



function pard=pardef
pard.plugininfo.name='Cluster with MIiSR';
end
