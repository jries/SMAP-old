classdef PlotLocsPreview<interfaces.WorkflowModule;
    properties
%          pixelsize      
    end
    methods
       function obj=PlotLocsPreview(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
%              obj.setInputChannels(1,'frame');
       end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
            
        end
        function pard=guidef(obj)
            pard.plugininfo.type='WorkflowModule'; 
        end
        function output=run(obj,data,p)
            output=[];
            if obj.getPar('loc_preview') 
                locs=data.data;
                if ~isempty(locs)&&~isempty(locs.xpix)
                    ax=findobj(obj.getPar('loc_outputfig').Children,'Type','Axes');
                    ax.NextPlot='add';
                    plot(locs.xpix,locs.ypix,'ko','Parent',ax,'MarkerSize',12)
                end
            end 


        end

    end
end
