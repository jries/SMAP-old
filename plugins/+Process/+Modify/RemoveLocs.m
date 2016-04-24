classdef RemoveLocs<interfaces.DialogProcessor
    methods
        function obj=RemoveLocs(varargin)   
            obj@interfaces.DialogProcessor(varargin{:});  
        end
        
        function out=run(obj,p)
            obj.setPar('undoModule','RemoveLocs');
            notify(obj.P,'backup4undo');
            [~,indroi]=obj.locData.getloc('xnm','position','roi');
            obj.locData.removelocs(indroi);
            obj.locData.regroup;   
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function updateGui(obj,event,data)
            if ~isempty(obj.locData.loc)
            fn=fieldnames(obj.locData.loc);
            obj.guihandles.fieldselect.String=fn;
            end
        end       

    end
end




function pard=guidef
pard.textb.object=struct('String','remove locs in ROI','Style','text');
pard.textb.position=[2,1];
pard.textb.Width=3;
pard.textb.object.TooltipString='';

end