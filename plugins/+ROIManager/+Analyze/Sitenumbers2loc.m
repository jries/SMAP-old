classdef Sitenumbers2loc<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=Sitenumbers2loc(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
%             obj.inputParameters={'se_layerson'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)  
            out=[];
            sites=obj.SE.sites;
            
            ld=obj.locData.loc;
            sitenumbers=0*ld.xnm;
            cellnumbers=0*ld.xnm;
            for k=1:length(sites)
                [~,ind]=obj.locData.getloc('xnm','Position',sites(k));
                sitenumbers(ind)=sites(k).ID;
                cellnumbers(ind)=sites(k).info.cell;
            end
            obj.locData.addloc('sitenumbers',sitenumbers);
            obj.locData.addloc('cellnumbers',cellnumbers);
            obj.locData.regroup;
          
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.t1.object=struct('String','Write site and cell numbers to locData.loc','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=4;


pard.plugininfo.type='ROI_Analyze';


end