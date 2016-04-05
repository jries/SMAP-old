classdef ShowHistory<interfaces.DialogProcessor
    methods
        function obj=ShowHistory(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;
        end
        
        function out=run(obj,p)
            phist.p=obj.locData.history{1};
            txt=struct2txt(phist,'p');
            listdlg('ListString',txt,'ListSize',[800,800]);
            
            out=txt;
           
        end
        function pard=pardef(obj)
            pard=pardef;
        end
        

    end
end




function pard=pardef
pard.plugininfo.name='Show History';

end