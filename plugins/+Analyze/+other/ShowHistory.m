classdef ShowHistory<interfaces.DialogProcessor
    methods
        function obj=ShowHistory(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;
        end
        
        function out=run(obj,p)
            hist=obj.locData.history;
            texta={};
            for k=1:length(hist)
                
                phist=hist{k};
                txt=struct2txt(phist,'');
                if isfield(phist,'name')
                    name=phist.name;
                else 
                    name='';
                end
                texta{end+1}=['Module' num2str(k) ': ' name];
                texta(end+1:end+length(txt))=txt;
                texta{end+1}='';
            end
            
            listdlg('ListString',texta,'ListSize',[800,800]);
            out=texta;
           
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        

    end
end




function pard=guidef
pard.plugininfo.name='Show History';
pard.plugininfo.type='ProcessorPlugin';
end