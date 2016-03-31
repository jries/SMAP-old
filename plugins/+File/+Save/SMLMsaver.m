classdef SMLMsaver<interfaces.DialogProcessor
    methods
        function obj=SMLMsaver(varargin)  
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'filelist_long','mainfile','mainGui','numberOfLayers','sr_layerson'};
        end
        
        function out=save(obj,p)
            obj.status('save sml file')
            fn=p.filelist_long.selection;
            ind=strfind(fn,'_sml');
            if isempty(ind)
                ind=strfind(fn,'_fitpos');
            end
            if isempty(ind)
                ind=length(fn)-3;
            end
            of=[fn(1:ind-1) '_sml.mat'];
            
            [f,path]=uiputfile(of);
            if f
                if isempty(strfind(f,'_sml'))
                    f(end-3:end)=[];
                    f=[f '_sml.mat'];
                end   
            par=obj.getAllParameters;
            
            par.saveroi=par.savevisible;
            savesml(obj.locData,[path f],par)
            
            end
            obj.status('save done')
          
        end
        function pard=pardef(obj)
            pard=pardef;
        end
        function run(obj,p)
            obj.save(p)
        end
    end
end




function pard=pardef
pard.savevisible.object=struct('Style','checkbox','Visible','on','String','only save visible','Value',0);
pard.savevisible.position=[1,1];
pard.savevisible.Width=4;
end