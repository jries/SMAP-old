classdef saver_csv<interfaces.DialogProcessor
    methods
        function obj=saver_csv(varargin)  
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'filelist_long','mainfile','mainGui','numberOfLayers','sr_layerson'};
        end
        
        function out=save(obj,p)
            obj.status('save csv file')
            fn=p.filelist_long.selection;
            [path,file,ext]=fileparts(fn);
            of=[path filesep file  '.csv'];
            
            [f,path]=uiputfile(of);
            if f
               
            par=obj.getAllParameters;
            
            par.saveroi=par.savevisible;
            saveLocalizationsCSV(obj.locData,[path f],par.saveroi,par.numberOfLayers,par.sr_layerson);
            obj.status('save done')
            end
            
           
          
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function run(obj,p)
            obj.save(p)
        end
    end
end




function pard=guidef
% pard.plugininfo={'csv saver'};
pard.savevisible.object=struct('Style','checkbox','Visible','on','String','only save visible','Value',1);
pard.savevisible.position=[1,1];
pard.savevisible.Width=4;

pard.plugininfo.type='SaverPlugin';
end