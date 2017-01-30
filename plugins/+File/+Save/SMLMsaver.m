classdef SMLMsaver<interfaces.DialogProcessor
    methods
        function obj=SMLMsaver(varargin)  
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'filelist_long','mainfile','mainGui','numberOfLayers','sr_layerson'};
        end
        
        function out=save(obj,p)
            obj.status('save sml file')
            lastfile=obj.getPar('lastSMLFile');
            
            if isempty(lastfile)
                fn=p.filelist_long.selection;
                ind=strfind(fn,'_sml');
                if isempty(ind)
                    ind=strfind(fn,'_fitpos');
                end
                if isempty(ind)
                    ind=length(fn)-3;
                end
                of=[fn(1:ind-1) '_sml.mat'];
            else
                of=lastfile;
            end
              
            
            [f,path]=uiputfile(of);
            if f
                if isempty(strfind(f,'_sml'))
                    f(end-3:end)=[];
                    f=[f '_sml.mat'];
                end   
            par=obj.getAllParameters;
            
            par.saveroi=par.savevisible;
%             savelocData=obj.locdata.copy;
            savesml(obj.locData,[path f],par)
            
            end
            obj.status('save done')
          
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
pard.savevisible.object=struct('Style','checkbox','Visible','on','String','only save visible','Value',0);
pard.savevisible.position=[1,1];
pard.savevisible.Width=4;
pard.savevisible.object.TooltipString='save only those filtered localizations that have been used to render the image';

pard.savefile.object=struct('Style','checkbox','Visible','on','String','Only: ','Value',0);
pard.savefile.position=[2,1];
pard.savefile.Width=1;
pard.savefile.object.TooltipString='save only selected file';

pard.dataselect.object=struct('Style','popupmenu','Visible','on','String',{{'empty'}});
pard.dataselect.position=[2,2];
pard.dataselect.Width=3;
pard.dataselect.object.TooltipString='save only selected file';

pard.syncParameters={{'filelist_short','dataselect',{'String'}}};
pard.plugininfo.type='SaverPlugin';

end