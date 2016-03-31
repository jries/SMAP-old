classdef ApplyTransform<interfaces.DialogProcessor
    methods
        function obj=ApplyTransform(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'cam_pixelsize_nm'};
        end
        
        function out=run(obj,p) 
            obj.setPar('undoModule','ApplyTransform');
            notify(obj.P,'backup4undo');
            load(p.Tfile)
            file=obj.locData.files.file(p.dataselect.Value);
            
            if p.transformwhat.Value==1||p.transformwhat.Value==3
            tiffs=file.tif;
            for k=1:length(tiffs)
                      
                if ~isempty(tiffs(k).info)&&isempty(strfind(tiffs(k).info.name,'_T.tif'))
                    p.roitiff=tiffs(k).info.roi;  
                    tiffn=tiffs(k);
                    tiffn.image=apply_transform_image(tiffs(k).image,transformation,p);
                    tiffn.info.name=strrep(tiffs(k).info.name, '.tif' ,'_T.tif');
                    obj.locData.files.file(p.dataselect.Value).tif(end+1)=tiffn;
                end
            end
            end
             if p.transformwhat.Value==1||p.transformwhat.Value==2
            loc=apply_transform_locs(obj.locData.loc,transformation,file,p);
            obj.locData.loc=copyfields(obj.locData.loc,loc);
            obj.locData.regroup;
             end
        end
        function pard=pardef(obj)
            pard=pardef(obj);
        end
        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
        end
        function initGui(obj)
            obj.addSynchronization('transformationfile',obj.guihandles.Tfile,'String');
        end
        function loadbutton(obj,a,b)
            fn=obj.guihandles.Tfile.String;
            [f,path]=uigetfile(fn,'Select transformation file _T.mat');
            if f
                obj.guihandles.Tfile.String=[path f];
                obj.setPar('transformationfile',[path f]);
            end      
        end  
    end
end


function pard=pardef(obj)
pard.texta.object=struct('String','dataset','Style','text');
pard.texta.position=[1,1];

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[2,1];

pard.datapart.object=struct('Style','popupmenu','String','all|reference|target');
pard.datapart.position=[3,1];

pard.setchannel.object=struct('Style','checkbox','String','set channel (ref=1,target=2)', 'Value',0);
pard.setchannel.position=[4,1];
pard.setchannel.Width=3;

pard.transformwhat.object=struct('Style','popupmenu','String','all|locs|tiffs');
pard.transformwhat.position=[5,1];

pard.Tfile.object=struct('Style','edit','String','settings/temp/temp_T.mat');
pard.Tfile.position=[8,1];
pard.Tfile.Width=3;

pard.loadbutton.object=struct('Style','pushbutton','String','load T','Callback',@obj.loadbutton);
pard.loadbutton.position=[8,4];

pard.syncParameters={{'filelist_short','dataselect',{'String'}}};

end