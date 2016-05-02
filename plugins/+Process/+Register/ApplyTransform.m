classdef ApplyTransform<interfaces.DialogProcessor
    % ApplyTransform applies transformation to localizations or associated
    % Tif images
    methods
        function obj=ApplyTransform(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'cam_pixelsize_nm'};
            obj.showresults=false;
            obj.history=true;
        end
        
        function out=run(obj,p) 
            out=[];
            obj.setPar('undoModule','ApplyTransform');
            notify(obj.P,'backup4undo');
            load(p.Tfile)
            if ~exist('transformation','var')
                out.error='selected transformation file does not have a valid transformation';
                return
            end
                
            if p.allfiles
                filenumbers=1:length(obj.locData.files.file);
            else
                filenumbers=p.dataselect.Value;
            end
            
            for f=1:length(filenumbers)
                file=obj.locData.files.file(filenumbers(f));

                if p.transformwhat.Value==1||p.transformwhat.Value==3
                    tiffs=file.tif;
                    for k=1:length(tiffs)

                        if ~isempty(tiffs(k).info)&&isempty(strfind(tiffs(k).info.name,'_T.tif'))
                            p.roitiff=tiffs(k).info.roi;  
                            tiffn=tiffs(k);
                            tiffn.image=apply_transform_image(tiffs(k).image,transformation,p);
                            tiffn.info.name=strrep(tiffs(k).info.name, '.tif' ,'_T.tif');
                            obj.locData.files.file(filenumbers(f)).tif(end+1)=tiffn;
                        end
                    end
                end
                 if p.transformwhat.Value==1||p.transformwhat.Value==2
                    loc=apply_transform_locs(obj.locData.loc,transformation,file,p);
                    obj.locData.loc=copyfields(obj.locData.loc,loc);

                 end
            end
            if p.transformwhat.Value==1||p.transformwhat.Value==2
                obj.locData.regroup;
            end
        end
        function pard=guidef(obj)
            pard=guidef(obj);
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


function pard=guidef(obj)
pard.texta.object=struct('String','dataset:','Style','text');
pard.texta.position=[1,1];

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[1,1.7];
pard.dataselect.Width=2;

pard.allfiles.object=struct('Style','checkbox','String','transform all files','Value',0);
pard.allfiles.position=[1,3.7];
pard.allfiles.Width=1.3;

pard.textb.object=struct('String','transform:','Style','text');
pard.textb.position=[2,1];

pard.datapart.object=struct('Style','popupmenu','String',{{'all (T->R)','all (R->T)','reference','target'}});
pard.datapart.position=[2,1.7];

pard.textc.object=struct('String','locs/tifs:','Style','text');
pard.textc.position=[3,1];

pard.setchannel.object=struct('Style','checkbox','String','set channel (ref=1,target=2)', 'Value',0);
pard.setchannel.position=[3,3];
pard.setchannel.Width=2;

pard.transformwhat.object=struct('Style','popupmenu','String',{{'both','locs','tiffs'}});
pard.transformwhat.position=[3,1.7];

pard.Tfile.object=struct('Style','edit','String','settings/temp/temp_T.mat');
pard.Tfile.position=[7,1];
pard.Tfile.Width=3;

pard.loadbutton.object=struct('Style','pushbutton','String','load T','Callback',@obj.loadbutton);
pard.loadbutton.position=[7,4];

pard.syncParameters={{'filelist_short','dataselect',{'String'}}};

pard.plugininfo.description='ApplyTransform applies transformation to localizations or associated Tif images';
pard.plugininfo.name='Apply Transformation';
end