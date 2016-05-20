classdef TifLoader<interfaces.WorkflowModule
    properties     
        imloader
        framestop
        framestart
        numberOfFrames
        timerfitstart
    end
    methods
        function obj=TifLoader(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1;
            obj.isstartmodule=true;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.inputParameters={'loc_subtractbg','loc_blocksize_frames'};            
            obj.guihandles.loadtifbutton.Callback={@loadtif_callback,obj};
            obj.addSynchronization('filelist_localize',obj.guihandles.tiffile,'String',{@loadtif_ext,obj});
        end
        function prerun(obj,p)
%             p=obj.getAllParameters;
            
            if ~exist(p.tiffile,'file')
                obj.status('TifLoader: localization file not found')
                 error('TifLoader: localization file not found')
            elseif p.locdata_empty
                obj.locData.clear;
            end
            obj.imloader=imageLoader(p.tiffile);     
            obj.imloader.onlineAnalysis=p.onlineanalysis;
            obj.imloader.waittime=p.onlineanalysiswaittime;
            obj.imloader.setImageNumber(p.framestart-1);
            obj.numberOfFrames=obj.imloader.info.numberOfFrames;
%             obj.imloader.currentImageNumber=p.framestart-1;
            if p.onlineanalysis
                obj.framestop=inf;
            else
                obj.framestop=p.framestop;
            end
            obj.framestart=p.framestart;
            if obj.getPar('loc_preview')
                previewframe=obj.getPar('loc_previewframe');               
                if p.loc_subtractbg
                dt=p.loc_blocksize_frames;
                frameload=max(1,previewframe-floor(dt/2));
                obj.framestart=frameload;
                obj.imloader.setImageNumber(frameload-1);
%                 obj.imloader.currentImageNumber=frameload-1;
                obj.framestop=frameload+dt-1;
                else
                    obj.imloader.setImageNumber(previewframe-1);
%                     obj.imloader.currentImageNumber=previewframe-1;
                    obj.framestop=previewframe;
                end               
            end
            obj.timerfitstart=tic;            
        end
        function run(obj,data,p)
            global SMAP_stopnow
            
            if nargin>1&& ~isempty(data)&&~isempty(data.data) %optional input channel
                file=data.data;
                obj.addFile(file)
            end
            id=1;
%             disp('run loader')
            imloader=obj.imloader;
            image=imloader.readNext;      
            
            while ~isempty(image)&&imloader.currentImageNumber<=obj.framestop&&~SMAP_stopnow
                
                datout=interfaces.WorkflowData;
%                 datout.set(image);
                datout.data=image;

                datout.frame=imloader.currentImageNumber;
                datout.ID=id;
                id=id+1;
                obj.output(datout)
                image=imloader.readNext;
                
                
                %display
                if mod(datout.frame,10)==0
                    obj.setPar('loc_currentframe',struct('frame',datout.frame,'image',image));
                    if p.onlineanalysis
                        
                        totalf=imloader.info.numberOfFrames-obj.framestart;
                        elapsed=toc(obj.timerfitstart);
                        totaltime=elapsed/(datout.frame-obj.framestart+1)*totalf;
                    else
                        totalf=min(obj.numberOfFrames,obj.framestop)-obj.framestart;
                        elapsed=toc(obj.timerfitstart);
                        totaltime=elapsed/(datout.frame-obj.framestart+1)*totalf;
                    end
                    statuss=['frame ' int2str(datout.frame-obj.framestart) ' of ' int2str(totalf) ...
                        '. Time: ' num2str(elapsed,'%4.0f') ' of ' num2str(totaltime,'%4.0f') 's'];
                    obj.status(statuss);
                end
            end
            dateof=interfaces.WorkflowData;
            dateof.frame=imloader.currentImageNumber+1;
            dateof.ID=id;
            dateof.eof=true;
            obj.output(dateof)
            disp('fitting done')
        end
        function addFile(obj,file)
            
            obj.imloader=imageLoader(file);
            if ~isempty(obj.imloader.info.metafile)
             obj.setPar('loc_metadatafile',obj.imloader.info.metafile);
            end
            obj.setPar('loc_fileinfo',obj.imloader.info);
% 
            obj.guihandles.tiffile.String=obj.imloader.file;
%              obj.setPar('loc_newfile',true);
             p=obj.getAllParameters;
             if p.onlineanalysis
                 obj.guihandles.framestop.String='inf';
             else
                obj.guihandles.framestop.String=int2str(obj.imloader.info.numberOfFrames);
             end
% 
        end
        function loadtif_ext(obj)
            p=obj.getAllParameters;
            obj.addFile(p.tiffile)
        end
    end
end

function loadtif_callback(a,b,obj)
p=obj.getGuiParameters;
[f,path]=uigetfile([fileparts(p.tiffile) filesep '*.tif']);
if f
    obj.addFile([path f]);
end  
end

function pard=guidef
pard.text.object=struct('Style','text','String','Image source file:');
pard.text.position=[1,1];
pard.text.Width=1.5;


pard.loadtifbutton.object=struct('Style','pushbutton','String','load images','Visible','on');
pard.loadtifbutton.position=[3,1];
pard.loadtifbutton.TooltipString=sprintf('Open raw camera image tif files. \n Either single images in directory. \n Or multi-image Tiff stacks');

pard.tiffile.object=struct('Style','edit','String',' ','HorizontalAlignment','right');
pard.tiffile.position=[2,1];
pard.tiffile.Width=4;

pard.onlineanalysis.object=struct('Style','checkbox','String','Online analysis. Waittime (s):','Value',0);
pard.onlineanalysis.position=[3,2.25];
pard.onlineanalysis.Width=1.75;
pard.onlineanalysis.TooltipString='Fit during acquisition. If checked, max frames is ignored. Waits until no more images are written to file.';

pard.onlineanalysiswaittime.object=struct('Style','edit','String','5');
pard.onlineanalysiswaittime.position=[3,4];
pard.onlineanalysiswaittime.Width=0.5;
pard.onlineanalysiswaittime.TooltipString='Waiting time for online analysis.';


pard.textf.object=struct('Style','text','String','Frame range');
pard.textf.position=[4.2,1.25];
pard.textf.Width=0.75;
pard.framestart.object=struct('Style','edit','String','1');
pard.framestart.position=[4.2,2];
pard.framestart.Width=0.7;

pard.framestop.object=struct('Style','edit','String','1000000');
pard.framestop.position=[4.2,2.7];
pard.framestop.Width=0.7;

pard.locdata_empty.object=struct('Style','checkbox','String','Empty localizations','Value',1);
pard.locdata_empty.position=[4.2,3.5];
pard.locdata_empty.Width=1.5;
pard.locdata_empty.TooltipString=sprintf('Empty localization data before fitting. \n Important if post-processing (eg drift correction) is perfromed as part of workflow');
pard.plugininfo.type='WorkflowModule'; 
end