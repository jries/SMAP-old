classdef TifLoader<interfaces.WorkflowModule
    properties     
        imloader
        framestop
        framestart
        numberOfFrames
        timerfitstart
        edgesize
    end
    methods
        function obj=TifLoader(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1;
            obj.isstartmodule=true;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.inputParameters={'loc_subtractbg','loc_blocksize_frames'};            
            obj.guihandles.loadtifbutton.Callback={@loadtif_callback,obj};
            obj.addSynchronization('filelist_localize',obj.guihandles.tiffile,'String',{@loadtif_ext,obj});
        end
        function prerun(obj,p)
            if ~exist(p.tiffile,'file')
                obj.status('TifLoader: localization file not found')
                 error('TifLoader: localization file not found')
            end
%             obj.imloader=imageLoader(p.tiffile);     
            obj.imloader=imageloaderAll(p.tiffile);     
            obj.imloader.onlineAnalysis=p.onlineanalysis;
            obj.imloader.waittime=p.onlineanalysiswaittime;
            obj.imloader.setImageNumber(p.framestart-1);
            obj.numberOfFrames=obj.imloader.metadata.numberOfFrames;
%             obj.numberOfFrames=obj.imloader.info.numberOfFrames;

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
            tall=0;
            while ~isempty(image)&&imloader.currentImageNumber<=obj.framestop&&~SMAP_stopnow
                
                datout=interfaces.WorkflowData;
%                 datout.set(image);


                if p.padedges
                    dr=obj.edgesize;
                    imgo=zeros(size(image)+2*dr,'like',image);
                    bg=myquantilefast(image(:),0.2,1000);
                    imgo=imgo+bg;
                    imgo(dr+1:end-dr,dr+1:end-dr)=image;
                    image=imgo;
                end
                datout.data=image;

                datout.frame=imloader.currentImageNumber;
                datout.ID=id;
                id=id+1;
                obj.output(datout)
                th=tic;
                image=imloader.readNext;
      
                
                
                %display
                if mod(datout.frame,10)==0
                    obj.setPar('loc_currentframe',struct('frame',datout.frame,'image',image));
                    if p.onlineanalysis
                        
                        totalf=imloader.metadata.numberOfFrames-obj.framestart;
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
                tall=tall+toc(th);
            end
            
            obj.setPar('tiffloader_loadingtime',tall);
            dateof=interfaces.WorkflowData;
            dateof.frame=imloader.currentImageNumber+1;
            dateof.ID=id;
            dateof.eof=true;
            obj.output(dateof)
            disp('fitting done')
        end
        function addFile(obj,file)
%         if 1
%         else
                
            obj.imloader=imageloaderAll(file);
%             if ~isempty(obj.imloader.info.metafile)
            if ~isempty(obj.imloader.metadata.allmetadata) 
             obj.setPar('loc_metadatafile',obj.imloader.metadata.allmetadata.metafile);
            end
%             end
% obj.imloader.metadata
            p=obj.getGuiParameters;
            fileinf=obj.imloader.metadata;
            if p.padedges
%                 locsettings=obj.getPar('loc_cameraSettings');
                roisize=13;
                dr=ceil((roisize-1)/2);
                fileinf.roi(1:2)=fileinf.roi(1:2)-dr;
                fileinf.roi(3:4)=fileinf.roi(3:4)+2*dr;
%                 fileinf.Width=fileinf.Width+2*dr;
%                 fileinf.Height=fileinf.Height+2*dr;
%                 obj.setPar('loc_cameraSettings',locsettings);
                obj.edgesize=dr;
            end
            obj.setPar('loc_fileinfo',fileinf);
% 
            obj.guihandles.tiffile.String=obj.imloader.file;
%              obj.setPar('loc_newfile',true);
             p=obj.getAllParameters;
             if p.onlineanalysis
                 obj.guihandles.framestop.String='inf';
             else
                obj.guihandles.framestop.String=int2str(obj.imloader.metadata.numberOfFrames);
             end
%         end
% 
        end
        function loadtif_ext(obj)
            p=obj.getAllParameters;
            obj.addFile(p.tiffile)
        end
        function loadedges(obj,a,b)
            tiffile=obj.getSingleGuiParameter('tiffile');
            if exist(tiffile,'file')
                obj.addFile(tiffile);
            end
        end
    end
end

function loadtif_callback(a,b,obj)
p=obj.getGuiParameters;
fe=bfGetFileExtensions;
[f,path]=uigetfile(fe,'select camera images',[fileparts(p.tiffile) filesep '*.tif']);
if f
    obj.addFile([path f]);
end  
end

function pard=guidef(obj)
pard.text.object=struct('Style','text','String','Image source file:');
pard.text.position=[1,1];
pard.text.Width=1.5;
pard.text.Optional=true;

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
pard.onlineanalysis.Optional=true;
pard.onlineanalysiswaittime.object=struct('Style','edit','String','5');
pard.onlineanalysiswaittime.position=[3,4];
pard.onlineanalysiswaittime.Width=0.5;
pard.onlineanalysiswaittime.TooltipString='Waiting time for online analysis.';
pard.onlineanalysiswaittime.Optional=true;

pard.textf.object=struct('Style','text','String','Frame range');
pard.textf.position=[4.2,1.25];
pard.textf.Width=0.75;
pard.textf.Optional=true;
pard.framestart.object=struct('Style','edit','String','1');
pard.framestart.position=[4.2,2];
pard.framestart.Width=0.7;
pard.framestart.Optional=true;
pard.framestop.object=struct('Style','edit','String','1000000');
pard.framestop.position=[4.2,2.7];
pard.framestop.Width=0.7;
pard.framestop.Optional=true;

pard.padedges.object=struct('Style','checkbox','String','Pad edges','Value',0,'Callback',@obj.loadedges);
pard.padedges.position=[4.2,3.5];
pard.padedges.Width=1;
pard.padedges.Optional=true;
pard.padedges.TooltipString='Pad edges with minimum values to allow detection of localizations close to edges. Usually not necessary.';
% pard.locdata_empty.object=struct('Style','checkbox','String','Empty localizations','Value',1);
% pard.locdata_empty.position=[4.2,3.5];
% pard.locdata_empty.Width=1.5;
% pard.locdata_empty.TooltipString=sprintf('Empty localization data before fitting. \n Important if post-processing (eg drift correction) is perfromed as part of workflow');
pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Loads tiff files. Can load single tiff files or tiff stacks, also while they are written. It also tries to locate the metadata.txt file from micromanger and passes it on to the camera converter.';
end