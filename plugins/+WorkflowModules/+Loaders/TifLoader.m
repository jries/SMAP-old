classdef TifLoader<interfaces.WorkflowModule
    properties     
        loaders
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
             obj.loaders={'auto',@imageloaderAll;'MMStack',@imageloaderMM;'MMSingle',@imageloaderMMsingle;'OME',@imageloaderOME;'simple Tif',@imageloaderTifSimple};
            initGui@interfaces.WorkflowModule(obj);
            obj.inputParameters={'loc_subtractbg','loc_blocksize_frames'};            
            obj.guihandles.loadtifbutton.Callback={@loadtif_callback,obj};
            obj.addSynchronization('filelist_localize',obj.guihandles.tiffile,'String',{@loadtif_ext,obj});
            obj.guihandles.loaderclass.String=obj.loaders(:,1);
        end
        function prerun(obj,p)
            if ~exist(p.tiffile,'file')
                obj.status('TifLoader: localization file not found')
                 error('TifLoader: localization file not found')
            end
%             obj.imloader=imageLoader(p.tiffile);  
            if isempty(obj.imloader)
                obj.addFile(p.tiffile)
            else
                try 
                    img=(obj.imloader.getimage(1));
                    if isemtpy(img)
                        obj.addFile(p.tiffile)
                    end
                catch
                    disp('reload file in image loader')
                    
                end
                
%                 obj.imloader=imageloaderAll(p.tiffile,obj.getPar('loc_fileinfo')); 
            end
%             try 
%                 img=obj.imloader.getimage(1);
%             catch err
%                 err
%                  obj.imloader=imageloaderAll(p.tiffile,obj.getPar('loc_fileinfo')); 
%             end
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
            image=imloader.readNext ; 
%             imloader.currentImageNumber
            tall=0;
            tfitall=0;
            
            % parallel
            if p.parallelload
                parp=gcp;
            end
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
                if p.mirrorimage
                    image=image(:,end:-1:1);
                end
                datout.data=image;

                datout.frame=imloader.currentImageNumber;
                datout.ID=id;
                id=id+1;
                th=tic;
                obj.output(datout)
                tfitall=tfitall+toc(th);
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
                    
                    numlocs=obj.getPar('fittedLocs');
                    sep='''';
                    statuss=['frame ' separatethousands(datout.frame-obj.framestart,sep) ' of ' separatethousands(totalf,sep) ...
                        '. Time: ' separatethousands(elapsed,sep,0) ' of ' separatethousands(totaltime,sep,0) sprintf('s.')...
                        separatethousands(numlocs,sep) ' locs, ' separatethousands(numlocs/elapsed,sep,0) ' locs/s.'];
                    obj.status(statuss);
                end
                tall=tall+toc(th);
            end
            
            obj.setPar('tiffloader_loadingtime',tall);
            obj.setPar('tiffloader_fittime',tfitall);
            dateof=interfaces.WorkflowData;
            dateof.frame=imloader.currentImageNumber+1;
            dateof.ID=id;
            dateof.eof=true;
            obj.output(dateof)
            if ~obj.getPar('loc_preview') %if real fitting: close
                obj.imloader.close;
%                 obj.imloader=[];
                disp('close imloader')
            end
            disp('fitting done')
        end
        function addFile(obj,file)
%         if 1
%         else
            try
                obj.imloader.close;
            catch err
%                 err
            end
             lc=obj.getSingleGuiParameter('loaderclass');
            loader=obj.loaders{lc.Value,2};
            obj.imloader=loader(file,obj.getPar('loc_fileinfo'),obj.P);
            
%             obj.imloader=imageloaderAll(file,obj.getPar('loc_fileinfo'),obj.P);
%             if ~isempty(obj.imloader.info.metafile)
            if ~isempty(obj.imloader.metadata.allmetadata)&&isfield(obj.imloader.metadata.allmetadata,'metafile')
             obj.setPar('loc_metadatafile',obj.imloader.metadata.allmetadata.metafile);
            end
            allmd=obj.imloader.metadata;
            warnmissingmeta(allmd);
%             end
% obj.imloader.metadata
            p=obj.getGuiParameters;
            fileinf=obj.imloader.metadata;
            if p.padedges
%                 locsettings=obj.getPar('loc_cameraSettings');
                dr=p.padedgesdr;
                %dr=ceil((roisize-1)/2);
                fileinf.roi(1:2)=fileinf.roi(1:2)-dr;
                fileinf.roi(3:4)=fileinf.roi(3:4)+2*dr;
%                 fileinf.Width=fileinf.Width+2*dr;
%                 fileinf.Height=fileinf.Height+2*dr;
%                 obj.setPar('loc_cameraSettings',locsettings);
                obj.edgesize=dr;
            end
            obj.setPar('loc_fileinfo',fileinf);
            obj.setPar('loc_filename',file);
            
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

function warnmissingmeta(md)
expected={'emgain','conversion','offset','EMon','cam_pixelsize_um'};
missing=setdiff(expected,fieldnames(md));
if isempty(missing)
    return
end
str={'essential metadata missing (set manually in camera settings): ', missing{:}};
warndlg(str)
end

function changeloader(a,b,obj)
file=obj.getSingleGuiParameter('tiffile');
try
addFile(obj,file)
catch err
    obj.status('error in loading file. Choose different tiff loader')
    err
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
pard.tiffile.Width=3;

pard.loaderclass.object=struct('Style','popupmenu','String','auto','Callback',{{@changeloader,obj}});
pard.loaderclass.position=[2,4];
pard.loaderclass.Width=1;
pard.loaderclass.Optional=true;

pard.onlineanalysis.object=struct('Style','checkbox','String','Online analysis. Waittime (s):','Value',0);
pard.onlineanalysis.position=[3,2];
pard.onlineanalysis.Width=1.75;
pard.onlineanalysis.TooltipString='Fit during acquisition. If checked, max frames is ignored. Waits until no more images are written to file.';
pard.onlineanalysis.Optional=true;
pard.onlineanalysiswaittime.object=struct('Style','edit','String','5');
pard.onlineanalysiswaittime.position=[3,3.75];
pard.onlineanalysiswaittime.Width=0.5;
pard.onlineanalysiswaittime.TooltipString='Waiting time for online analysis.';
pard.onlineanalysiswaittime.Optional=true;

pard.parallelload.object=struct('Style','checkbox','String','parallel','Value',0);
pard.parallelload.position=[3,4.25];
pard.parallelload.Width=0.75;
pard.parallelload.TooltipString='Parallel load and process. Makes sense only for very slow loading process.';
pard.parallelload.Optional=true;


pard.textf.object=struct('Style','text','String','Frame range');
pard.textf.position=[4.2,1.25];
pard.textf.Width=0.75;
pard.textf.Optional=true;
pard.framestart.object=struct('Style','edit','String','1');
pard.framestart.position=[4.2,2];
pard.framestart.Width=0.5;
pard.framestart.Optional=true;
pard.framestop.object=struct('Style','edit','String','1000000');
pard.framestop.position=[4.2,2.5];
pard.framestop.Width=0.5;
pard.framestop.Optional=true;

pard.padedges.object=struct('Style','checkbox','String','Pad edges','Value',0,'Callback',@obj.loadedges);
pard.padedges.position=[4.2,3];
pard.padedges.Width=1;
pard.padedges.Optional=true;
pard.padedges.TooltipString='Pad edges with minimum values to allow detection of localizations close to edges. Usually not necessary.';

pard.padedgesdr.object=struct('Style','edit','String','9','Callback',@obj.loadedges);
pard.padedgesdr.position=[4.2,4];
pard.padedgesdr.Width=.3;
pard.padedgesdr.Optional=true;
pard.padedgesdr.TooltipString='Pad edges with minimum values to allow detection of localizations close to edges. Usually not necessary.';

pard.mirrorimage.object=struct('Style','checkbox','String','mirror');
pard.mirrorimage.position=[4.2,4.3];
pard.mirrorimage.Width=.7;
pard.mirrorimage.Optional=true;
pard.mirrorimage.TooltipString='Mirror image. Useful for conventional vs EM gain, that swaps image.';

% pard.locdata_empty.object=struct('Style','checkbox','String','Empty localizations','Value',1);
% pard.locdata_empty.position=[4.2,3.5];
% pard.locdata_empty.Width=1.5;
% pard.locdata_empty.TooltipString=sprintf('Empty localization data before fitting. \n Important if post-processing (eg drift correction) is perfromed as part of workflow');
pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Loads tiff files. Can load single tiff files or tiff stacks, also while they are written. It also tries to locate the metadata.txt file from micromanger and passes it on to the camera converter.';
end