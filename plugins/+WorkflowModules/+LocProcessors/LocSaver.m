classdef LocSaver<interfaces.WorkflowModule;
    properties
        timer
        filenumber
        fileinfo
        locDatatemp;
         deltaframes;
         index;
        numsaved
        frames;
        saveframes=20;
        
    end
    methods
       function obj=LocSaver(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
        end
        function pard=guidef(obj)
            pard.plugininfo.type='WorkflowModule'; 
            pard.plugininfo.description='Saves the fitted localizations as a SMAP *.sml file. When fitting via a network, fitting a local copy which is then moved to the destination can be faster.';
            pard.savelocal.object=struct('Style','checkbox','String','save local and copy','Value',0);
            pard.savelocal.object.TooltipString='Select this if you fit via a network and the saving of the localizations is very long (stauts bar stops for a long time at last frames).';
            pard.savelocal.position=[1,1];
            pard.savelocal.Width=2;
            pard.savelocal.Optional=true;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
            if obj.getPar('loc_preview')
                return
            end
            numberofframes=obj.getPar('loc_fileinfo').numberOfFrames;
            obj.deltaframes=floor(numberofframes/obj.saveframes);
             obj.index=round(obj.deltaframes/2);
            obj.numsaved=0;
%             obj.frames=struct('image',[],'frame',[]);
            obj.frames=[];
            obj.locDatatemp=interfaces.LocalizationData;
            obj.filenumber=1;
             obj.locDatatemp.files.file=locsaveFileinfo(obj);  
             p=obj.parent.getGuiParameters(true,true);
             p.name=obj.parent.pluginpath;
             
             obj.locDatatemp.addhistory(p);
             
             obj.fileinfo=obj.getPar('loc_cameraSettings');    %not used?  
             obj.timer=tic;
        end
        function output=run(obj,data,p)
            output=[];
            if obj.getPar('loc_preview')
                return
            end
            locs=data.data;%get;
            
            if ~isempty(locs)&&~isempty(locs.frame)
                maxfitdist=5;
                indin=abs(locs.xpix-locs.peakfindx)<maxfitdist & abs(locs.ypix-locs.peakfindy)<maxfitdist;
                locdat=interfaces.LocalizationData;
                locdat.loc=fitloc2locdata(obj,locs,indin);
                obj.locDatatemp.addLocData(locdat);
                    if locs.frame(end)>obj.index && obj.numsaved<obj.saveframes
                        obj.numsaved=obj.numsaved+1;
                        obj.index=obj.index+obj.deltaframes;
                        if isempty(obj.frames)
                            obj.frames=obj.getPar('loc_currentframe');
                        else
                            obj.frames(obj.numsaved)=obj.getPar('loc_currentframe');
                        end
                    end
            end
            
            if data.eof %save locs

                filenameold=obj.locDatatemp.files.file(1).name;
                filename=filenameold;
                ind=2;
                while exist(filename,'file')
                    filename=[filenameold(1:end-7) num2str(ind) '_sml.mat'];
                    ind=ind+1;
                end
                mainfile=filename;
                if p.savelocal
                    filenameremote=filename;
                    filename=[pwd filesep 'temp.sml'];
                    mainfile=filename;
                end
                fitpar=obj.parent.getGuiParameters(true).children;
                fitpar.fittime=toc(obj.timer);
                fitpar.loadtifftime=obj.getPar('tiffloader_loadingtime');
                try
                disp([num2str(length(obj.locDatatemp.loc.xnm)) ' localizations in ' num2str(fitpar.fittime) ' seconds.']);
                catch err
                    err
                end
                obj.locDatatemp.files.file.raw=obj.frames;
                obj.locDatatemp.savelocs(filename,[],struct('fitparameters',fitpar));
                
                if p.savelocal
                    movefile(filename,filenameremote);
                end
%               write to main GUI
%                 obj.locData.clear;
                obj.locData.setLocData(obj.locDatatemp);

                initGuiAfterLoad(obj);
                obj.setPar('mainfile',mainfile);
                [path,file]=fileparts(filename);
                try
                imageout=makeSRimge(obj.locDatatemp);
                options.comp='jpeg';
                options.color=true;
                s=size(imageout);
                sr=ceil(s/16)*16;
                imageout(sr(1),sr(2),1)=0;
                saveastiff(uint16(imageout/max(imageout(:))*(2^16-1)),[path filesep file '.tif'],options)
                catch err
                    err
                end
                output=data;
            end
            
        end

    end
end

function imout=makeSRimge(locDatatemp)
channelfile='settings/workflows/FitTif_Channelsettings.mat';
pall=load(channelfile);
p=pall.globalParameters;
p.sr_pixrec=20;
p.layer=1;
% p.gaussfac=0.4;
minx=min(locDatatemp.loc.xnm);maxx=max(locDatatemp.loc.xnm);
miny=min(locDatatemp.loc.ynm);maxy=max(locDatatemp.loc.ynm);
p.sr_pos=[(minx+maxx)/2 (miny+maxy)/2];
p.sr_size=[(maxx-minx)/2 (maxy-miny)/2];
p.sr_layerson=1;
rawimage=renderSMAP(locDatatemp.loc,p);
layers.images.finalImages=drawerSMAP(rawimage,p);
imoutt=displayerSMAP(layers,p);
imout=imoutt.image;
% imout=TotalRender(locDatatemp,pall.defaultChannelParameters);
end

