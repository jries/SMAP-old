classdef LocSaver<interfaces.WorkflowModule;
    properties
        filenumber
        fileinfo
        locDatatemp;
        
    end
    methods
       function obj=LocSaver(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
        end
        function pard=guidef(obj)
            pard.plugininfo.type='WorkflowModule'; 
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
            if obj.getPar('loc_preview')
                return
            end
            obj.locDatatemp=interfaces.LocalizationData;
            obj.filenumber=1;
             obj.locDatatemp.files.file=locsaveFileinfo(obj);  
             p=obj.parent.getGuiParameters(true,true);
             p.name=obj.parent.pluginpath;
             
             obj.locDatatemp.addhistory(p);
             
             obj.fileinfo=obj.getPar('loc_cameraSettings');    %not used?       
        end
        function output=run(obj,data,p)
            output=[];
            if obj.getPar('loc_preview')
                return
            end
            locs=data.data;%get;
            
            if ~isempty(locs)
                maxfitdist=5;
                indin=abs(locs.xpix-locs.peakfindx)<maxfitdist & abs(locs.ypix-locs.peakfindy)<maxfitdist;
                locdat=interfaces.LocalizationData;
                locdat.loc=fitloc2locdata(obj,locs,indin);
                obj.locDatatemp.addLocData(locdat);
            end
            
            if data.eof %save locs

                filenameold=obj.locDatatemp.files.file(1).name;
                filename=filenameold;
                ind=2;
                while exist(filename,'file')
                    filename=[filenameold(1:end-7) num2str(ind) '_sml.mat'];
                    ind=ind+1;
                end
                fitpar=obj.parent.getGuiParameters(true).children;
                obj.locDatatemp.savelocs(filename,[],struct('fitparameters',fitpar));
                if isempty(obj.locData.loc)
                    obj.locData.addLocData(obj.locDatatemp);
                end
                [path,file]=fileparts(filename);
                imageout=makeSRimge(obj.locDatatemp);
                options.comp='jpeg';
                options.color=true;
                s=size(imageout);
                sr=ceil(s/16)*16;
                imageout(sr(1),sr(2),1)=0;
                saveastiff(uint8(imageout/max(imageout(:))*255),[path filesep file '.tif'],options)
                output=data;
            end
            
        end

    end
end

function imout=makeSRimge(locDatatemp)
channelfile='settings/FitTif_Channelsettings.mat';
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

