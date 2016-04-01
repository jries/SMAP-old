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
%              obj.setInputChannels(1,'frame');
        end
        function pard=pardef(obj)
            pard=[];
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
            if obj.getPar('loc_preview')
                return
            end
            obj.locDatatemp=interfaces.LocalizationData;
            obj.filenumber=1;%obj.locData.files.filenumberEnd+1;
%             obj.locData.files.filenumberEnd=1;%obj.filenumber;
%             obj.locData.loc=struct('frame',1,'filenumber',1,'xnm',0,'ynm',0,'channel',-1);
%             fn=fieldnames(obj.locData.loc);
%             obj.setPar('locFields',fn);
%             obj.locData.empty;
            tifs=struct('image',[],'info',[]);
            finfo=obj.getPar('loc_fileinfo');
            sfile=finfo.basefile;
            filename=[sfile '_sml.mat'];  
            infost=obj.getPar('cameraSettings');
            infost=copyfields(infost,finfo);
            filestruct=struct('info',infost,'average',[],'name',filename,...
                'number',obj.filenumber,...
                'numberOfTif',0,'tif',tifs);
%             filestruct.info=copyfields(filestruct.info,finfo);^
             if ~isfield(filestruct.info,'roi')||isempty(filestruct.info.roi)
                 filestruct.info.roi=[0 0 filestruct.info.Width filestruct.info.Height];
             end
             obj.locDatatemp.files.file=filestruct;  
                

%             obj.localtimervalue=tic;  
             obj.fileinfo=obj.getPar('cameraSettings');
%             obj.setPar('currentfileinfo',obj.fileinfo);
% 
%             obj.setPar('filelist_long',filelist,'String')
%             obj.setPar('filelist_short',filelists,'String')
            
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
            end
            
        end

    end
end

function locdat=fitloc2locdata(obj,locs,indin)
fieldsremove={'xerrpix','yerrpix','PSFxpix','PSFypix','xpix','ypix'};
fn=fieldnames(locs);
keepfields=setdiff(fn,fieldsremove);

for k=1:length(keepfields)
    locdat.(keepfields{k})=locs.(keepfields{k})(indin);
end

pixelsize=obj.fileinfo.pixsize*1000;
if isfield(obj.fileinfo,'roi')&&~isempty(obj.fileinfo.roi)
roi=obj.fileinfo.roi;
else
    roi=zeros(2);
end
% pixelsize=obj.globpar.parameters.cameraSettings.pixsize*1000;
% fields={'frame','phot','bg','bgerr','photerr'};
% locdat=copyfields([],locs,fields);
% locdat.frame=locs.frame(indin);
% locdat.phot=locs.phot(indin);
% locdat.bg=locs.bg(indin);
% if isfield(locs,'photerr')
% locdat.photerr=locs.photerr(indin);
% end
% if isfield(locs,'bgerr')
% locdat.bgerr=locs.bgerr(indin);
% end

locdat.xnm=(locs.xpix(indin)+roi(1))*pixelsize;
locdat.ynm=(locs.ypix(indin)+roi(2))*pixelsize;

if isfield(locs,'xerrpix')
locdat.xerr=locs.xerrpix(indin)*pixelsize;
else
    locdat.xerr=locdat.xnm*0+1;
end

if isfield(locs,'yerrpix')
locdat.yerr=locs.yerrpix(indin)*pixelsize;
else
    locdat.yerr=locdat.xerr;
end
if isfield(locs,'PSFxpix')
locdat.PSFxnm=locs.PSFxpix(indin)*pixelsize;
else
    locdat.PSFxnm=0*locdat.xnm+100;
end
if isfield(locs,'PSFypix')
    locdat.PSFynm=locs.PSFypix(indin)*pixelsize;
else
    locdat.PSFynm=locdat.PSFxnm;
end

% if isfield(locs,'gradient3Dellipticity')
% locdat.gradient3Dellipticity=locs.gradient3Dellipticity(indin);
% end


locdat.locprecnm=sqrt((locdat.xerr.^2+locdat.yerr.^2)/2);
locdat.filenumber=uint8(0*locdat.xnm+obj.filenumber);
locdat.channel=0*locdat.xnm;
end

