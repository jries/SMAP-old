classdef OnlineReconstruction<interfaces.WorkflowModule;
    properties
        localupdatetime=10
        localtimervalue
        filenumber
        fileinfo
        locDatatemp;
        
    end
    methods
       function obj=OnlineReconstruction(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
%              obj.setInputChannels(1,'frame');
        end
        function pard=pardef(obj)
            pard=pardef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
            if obj.getPar('loc_preview')
                return
            end
            obj.localupdatetime=obj.getGuiParameters.loc_updatetime;
            obj.filenumber=1;%obj.locData.files.filenumberEnd+1;
            obj.locData.files.filenumberEnd=1;%obj.filenumber;
            obj.locData.loc=struct('frame',1,'filenumber',1,'xnm',0,'ynm',0,'channel',-1);
            fn=fieldnames(obj.locData.loc);
            obj.setPar('locFields',fn);
            obj.locData.empty;
            tifs=struct('image',[],'info',[]);
            finfo=obj.getPar('loc_fileinfo');
            sfile=finfo.basefile;
            filename=[sfile '_sml.mat'];           
            filestruct=struct('info',obj.getPar('cameraSettings'),'average',[],'name',filename,...
                'number',obj.filenumber,...
                'numberOfTif',0,'tif',tifs);
                obj.locData.files.file=filestruct;  
                
            filelist={obj.locData.files.file(:).name};
            for k=1:length(filelist)
                [~,filelists{k}]=fileparts(filelist{k});
            end
            obj.localtimervalue=tic;  
            obj.fileinfo=obj.getPar('cameraSettings');
            obj.setPar('currentfileinfo',obj.fileinfo);

            obj.setPar('filelist_long',filelist,'String')
            obj.setPar('filelist_short',filelists,'String')
            obj.locDatatemp=interfaces.LocalizationData;
        end
        function output=run(obj,data,p)
            output=[];
            if obj.getPar('loc_preview')||~p.update_check
                return
            end
            locs=data.data;%get;
            
            if ~isempty(locs)
                maxfitdist=5;
                indin=abs(locs.xpix-locs.peakfindx)<maxfitdist & abs(locs.ypix-locs.peakfindy)<maxfitdist;
                locdat=interfaces.LocalizationData;
                locdat.loc=fitloc2locdata(obj,locs,indin);
                obj.locDatatemp.addLocData(locdat);

                if toc(obj.localtimervalue)>obj.localupdatetime||data.eof 
                    obj.locData.addLocData(obj.locDatatemp); %Careful: does this add it many time? need to intialize obj.locDatatemp?
                   initGuiAfterLoad(obj);  
%                     obj.locData.regroup;
%                     obj.setPar('locFields',fieldnames(obj.locData.loc));
                    notify(obj.P,'sr_render')
                    drawnow
                    obj.localtimervalue=tic;
                    obj.locDatatemp.empty; %??? new

                end   
            
           
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

function pard=pardef
pard.loc_updatetime.object=struct('Style','edit','String','30');
pard.loc_updatetime.position=[1,2.5];
pard.loc_updatetime.Width=0.5;

pard.update_check.object=struct('Style','checkbox','String','Render update (s):','Value',1);
pard.update_check.position=[1,1];
pard.update_check.Width=1.5;

      
end