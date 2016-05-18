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
        function pard=guidef(obj)
            pard=guidef;
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
            obj.locData.clear;
            obj.locData.files.file=locsaveFileinfo(obj);
                
            filelist={obj.locData.files.file(:).name};
            for k=1:length(filelist)
                [~,filelists{k}]=fileparts(filelist{k});
            end
            
            obj.fileinfo=obj.getPar('loc_cameraSettings');
            obj.setPar('currentfileinfo',obj.fileinfo);
            obj.setPar('filelist_long',filelist,'String')
            obj.setPar('filelist_short',filelists,'String')
            pos=(obj.fileinfo.roi(1:2)+obj.fileinfo.roi(3:4)/2)*obj.fileinfo.pixsize*1000;
            srs=obj.fileinfo.roi(3:4)*obj.fileinfo.pixsize*1000;
            pixels=max(srs./obj.getPar('sr_sizeRecPix'));
            obj.setPar('sr_pos',pos);
             obj.setPar('sr_pixrec',pixels);
            obj.locDatatemp=interfaces.LocalizationData;
             p=obj.parent.getGuiParameters(true,true);
             p.name=obj.parent.pluginpath;
             obj.locData.addhistory(p);
              obj.localtimervalue=tic;
        end
        function output=run(obj,data,p)
            output=[];
            if obj.getPar('loc_preview')||~obj.getSingleGuiParameter('update_check')
                return
            end
            locs=data.data;%get;
            
            if ~isempty(locs)&& ~data.eof
                maxfitdist=5;
                indin=abs(locs.xpix-locs.peakfindx)<maxfitdist & abs(locs.ypix-locs.peakfindy)<maxfitdist;
                locdat=interfaces.LocalizationData;
                locdat.loc=fitloc2locdata(obj,locs,indin);
                obj.locDatatemp.addLocData(locdat);

                if toc(obj.localtimervalue)>obj.getSingleGuiParameter('loc_updatetime')||data.eof 
                    obj.locData.addLocData(obj.locDatatemp); %Careful: does this add it many time? need to intialize obj.locDatatemp?
                   initGuiAfterLoad(obj);  %resets the view always! 
                    notify(obj.P,'sr_render')
                    drawnow
                    obj.localtimervalue=tic;
                    obj.locDatatemp.clear; %??? new
                end                    
            end           
        end
    end
end


function pard=guidef
pard.loc_updatetime.object=struct('Style','edit','String','30');
pard.loc_updatetime.position=[1,2.5];
pard.loc_updatetime.Width=0.5;
pard.loc_updatetime.TooltipString=sprintf('Render fitted data every XX seconds.');

pard.update_check.object=struct('Style','checkbox','String','Render update (s):','Value',1);
pard.update_check.position=[1,1];
pard.update_check.Width=1.5;
pard.update_check.TooltipString=sprintf('If checked, the fitted localizations are directly rendered and can be analyzed during the fit.');
 
pard.plugininfo.type='WorkflowModule'; 
end