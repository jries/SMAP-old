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
              obj.run(); %clear
        end
        function output=run(obj,data,p)
            persistent templocs numlocs
            if nargin<2
                templocs=[];numlocs=[];
                return
            end
            
            output=[];
            if obj.getPar('loc_preview')||~obj.getSingleGuiParameter('update_check')
                return
            end
            locs=data.data;%get;
            
            if ~isempty(locs)&&~isempty(locs.frame)
                maxfitdist=5;
                indin=abs(locs.xpix-locs.peakfindx)<maxfitdist & abs(locs.ypix-locs.peakfindy)<maxfitdist;
               
%                 locdat=interfaces.LocalizationData;
%                 locdat.loc=fitloc2locdata(obj,locs,indin);
%                 obj.locDatatemp.addLocData(locdat);
                fn=fieldnames(locs);
                if isempty(templocs)
                    for k=1:length(fn)
                        templocs.(fn{k})=locs.(fn{k})(indin);
                    end
                    numlocs=length(templocs.(fn{1}));
                else
                    newlocs=length(locs.(fn{1}));
                    if numlocs+newlocs>length(templocs.(fn{1}))
                        newlen=min(1000,2*(numlocs+length(locs.(fn{1}))));
                        for k=1:length(fn)
                            templocs.(fn{k})(newlen)=0;
                        end
                    end
                    for k=1:length(fn)
                        templocs.(fn{k})(numlocs+1:numlocs+sum(indin))=locs.(fn{k})(indin);
                    end
                    numlocs=numlocs+sum(indin);
                end
                if toc(obj.localtimervalue)>obj.getSingleGuiParameter('loc_updatetime')||data.eof 
                    locdat=interfaces.LocalizationData;
                    locdat.loc=fitloc2locdata(obj,templocs,1:numlocs);
%                     obj.locDatatemp.addLocData(locdat);
                    templocs=[];
                    numlocs=[];
                    
                    obj.locData.addLocData(locdat); %Careful: does this add it many time? need to intialize obj.locDatatemp?
                    initGuiAfterLoad(obj,false);  %resets the view always! 
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
pard.plugininfo.description='Passes on the fitted localizations to the rendering engine of SMAP. This happens at user-defined time intervals.';
end