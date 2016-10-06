classdef WorkflowFitter<interfaces.WorkflowModule
    properties
%         imagestack
%         stackinfo
        stackind
        numberInBlock=0;
        newID=1;
        fittedlocs=0;
    end
    methods

        function obj= WorkflowFitter(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
        end
        function pard=guidef(obj)
            pard=guidef;
        end

        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.inputParameters={'loc_ROIsize','loc_cameraSettings'};
            obj.stackind=0;

        end
        function fitinit(obj) %dummy function
        end
        function prerun(obj,p)
            global fitterstackinfo fitterimagestack fitterbgstack
%             obj.numberInBlock=1; %round(5500*100/roisize^2);
            obj.fitinit;
            roisize=obj.getPar('loc_ROIsize');
%             obj.numberInBlock=round(5500*100/roisize^2);
            
%             disp(['number in block: ' num2str(obj.numberInBlock)]);
            obj.stackind=0;
            obj.fittedlocs=0;
            fitterimagestack=zeros(roisize,roisize,obj.numberInBlock,'single');
            if obj.inputChannels==2
                fitterbgstack=fitterimagestack;
            else
                fitterbgstack=[];
            end
%             infos=struct('x',0,'y',0,'frame',0);
%             fitterstackinfo=infos;
            fitterstackinfo.x=zeros(obj.numberInBlock,1,'single');
            fitterstackinfo.y=zeros(obj.numberInBlock,1,'single');
            fitterstackinfo.frame=zeros(obj.numberInBlock,1,'double');
            obj.fitinit;
            
            
            
%             obj.setInputChannels(2,'frame');
%             p=obj.getGuiParameters.par;
%             obj.loc_ROIsize=p.loc_ROIsize;
           
        end
        function out=run(obj,data,p)  
            
            global fitterimagestack fitterstackinfo fitterbgstack 
            persistent reporttimer
            out=[];
            passbg=(obj.inputChannels==2);%~isempty(fitterbgstack);
            fninfo=fieldnames(fitterstackinfo);
            if ~iscell(data) 
                dstruc=data.data;
                eof=data.eof;
            else
                dstruc=data{1}.data;%get;
                eof=data{1}.eof;
                
            end
            
            if isempty(reporttimer)
                reporttimer=tic;
            end
            if ~isempty(dstruc)&&~isempty(dstruc.img)       
                 imgstack=dstruc.img;
                 
                 stackinf=dstruc.info;
                 s=size(imgstack);
                 if length(s)==2
                     s(3)=1;
                 end                 
                 
                 obj.fittedlocs=s(3)+obj.fittedlocs;
                 if toc(reporttimer)>2
                     obj.setPar('fittedLocs',obj.fittedlocs);
                     reporttimer=tic;
                 end
                
                if passbg
                    bgstack=data{2}.data.img;
                end
                
                numberInBlockh=obj.numberInBlock;
                if numberInBlockh>1 %make blocks
                 stackindh=obj.stackind;%pointer to last element
                 stackindh=stackindh+1; %new pointer
                 
                 %avoid loop
                 imagesleft=s(3);
                 startinstack=1;
                 
                 newstackind=stackindh+imagesleft-1;
                 while newstackind>numberInBlockh
                     imagestowrite=numberInBlockh-stackindh+1;
                     
                     %images
                     fitterimagestack(:,:,stackindh:end)=imgstack(:,:,startinstack:startinstack+imagestowrite-1);
                     if passbg %background
                        fitterbgstack(:,:,stackindh:end)=bgstack(:,:,startinstack:startinstack+imagestowrite-1);
                     end
                     for fsi=1:length(fninfo) %stackinfor
                        fitterstackinfo.(fninfo{fsi})(stackindh:end)=stackinf.(fninfo{fsi})(startinstack:startinstack+imagestowrite-1);
                     end
                     
                     %do fitting
                     if passbg
                        locs=obj.fit(fitterimagestack,fitterbgstack,fitterstackinfo);
                     else
                         locs=obj.fit(fitterimagestack,fitterstackinfo);
                     end
                     outputlocs(obj,locs,fitterstackinfo,obj.newID,eof);
                     obj.newID=obj.newID+1;
                     
                     imagesleft=imagesleft-imagestowrite;                   
                     startinstack=startinstack+imagestowrite;         
                     stackindh=1;
                     newstackind=imagesleft;
                 end
                 fitterimagestack(:,:,stackindh:stackindh+imagesleft-1)=imgstack(:,:,startinstack:end);
                 if passbg
                    fitterbgstack(:,:,stackindh:stackindh+imagesleft-1)=bgstack(:,:,startinstack:end);
                 end
                 for fsi=1:length(fninfo) %stackinfor
                        fitterstackinfo.(fninfo{fsi})(stackindh:stackindh+imagesleft-1)=stackinf.(fninfo{fsi})(startinstack:end);
                 end
                 obj.stackind=stackindh+imagesleft-1;
                 
                else %go framewise
                     if passbg
                        locs=obj.fit(imgstack,bgstack,stackinf);
                     else
                        locs=obj.fit(imgstack,stackinf);
                     end
                     outputlocs(obj,locs,stackinf,obj.newID,eof);
                     obj.newID=obj.newID+1;
                end
            end
            
            
            if eof
                if obj.numberInBlock>1
                    for fsi=1:length(fninfo) %stackinfor
                        fitterstackinfo.(fninfo{fsi})=fitterstackinfo.(fninfo{fsi})(1:obj.stackind);
                    end
                    if passbg
                        locs=obj.fit(fitterimagestack(:,:,1:obj.stackind),fitterbgstack(:,:,1:obj.stackind),fitterstackinfo);
                    else
                        locs=obj.fit(fitterimagestack(:,:,1:obj.stackind),fitterstackinfo);
                    end


                    outputlocs(obj,locs,fitterstackinfo,obj.newID,true);
                else
                    outputlocs(obj,[],[],obj.newID,true);
                end
            end

        end
       
    end
end

function outputlocs(obj,locs,stackinfo,tag,eof)
dato=interfaces.WorkflowData;
dato.ID=tag;
dato.eof=eof;
dato.data=(locs);
obj.output(dato);

end

function pard=guidef

pard=[];

end