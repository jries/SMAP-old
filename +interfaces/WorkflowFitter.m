classdef WorkflowFitter<interfaces.WorkflowModule
    properties
%         imagestack
%         stackinfo
        stackind
        numberInBlock=2500;
        newID=1;
    end
    methods

        function obj= WorkflowFitter(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
        end
        function pard=pardef(obj)
            pard=pardef;
        end

        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.inputParameters={'ROI_size','EMexcessNoise'};
            obj.stackind=0;

        end
        function fitinit(obj) %dummy function
        end
        function prerun(obj,p)
            global fitterstackinfo fitterimagestack fitterbgstack
            obj.stackind=0;
            roisize=obj.getPar('ROI_size');
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
%             obj.ROI_size=p.ROI_size;
           
        end
        function out=run(obj,data,p)  
            
            global fitterimagestack fitterstackinfo fitterbgstack
            out=[];
            passbg=~isempty(fitterbgstack);
            fninfo=fieldnames(fitterstackinfo);
            if ~iscell(data) 
                dstruc=data.data;
                eof=data.eof;
            else
                dstruc=data{1}.data;%get;
                eof=data{1}.eof;
                
            end
            if ~isempty(dstruc)&&~isempty(dstruc.img)       
                 imgstack=dstruc.img;
                 
                 stackinf=dstruc.info;
                 s=size(imgstack);
                 if length(s)==2
                     s(3)=1;
                 end                 
                 
                
                if passbg
                    bgstack=data{2}.data.img;
                end

                
                 stackindh=obj.stackind;%pointer to last element
                 stackindh=stackindh+1; %new pointer
                 numberInBlockh=obj.numberInBlock;
                 %avoid loop
                 imagesleft=s(3);
                 startinstack=1;
                 
                 newstackind=obj.stackind+imagesleft;
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
            end
            if eof
                for fsi=1:length(fninfo) %stackinfor
                    fitterstackinfo.(fninfo{fsi})=fitterstackinfo.(fninfo{fsi})(1:obj.stackind);
                end
                if passbg
                    locs=obj.fit(fitterimagestack(:,:,1:obj.stackind),fitterbgstack(:,:,1:obj.stackind),fitterstackinfo);
                else
                    locs=obj.fit(fitterimagestack(:,:,1:obj.stackind),fitterstackinfo);
                end

                
                outputlocs(obj,locs,fitterstackinfo,obj.newID,true);
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

function pard=pardef

pard=[];

end