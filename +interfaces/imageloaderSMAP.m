classdef imageloaderSMAP<interfaces.GuiParameterInterface
    %imageloaderSMAP Superclass for image loader classes
    
    properties
        metadata=interfaces.metadataSMAP;
        file='';
        onlineAnalysis=false;
        waittime=5;
        currentImageNumber;
        allmetadatatags;
        calibrationFile='settings/cameras.mat';
        
        ismultichannel=false;
        multiloader={};
        multiloadermetadata;
    end
    methods
       function obj=imageloaderSMAP(varargin)
           obj.metadata=interfaces.metadataSMAP;
           if nargin>2 && ~isempty(varargin{3})
               if isa(varargin{3},'interfaces.ParameterData')
                    obj.P=varargin{3};
                    obj.calibrationFile=obj.getGlobalSetting('cameraSettingsFile');
               else
                   obj.calibrationFile=varargin{3};
               end
           end
           if nargin>1 && ~isempty(varargin{2})
                obj.updatemetadata(varargin{2});
                obj.multiloadermetadata=varargin{2};
           end

            if nargin>0 && ~isempty(varargin{1})
                obj.open(varargin{1});
            end

        end
    end
    methods (Abstract=true)
        openi(obj,filename);
        image=getimagei(obj,frame);
        allmd=getmetadatatagsi(obj)
    end
    methods
        function image=getimage(obj,frame)
            if obj.ismultichannel
                for k=length(obj.multiloader):-1:1
                    image{k}=obj.multiloader{k}.getimagei(frame);
                end
            else
                image=obj.getimagei(frame);
            end
        end
        function image=readNext(obj)
            obj.currentImageNumber=obj.currentImageNumber+1;
            image=obj.getimageonline(obj.currentImageNumber);
        end
        
        function image=getimageonline(obj,number)
            image=obj.getimage(number);
            if isempty(image)&&obj.onlineAnalysis 
                    disp('wait')
                    pause(obj.waittime*2)
                    image=obj.getimage(number);
            end
       
        end
        
        function setImageNumber(obj,number)
            obj.currentImageNumber=number;
        end
        function images=getmanyimages(obj,numbers,format)
            loadfun=@obj.getimageonline;
            if nargin<3
                format='cell';
            end
            if nargin<2|| isempty(numbers)
                numbers=1:obj.metadata.numberOfFrames;
            end
            switch format
                case 'cell'
                    for k=length(numbers):-1:1
                        images{k}=loadfun(numbers(k));
                        if isempty(images{k})
                            loadfun=@obj.getimage;
                        end
                    end
                case 'mat'
                    for k=length(numbers):-1:1
                        imh=loadfun(numbers(k));
                        if ~isempty(imh)
                        images(:,:,k)=imh;
                        else
                            loadfun=@obj.getimage;
                        end
                    end
            end
            
        end
        function file=getbasefile(obj)
            file=obj.file;
        end
        function updatemetadata(obj, md)
            if isempty(md)
                return;
            end
            fn=fieldnames(md);fn2=properties(obj.metadata);
           fna=intersect(fn,fn2);
           obj.metadata=copyfields(obj.metadata,md,fna);
        end
        
        function val=gettag(obj,tag)
            if isempty(obj.allmetadatatags)
                obj.allmetadatatags=obj.getmetadatatags;
            end
            val=[];
            if ~isempty(obj.allmetadatatags)
                ind=find(strcmp(obj.allmetadatatags(:,1),tag),1,'first');
                if ~isempty(ind)
                    val=obj.allmetadatatags{ind,2};
                end
                
            end
          
        end
        
        function metao=getmetadata(obj)
            getmetadatacam(obj);
             obj.metadata.basefile=obj.getbasefile;
             metao=obj.metadata;   
        end
        function metao=getmetadatacam(obj)
            try
                camfile=obj.getGlobalSetting('cameraSettingsFile');
            catch err
                camfile=obj.calibrationFile;
                display('could not find camera file in global settings. Using default file.')
%               ererwe
              
            end
            try
                usedef=obj.getPar('useDefaultCam');
            catch
                usedef=true;
            end
            md=getCameraCalibration(obj,[],usedef,camfile);
            if isempty(md)
                metao=[];
                return
            end
            fn=fieldnames(md);
            for k=1:length(fn)
                if ~isempty(md.(fn{k}))&&isprop(obj.metadata,fn{k})
                    obj.metadata.(fn{k})=md.(fn{k});
                    obj.metadata.assigned.(fn{k})=true;
                end
            end
            obj.metadata.allmetadata=obj.allmetadatatags;
            obj.metadata.imagefile=obj.file;
%             if isempty(obj.metadata.basefile)
               
%             end
        metao=obj.metadata;
        end
        function open(obj,file)
            if iscell(file) %multiple channels in multiple files
                obj.ismultichannel=true;
                for k=1:length(file)
                    obj.multiloader{k}= obj.newimageloader(file{k});
                end
                obj.metadata=obj.multiloader{1}.metadata; %copy metadata and infor from first
            else
                obj.openi(file)
            end
            
        end
        function close(obj)
            display(['close not implemented in ' class(obj)])
        end
        function allmd=getmetadatatags(obj)
            if obj.ismultichannel
                for k=1:length(obj.multiloader)
                    allmd=obj.multiloader{k}.getmetadatatags;
                end
                obj.metadata=obj.multiloader{1}.metadata;
            else
                allmd=obj.getmetadatatagsi;
            end
        end
        
        function il=newimageloader(obj,file)
            il=eval(class(obj));
            il.P=obj.P;
            il.calibrationFile=obj.calibrationFile;
            il.onlineAnalysis=obj.onlineAnalysis;
            il.waittime=obj.waittime;
            if ~isempty(obj.multiloadermetadata) 
                il.updatemetadata(obj.multiloadermetadata)
            end
            il.open(file);
        end
    end
    
end

