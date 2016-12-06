classdef imageloaderSMAP<interfaces.GuiParameterInterface
    %imageloaderSMAP Superclass for image loader classes
    
    properties
        metadata=interfaces.metadataSMAP;
        file='';
        onlineAnalysis=false;
        waittime=5;
        currentImageNumber;
        allmetadatatags;
    end
    methods
       function obj=imageloaderSMAP(varargin)
           obj.metadata=interfaces.metadataSMAP;
           if nargin>2
               obj.P=varargin{3};
           end
           if nargin>1
               if ~isempty(varargin{2})
                obj.updatemetadata(varargin{2});
               end
           end
%            obj.getPar
            if nargin>0
                obj.open(varargin{1});
            end

        end
    end
    methods (Abstract=true)
        open(obj,filename);
%         md=getmetadata(obj);
        image=getimage(obj,frame);
        allmd=getmetadatatags(obj)
    end
    methods
        function image=readNext(obj)
            obj.currentImageNumber=obj.currentImageNumber+1;
            image=obj.getimage(obj.currentImageNumber);
        end
        function setImageNumber(obj,number)
            obj.currentImageNumber=number;
        end
        function images=getmanyimages(obj,numbers,format)
            if nargin<3
                format='cell';
            end
            if isempty(numbers)
                numbers=1:obj.metadata.numberOfFrames;
            end
            switch format
                case 'cell'
                    for k=length(numbers):-1:1
                        images{k}=obj.getimage(numbers(k));
                    end
                case 'mat'
                    for k=length(numbers):-1:1
                        images(:,:,k)=obj.getimage(numbers(k));
                    end
            end
            
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
            ind=find(strcmp(obj.allmetadatatags(:,1),tag),1,'first');
            if ~isempty(ind)
                val=obj.allmetadatatags{ind,2};
            else
                val=[];
            end
          
        end
        
        function metao=getmetadata(obj)
            metao=getmetadatacam(obj);
        end
        function metao=getmetadatacam(obj)
            md=getCameraCalibration(obj);
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
        metao=obj.metadata;
        end
    end
    
end

