classdef imageloaderSMAP<handle
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
           if nargin>1
                obj.updatemetadata(varargin{2});
           end
%            obj.getPar
            if nargin>0
                obj.open(varargin{1});
            end

        end
    end
    methods (Abstract=true)
        open(obj,filename);
        md=getmetadata(obj);
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
    end
    
end

