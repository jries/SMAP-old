classdef imageloaderSMAP<handle
    %imageloaderSMAP Superclass for image loader classes
    
    properties
        metadata=interfaces.metadataSMAP;
        file='';
        onlineAnalysis=false;
        waittime=5;
        currentImageNumber;
    end
    methods
       function obj=imageloaderSMAP(varargin)
           obj.metadata=interfaces.metadataSMAP;
            if nargin>0
                obj.open(varargin{1});
            end
        end
    end
    methods (Abstract=true)
        open(obj,filename);
        md=getmetadata(obj);
        image=getimage(obj,frame);
    end
    methods
        function image=readNext(obj)
            obj.currentImageNumber=obj.currentImageNumber+1;
            image=obj.getimage(obj.currentImageNumber);
        end
        function setImageNumber(obj,number)
            obj.currentImageNumber=number;
        end
        function images=getmanyimages(obj,numbers)
            for k=length(numbers):-1:1
                images{k}=obj.getimage(numbers(k));
            end
            
        end
    end
    
end

