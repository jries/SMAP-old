classdef imageloaderSMAP<handle
    %imageloaderSMAP Superclass for image loader classes
    
    properties
        metadata=interfaces.metadataSMAP;
        file='';
        onlineAnalysis=false;
        waittime=5;
    end
    methods
       function obj=imageloaderSMAP(varargin)
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
    end
    
end

