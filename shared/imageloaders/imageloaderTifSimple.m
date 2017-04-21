classdef imageloaderTifSimple<interfaces.imageloaderSMAP
    %imageloaderMM image loader for micromanager  tiff stack files
    %   Detailed explanation goes here
    
    properties
%         reader
    end
    
    methods
        function obj=imageloaderTifSimple(varargin)
            obj@interfaces.imageloaderSMAP(varargin{:});
        end
        function open(obj,file)
        
%             obj.reader = javaObjectEDT('org.micromanager.acquisition.TaggedImageStorageMultipageTiff',fileparts(file), false, [], false, false, true);
            obj.file=file;
%             obj.reader=bfGetReader(file);
            md=obj.getmetadata;
%             [p,f]=fileparts(file);
            obj.metadata.basefile=file;
            
        end
        function image=getimage(obj,frame)
            image=imread(obj.file,'Index',frame);
        end
        
        function close(obj)
%             obj.reader.close
%             clear(obj.reader)
        end
        
        function image=getimageonline(obj,number)
            image=obj.getimage(number);
            if isempty(image)&&obj.onlineAnalysis 
                    disp('wait')
%                     obj.reader.close;
%                     delete(obj.reader)
                    pause(obj.waittime*2)
%                     obj.reader = javaObjectEDT('org.micromanager.acquisition.TaggedImageStorageMultipageTiff',fileparts(obj.file), false, [], false, false, true);
                    image=obj.getimage(number);
            end
        end
        
        function allmd=getmetadatatags(obj)
            allmd={'Format','SimpleTif'};
            obj.allmetadatatags=allmd;
                
        
        end
        
    end
    
end

