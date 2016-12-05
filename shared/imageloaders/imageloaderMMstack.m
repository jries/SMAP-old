classdef imageloaderMMstack<interfaces.imageloaderSMAP
    %imageloaderMM image loader for micromanager single tiff files
    %   Detailed explanation goes here
    
    properties
        calfile='settings/CameraCalibration.xls';
        stack
    end
    
    methods
        function obj=imageloaderMMstack(varargin)
            obj@interfaces.imageloaderSMAP(varargin{:});
        end
        function open(obj,file)
            obj.file=file;
            md=obj.getmetadata;
            obj.stack.currentfile=1;
            obj.updatestack(md.allmetadata);
            warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning')
            obj.stack.tiffh=Tiff(obj.file,'r');
        end
%         function mdo=getmetadata(obj)
%             mdo=getmetadataMM(obj);
%              
%         end
        function image=getimage(obj,frame)
            image=readstack(obj,frame);
        end
        function updatestack(obj,info)
            
            info=getimageinfo(obj.file);
             numframes=[info.allfiles(1:end-1).numberOfFrames];
            obj.stack.imageoffset=cumsum([0 numframes]);
            obj.stack.lastimage=cumsum([info.allfiles(1:end).numberOfFrames]);
            obj.stack.files={info.allfiles(:).name};
        end
        function [imagenumber,filenumberc,filenumber]=getstacknumber(obj,frame)
            filenumber=find(frame<=obj.stack.lastimage,1,'first');
            filenumberc=filenumber;
            if isempty(filenumber)
                filenumberc=length(obj.stack.lastimage);
            end
            imagenumber=frame-obj.stack.imageoffset(filenumberc);
        end
        function allmd=getmetadatatags(obj)
%             metafile=[fileparts(obj.file) filesep 'metadata.txt'];
            allmd=getmetadataMMnew(obj.file);
        end
    end
    
end

function image=readstack(obj,imagenumber,recursions)
    if nargin<3
        recursions=0;
    end
    maxrecursions=1;
    image=[]; %default, if nothing can be read: end.
%     image=imagenumber;
    [imagenumber,filenumber]=obj.getstacknumber(imagenumber);
    if filenumber~=obj.stack.currentfile
        obj.stack.tiffh.close;
        newfile=obj.stack.files{filenumber};
%         obj.info=getimageinfo(newfile);
        obj.updatestack(getimageinfo(newfile));
        obj.stack.tiffh=Tiff(newfile,'r');
        obj.stack.currentfile=filenumber;
    end
%     if ~isvalid(obj.stack.tiffh)
%         obj.stack.tiffh=Tiff(obj.stack.files{obj.stack.currentfile},'r');  
%     end
    try
        th=obj.stack.tiffh;
        th.setDirectory(imagenumber);
        image=th.read;                    
    catch err

        if obj.onlineAnalysis
            if recursions<maxrecursions
                pause(obj.waittime*2);
                obj.stack.currentfile=-1; %update stuff
                image=readstack(obj,imagenumber,recursions+1);
            else
                disp('after waiting no more files')
            end
        else

            disp('image out of range');
            image=obj.stack.tiffh;
        end
    end
end