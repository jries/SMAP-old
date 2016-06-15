classdef imageloaderOME<interfaces.imageloaderSMAP
    %imageloaderMM image loader for micromanager single tiff files
    %   Detailed explanation goes here
    
    properties
        calfile='settings/CameraCalibration.xls';
        reader
    end
    
    methods
        function obj=imageloaderOME(varargin)
            obj@interfaces.imageloaderSMAP(varargin{:});
        end
        function open(obj,file)
            obj.file=file;
            obj.reader=bfGetReader(file);
            md=obj.getmetadata;
            
        end
        function mdo=getmetadata(obj)
            mdo=getmetadataome(obj);
             
        end
        function image=getimage(obj,frame)
            image=readstack(obj,frame);
        end
        function updatestack(obj,info)
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
    end
    
end
function meta=getmetadataome(obj)
reader=obj.reader;
omemeta=reader.getMetadataStore;
globalmeta=reader.getGlobalMetadata;
[~,~,ext]=fileparts(obj.file);
switch ext
    case '.nd2' %Nikon
        meta=getMetaNd2(globalmeta);
end
m2.pixsize=double(omemeta.getPixelsPhysicalSizeX(0).value());
obj.metadata=copyfields(obj.metadata,meta);
obj.metadata=copyfields(obj.metadata,m2);
end

function meta=getMetaNd2(m)
searchstr={'exposure','Exposure:','emgain','Multiplier:'};
    str=m.get('sSpecSettings'); %
    for k=1:2:length(searchstr);
        ind=strfind(str,searchstr{k+1})+length(searchstr{k+1});
        ind2=ind+1;
        while str(ind2)>='0'&& str(ind2)<='9'
            ind2=ind2+1;
        end
        meta.(searchstr{k})=double(str2num(str(ind:ind2-1)));
    end
    meta.EMon=str2num(m.get('EnableGainMultiplier'));
    
end




function image=readstack(obj,imagenumber,recursions)
    if nargin<3
        recursions=0;
    end
    maxrecursions=1;
    image=[]; %default, if nothing can be read: end.

    [imagenumber,filenumber]=obj.getstacknumber(imagenumber);
    if filenumber~=obj.stack.currentfile
        obj.stack.tiffh.close;
        newfile=obj.stack.files{filenumber};
%         obj.info=getimageinfo(newfile);
        obj.updatestack(getimageinfo(newfile));
        obj.stack.tiffh=Tiff(newfile,'r');
        obj.stack.currentfile=filenumber;
    end
    try
        th=obj.stack.tiffh;
        th.setDirectory(imagenumber);
        image=th.read;                    
    catch

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
        end
    end
end