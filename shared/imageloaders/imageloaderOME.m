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

    end
    
end
function metao=getmetadataome(obj)
reader=obj.reader;
omemeta=reader.getMetadataStore;
globalmeta=reader.getGlobalMetadata;
[ph,fh,ext]=fileparts(obj.file);
fn={};
switch ext
    case '.nd2' %Nikon
        meta=getMetaNd2(globalmeta);
        fn=fieldnames(meta);
    otherwise
        meta=[];
end
try
m2.pixsize=double(omemeta.getPixelsPhysicalSizeX(0).value());
fn=[fn fieldnames(m2)];
catch
    m2=[];
end
obj.metadata=copyfields(obj.metadata,meta);
obj.metadata=copyfields(obj.metadata,m2);
obj.metadata.Width=double(omemeta.getPixelsSizeX(0).getValue());
obj.metadata.Height=double(omemeta.getPixelsSizeY(0).getValue());
obj.metadata.numberOfFrames=max(double(omemeta.getPixelsSizeT(0).getValue()),double(omemeta.getPixelsSizeZ(0).getValue()));
obj.metadata.basefile=[ph filesep fh];

fn=[fn {'Width','Height','numberOfFrames','basefile'}];

for k=1:length(fn)
    obj.metadata.assigned.(fn{k})=true;
end
metao=obj.metadata;
end

function meta=getMetaNd2(m)

%CameraUniqueName, Conversion Gain, Vertical Shift Speed, Readout Speed
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




function image=readstack(obj,imagenumber)
   if imagenumber<=obj.reader.getImageCount()
       image=bfGetPlane(obj.reader,imagenumber);
   else
       image=[];
   end
end