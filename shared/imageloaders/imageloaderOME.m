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

[ph,fh,ext]=fileparts(obj.file);
fn={};
switch ext
    case '.nd2' %Nikon
        meta=getMetaNd2(reader);
        fn=fieldnames(meta);
    case '.lif'
        meta=getMetaLif(reader);
        fn=fieldnames(meta);
        %determine series
        seri=reader.getSeriesCount;
        for k=seri:-1:1
            reader.setSeries(k-1);
            numim(k)=reader.getImageCount;
        end
        [~,largseries]=max(numim);
        reader.setSeries(largseries-1);
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
obj.metadata.roi=[0 0 obj.metadata.Width obj.metadata.Height];

fn=[makehorz(fn) makehorz({'Width','Height','numberOfFrames','basefile'})];

for k=1:length(fn)
    obj.metadata.assigned.(fn{k})=true;
end
metao=obj.metadata;
end

function md=getMetaLif(reader)

cm=reader.getCoreMetadataList;
cm1=cm.get(1);
sm=cm1.seriesMetadata;
md.emgain=str2double(sm.get('ProcessingHistory|ATLCameraSettingDefinition|EMGainValue'));
md.conversion=str2double(sm.get('ProcessingHistory|ATLCameraSettingDefinition|GainValue'));
md.EMon=str2double(sm.get('ProcessingHistory|ATLCameraSettingDefinition|CanDoEMGain'));
md.exposure=1000*str2double(sm.get('ProcessingHistory|ATLCameraSettingDefinition|WideFieldChannelConfigurator|SameExposureTime'));
md.timediff=md.exposure;

k=sm.keys;
ind=1;
while k.hasNext
    kh=k.nextElement;
    allmd{ind}=([kh ',' sm.get(kh)]);
    ind=ind+1;

end
md.allmetadata.omeLif=allmd;
end

function meta=getMetaNd2(reader)


m=reader.getGlobalMetadata;


    meta.EMon=str2num(m.get('EnableGainMultiplier'));
    meta.exposure= str2num(m.get('Exposure'));
    meta.emgain=str2num(m.get('GainMultiplier'));
    meta.conversion=str2num(m.get('ConversionGain'));
    meta.timediff=meta.exposure;
end




function image=readstack(obj,imagenumber)
   if imagenumber<=obj.reader.getImageCount()
       image=bfGetPlane(obj.reader,imagenumber);
   else
       image=[];
   end
end