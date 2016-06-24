classdef imageloaderMMsingle<interfaces.imageloaderSMAP
    %imageloaderMM image loader for micromanager single tiff files
    %   Detailed explanation goes here
    
    properties
        calfile='settings/CameraCalibration.xls';
        separate
        separatefiles
    end
    
    methods
        function obj=imageloaderMMsingle(varargin)
            obj@interfaces.imageloaderSMAP(varargin{:});
        end
        function open(obj,file)
            obj.file=file;
            md=obj.getmetadata;
            obj.separatefiles=md.allmetadata.files;
            obj.separate.numfiles=md.allmetadata.frames;
            obj.separate.path=md.allmetadata.path;
%             obj.currentImageNumber=0;
            obj.separate.fmt_s=imformats('tif');
        end
        function mdo=getmetadata(obj)
            mdo=getmetadataMM(obj);
             
        end
        function image=getimage(obj,frame)
            image=readseparate(obj,frame);
        end
    end
    
end

function image=readseparate(obj,number)
separate=obj.separate;
% lenfiles=length(separate.files);
lenfiles=separate.numfiles;
if lenfiles<number 
    if obj.onlineAnalysis %ask for image that is not in list
        lastfile= obj.separatefiles{lenfiles};
        thisname= generateFileName(lastfile,lenfiles,obj.metadata.allmetadata.numberNameRange,number);
        thisfile=[separate.path filesep thisname];
        if ~exist(thisfile,'file')
            disp('wait')
            pause(obj.waittime*2)
        end
        if ~exist(thisfile,'file')
            image=[];
            return
        else
            if number>lenfiles
                obj.separatefiles{number+1000}='';
            end
            obj.separatefiles{number}=thisname;
            obj.separate.numfiles=max(lenfiles,number);
            obj.metadata.numberOfFrames=max(lenfiles,number);
        end 
    else
        image=[];
        return
    end
else
    thisfile=[separate.path filesep obj.separatefiles{number}];
end
image=myimread(thisfile,separate.fmt_s);     
end

function newfile=generateFileName(oldfile, oldfilenumber,indfbar,newfilenumber)
oldfilenamenumber=str2double(oldfile(indfbar(1):indfbar(2)));
newfilenamenumber=oldfilenamenumber-oldfilenumber+newfilenumber;
lenfield=indfbar(2)-indfbar(1)+1;
newfilestr=num2str(newfilenamenumber,['%0' num2str(lenfield) 'i']);
newfile=[oldfile(1:indfbar(1)-1) newfilestr oldfile(indfbar(2)+1:end)];
end