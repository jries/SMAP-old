classdef Loader_tif<interfaces.DialogProcessor
    methods
        function obj=Loader_tif(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'mainGui','filelist_short'};
        end
        
        function out=load(obj,p,file,mode)
            if nargin<4
                mode=getfilemode(file);
            end
            p=obj.getAllParameters;
            loadfile(obj,p,file,mode);
        end
        
        function pard=pardef(obj)
            pard=pardef;
        end
        function run(obj,p)
            [f,p]=uigetfile(obj.info.extensions);
            obj.load(p,[p f]);
            initGuiAfterLoad(obj);
        end
        function empty(obj,file,isadd)
                obj.locData.empty('filter');
        end        

    end
end




function pard=pardef
info.name='tif loader';
info.extensions={'*.tif';'*.*'};
info.dialogtitle='select any Tif file';
pard.plugininfo=info;
end

function loadfile(obj,p,file,mode)  
if ~strcmp(mode,'tif')
   disp('no recognized Tiff file')
   return
end

if length(obj.locData.files.file)>1
    s=p.filelist_short.String;
    f=listdlg('ListString',s);
else
    f=1;
end
obj.locData.files.file(f).numberOfTif=obj.locData.files.file(f).numberOfTif+1;
imout=gettif(file);
imout.info.pixsize=obj.locData.files.file(f).info.pixsize;
obj.locData.files.file(f).tif(obj.locData.files.file(f).numberOfTif)=imout;

    
end
        

function imout=gettif(file)
imout.image=imread(file); 
sim=size(imout.image);
imout.info.Width=sim(1);
imout.info.Height=sim(2);
imout.info.roi=getRoiTif(file);
imout.info.name=file;
end
