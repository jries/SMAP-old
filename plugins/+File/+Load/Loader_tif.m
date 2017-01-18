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
        
        function pard=guidef(obj)
            pard=guidef;
        end
        function run(obj,p)
            [f,p]=uigetfile(obj.info.extensions);
            obj.load(p,[p f]);
            initGuiAfterLoad(obj);
        end
        function clear(obj,file,isadd)
                obj.locData.clear('filter');
        end        

    end
end




function pard=guidef
info.name='tif loader';
info.extensions={'*.tif';'*.*'};
info.dialogtitle='select any Tif file';
pard.plugininfo=info;
pard.plugininfo.type='LoaderPlugin';
end

function loadfile(obj,p,file,mode)  
if ~strcmp(mode,'tif')
   disp('no recognized Tiff file')
   return
end

if obj.locData.files.filenumberEnd==0
    [p,f]=fileparts(file);
    filename=[p filesep f];
   obj.locData.addfile(filename) 
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

  initGuiAfterLoad(obj)    
end
        


